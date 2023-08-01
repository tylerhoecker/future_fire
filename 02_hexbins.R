# ------------------------------------------------------------------------------
# Description: Divides climate space of forest into a tessellation of hexagons, 
# with infomation about urn severity, productivity, fire occurrence, and climate 
# summarized within each.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load packages
# ------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(hexbin)

# ------------------------------------------------------------------------------
# Read in previously saved dataframe of fire, veg, climate dataframe
# ------------------------------------------------------------------------------
pyroclimate_input_df <- readRDS('data/process/pyroclimate_input_df.Rdata')
forest_clim_df <- readRDS('data/process/forest_clim_df.Rdata') %>% 
  as_tibble()

# ------------------------------------------------------------------------------
# Hexbins: delineate tesselation of hexbins based on climate axes
# ------------------------------------------------------------------------------
# Set minimum statistical support for pyromes
# Minimum number of observations (ie Landsat pixels)
sample_min <- 1000
# Minimum number of unique fire events observed
min_fires <- 10

# Number of bins in x direction - this is 'arbitary' but optimizes spatial footprint
# of climate bins and movement between bins under future scenario
n_bins <- 30

# Use same range of bins for burned and unburned areas
# Must span largest possible range of areas that will be considered
# For some reason all forest doesn't span a larger range that the burned points
# Probably due to finer resolution of burned points landing on additional TerraClimate cells
xrange <- range(c(forest_clim_df$def_hist,forest_clim_df$def_2C)) #,pyroclimate_input_df$def_hist,pyroclimate_input_df$def_2C
yrange <- range(c(forest_clim_df$aet_hist,forest_clim_df$aet_2C)) # ,pyroclimate_input_df$aet_hist,pyroclimate_input_df$aet_2C

# Size of bins in mm def
(xrange[2]-xrange[1])/n_bins
# Size of bins in mm aet
(yrange[2]-yrange[1])/n_bins

# Custom wrapper function to do this repeatedly, below
hex_custom <- function(x, y){
  hex <- hexbin(x, y, xbins= n_bins, xbnds = xrange, ybnds = yrange, IDs = TRUE)
  hex_df <- data.frame(hcell2xy(hex),  hexID = hex@cell, counts = hex@count)
  return(list('hex' = hex,
              'hex_df' = hex_df))
}

# Create for forest... just for plotting
hexdf_forest <- hex_custom(forest_clim_df$def_hist, forest_clim_df$aet_hist)[['hex_df']]

# Same for the burned points
hexdf_burned <-  hex_custom(pyroclimate_input_df$def_hist, pyroclimate_input_df$aet_hist)[['hex_df']] 

# Create same hexbins in the complete veg-fire df... 
pyrome_df <- pyroclimate_input_df %>% 
  mutate(cell_hist = hex_custom(def_hist, aet_hist)[['hex']]@cID,
         cell_2C = hex_custom(def_2C, aet_2C)[['hex']]@cID) %>% 
  # Identify whether or not points move groups between periods
  mutate(change = ifelse(cell_2C == cell_hist, 0, 1)) %>% 
  group_by(cell_hist) %>%
  mutate(sample_n = n()) %>% 
  mutate(fires = length(unique(fire.id))) %>% 
  ungroup() 

# Investigate the representativeness of the pyroclimates
pyrome_df %>% 
  group_by(cell_hist) %>% 
  summarize(fires = length(unique(fire.id))) %>% 
  summary()

n_fires <- pyrome_df %>% 
  group_by(cell_hist) %>% 
  summarize(fires = mean(fires))

pyrome_df <- pyrome_df %>% 
  # Some climate bins are poorly represented... remove these 
  # filter(sample_n >= sample_min) 
  filter(fires >= min_fires)

unique_frp <- pyrome_df %>%
  distinct(cell_hist,frp_pt,log_frp_pt)

# Same procedure for all forested points 
forest_pyrome_df <- forest_clim_df %>% 
  mutate(cell_hist = hex_custom(def_hist, aet_hist)[['hex']]@cID,
         cell_2C = hex_custom(def_2C, aet_2C)[['hex']]@cID) 

saveRDS(forest_pyrome_df, 'data/process/forest_pyrome_df.Rdata')

pyrome_df <- left_join(pyrome_df, unique_frp)
saveRDS(pyrome_df, 'data/process/pyrome_df.Rdata')

# Summarize NDVI and CBI by climate bin
pyrome_bins <- pyrome_df %>% 
  group_by(cell_hist) %>% 
  summarize(ndvi_bin = mean(ndvi),
            cbi_bin = mean(cbi),
            frp_bin = mean(frp_pt))

# Update this hexbin df used for plotting to reflect those filtered out
hexdf_burned_n <- hexdf_burned %>% 
  filter(hexID %in% pyrome_df$cell_hist) %>% 
  left_join(., pyrome_bins, by = c('hexID' = 'cell_hist')) %>% 
  left_join(., n_fires, by = c('hexID' = 'cell_hist')) 

# ------------------------------------------------------------------------------
# Plot climate space that is represented by these various geographies 
# (all forest, burned forest, burned forest with > min_sample points)
# ------------------------------------------------------------------------------
eg_data <- pyrome_df %>% 
  group_by(cell_hist) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  slice_sample(n = 35)

# This is the climate space that is represented...
ggplot() +
  geom_blank(data = hexdf_forest, aes(x = x, y = y)) +
  geom_hex(data = hexdf_forest, aes(x = x, y = y),
           color = 'grey50', fill = '#585858', stat = 'identity', size = .1) +
  geom_hex(data = hexdf_burned, aes(x = x, y = y),
           color = 'grey50', fill = '#ff0000', stat = 'identity', size = .1) + #alpha = 0.7,
  # geom_hex(data = hexdf_burned, aes(x = x, y = y, fill = test529),
  #          color = 'black', alpha = 1, stat = 'identity') +
  geom_hex(data = hexdf_burned_n, aes(x = x, y = y, fill = fires),
           size = .1, color = 'black', stat = 'identity') +
  #geom_text(data = hexdf_burned_n, aes(x = x, y = y)) +
  # geom_hex(data = hexdf_nw, aes(x = x, y = y),
  #          color = 'darkgreen', fill = 'transparent', alpha = 0.3, stat = 'identity') +
  # geom_point(data = eg_data,
  #              aes(x = def_hist, y = aet_hist), color = 'black', alpha = 0.7) +
  geom_curve(data = eg_data,
               aes(x = def_hist, y = aet_hist, xend = def_2C, yend = aet_2C),
               arrow = arrow(length=unit(0.15,"cm")), 
             color = 'white', curvature = -0.5, size = 0.4) +
  # ggrepel::geom_label_repel(data = eg_data, aes(x = def_hist, y = aet_hist, label = region), 
  #                           size = 2, box.padding = 0.05) +
  scale_fill_gradient('Fires observed', low = '#d70000', high = '#4c0000',
                      limits = c(0,150)) +
  # scale_fill_gradient2('Climate Hexbin',
  #                      low = '#8c510a', mid = '#f6e8c3', high = '#01665e',
  #                      midpoint = 500) +
  labs(y = 'Actual evapotranspiration (mm)', x = 'Climatic water deficit (mm)') +
  theme_bw(base_size = 8) +
  theme(legend.position = 'none',
        legend.background = element_blank(),
        axis.text = element_text(color = 'black'))

ggsave('climate_hex_burned.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs/',
       height = 3.1, width = 3.1, dpi = 500)

# ------------------------------------------------------------------------------
# Plot NDVI, CBI and FRP in climate space
# ------------------------------------------------------------------------------
# FRP
frp_hex <- 
  ggplot() +
  geom_blank(data = hexdf_forest, aes(x = x, y = y)) +
  geom_hex(data = hexdf_forest, aes(x = x, y = y),
           color = 'grey50', fill = '#585858', stat = 'identity') +
  # geom_hex(data = hexdf_burned, aes(x = x, y = y),
  #          color = 'grey50', fill = 'black', stat = 'identity') + #alpha = 0.7,
  geom_hex(data = hexdf_burned_n, aes(x = x, y = y, fill = frp_bin ),
           color = 'black', stat = 'identity') +
  geom_curve(data = eg_data,
             aes(x = def_hist, y = aet_hist, xend = def_2C, yend = aet_2C),
             arrow = arrow(length=unit(0.15,"cm")), 
             color = 'black', curvature = -0.5, size = 1.5) +
  scale_fill_distiller('Fire rotation \nperiod (yrs)', palette = 'RdPu', direction =1,
                       limits = c(50,250), oob = scales::squish) +
  labs(y = 'Actual evapotranspiration (mm)', x = 'Climatic water deficit (mm)') +
  theme_bw(base_size = 24) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = c(0.82,0.83),
        legend.background = element_blank(),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 24)) 
frp_hex

# CBI
cbi_hex <- 
  ggplot() +
  geom_blank(data = hexdf_forest, aes(x = x, y = y)) +
  geom_hex(data = hexdf_forest, aes(x = x, y = y),
           color = 'grey50', fill = '#585858', stat = 'identity') +
  # geom_hex(data = hexdf_burned, aes(x = x, y = y),
  #          color = 'grey50', fill = 'black', stat = 'identity') + #alpha = 0.7,
  geom_hex(data = hexdf_burned_n, aes(x = x, y = y, fill = cbi_bin),
           color = 'black', stat = 'identity') +
  geom_curve(data = eg_data,
             aes(x = def_hist, y = aet_hist, xend = def_2C, yend = aet_2C),
             arrow = arrow(length=unit(0.15,"cm")), 
             color = 'black', curvature = -0.5, size = 1.5) +
  scale_fill_distiller('Burn severity \n(CBI)', palette = 'OrRd', direction = 1,
                       limits = c(1.3,2), oob = scales::squish) +
  labs(y = 'Actual evapotranspiration (mm)', x = 'Climatic water deficit (mm)') +
  theme_bw(base_size = 24) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = c(0.82,0.83),
        legend.background = element_blank(),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 24)) 

# NDVI
ndvi_hex <- 
  ggplot() +
  geom_blank(data = hexdf_forest, aes(x = x, y = y)) +
  geom_hex(data = hexdf_forest, aes(x = x, y = y),
           color = 'grey50', fill = '#585858', stat = 'identity') +
  # geom_hex(data = hexdf_burned, aes(x = x, y = y),
  #          color = 'grey50', fill = 'black', stat = 'identity') + #alpha = 0.7,
  geom_hex(data = hexdf_burned_n, aes(x = x, y = y, fill = ndvi_bin),
           color = 'black', stat = 'identity') +
  geom_curve(data = eg_data,
             aes(x = def_hist, y = aet_hist, xend = def_2C, yend = aet_2C),
             arrow = arrow(length=unit(0.15,"cm")), 
             color = 'black', curvature = -0.5, size = 1.5) +
  scale_fill_distiller('Productivity \n(NDVI)', palette = 'YlGn', direction = 1,
                       limits = c(0.45,0.80), oob = scales::squish) +
  labs(y = 'Actual evapotranspiration (mm)', x = 'Climatic water deficit (mm)') +
  theme_bw(base_size = 24) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = c(0.82,0.83),
        legend.background = element_blank(),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 24)) 


library(cowplot)
panel_hex <- plot_grid(frp_hex, ndvi_hex, cbi_hex, nrow = 2, ncol = 2)
panel_hex
ggsave('three_panel_hex.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs/',
       height = 18, width = 18, dpi = 500)

# ------------------------------------------------------------------------------
# Mapping
# ------------------------------------------------------------------------------
# Plot on map
# State boundaries for reference
boundary <- read_sf("C:\\Users\\hoecker\\Work\\GIS\\cb_2018_us_state_20m\\cb_2018_us_state_20m.shp") %>%
  st_transform(., crs = 4326)
# 
# # Remake spatial version of veg-fire df now that it includes hexbin info
pyrome_sf <- st_as_sf(pyrome_df, coords = c('x', 'y'), crs = 4326)

plot_df <- pyrome_df %>% 
  slice_sample(n = 100000)

plot_sf <- st_as_sf(plot_df, coords = c('x', 'y'), crs = 4326)
eg_sf <- st_as_sf(eg_data, coords = c('x', 'y'), crs = 4326)

tm_shape(boundary, bbox = pyrome_sf) +
  tm_borders() +
  tm_fill(col = 'grey10') +
  tm_shape(plot_sf) +
  tm_dots(title = 'Climate Bin',
          col = 'cell_hist',
          size = 0.003,
          palette = c('#8c510a','#f6e8c3','#01665e'),
          style = 'cont',
          shape = 15) +
  tm_shape(eg_sf) +
  tm_dots(size = 0.2, col = 'grey90', border.col = 'red', shape = 21) +
  tm_text(text = 'region', col = 'black', auto.placement=F, just=c("left", "top"), 
          bg.color = 'white', bg.alpha = 0.6) +
  tm_layout(legend.outside = TRUE)

# Mysteriously starting throwing an error...
# tmap_save(filename = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/climate_map_egpts.png',
#           width = 6, height = 10, dpi = 500)

# ------------------------------------------------------------------------------
# Plotting the 'pyrome' space
# ------------------------------------------------------------------------------
pyrome_plot <- pyrome_df %>%
  filter(cell_hist %in% eg_data$cell_hist)

eg_labs <- paste0(eg_data$region, ' (', eg_data$cell_hist, ')')
names(eg_labs) <- eg_data$cell_hist

ggplot(pyrome_plot) +
  geom_hex(aes(x = ndvi, y = cbi, fill = ..ndensity..), color = 'black', bins = 15) +
  facet_wrap(~cell_hist, labeller = labeller(cell_hist = eg_labs)) +
  scale_fill_viridis_c(option = "plasma") +
  labs(x = 'Prefire NDVI', y = 'Modeled Composite Burn Index (CBI)') +
  theme_bw()

ggsave('cbi-ndvi_egpts.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
       height = 8, width = 9, unit = 'in', dpi = 300)

ggplot(pyrome_plot) +
  geom_hex(aes(x = fri, y = cbi, fill = ..ndensity..), color = 'black', bins = 15) +
  facet_wrap(~cell_hist, labeller = labeller(cell_hist = eg_labs)) +
  scale_fill_viridis_c(option = "plasma") +
  xlim(c(0,300)) +
  labs(x = 'Fire Rotaion Period (yrs)', y = 'Modeled Composite Burn Index (CBI)') +
  theme_bw()

ggsave('frp-cbi_egpts.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
       height = 8, width = 9, unit = 'in', dpi = 300)

ggplot(pyrome_plot) +
  geom_hex(aes(x = frs, y = cbi, fill = ..ndensity..), color = 'black', bins = 15) +
  facet_wrap(~cell_hist, labeller = labeller(cell_hist = eg_labs)) +
  scale_fill_viridis_c(option = "plasma") +
  #xlim(c(0,300)) +
  labs(x = 'Fire Resistance Score', y = 'Modeled Composite Burn Index (CBI)') +
  theme_bw()

ggplot(pyrome_plot) +
  geom_hex(aes(x = ndvi, y = fri, fill = ..ndensity..), color = 'black', bins = 15) +
  facet_wrap(~cell_hist, labeller = labeller(cell_hist = eg_labs)) +
  scale_fill_viridis_c(option = "plasma") +
  labs(x = 'Prefire NDVI', y = 'Fire Rotaion Period (yrs)') +
  theme_bw()

ggsave('ndvi-frp_egpts.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
       height = 8, width = 9, unit = 'in', dpi = 300)
# 
# ------------------------------------------------------------------------------
# Examples...
# ------------------------------------------------------------------------------
# Randomly select a sample point from each region
# set.seed(1234)
# eg_data <- pyrome_df %>% 
#   group_by(region) %>% 
#   slice_sample(n = 1)
# 
# hexdf_burned <- hexdf_burned %>%
#   mutate(test529 = ifelse(hexID == 529, 50, 1),
#          test529 = ifelse(hexID == 562, 100, test529))
# 
# test529_hex <- filter(hexdf_forest, hexID == 529)


# Plot the reference and future pyromes for the example points
eg_list <- eg_data %>%
  split(seq(nrow(.)))


eg_hex <- function(x, y){
  
  hex <- hexbin(x, y, xbins = 40, shape = 0.33, IDs = TRUE) #, xbnds = c(0,1), ybnds = c(0,3), 
  hex_df <- data.frame(hcell2xy(hex),  hexID = hex@cell, counts = hex@count)
  return(list('hex' = hex,
              'hex_df' = hex_df))
}

#eg = eg_list[[2]]

example_plots <- function(eg){

  print(eg$region)

  ref_df <- pyrome_df %>%
    filter(cell_hist == eg$cell_hist)
  
  fut_df <- pyrome_df %>%
    filter(cell_hist == eg$cell_2C)

  # In pyrome space
  # xdim = 'ndvi'
  # ydim = 'cbi'
  space_plot <- function(period_df, xdim, ydim){
    
    hex <- eg_hex(period_df[[xdim]], period_df[[ydim]])[['hex_df']]

    ggplot() +
      geom_hex(data = hex, aes(x = x, y = y, fill = counts), color = 'black', stat = 'identity') +
      stat_ellipse(aes(x = ref_df[[xdim]], y = ref_df[[ydim]]),
                   type = 'norm', level = 0.25,
                   color = 'white') +
      stat_ellipse(aes(x = ref_df[[xdim]], y = ref_df[[ydim]]),
                   type = 'norm', level = 0.50,
                   color = 'white') +
      stat_ellipse(aes(x = fut_df[[xdim]], y = fut_df[[ydim]]),
                   type = 'norm', level = 0.25, linetype = 2,
                   color = 'white') +
      stat_ellipse(aes(x = fut_df[[xdim]], y = fut_df[[ydim]]),
                   type = 'norm', level = 0.50, linetype = 2,
                   color = 'white') +
      scale_fill_viridis_c(option = "plasma", guide = 'none') +
      labs(x = toupper(xdim), y = toupper(ydim)) +
      #coord_cartesian(xlim = c(0.35,0.9), ylim = c(0,3)) +
      theme_bw() 
  }
  
  ref_cbi_fri <- space_plot(ref_df, 'cbi','frp')
  ref_ndvi_cbi <- space_plot(ref_df, 'ndvi','cbi')
  ref_ndvi_fri <- space_plot(ref_df,'ndvi','frp')
  fut_cbi_fri <- space_plot(fut_df, 'cbi','frp')
  fut_ndvi_cbi <- space_plot(fut_df, 'ndvi','cbi')
  fut_ndvi_fri <- space_plot(fut_df,'ndvi','frp')
  
  print(cowplot::plot_grid(ref_cbi_fri, ref_ndvi_fri, ref_ndvi_cbi,
                           fut_cbi_fri, fut_ndvi_fri, fut_ndvi_cbi,
                           nrow = 2, labels = eg$region))

  ggsave(paste0('ref_fut_eg', eg$region, '.png'),
                path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
                height = 5, width = 8, unit = 'in', dpi = 300)
}

walk(eg_list, example_plots)

# ------------------------------------------------------------------------------
# END
# ------------------------------------------------------------------------------


