library(tidyverse)
library(sf)
library(hexbin)
library(tmap)

# ------------------------------------------------------------------------------
# Read in previously saved dataframe of fire, veg, climate dataframe
# ------------------------------------------------------------------------------
fire_veg_clim_df <- readRDS('data/fire_veg_clim_df.Rdata')
forest_clim_df <- readRDS('data/forest_clim_df.Rdata')

# ------------------------------------------------------------------------------
# Hexbins: delineate tesselation of hexbins based on climate axes
# ------------------------------------------------------------------------------
# Sample size - used later here for convenience - equal support among pyromes
sample_min <- 1000

# Number of bins in x direction
n_bins <- 30

# Use same range of bins for burned and unburned areas
# Must span largest possible range of areas that will be considered
xrange <- range(c(forest_clim_df$def_hist,forest_clim_df$def_2C))
yrange <- range(c(forest_clim_df$aet_hist,forest_clim_df$aet_2C))

# Custom wrapper function to do this repeatedly, below
hex_custom <- function(x, y){
  hex <- hexbin(x, y, xbins= n_bins, xbnds = xrange, ybnds = yrange, IDs = TRUE)
  hex_df <- data.frame(hcell2xy(hex),  hexID = hex@cell, counts = hex@count)
  return(list('hex' = hex,
              'hex_df' = hex_df))
}

# Create for forest
hexdf_forest <- hex_custom(forest_clim_df$def_hist, forest_clim_df$aet_hist)[['hex_df']]

# Delete...
# Create for nw forest 
#hexdf_nw <- hex_custom(nw_clim_df$def_hist, nw_clim_df$aet_hist)[['hex_df']] 

# Same for the burned points
hexdf_burned <-  hex_custom(fire_veg_clim_df$def_hist, fire_veg_clim_df$aet_hist)[['hex_df']] 
# Filter out bins with less than minimum sample size
hexdf_burned_n <- filter(hexdf_burned, counts > sample_min)
  

# Create same hexbins in the complete veg-fire df... 
pyrome_df <- fire_veg_clim_df %>% 
  mutate(cell_hist = hex_custom(def_hist, aet_hist)[['hex']]@cID,
         cell_2C = hex_custom(def_2C, aet_2C)[['hex']]@cID) %>% 
  # Identify whether or not points move groups between periods
  mutate(change = ifelse(cell_2C == cell_hist, 0, 1)) %>% 
  group_by(cell_hist) %>%
  mutate(sample_n = n()) %>% 
  ungroup() %>% 
  # Some climate bins are poorly represented... remove these 
  filter(sample_n >= sample_min) 

# Same procedure for all forested points in the NW
forest_pyrome_df <- forest_clim_df %>% 
  mutate(cell_hist = hex_custom(def_hist, aet_hist)[['hex']]@cID,
         cell_2C = hex_custom(def_2C, aet_2C)[['hex']]@cID) %>% 
  # Identify whether or not points move groups between periods
  mutate(change = ifelse(cell_2C == cell_hist, 0, 1)) 

# Randomly select a sample point from each region
set.seed(1234)
eg_data <- pyrome_df %>% 
  group_by(region) %>% 
  slice_sample(n = 1)

# ------------------------------------------------------------------------------
# Save data that will be needed in part 03
# ------------------------------------------------------------------------------
saveRDS(pyrome_df, 'data/pyrome/pyrome_df.Rdata')
saveRDS(forest_pyrome_df, 'data/pyrome/forest_pyrome_df.Rdata')


# ------------------------------------------------------------------------------
# Plot climate space that is represented by these various geographies 
# (all forest, burned forest, burned forest with > min_sample points)
# ------------------------------------------------------------------------------
# This is the climate space that is represented...
ggplot() +
  geom_hex(data = hexdf_forest, aes(x = x, y = y),
           color = 'black', fill = 'grey80', stat = 'identity') +
  geom_hex(data = hexdf_burned, aes(x = x, y = y),
           color = 'black', fill = 'grey30', stat = 'identity') + #alpha = 0.7, 
  geom_hex(data = hexdf_burned_n, aes(x = x, y = y, fill = hexID),
           color = 'black', alpha = 1, stat = 'identity') +
  # geom_hex(data = hexdf_nw, aes(x = x, y = y),
  #          color = 'darkgreen', fill = 'transparent', alpha = 0.3, stat = 'identity') +
  geom_point(data = eg_data,
               aes(x = def_hist, y = aet_hist)) +
  geom_curve(data = eg_data,
               aes(x = def_hist, y = aet_hist, xend = def_2C, yend = aet_2C),
               arrow = arrow(length=unit(0.20,"cm")), curvature = -0.5) +
  ggrepel::geom_label_repel(data = eg_data, aes(x = def_hist, y = aet_hist, label = region), 
                            size = 2, box.padding = 0.05) +
  scale_fill_gradient2('Climate Hexbin', 
                       low = '#8c510a', mid = '#f6e8c3', high = '#01665e', 
                       midpoint = 500) +
  labs(y = 'AET', x = 'DEF') +
  theme_bw() +
  theme(legend.position = c(0.86,0.82),
        legend.background = element_blank())

ggsave('whole_climate_hex.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
       height = 6, width = 7, dpi = 300)

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
  theme_bw()

ggsave('cbi-ndvi_egpts.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
       height = 8, width = 9, unit = 'in', dpi = 300)

ggplot(pyrome_plot) +
  geom_hex(aes(x = frp, y = cbi, fill = ..ndensity..), color = 'black', bins = 15) +
  facet_wrap(~cell_hist, labeller = labeller(cell_hist = eg_labs)) +
  scale_fill_viridis_c(option = "plasma") +
  theme_bw()

ggsave('frp-cbi_egpts.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
       height = 8, width = 9, unit = 'in', dpi = 300)

ggplot(pyrome_plot) +
  geom_hex(aes(x = ndvi, y = frp, fill = ..ndensity..), color = 'black', bins = 15) +
  facet_wrap(~cell_hist, labeller = labeller(cell_hist = eg_labs)) +
  scale_fill_viridis_c(option = "plasma") +
  theme_bw()

ggsave('ndvi-frp_egpts.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
       height = 8, width = 9, unit = 'in', dpi = 300)
# 
# ------------------------------------------------------------------------------
# An example...
# ------------------------------------------------------------------------------

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
                height = 5, width = 12, unit = 'in', dpi = 300)
}

walk(eg_list, example_plots)