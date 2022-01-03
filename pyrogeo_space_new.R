library(tidyverse)
library(raster)
library(sf)
library(hexbin)
library(tmap)
library(EMDomics)
library(emdist)

# ------------------------------------------------------------------------------
# Read in previously saved dataframe of fire, veg, climate dataframe
# ------------------------------------------------------------------------------
fire_veg_clim_df <- readRDS('data/fire_veg_clim_df.Rdata')
forest_climate <- readRDS('data/forest_climate.Rdata')

# ------------------------------------------------------------------------------
# Hexbins: delineate tesselation of hexbins based on climate axes
# ------------------------------------------------------------------------------
# Sample size - used later here for convenience - equal support among pyromes
sample_min <- 5000


# Number of bins in x direction
n_bins <- 30

# Use same range of bins for burned and unburned areas
xrange <- range(c(forest_climate$def_hist,forest_climate$def_2C))
yrange <- range(c(forest_climate$aet_hist,forest_climate$aet_2C))

# Custom wrapper function to do this repeatedly, below
hex_custom <- function(x, y){
  hex <- hexbin(x, y, xbins= n_bins, xbnds = xrange, ybnds = yrange, IDs = TRUE)
  hex_df <- data.frame(hcell2xy(hex),  hexID = hex@cell, counts = hex@count)
  return(list('hex' = hex,
              'hex_df' = hex_df))
}

# Create for forest
hexdf_forest <- hex_custom(forest_climate$def_hist, forest_climate$aet_hist)[['hex_df']]

# Same for the burned points
hexdf_burned <-  hex_custom(fire_veg_clim_df$def_hist, fire_veg_clim_df$aet_hist)[['hex_df']] 
hexdf_500 <- filter(hexdf_burned, counts > 500)
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
  filter(sample_n >= 500) 

# Some climate bins are poorly represented... remove these 
# This is written into EMD function, may not longer be necessary but leaving there as well for now 
pyrome_equal_df <- pyrome_df %>% 
  filter(sample_n > sample_min) 

length(unique(pyrome_equal_df$cell_hist))

# %>%
#   # Then, take a random sample of pyrome_n points from the remaining cells
#   slice_sample(n = pyrome_n) %>% 
#   ungroup()


# Randomly select a sample point from each region
set.seed(1234)
eg_data <- pyrome_df %>% 
  group_by(region) %>% 
  slice_sample(n = 1)

#This is the climate space that is represented...
ggplot() +
  geom_hex(data = hexdf_forest, aes(x = x, y = y),
           color = 'black', fill = 'grey70', stat = 'identity') +
  geom_hex(data = hexdf_burned, aes(x = x, y = y),
           color = 'black', fill = 'grey30', alpha = 0.7, stat = 'identity') +
  geom_hex(data = hexdf_500, aes(x = x, y = y),
           color = 'black', fill = 'grey5', alpha = 0.7, stat = 'identity') +
  geom_hex(data = hexdf_burned_n, aes(x = x, y = y, fill = hexID),
           color = 'black', alpha = 1, stat = 'identity') +
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
pyrome_plot <- pyrome_equal_df %>%
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
# 
# ------------------------------------------------------------------------------
# An example...
# ------------------------------------------------------------------------------

# Plot the reference and future pyromes for the example points
eg_list <- eg_data %>%
  split(seq(nrow(.)))


eg_hex <- function(x, y){
  hex <- hexbin(x, y, xbins= 40, xbnds = c(0,1), ybnds = c(0,3), shape = 0.33, IDs = TRUE)
  hex_df <- data.frame(hcell2xy(hex),  hexID = hex@cell, counts = hex@count)
  return(list('hex' = hex,
              'hex_df' = hex_df))
}

#eg = eg_list[[1]]

example_plots <- function(eg){

  print(eg$region)

  ref_df <- pyrome_equal_df %>%
    filter(cell_hist == eg$cell_hist)

  ref_hex <- eg_hex(ref_df$ndvi, ref_df$cbi)[['hex_df']]

  fut_df <- pyrome_equal_df %>%
    filter(cell_hist == eg$cell_2C)

  fut_hex <- eg_hex(fut_df$ndvi, fut_df$cbi)[['hex_df']]

  # In pyrome space
  ref_space <-
    ggplot(ref_hex) +
    geom_hex(aes(x = x, y = y, fill = counts), color = 'black', stat = 'identity') +
    geom_segment(data = eg, aes(x = ndvi, y = cbi,
                                xend = mean(fut_df$ndvi), yend = mean(fut_df$cbi)),
                 size = 0.7, color = 'white') +
    geom_segment(data = ref_df, aes(x = mean(ndvi), y = mean(cbi),
                                    xend = mean(fut_df$ndvi), yend = mean(fut_df$cbi)),
                 size = 0.7, color = 'white') +
    geom_point(data = eg, aes(x = ndvi, y = cbi),
               shape = 21, fill = 'white', color = 'black', size = 3, stroke = 1) +
    stat_ellipse(data = ref_df, aes(x = ndvi, y = cbi),
                 type = 'norm', level = 0.25,
                 color = 'white') +
    stat_ellipse(data = ref_df, aes(x = ndvi, y = cbi),
                 type = 'norm', level = 0.50,
                 color = 'white') +
    stat_ellipse(data = fut_df, aes(x = ndvi, y = cbi),
                 type = 'norm', level = 0.25, linetype = 2,
                 color = 'white') +
    stat_ellipse(data = fut_df, aes(x = ndvi, y = cbi),
                 type = 'norm', level = 0.50, linetype = 2,
                 color = 'white') +
    geom_point(data = fut_df, aes(x = mean(ndvi), y = mean(cbi)),
               color = 'white', shape = 3, size = 2) +
    geom_point(data = ref_df, aes(x = mean(ndvi), y = mean(cbi)),
               color = 'white', shape = 3, size = 3) +

    scale_fill_viridis_c(option = "plasma", guide = 'none') +
    coord_cartesian(xlim = c(0.35,0.9), ylim = c(0,3)) +
    theme_bw() +
    labs(x = 'Productivity (NDVI)', y = 'Fire severity (CBI)', title = 'Reference climate')

  future_space <-
    ggplot(fut_hex) +
    geom_hex(aes(x = x, y = y, fill = counts), color = 'black', stat = 'identity') +
    geom_segment(data = eg, aes(x = ndvi, y = cbi,
                                xend = mean(fut_df$ndvi), yend = mean(fut_df$cbi)),
                 size = 0.7, color = 'white') +
    geom_segment(data = ref_df, aes(x = mean(ndvi), y = mean(cbi),
                                    xend = mean(fut_df$ndvi), yend = mean(fut_df$cbi)),
                 size = 0.7, color = 'white') +
    geom_point(data = eg, aes(x = ndvi, y = cbi),
               shape = 21, fill = 'white', color = 'black', size = 3, stroke = 1) +
    stat_ellipse(data = ref_df, aes(x = ndvi, y = cbi),
                 type = 'norm', level = 0.25,
                 color = 'white') +
    stat_ellipse(data = ref_df, aes(x = ndvi, y = cbi),
                 type = 'norm', level = 0.50,
                 color = 'white') +
    stat_ellipse(data = fut_df, aes(x = ndvi, y = cbi),
                 type = 'norm', level = 0.25, linetype = 2,
                 color = 'white') +
    stat_ellipse(data = fut_df, aes(x = ndvi, y = cbi),
                 type = 'norm', level = 0.50, linetype = 2,
                    color = 'white') +
    geom_point(data = ref_df, aes(x = mean(ndvi), y = mean(cbi)),
               color = 'white', shape = 3, size = 3) +
    geom_point(data = fut_df, aes(x = mean(ndvi), y = mean(cbi)),
               color = 'white', shape = 3, size = 3) +
    scale_fill_viridis_c(option = "plasma", guide = 'none') +
    coord_cartesian(xlim = c(0.35,0.9), ylim = c(0,3)) +
    theme_bw() +
    labs(x = 'Productivity (NDVI)', y = '', title = 'Future climate')


  print(cowplot::plot_grid(ref_space, future_space,
                           nrow = 1, labels = eg$region))

  ggsave(paste0('ref_fut_eg', eg$region, '.png'),
                path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
                height = 5, width = 12, unit = 'in', dpi = 300)
}

walk(eg_list, example_plots)


# ------------------------------------------------------------------------------
# Measure departure distance: Earth Mover Distance / Wasserstein Metric
# ------------------------------------------------------------------------------
rescale01 <- function(x){(x-min(x))/(max(x)-min(x))}

# An example illustrating effect of differences in sample size
test_samples <- function(eg_pt, ref_n = 5000, fut_n = 50, var = 'cbi', main_df){
  
  prep_fun <- function(time){
    if(time == 'ref'){
      cell = 'cell_hist'
      n = ref_n 
    } else {
      cell = 'cell_2C'
      n = fut_n
    }
    dat_re <- main_df %>% 
      filter(cell_hist == eg_pt[[cell]]) %>% 
      dplyr::select(ndvi,cbi) %>% 
      slice_sample(n = n) %>% 
      pull(var) %>% 
      rescale01() 
    dat_dens <- density(dat_re, bw = 0.05) 
    dat_mat <- cbind(dat_dens$x, dat_dens$y)
    return(list('dat_re' = dat_re, 'dat_mat' = dat_mat))
  }
  
  ref_dat <- prep_fun(time = 'ref')
  fut_dat <- prep_fun(time = '2C')
  
  emd_result <- round(emdist::emd(ref_dat$dat_mat,fut_dat$dat_mat, max.iter = 5000), 3)
  
  if(ref_n == fut_n){sam_lab = 'equal'}else{sam_lab = 'unequal'}
  
  png(filename = paste0('C:/Users/hoecker/Work/Postdoc/Future_Fire/Sample_Figs/',
                        eg_pt$region,'_',sam_lab,'_samp_hist.png'),
      height = 5, width = 8, units = 'in', res = 200)
  par(mfrow = c(1,2))
  hist(ref_dat$dat_re, breaks = 20, xlab = 'Rescaled variable', 
       main = paste(eg_pt$region, 'EMD = ', emd_result))
  hist(fut_dat$dat_re, col = 'grey10', add = T, breaks = 20)
  hist(ref_dat$dat_re, freq = F, breaks = 20, xlab = 'Rescaled variable', main = '')
  hist(fut_dat$dat_re, col = 'grey10', add = T, , freq = F, breaks = 20)
  dev.off()
  
  return(data.frame('region' = eg_pt$region,
                    'emd' = emd_result))
}

unequal_df <- map_df(eg_list, test_samples, ref_n = 5000, fut_n = 500, var = 'cbi', main_df = pyrome_df) %>% 
  rename(emd_unequal = emd)
equal_df <- map_df(eg_list, test_samples, ref_n = 5000, fut_n = 5000, var = 'cbi', main_df = pyrome_df) %>% 
  rename(emd_equal = emd)

comparison_df <- full_join(unequal_df, equal_df)


ggplot(comparison_df, aes(x = emd_unequal, y = emd_equal)) +
  geom_point() +
  geom_smooth(method = 'lm')

# Functions to compare the reference and +2C pyrome space for every point

# First, go through all reference-future climate cell pairs and trim down to 
# a list of unique pairs
cell_pairs <- pyrome_df %>% 
  dplyr::select(cell_hist, cell_2C) %>% 
  distinct() %>% 
  split(seq(nrow(.)))

pair = cell_pairs[[5]]

calc_emds <- function(pair, index){
  
  pair = c(as.numeric(pair[1,2]), as.numeric(pair[1,1]))
  
  print(paste('Cell pair:',index,'-',pair[1],pair[2]))
  
  # Skip pairs where reference and future are the same cell
  if(pair[1] == pair[2]){
    return(
      data.frame(
        'cell_hist'= pair[1],
        'cell_2C'  = pair[2],
        'cbi_emd'  = 0,
        'ndvi_emd' = 0,
        'cbi_dir'  = 0,
        'ndvi_dir' = 0 
      )
    )
  }
  
  # Create a dataframe of reference and future pyrome for that point
  ref_df <- filter(pyrome_equal_df, cell_hist == pair[1]) 
  fut_df <-  filter(pyrome_equal_df, cell_hist == pair[2]) 
  
  # Not enough points in the reference climate bin
  if(nrow(ref_df) == 0){
    return(
      data.frame(
        'cell_hist'= pair[1],
        'cell_2C'  = pair[2],
        'cbi_emd'  = NA,
        'ndvi_emd' = NA,
        'cbi_dir'  = NA,
        'ndvi_dir' = NA,
        'novel'    = '0 reference'
        #'fri_emd'  = NA
      )
    )
  }
  
  # Not enough points in the future climate bin
  if(nrow(fut_df) == 0){
    return(
      data.frame(
        'cell_hist'= pair[1],
        'cell_2C'  = pair[2],
        'cbi_emd'  = NA,
        'ndvi_emd' = NA,
        'cbi_dir'  = NA,
        'ndvi_dir' = NA,
        'novel'    = '0 future'
        #'fri_emd'  = NA
      )
    )
  }
  
 # Calculate the Earth Mover's Distance metric 
  prep_mat_fun <- function(df, var){
    dens <- df %>% 
      dplyr::select(var) %>% 
      pull(var) %>% 
      rescale01() %>% 
      density(bw = 0.05) 
    mat <- cbind(dens$x, dens$y)
    return(mat)
  }
  
  
  



  
  return(
    data.frame(
      'cell_hist'= pair[1],
      'cell_2C'  = pair[2],
      'cbi_emd'  = results$emd['cbi','emd'],
      'ndvi_emd' = results$emd['ndvi','emd'],
      'cbi_dir'  = mean(fut_df$cbi) - mean(ref_df$cbi),
      'ndvi_dir' = mean(fut_df$ndvi) - mean(ref_df$ndvi),
      'novel'    = 'No'
    )
  )
  
}

emd_df <- map2_df(cell_pairs, seq_along(cell_pairs), calc_emds)


pyrome_list <- pyrome_equal_df %>% 
  split(seq(nrow(.)))

calc_pt_dirs <- function(pt_df){
  
  print(pt_df$ID)
  
  fut_df <-  filter(pyrome_equal_df, cell_hist == pt_df$cell_2C) 
  
  # Skip pairs where reference and future are the same cell
  if(pt_df$cell_hist == pt_df$cell_2C){
    return(
      data.frame(
        'ID' = pt_df$ID,
        'cell_hist'= pt_df$cell_hist,
        'cell_2C'  = pt_df$cell_2C,
        'cbi_pt_dir'  = 0,
        'ndvi_pt_dir' = 0)
    )
  }
  
  # Not enough points in the reference or future climate bins
  if(nrow(fut_df) == 0){
    return(
      data.frame(
        'ID' = pt_df$ID,
        'cell_hist'= pt_df$cell_hist,
        'cell_2C'  = pt_df$cell_2C,
        'cbi_pt_dir'  = NA,
        'ndvi_pt_dir' = NA)
    )
  }
  
  return(
    data.frame(
      'ID' = pt_df$ID,
      'cell_hist'= pt_df$cell_hist,
      'cell_2C'  = pt_df$cell_2C,
      'cbi_pt_dir'  = mean(fut_df$cbi) - pt_df$cbi,
      'ndvi_pt_dir' = mean(fut_df$ndvi) - pt_df$ndvi)
  )
  
  
}

pt_dirs_df <- map_df(pyrome_list, calc_pt_dirs)

# Not sure why these columns made it in here...
emd_pt_df <- full_join(emd_df, pt_dirs_df)
pyrome_emd <- full_join(pyrome_sf, emd_pt_df)
saveRDS(pyrome_emd, 'C:/Users/hoecker/GitHub/future_fire/data/pyrome/pyrome_emd_equalN.Rdata')
st_write(pyrome_emd, 'C:/Users/hoecker/GitHub/future_fire/data/pyrome/pyrome_emd_equalN.shp', 
         append = FALSE)

























# OLD
# Comparing sample sizes
pyrome_emd_min_1000 <- readRDS('C:/Users/hoecker/GitHub/future_fire/data/pyrome/pyrome_emd_min_1000.Rdata')
pyrome_emd_min_10 <- readRDS('C:/Users/hoecker/GitHub/future_fire/data/pyrome/pyrome_emd_min_10.Rdata')

pyrome_emd_1000 <- pyrome_emd_min_1000 %>% 
  st_drop_geometry()

pyrome_emd_10 <- pyrome_emd_min_10 %>% 
  st_drop_geometry()

cbi_1000 <- 
ggplot(pyrome_emd_1000) +
  geom_histogram(aes(x = cbi_emd, y = ..count../1000), 
                 binwidth = 0.1, fill = 'black') +
  coord_cartesian(xlim = c(0,5)) +
  theme_bw() +
  ggtitle('Min. n = 1000')

cbi_10 <- 
  ggplot(pyrome_emd_10) +
  geom_histogram(aes(x = cbi_emd, y = ..count../1000), 
                 binwidth = 0.1, fill = 'black') +
  coord_cartesian(xlim = c(0,5)) +
  theme_bw() +
  ggtitle('Min. n = 10')

cowplot::plot_grid(cbi_1000, cbi_10)

fri_hist <- 
ggplot(pyrome_emd_df) +
  geom_histogram(aes(x = fri_emd, y = ..density..), fill = 'black') +
  #facet_wrap(~region) +
  coord_cartesian(xlim = c(0,6)) +
  theme_bw()

ndvi_hist <- 
ggplot(pyrome_emd_df) +
  geom_histogram(aes(x = ndvi_emd, y = ..density..), fill = 'black') +
  #facet_wrap(~region) +
  coord_cartesian(xlim = c(0,6)) +
  theme_bw()

cowplot::plot_grid(cbi_hist, fri_hist, ndvi_hist, nrow = 1)



plot_test <- pyrome_emd %>% 
  group_by(region) %>% 
  slice_head(prop = .01) 

tm_shape(boundary, bbox = plot_test) +
  tm_borders() +
  tm_fill(col = 'grey10') +
  tm_shape(plot_test) +
  tm_dots(title = 'Departure (EMD)',
          col = 'cbi_dir',
          size = 0.001, 
          palette = "plasma",
          midpoint = 0,
          #breaks = c(0,0.25,0.5,1,1.5,2,2.5,3,4),
          colorNULL = "#F0F921FF",
          style = 'cont',
          shape = 15) +
  tm_layout(legend.position = c('left','bottom'))

tmap_save(filename = 'depature_map.png', width = 7, height = 10, unit = 'in', dpi = 600)


ggplot(pyrome_emd) +
  geom_point(aes(x = x, y = y, color = cbi_emd))























# Old version of this:
# Create same hexbins in the complete veg-fire df... 
pyrome_df <- fire_veg_clim_df %>% 
  mutate(cell = hexbin(x = def, y = aet, 
                       xbins= n_bins, 
                       xbnds = range(forest_climate$def), ybnds = range(forest_climate$aet),
                       IDs = TRUE)@cID) %>% 
  pivot_wider(names_from = period, values_from = c(aet, def, cell))
# Do stuff here to identify whether or not points move groups between periods
group_by(ID) %>% 
  mutate(change = as.factor(ifelse(length(unique(cell)) > 1, 1, 0))) %>% 
  ungroup() %>%
  # Some climate bins are poorly represented... remove these 
  # This is written into EMD function, may not longer be necessary but leaving there as well for now 
  group_by(cell) %>%
  mutate(sample = n()) %>%
  filter(sample > 1000) %>%
  # Then, take a random sample of 1000 points from the remaining cells
  slice_sample(n = 1000) %>% 
  ungroup()

fire_veg_clim_df <- cbind(fire_veg_df, climate_pts) %>% 
  #sample_n(100) %>% 
  filter(aet_hist > 0) %>% 
  ungroup() %>% 
  pivot_longer(cols = c(aet_2C, aet_hist, def_2C, def_hist), 
               names_to = 'temp_name', values_to = 'value') %>% 
  #dplyr::select(-ID) %>% 
  separate(temp_name, into = c('variable','period'), sep = '_') %>% 
  pivot_wider(names_from = variable, values_from = value) 

# Do the same for all forest, to compare the space with burned forest
# A shapefile of 95,000 unburned forest points I created at some point... 
forest_sf <- st_read('data/forest/unburned.shp')

# Extract climate data
forest_climate <- extract(climate_stack, forest_sf, df = T) %>% 
  filter(aet_hist > 0) %>% 
  ungroup() %>% 
  pivot_longer(cols = c(aet_2C, aet_hist, def_2C, def_hist), 
               names_to = 'temp_name', values_to = 'value') %>% 
  #dplyr::select(-ID) %>% 
  separate(temp_name, into = c('variable','period'), sep = '_') %>% 
  pivot_wider(names_from = variable, values_from = value) 


# ggplot(plot_df, aes(x = def, y = aet)) +
#   geom_hex(data = hexdf_forest, aes(x = x, y = y), 
#            color = 'black', fill = 'grey70', stat = 'identity') +
#   geom_hex(data = hexdf_burned, aes(x = x, y = y, fill = hexID), color = 'black', alpha = 0.6, stat = 'identity') +
#   geom_point(aes(alpha = change), size = 1) +
#   geom_line(data = filter(plot_df, change == 0),
#             aes(group = ID, alpha = change), size = 0.75) +
#   geom_line(data = filter(plot_df, change == 1),
#             aes(group = ID, alpha = change), size = 0.75,
#             arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "open")) +
#   scale_fill_gradient2('Climate Hexbin', low = '#8c510a', mid = '#f6e8c3', high = '#01665e', midpoint = 250) +
#   #scale_color_manual('+2C group shift?', values = c('white','black'), labels = c('No','Yes')) +
#   scale_alpha_manual('+2C group shift?', values = c(0.3,0.8), labels = c('No','Yes')) +
#   labs(y = 'AET', x = 'DEF') +
#   theme_bw()
# 
# ggsave('hex_move_eg.png',path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
#        height = 6, width = 8, dpi = 300)





# for_template <- raster('data/forest/forest_na.tif')
# pyrome_emd_rast <- rasterize(pyrome_emd, for_template, field = 'cbi_emd')


# ------------------------------------------------------------------------------
# For later, incorporate/sample unburned points and their NDVI
# ------------------------------------------------------------------------------

# Read in burned forest raster
# This was created in QGIS using the 'GDAL Rasterize (Overwrite with Fixed Value)' tool 
# This fixed value was 0, such that the forest mask was updated to exclude burned areas
# Original forest mask unchanged
unburned_rast <- raster('data/forest/unburned_na.tif')

#unburned_rast <- reclassify(unburned_rast, cbind(0, NA))
#writeRaster(unburned_rast, 'data/forest/unburned_na.tif')

# Pick random sample equal in size to burned cells from non-burned cells
unburned_samp <- sampleRandom(x = unburned_rast, 
                              size = length(fire_veg_df[['x']]), # Update eventually with the full thing?
                              na.rm = TRUE,
                              xy = TRUE)

# plot(unburned_rast)
# points(unburned_samp)

# Make unburned sample a spatial sf
unburned_sf <- st_as_sf(as.data.frame(unburned_samp), 
                        coords = c('x', 'y'), 
                        crs = 4326)

# Get region for unburned points
ecoreg_sf <- st_read('data/ecoregions/ecoregions_edc.shp') %>% 
  st_transform(., crs = 4326)

unburned_sf <- unburned_sf %>% 
  st_join(., dplyr::select(ecoreg_sf, Region_ID)) %>% 
  filter(!is.na(Region_ID))

# Get NDVI for unburned points
ndvi <- list.files('data/ndvi/', full.names = T) %>% 
  map(., raster) 

ndvi_df <- map(ndvi, extract, unburned_sf, df = T, sp = T) %>% 
  map_df(., as.data.frame) %>% 
  pivot_longer(cols = starts_with('NDVI_'), names_to = 'name', values_to = 'ndvi') %>% 
  filter(!is.na(ndvi)) %>% 
  mutate(region = paste0('r', Region_ID),
         year = 2019,
         cbi = 0) %>% 
  dplyr::select(region, x = coords.x1, y = coords.x2, year, cbi, ndvi)

burn_unburn_df <- rbind(fire_veg_df, ndvi_df)


# ------------------------------------------------------------------------------
# Calculate FRP at the ecoregion level
# ------------------------------------------------------------------------------
# Calculate FRP for each ecoregion
ecoreg_area <- ecoreg_sf %>% 
  mutate(area = st_area(.),
         region = paste0('r', Region_ID)) %>% 
  dplyr::select(region, area) %>% 
  st_set_geometry(NULL)

burned_area <- fire_veg_df %>% 
  group_by(region) %>% 
  # Each pixel is 30x30 m * the proportion of the dataframe that was subsampled * the subsample that Sean provided
  summarise(burned_area = (n()*(30*30)*(1/subprop)*20) )  

# Approximate fire return interval, as the inverse of FRP 
fri_df <- full_join(ecoreg_area, burned_area) %>% 
  mutate(burned_area = ifelse(is.na(burned_area), 0, burned_area),
         fri = 1/(as.numeric(burned_area/area)/35),
         fri = ifelse(fri == 'Inf', 5000, fri)) %>% 
  dplyr::select(region, fri)

complete_df <- full_join(burn_unburn_df, fri_df)

# forest_rast <- reclassify(forest_rast, cbind(0,NA))
# writeRaster(forest_rast, 'data/forest/forest_na.tif')











# Old/Extra
# Loop through each ecoregion...
ecoreg_sf <- st_read('data/ecoregions/ecoregions_edc.shp') %>% 
  st_transform(., crs = 4326)

# Temp testing
ecoreg_i <- ecoreg_sf[6,] 

# FUNCTION WILL START HERE...
# Mask/crop the forest_rast to the ecoregion
forest_er <- crop(forest_rast, ecoreg_i)


# Calculate proportion of the ecoregion that burned out of total
length(burned_cells)/
  length(fire_veg_df_full[['x']])*100/length(forest_sf)

# Add FRP column to fire dataframe

# Mask the forest-ecoregion area to only the unburned area (mask out burned cells)
burned_cells <- cellFromXY(forest_er, cbind(fire_veg_df_full[,c('x','y')]))

# Randomly sample points from the unburned area
burned_er <- forest_er
burned_er[burned_cells] <- NA



# Crop NDVI to the ecoregion

# Extract NDVI for the random sample of points

# Save the unburned data:
# region, x, y, fire.id = NA, fire.year = 19, cbi = 0, ndvi = ndvi

# Join with fire dataframe




# This is the proportion of forested area that has burned... assuming Sean gave me a 1% subset of the data
length(fire_veg_df_full[['x']])*100/length(forest_sf)



plot(forest_sf)

# Need ecoregion shapefiles and area
ecoreg_sf <- st_read('data/ecoregions/ecoregions_edc.shp')


#rm(fire_veg_df_full)


# Mask out burned cells
# burned_cells <- cellFromXY(forest_rast, cbind(fire_veg_df_full[,c('x','y')]))
# unburned_rast <- forest_rast
# unburned_rast[burned_cells] <- NA # This is pretty slow...

# Read in burned areas perimeters
mtbs_perims <- st_read("C:\\Users\\hoecker\\Work\\GIS\\fire_data\\mtbs_1984_2018\\west_mtbs_perims.shp") %>% 
  st_transform(crs = 4326)

# Rasterize burned areas to the same grid as the forest mask, make the values 0
burned_forest <- fasterize::fasterize(mtbs_perims, forest_rast)









# Pick some representative bins from each region to explore. Don't have future climate data, 
# so these are just a hypothetical selection, ie, bin 88 historical and 33 future
hist_bin <- 19
future_bin <- 8

fire_veg_test <- fire_veg_df %>% 
  filter(cell == hist_bin | cell == future_bin) %>% 
  filter(cbi > 0) %>% 
  group_by(cell) %>% 
  slice_sample(., n = 10000)

ggplot(fire_veg_test, aes(x = ndvi, y = cbi, color = as.factor(cell), fill = as.factor(cell))) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(group = as.factor(cell))) +
  scale_color_manual('Climate',
                     values = c('#8c510a','#01665e'), 
                     labels = c('Future (8)','Reference (19)')) +
  scale_fill_manual('Climate',
                     values = c('#8c510a','#01665e'), 
                     labels = c('Future (8)','Reference (19)')) +
  theme_bw() +
  labs(x = 'NDVI', y = 'CBI')



