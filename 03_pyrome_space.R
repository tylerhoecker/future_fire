library(tidyverse)
library(emdist)

# ------------------------------------------------------------------------------
# Load dataframes created by parts 01 and 02
# ------------------------------------------------------------------------------
pyrome_df <- readRDS('data/pyrome/pyrome_df.Rdata') %>% 
  # Filter out points where FRP was 0, this is due to the resampling/coarsening 
  # burned forest area... they are 0.01% of the data: sum(is.na(pyrome_df$frp))/length(pyrome_df$frp)
  filter(!is.na(frp))
forest_pyrome_df <- readRDS('data/pyrome/forest_pyrome_df.Rdata')

# ------------------------------------------------------------------------------
# Measure departure distance: Earth Mover Distance / Wasserstein Metric
# ------------------------------------------------------------------------------
# Use this function to rescale data from 0-1, so comparisons among dimensions make sense
rescale01 <- function(x){(x-min(x))/(max(x)-min(x))}

# Create a list of unique reference-future climate cell pairs (its 1723)
cell_pairs <- forest_pyrome_df %>% 
  dplyr::select(cell_hist, cell_2C) %>% 
  distinct() %>% 
  split(seq(nrow(.)))

# Add columns for all variables that are rescaled from 0-100
pyrome_df <- pyrome_df %>% 
  mutate(across(c(cbi, ndvi, frp), ~ rescale01(.x)*100, .names = 'z_{.col}'))

# Testing: index = 1; pair = cell_pairs[[index]]

calc_emds <- function(pair, index){
  
  pair = c(as.numeric(pair[1,1]),as.numeric(pair[1,2]))
  
  print(paste0('Cell pair ',index,': ',pair[1],',',pair[2]))
  
  # Skip pairs where reference and future are the same cell
  if(pair[1] == pair[2]){
    return(
      data.frame(
        'cell_hist'= pair[1],
        'cell_2C'  = pair[2],
        'cbi_emd'  = 0,
        'ndvi_emd' = 0,
        'frp_emd'  = 0,
        'cbi_dir'  = 0,
        'ndvi_dir' = 0,
        'frp_dir' = 0,
        'case' = 'No change'
      )
    )
  }
  
  # Create a dataframe of reference and future pyrome for that pairing
  ref_df <- filter(pyrome_df, cell_hist == pair[1]) 
  fut_df <-  filter(pyrome_df, cell_hist == pair[2]) 
  
  
  # Not enough points in the reference climate bin
  if(nrow(ref_df) == 0){
    return(
      data.frame(
        'cell_hist'= pair[1],
        'cell_2C'  = pair[2],
        'cbi_emd'  = NA,
        'ndvi_emd' = NA,
        'frp_emd'  = NA,
        'cbi_dir'  = NA,
        'ndvi_dir' = NA,
        'frp_dir'  = NA,
        'case'     = '0 reference'
      )
    )
  }
  
  # Not enough points in the future climate bin
  if(nrow(fut_df) == 0 & nrow(ref_df) != 0){
    return(
      data.frame(
        'cell_hist'= pair[1],
        'cell_2C'  = pair[2],
        'cbi_emd'  = NA,
        'ndvi_emd' = NA,
        'frp_emd'  = NA,
        'cbi_dir'  = NA,
        'ndvi_dir' = NA,
        'frp_dir'  = NA,
        'case'     = '0 future'
      )
    )
  }
  
  # Prepare the data to calculate EMD 
  prep_fun <- function(df,var){
    dat_re <- df %>% 
      pull(var) 
    dat_dens <- density(dat_re, bw = 0.05) 
    dat_mat <- cbind(dat_dens$x, dat_dens$y)
    return(list('dat_re' = dat_re, 'dat_mat' = dat_mat))
  }
  
  # CBI EMD
  ref_dat_cbi <- prep_fun(ref_df,'z_cbi')
  fut_dat_cbi <- prep_fun(fut_df, 'z_cbi')
  emd_cbi <- round(emdist::emd(ref_dat_cbi$dat_mat, fut_dat_cbi$dat_mat, 
                               max.iter = 10000), 5)
  
  # NDVI EMD
  ref_dat_ndvi <- prep_fun(ref_df,'z_ndvi')
  fut_dat_ndvi <- prep_fun(fut_df, 'z_ndvi')
  emd_ndvi <- round(emdist::emd(ref_dat_ndvi$dat_mat, fut_dat_ndvi$dat_mat, 
                                max.iter = 10000), 5)
  
  # FRP EMD
  ref_dat_frp <- prep_fun(ref_df,'z_frp')
  fut_dat_frp <- prep_fun(fut_df, 'z_frp')
  emd_frp <- round(emdist::emd(ref_dat_frp$dat_mat, fut_dat_frp$dat_mat, 
                               max.iter = 10000), 5)
  
  return(
    data.frame(
      'cell_hist'= pair[1],
      'cell_2C'  = pair[2],
      'cbi_emd'  = emd_cbi,
      'ndvi_emd' = emd_ndvi,
      'frp_emd'  = emd_frp,
      'cbi_dir'  = round(mean(fut_df$cbi) - mean(ref_df$cbi), 2),
      'ndvi_dir' = round(mean(fut_df$ndvi) - mean(ref_df$ndvi), 2),
      'frp_dir' =  round(mean(fut_df$frp) - mean(ref_df$frp), 2),
      'case'    = 'normal'
    )
  )
}

# Run the function over all pairs
emd_df <- map2_df(cell_pairs, seq_along(cell_pairs), calc_emds)
saveRDS(emd_df, 'data/emd/emd_df_2022_01_13.Rdata')

# Join to complete list of points (the west-wide sample)
pyrome_emd_df <- full_join(forest_pyrome_df, emd_df)

# ------------------------------------------------------------------------------
# Rasterize output 
# ------------------------------------------------------------------------------
# Convert these to a spatial objects (sf)
# Convert to spatial object (sf)
pyrome_emd_sf <- sf::st_as_sf(pyrome_emd_df, coords = c('x', 'y'), crs = 4326)

# Write out as shapefile
sf::st_write(pyrome_emd_sf, 'C:/Users/hoecker/GitHub/future_fire/data/pyrome/emd_Jan7.shp',
             append = FALSE)

# Read in raster of forest area
forest_rast <- raster('data/forest/forest_30m.tif')

# Then aggregate the fine forest_rast to 4 km (1/24th deg.) resolution
nw_coarse <- aggregate(forest_rast, fact = c((1/24)/xres(nw_forest),(1/24)/yres(nw_forest)), 
                       fun = mean, na.rm = T)
nw_coarse  <- reclassify(forest_rast, matrix(c(0, 0.1, NA,  0.1, 1, 1), ncol=3, byrow=TRUE))

# Rasterize and write out the TIFFs

list('cbi_emd','ndvi_emd','frp_emd','cbi_dir','ndvi_dir','frp_dir') %>%  
  walk(function(x) {
    print(paste0('Rasterizing ',x,'...'))
    r <- terra::rasterize(nw_emd_sf, nw_coarse, 
                          field = x, fun = max, na.rm = T)
    print(paste0('Writing data/pyrome/',x,'.tiff'))
    writeRaster(r, filename = paste0('data/pyrome/nw_',x,'.tiff'))
    return(r)
  }) 

# ------------------------------------------------------------------------------
# END
# ------------------------------------------------------------------------------

pyrome_emd_sf <- st_as_sf(pyrome_emd_df, coords = c('x', 'y'), crs = 4326)
nw_emd_sf <- st_as_sf(nw_emd_df, coords = c('x', 'y'), crs = 4326)
# Use a rasterized version of MTBS burned area perimeters to rasterize the EMD results
burned_rast <- raster('data/forest/burned_forest.tif')
burned_rast <- reclassify(burned_rast, cbind(0, 1))
burned_rast <- reclassify(burned_rast, cbind(NA, 0))

# Aggregate this to a raster with coarser resolution. MACA data resolution is 1/24th of a degree.
# Note this is a different resolution and result than a similar step in part 01
coarse_rast <- aggregate(burned_rast, fact = c((1/24)/xres(burned_rast),(1/24)/yres(burned_rast)), 
                         fun = max, na.rm = T)

# Same as above, for NW. Note that this step is the same as in part 01, but it's pretty quick to rerun.
# Extact climate data for ALL unburned forest points (at coarsened resolution)
nw_region <- st_read("C:/Users/hoecker/Work/GIS/doi_12_unified_regions_20180801_shapefile/region9_boundary.shp")
nw_region <- st_transform(nw_region, crs = 4326)

# Create a 4 km raster for forest area in the NW
# First crop forest_rast to NW
nw_forest <- crop(forest_rast, nw_region) #as_Spatial() 

list('cbi_emd','ndvi_emd','frp_emd','cbi_dir','ndvi_dir','frp_dir') %>%  
  walk(function(x) {
    print(paste0('Rasterizing ',x,'...'))
    r <- terra::rasterize(pyrome_emd_sf, coarse_rast, 
                          field = x, fun = max, na.rm = T)
    print(paste0('Writing data/pyrome/',x,'.tiff'))
    writeRaster(r, filename = paste0('data/pyrome/',x,'.tiff'))
    return(r)
  }) 

# Join to complete list of points (the west-wide sample and the complete NW df)
pyrome_emd_df <- full_join(pyrome_df, emd_df)

# At some point, account for forest points in NW will hexbins that didn't burn since 1984
nw_pyrome_df$cell_hist[setdiff(nw_pyrome_df$cell_hist, pyrome_df$cell_hist)]
nw_emd_df <- left_join(nw_pyrome_df, emd_df)

# Save dataframes for quick access later
# saveRDS(pyrome_emd_df, 'data/pyrome/pyrome_emd_df.Rdata')






# Start with the NW region
nw_region <- st_read("C:\\Users\\hoecker\\Work\\GIS\\doi_12_unified_regions_20180801_shapefile\\region9_boundary.shp")
nw_region <- st_transform(nw_region, crs = 4326)

# Read in EMD data
pyrome_emd_df <- readRDS('pyrome_emd_df.Rdata')

# Convert this to a spatial object (sf)
pyrome_emd_sf <- st_as_sf(pyrome_emd_df, coords = c('x', 'y'), crs = 4326)
pyrome_emd_nw <- st_crop(pyrome_emd_sf, nw_region)

# pyrome_emd_sp <- SpatialPointsDataFrame(coords = cbind(pyrome_emd_df$x,pyrome_emd_df$y), data = pyrome_emd_df,
#                                         proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))



# NEXT ------------------------------------------------------------------------
pyrome_list <- pyrome_df %>% 
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







# An example illustrating effect of differences in sample size

# This was used to test the impact of uneven sampling:
# Main point is that uneven is OK, as long as there are more than ~500 points 
test_samples <- function(eg_pt, ref_n = 5000, fut_n = 50, var = 'cbi', main_df){
  
  # Prepare the data to calculate EMD 
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
















# OLD

# old example space plots
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



