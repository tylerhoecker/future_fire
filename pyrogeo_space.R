library(tidyverse)
library(raster)
library(sf)
library(hexbin)
library(tmap)

# Read in processed fire, veg data from Sean
fire_veg_df_full <- list('r19sev.pixel.data.through.2019.csv',
                    'r2sev.pixel.data.through.2019.csv') %>% 
  set_names(c('norock','norcal')) %>% 
  map_df(., read_csv, .id = 'region') 

fire_veg_df <- fire_veg_df_full 

# Read in AET and DEF TerraClim data for historical period
clim_hist <- list('aet.1981.2010.tif', 'def.1981.2010.tif') %>% 
  set_names(c('aet','def')) %>% 
  stack()

# Add climate data as columns to dataframe
# Make dataframe spatial first with coordinates
fire_veg_sf <- st_as_sf(fire_veg_df, coords = c('x', 'y'), crs = 4326)
# Extract values from climate rasters
fire_veg_df <- cbind(fire_veg_df, extract(clim_hist, fire_veg_sf))


# Number of bins in x direction
n_bins <- 5
#xrange <- c(0,200)
#yrange <- c(0,110)

# Create hexbins and save the hexbin info as it's own df
h <- hexbin(x = fire_veg_df$def, y = fire_veg_df$aet, 
            xbins= n_bins, 
            shape = 1, 
            IDs = TRUE)
hexdf <- data.frame(hcell2xy (h),  hexID = h@cell, counts = h@count)

# Create same hexbins in the complete veg-fire df
fire_veg_df <- fire_veg_df %>% 
  mutate(cell = hexbin(x = def, y = aet, 
                       xbins= n_bins, 
                       shape = 1, 
                       IDs = TRUE)@cID) 

# Create a smaller version for plotting etc
fire_veg_small <- fire_veg_df %>% 
  group_by(region) %>% 
  slice_sample(., n = 10000)

# Plot in bivariate climate space 
ggplot(fire_veg_test, aes(x = def, y = aet)) +
  geom_hex(data = hexdf, aes(x = x, y = y, fill = hexID), color = 'black', alpha = 0.8, stat = 'identity') +
  geom_text(data = hexdf, aes(x = x, y = y, label = hexID)) +
  geom_point(alpha = 0.1, size = 1) +
  scale_fill_gradient2('', low = '#8c510a', mid = '#f6e8c3', high = '#01665e', midpoint = 18) +
  scale_color_manual('+2C group shift?', values = c('white','black'), labels = c('No','Yes')) +
  labs(y = 'AET', x = 'DEF') +
  theme_bw() +
  coord_cartesian(xlim = c(150,900), ylim = c(200,750))

# Plot on map
# State boundaries for reference
boundary <- read_sf("/Users/tylerhoecker/Box/Work/PhD/GIS/cb_2018_us_state_20m/cb_2018_us_state_20m.shp") %>% 
  st_transform(., crs = 4326) 

# Remake spatial version of veg-fire df now that it includes hexbin info
fire_veg_sf <- st_as_sf(fire_veg_small, coords = c('x', 'y'), crs = 4326)
fire_veg_test_sf <- st_as_sf(fire_veg_test, coords = c('x', 'y'), crs = 4326)

tm_shape(boundary, bbox = fire_veg_sf) +
  tm_borders() +
  tm_fill(col = 'grey90') +
  tm_shape(fire_veg_test_sf) +
  tm_dots(title = 'Climate Bin',
          col = 'cell',
          #size = 0.01, 
          palette = c('#8c510a','#f6e8c3','#01665e'),
          style = 'cont',
          shape = 15) +
  tm_layout(legend.outside = TRUE)
# +
#   tm_text('cell', size = 0.5)

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



