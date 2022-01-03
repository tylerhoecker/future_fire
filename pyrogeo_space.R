library(tidyverse)
library(raster)
library(sf)
library(hexbin)
library(tmap)


# Read in processed fire, veg data from Sean
fire_veg_df_full <- list.files('data/fire/', full.names = T) %>% 
  set_names(gsub("sev\\..*","",list.files('data/fire/'))) %>% 
  map_df(., read_csv, .id = 'region') 
  # Do I need to filter out points that have burned multiple times?

# Subset this 
fire_veg_df <- fire_veg_df_full %>% 
  group_by(region) %>% 
  slice_sample(., prop = .01) 

# Read in forest mask
forest_rast <- raster('data/forest/fw.mask.dd.tif')
forest_rast <- projectRaster(forest_rast, crs = 4326)


unburned_samp <- sampleRandom(x = unburned_rast, 
                              size = length(fire_veg_df_full[['x']]),
                              na.rm = TRUE,
                              xy = TRUE)
plot(unburned_rast)


# Make unburned sample a spatial sf

# Get region for unburned points

# Get NDVI for unburned points



# Pick random sample equal in size to burned cells from non-burned cells





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

# Read in AET and DEF TerraClim data for historical period
climate <- list.files('data/terraclimate/', full.names = T) %>% 
  set_names(gsub("\\..*","",list.files('data/terraclimate/'))) %>%
  stack()

# Add climate data as columns to dataframe
# Make dataframe spatial first with coordinates
fire_veg_sf <- st_as_sf(fire_veg_df, coords = c('x', 'y'), crs = 4326)

# Crop this for faster processing
climate_crop <- crop(climate, fire_veg_sf)

# Extract values from climate rasters
climate_pts <- extract(climate_crop, fire_veg_sf, df = T)

# C-Bind to the fire-veg dataframe... this is slow!
fire_veg_clim_df <- cbind(fire_veg_df, climate_pts) %>% 
  #sample_n(100) %>% 
  filter(aet_hist > 0) %>% 
  ungroup() %>% 
  pivot_longer(cols = c(aet_2C, aet_hist, def_2C, def_hist), 
               names_to = 'temp_name', values_to = 'value') %>% 
  dplyr::select(-ID) %>% 
  separate(temp_name, into = c('variable','period'), sep = '_') %>% 
  pivot_wider(names_from = variable, values_from = value) 

# Number of bins in x direction
n_bins <- 20
#xrange <- c(0,200)
#yrange <- c(0,110)

# Create hexbins and save the hexbin info as it's own df
h <- hexbin(x = fire_veg_clim_df$def, y = fire_veg_clim_df$aet, 
            xbins= n_bins, 
            shape = 1, 
            IDs = TRUE)
hexdf <- data.frame(hcell2xy (h),  hexID = h@cell, counts = h@count)

# Create same hexbins in the complete veg-fire df... this is slow!
pyrome_df <- fire_veg_clim_df %>% 
  mutate(cell = hexbin(x = def, y = aet, 
                       xbins= n_bins, 
                       shape = 1, 
                       IDs = TRUE)@cID) %>% 
  # Do stuff here to identify whether or not points move groups between periods
  # Its slow so just do a subset
  group_by(x,y) %>% 
  mutate(ptID = cur_group_id()) %>% 
  group_by(ptID) %>% 
  mutate(change = as.factor(ifelse(length(unique(cell)) > 1, 1, 0)))


# How many points per group?
pyrome_df %>% 
  group_by(cell) %>% 
  tally() %>% 
  arrange(sample = desc(n))

# Create a(n even smaller) smaller version for plotting etc
plot_df <-  pyrome_df %>%
  #arrange(ptID) %>% 
  group_by(region) %>% 
  slice_head(n = 50) 


ggplot(plot_df, aes(x = def, y = aet)) +
  geom_hex(data = hexdf, aes(x = x, y = y, fill = hexID), color = 'black', alpha = 0.8, stat = 'identity') +
  geom_text(data = hexdf, aes(x = x, y = y, label = hexID)) +
  geom_point(aes(alpha = change), size = 1) +
  geom_line(data = filter(plot_df, change == 0),
            aes(group = ptID, alpha = change), size = 0.75) +
  geom_line(data = filter(plot_df, change == 1),
            aes(group = ptID, alpha = change), size = 0.75, 
            arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "open")) +
  scale_fill_gradient2('Climate Hexbin', low = '#8c510a', mid = '#f6e8c3', high = '#01665e', midpoint = 250) +
  #scale_color_manual('+2C group shift?', values = c('white','black'), labels = c('No','Yes')) +
  scale_alpha_manual('+2C group shift?', values = c(0.3,0.8), labels = c('No','Yes')) +
  labs(y = 'AET', x = 'DEF') +
  theme_bw() 

ggsave('hex_move_eg.png',path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
       height = 6, width = 8, dpi = 300)


# Plot on map
# State boundaries for reference
boundary <- read_sf("C:\\Users\\hoecker\\Work\\GIS\\cb_2018_us_state_20m\\cb_2018_us_state_20m.shp") %>% 
  st_transform(., crs = 4326) 

# Remake spatial version of veg-fire df now that it includes hexbin info
pyrome_sf <- st_as_sf(pyrome_df, coords = c('x', 'y'), crs = 4326)
#fire_veg_test_sf <- st_as_sf(fire_veg_test, coords = c('x', 'y'), crs = 4326)

tm_shape(boundary, bbox = pyrome_sf) +
  tm_borders() +
  tm_fill(col = 'grey10') +
  tm_shape(pyrome_sf) +
  tm_dots(title = 'Climate Bin',
          col = 'cell',
          size = 0.01, 
          palette = c('#8c510a','#f6e8c3','#01665e'),
          style = 'cont',
          shape = 15) +
  tm_layout(legend.outside = TRUE)








# Mask out burned cells
# burned_cells <- cellFromXY(forest_rast, cbind(fire_veg_df_full[,c('x','y')]))
# unburned_rast <- forest_rast
# unburned_rast[burned_cells] <- NA # This is pretty slow...










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



