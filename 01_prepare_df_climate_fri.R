# Prepares dataframe using data from Sean
# Extracts climate data for all points 
# Calculates FRP for all points, based on proportion of 4 km cell around each
#    point that burned in 35 year obs record

library(tidyverse)
library(raster)
library(terra)
library(sf)

# ------------------------------------------------------------------------------
# Fire and vegetation: Read in data on burned area from Sean; calculate FRP using
# information on ecoregion area; add in unburned samples
# ------------------------------------------------------------------------------
# Read in processed fire, veg data from Sean
fire_veg_df_full <- list.files('data/parks_data/', full.names = T) %>% 
  set_names(gsub("sev\\..*","",list.files('data/parks_data/'))) %>% 
  map_df(., read_csv, .id = 'region') 
# Do I need to filter out points that have burned multiple times?

# What proportion to subset?
subprop <- 1

# Subset this 
fire_veg_df <- fire_veg_df_full %>% 
  group_by(region) %>% 
  slice_sample(., prop = subprop) %>% 
  dplyr::select(region, x, y, year = fire.year, cbi, ndvi) %>% 
  filter(cbi > 0)

# Make dataframe spatial first with coordinates
fire_veg_sf <- st_as_sf(fire_veg_df, coords = c('x', 'y'), crs = 4326)

# Read in raster of burned areas - will be come master grid
burned_rast <- raster('data/forest/mtbs_rast.tif')
# Turn this grid into 0s (unburned) and 1s (burned), so that the mean of the aggregate cell is the 
# proportion of the larger cell that was burned
burned_rast <- reclassify(burned_rast, cbind(0, 1))
burned_rast <- reclassify(burned_rast, cbind(NA, 0))

# Read in raster of forest area... this is already set up as 0,1
# forest_rast <- raster('data/forest/forest.tif')
# Resample to identical grid to burned area raster
# forest_rast <- terra::resample(forest_rast, burned_rast, method = 'ngb')

# OR, after run once (it takes a while)
forest_rast <- raster('data/forest/forest_30m.tif')

# Need to create a version of the 30m forest mask that is NA,1 to mask the burned area
forest_na <- reclassify(forest_rast, cbind(0, NA))

# Now mask the burned area rast to the forest (it includes burned non-forest)
burned_forest_rast <- mask(burned_rast, forest_na)

# Aggregate to a 13.6 km grid (1/8th degree in both directions)
burn_frp_grid <- aggregate(burned_forest_rast, 
                           fact = c((1/8)/xres(burned_rast),(1/8)/yres(burned_rast)), 
                           fun = mean, na.rm = T)

# Aggregate this forested raster to a 13.6 km grid and calc. prop, as above
forest_frp_grid <- aggregate(forest_rast, 
                             fact = c((1/8)/xres(forest_rast),(1/8)/yres(forest_rast)), 
                      fun = mean, na.rm = T)


# ------------------------------------------------------------------------------
# Calculate fire return interval (fri) from burned area raster
# ------------------------------------------------------------------------------
# Extract values for each point of the fire-veg dataframe
#burned_area <- terra::area(burn_frp_grid)
frp_pts <- terra::extract(burn_frp_grid, fire_veg_sf, df = T)
for_pts <- terra::extract(forest_frp_grid, fire_veg_sf, df = T)

fire_veg_df <- fire_veg_df %>% 
  ungroup() %>% 
  mutate(frp = 35/(frp_pts$mtbs_rast/for_pts$forest_30m),
         frp = ifelse(frp == 'Inf', 500, frp),
         frp = ifelse(frp > 500, 500, frp)) #%>%  # 35 years divided by proportion burned = FRP, time required to burn study area
#dplyr::select(-prop, -area)

# Map frp 
plot_df <- fire_veg_df %>% 
  slice_sample(n = 100000)

plot_sf <- st_as_sf(plot_df, coords = c('x', 'y'), crs = 4326)

tm_shape(boundary, bbox = plot_sf) +
  tm_borders() +
  tm_fill(col = 'grey10') +
  tm_shape(plot_sf) +
  tm_dots(col = 'frp',
          size = 0.003,
          #palette = c('#8c510a','#f6e8c3','#01665e'),
          style = 'quantile',
          shape = 15) 

# tmap_save(filename = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/frp_map.png',
#           width = 4, height = 8, dpi = 500)

# ------------------------------------------------------------------------------
# Plot some general relationships 
# ------------------------------------------------------------------------------
# Explore CBI-NDVI relationship
ggplot() +
  geom_point(data = slice_sample(fire_veg_df, n = 5000), aes(x = ndvi, y = cbi),
             alpha = 0.3) +
  geom_smooth(data = slice_sample(fire_veg_df, n = 10000), aes(x = ndvi, y = cbi),
              size = 2) +
  theme_bw() +
  coord_cartesian(ylim = c(0.5, 2.5))

ggsave('ndvi-cbi_relationship.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
       height = 6, width = 8, dpi = 300)

ggplot() +
  geom_point(data = slice_sample(fire_veg_df, n = 5000), aes(x = frp, y = cbi),
             alpha = 0.3) +
  geom_smooth(data = slice_sample(fire_veg_df, n = 20000), aes(x = frp, y = cbi),
              size = 2) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 3))

ggsave('frp-cbi_relationship.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
       height = 6, width = 8, dpi = 300)



# ------------------------------------------------------------------------------
# Climate: Read in MACA climate data and add to dataframe with fire and veg data
# ------------------------------------------------------------------------------
# Read in AET and DEF TerraClim data for historical period
climate_stack <- list.files('data/terraclimate/', full.names = T) %>% 
  set_names(gsub("\\..*","",list.files('data/terraclimate/'))) %>%
  stack()

# Extract values from climate rasters
climate_fire_pts <- terra::extract(climate_stack, fire_veg_sf, df = T)

# C-Bind to the fire-veg dataframe... 
fire_veg_clim_df <- cbind(fire_veg_df, climate_fire_pts) %>% 
  filter(aet_hist > 0) 

saveRDS(fire_veg_clim_df, 'data/fire_veg_clim_df.Rdata')

# Do the same for all forest...
# All forest area, aggregated to ~4 km resolution
forest_4km <- terra::aggregate(forest_rast, 
                               fact = c((1/24)/xres(forest_rast),(1/24)/yres(forest_rast)), 
                               fun = mean, na.rm = T)
forest_4km  <- reclassify(forest_4km, matrix(c(0, 0.1, NA,  0.1, 1, 1), ncol=3, byrow=TRUE))
# Get centroids
forest_pts_sp <- rasterToPoints(forest_4km, spatial = T)
# Extract climate data
forest_climate <- terra::extract(climate_stack, forest_pts_sp, df = T) 

# Transform to dataframe and filter out non-forest cells, unrealistic AET vals
forest_clim_df <- forest_pts_sp %>% 
  as.data.frame() %>% 
  cbind(., forest_climate) %>% 
  filter(forest_30m == 1) %>% 
  dplyr::select(-forest_30m) %>% 
  filter(aet_hist > 0) 

saveRDS(forest_clim_df, 'data/forest_clim_df.Rdata')



# Old, save for now
# -----------------------

# A shapefile of 95,000 unburned forest points I created at some point... 
forest_sf <- st_read('data/forest/unburned.shp')

# Extract climate data
forest_climate <- terra::extract(climate_stack, forest_sf, df = T) %>% 
  filter(aet_hist > 0) 

saveRDS(forest_climate, 'data/forest_climate.Rdata')

# Extact climate data for ALL unburned forest points (at coarsened resolution)
nw_region <- st_read("C:/Users/hoecker/Work/GIS/doi_12_unified_regions_20180801_shapefile/region9_boundary.shp")
nw_region <- st_transform(nw_region, crs = 4326)

# Create a 4 km raster for forest area in the NW
# First crop forest_rast to NW
nw_forest <- crop(forest_rast, nw_region) #as_Spatial() 

# Then aggregate the fine forest_rast to 4 km (1/24th deg.) resolution
nw_coarse <- aggregate(nw_forest, fact = c((1/24)/xres(nw_forest),(1/24)/yres(nw_forest)), 
                         fun = mean, na.rm = T)
nw_coarse  <- reclassify(nw_coarse, matrix(c(0, 0.5, NA,  0.5, 1, 1), ncol=3, byrow=TRUE))
# Get centroids
nw_pts <- rasterToPoints(nw_coarse, spatial = T)
# Extract climate data
nw_climate <- terra::extract(climate_stack, nw_pts, df = T) 

# Save the centroids of this 4 km raster for later use (to map predictions across all points in PNW, including unburned)
nw_pts_sf <- st_as_sf(nw_pts)
saveRDS(nw_pts_sf, 'data/nw_pts_sf.Rdata')


nw_clim_df <- nw_pts %>% 
  as.data.frame() %>% 
  cbind(., nw_climate) %>% 
  filter(forest_30m == 1, ) %>% 
  dplyr::select(-forest_30m) %>% 
  filter(aet_hist > 0) 

saveRDS(nw_clim_df, 'data/nw_clim_df.Rdata')

