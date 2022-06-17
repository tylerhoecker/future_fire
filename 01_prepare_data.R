# Prepares dataframe of fire niche dimensions
# CBI (severity) and NDVI (productivity) from data from Sean Parks
# Calculates FRI from MTBS perimeter data
# Extracts fire resistance scores from Jens Stevens dataset
# Extracts TerraClimate climate data from TIFFs provided by SP

# ------------------------------------------------------------------------------
# Load packages
# ------------------------------------------------------------------------------
library(tidyverse)
library(raster)
library(terra)
library(sf)
library(tmap)

# ------------------------------------------------------------------------------
# Use dataframe which includes CBI and NDVI as the sample
# Eventually, the centroid of every burned 30x30m pixel in the western US
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
  dplyr::select(region, x, y, year = fire.year, fire.id, cbi, ndvi) %>% 
  filter(cbi > 0.1)

# Make dataframe spatial first with coordinates
fire_veg_sf <- st_as_sf(fire_veg_df, coords = c('x', 'y'), crs = 4326)

# ------------------------------------------------------------------------------
# Extract fire resistance scores
# ------------------------------------------------------------------------------
frs_rast <- rast('data/composition/frs_spatial.tif')
frs_rast <- terra::project(frs_rast, y = "epsg:4326")

frs_pts <- terra::extract(frs_rast, st_coordinates(fire_veg_sf))

# ------------------------------------------------------------------------------
# Create master grid (based on fire resistance grid)
# ------------------------------------------------------------------------------
# Stevens et al. FRS gridded dataset has the most constrained extent, so use 
# this to create a master grid for masking out proceeding layers
# Done --------------
# mast_rast <- terra::aggregate(frs_rast,
#                               fact = c((1/96)/xres(frs_rast),(1/96)/yres(frs_rast)), 
#                               fun = 'modal', na.rm = T) %>% 
#   classify(., matrix(c(0,1,1), ncol=3, byrow=TRUE)) 
# 
# writeRaster(mast_rast, 'data/process/mast_rast.tif')
## -------------------

mast_rast <- rast('data/process/mast_rast.tif')
# ------------------------------------------------------------------------------
# Calculate fire return interval (fri) from burned area raster
# ------------------------------------------------------------------------------
# Read in raster of burned areas. This could be created by rasterizing MTBS
# polygons to a raster grid
burned_rast <- rast('data/forest/mtbs_rast.tif')

# Resample burned area raster to identical grid as mast_rast
burned_rast <- terra::resample(burned_rast, mast_rast, method = 'near')

# Now mask the burned area rast to the forest (it includes burned non-forest)
burned_forest_rast <- mask(burned_rast, mast_rast)

writeRaster(burned_forest_rast, 'data/process/burned_forest_rast.tif')

# Turn this grid into 0s (unburned) and 1s (burned), so that the mean of the aggregate cell is the 
# proportion of the larger cell that was burned
burned_forest_rast <- classify(burned_forest_rast, cbind(0, 1))
burned_forest_rast <- classify(burned_forest_rast, cbind(NA, 0))

# Aggregate to a 13.6 km grid (1/8th degree in both directions)
burn_frp_grid <- terra::aggregate(burned_forest_rast, 
                                  fact = c((1/8)/xres(burned_forest_rast),
                                           (1/8)/yres(burned_forest_rast)), 
                                  fun = mean, na.rm = T)

# Aggregate the mast_rast grid (indicating forested area) to same grid
forest_frp_grid <-  mast_rast %>% 
  classify(., matrix(c(0,1,1), ncol=3, byrow=TRUE)) %>% 
  classify(., cbind(NA, 0)) %>% 
  terra::aggregate(.,
                   fact = c((1/8)/xres(frs_rast),(1/8)/yres(frs_rast)), 
                   fun = mean, na.rm = T)

# Extract values for each point of the fire-veg dataframe
frp_pts <- terra::extract(burn_frp_grid, st_coordinates(fire_veg_sf))
for_pts <- terra::extract(forest_frp_grid, st_coordinates(fire_veg_sf))

# Calculate fire-return interval as proportion of forested area that burned in 35 years
fri_pts <- 35/(frp_pts$mtbs_rast/for_pts$frs_spatial)

# ------------------------------------------------------------------------------
# Climate: Read in MACA climate data and extract from unburned and burned forested points
# ------------------------------------------------------------------------------
# Read in AET and DEF TerraClim data for historical period
climate_stack <- list.files('data/terraclimate/', full.names = T) %>% 
  set_names(gsub("\\..*","",list.files('data/terraclimate/'))) %>% 
  stack() %>% 
  terra::rast()

# Extract values from climate rasters
climate_fire_pts <- terra::extract(climate_stack, st_coordinates(fire_veg_sf))

# Build climate-only dataframe for all forest points (vs burned points in other df)
# Get centroids
forest_pts_sp <- as.points(mast_rast)
# Extract climate data
forest_climate <- terra::extract(climate_stack, forest_pts_sp, xy = TRUE) 

# Transform to dataframe and filter out non-forest cells, unrealistic AET vals
forest_clim_df <- forest_pts_sp %>% 
  as.data.frame() %>% 
  cbind(., forest_climate) %>% 
  filter(frs_spatial == 1) %>% 
  dplyr::select(-frs_spatial) %>% 
  filter(aet_hist > 0) 

saveRDS(forest_clim_df, 'data/process/forest_clim_df.Rdata')


# ------------------------------------------------------------------------------
# Add Climtae, FRI and FRS to the dataframe with burned area info
# ------------------------------------------------------------------------------
# C-Bind to the fire-veg dataframe... 
fire_veg_clim_df <- cbind(fire_veg_df, climate_fire_pts) 

pyroclimate_input_df <- fire_veg_clim_df %>% 
  ungroup() %>% 
  mutate(fri = fri_pts,
         frs = frs_pts$frs_spatial) %>% 
  filter(fri != 'Inf', !is.na(fri)) %>% # These are weird cases resulting from resampling, I think. It's not very many (~20k/9mil).
  mutate(log_fri = log(fri)) %>% 
  # The FRS grid uses a more conservative definition of forest, so filter out NAs from FRS
  filter(!is.na(frs)) 

# Previous method... 
# fri = ifelse(fri == 'Inf', 500, fri),
# fri = ifelse(fri > 500, 500, fri),

saveRDS(pyroclimate_input_df, 'data/process/pyroclimate_input_df.Rdata')


# Maps for checking 
plot_df <- fire_veg_df %>% 
  slice_sample(n = 100000)

plot_sf <- st_as_sf(plot_df, coords = c('x', 'y'), crs = 4326)
boundary <- read_sf("C:\\Users\\hoecker\\Work\\GIS\\cb_2018_us_state_20m\\cb_2018_us_state_20m.shp") %>%
  st_transform(., crs = 4326)

tm_shape(boundary, bbox = plot_sf) +
  tm_borders() +
  tm_fill(col = 'grey10') +
  tm_shape(plot_sf) +
  tm_dots(col = 'cbi',
          size = 0.003,
          #palette = c('#8c510a','#f6e8c3','#01665e'),
          style = 'quantile',
          shape = 15) 


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
  geom_point(data = slice_sample(fire_veg_df, n = 5000), aes(x = fri, y = cbi),
             alpha = 0.3) +
  geom_smooth(data = slice_sample(fire_veg_df, n = 20000), aes(x = fri, y = cbi),
              size = 2) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 3))

ggsave('frp-cbi_relationship.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
       height = 6, width = 8, dpi = 300)

ggplot() +
  geom_point(data = slice_sample(fire_veg_df, n = 5000), aes(x = fri, y = frs),
             alpha = 0.3) +
  geom_smooth(data = slice_sample(fire_veg_df, n = 20000), aes(x = fri, y = frs),
              size = 2, method = 'lm') +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1)) 

ggsave('frp-cbi_relationship.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
       height = 6, width = 8, dpi = 300)




# Old, save for now
# -----------------------
# 
# # A shapefile of 95,000 unburned forest points I created at some point... 
# forest_sf <- st_read('data/forest/unburned.shp')
# 
# # Extract climate data
# forest_climate <- terra::extract(climate_stack, forest_sf, df = T) %>% 
#   filter(aet_hist > 0) 
# 
# saveRDS(forest_climate, 'data/forest_climate.Rdata')
# 
# # Extact climate data for ALL unburned forest points (at coarsened resolution)
# nw_region <- st_read("C:/Users/hoecker/Work/GIS/doi_12_unified_regions_20180801_shapefile/region9_boundary.shp")
# nw_region <- st_transform(nw_region, crs = 4326)
# 
# # Create a 4 km raster for forest area in the NW
# # First crop forest_rast to NW
# nw_forest <- crop(forest_rast, nw_region) #as_Spatial() 
# 
# # Then aggregate the fine forest_rast to 4 km (1/24th deg.) resolution
# nw_coarse <- aggregate(nw_forest, fact = c((1/24)/xres(nw_forest),(1/24)/yres(nw_forest)), 
#                          fun = mean, na.rm = T)
# nw_coarse  <- reclassify(nw_coarse, matrix(c(0, 0.5, NA,  0.5, 1, 1), ncol=3, byrow=TRUE))
# # Get centroids
# nw_pts <- rasterToPoints(nw_coarse, spatial = T)
# # Extract climate data
# nw_climate <- terra::extract(climate_stack, nw_pts, df = T) 
# 
# # Save the centroids of this 4 km raster for later use (to map predictions across all points in PNW, including unburned)
# nw_pts_sf <- st_as_sf(nw_pts)
# saveRDS(nw_pts_sf, 'data/nw_pts_sf.Rdata')
# 
# 
# nw_clim_df <- nw_pts %>% 
#   as.data.frame() %>% 
#   cbind(., nw_climate) %>% 
#   filter(forest_30m == 1, ) %>% 
#   dplyr::select(-forest_30m) %>% 
#   filter(aet_hist > 0) 
# 
# saveRDS(nw_clim_df, 'data/nw_clim_df.Rdata')
# # Read in raster of forest area... this is already set up as 0,1
# # forest_rast <- raster('data/forest/forest.tif')
# # Resample to identical grid to burned area raster
# # forest_rast <- terra::resample(forest_rast, burned_rast, method = 'ngb')
# 
# # OR, after run once (it takes a while)
# forest_rast <- raster('data/forest/forest_30m.tif')
# 
# # Need to create a version of the 30m forest mask that is NA,1 to mask the burned area
# forest_na <- reclassify(forest_rast, cbind(0, NA))
# 


