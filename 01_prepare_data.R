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

# ------------------------------------------------------------------------------
# Use dataframe which includes CBI and NDVI as the sample
# Eventually, the centroid of every burned 30x30m pixel in the western US
# ------------------------------------------------------------------------------
# Read in processed fire, veg data from Sean
fire_veg_df_full <- list.files('data/parks_data/', full.names = T) %>% 
  set_names(gsub("sev\\..*","",list.files('data/parks_data/'))) %>% 
  map_df(., read_csv, .id = 'region') 

# Optionally, subset a proporiton of points
subprop <- 1

# Subset this 
fire_veg_df <- fire_veg_df_full %>% 
  group_by(region) %>% 
  slice_sample(., prop = subprop) %>% 
  dplyr::select(region, x, y, year = fire.year, fire.id, cbi, ndvi) %>% 
  filter(cbi > 0.1) %>% 
  ungroup()

# Make dataframe spatial first with coordinates
fire_veg_sf <- st_as_sf(fire_veg_df, coords = c('x', 'y'), crs = 4326)

# ------------------------------------------------------------------------------
# Create master grid (Stevens et al. fire-resistance traits product)
# ------------------------------------------------------------------------------
# Stevens et al. FRS gridded dataset has the most constrained extent, so use
# this to create a master grid for masking out proceeding layers
# This produces a course version of the Stevens layer 
# Done -------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Extract fire resistance scores
# ------------------------------------------------------------------------------
frs_rast <- rast('data/composition/frs_spatial.tif')
frs_rast <- terra::project(frs_rast, y = "epsg:4326")

mast_rast <- terra::aggregate(frs_rast,
                              fact = c((1/96)/xres(frs_rast),(1/96)/yres(frs_rast)),
                              fun = 'modal', na.rm = T) %>%
  classify(., matrix(c(0,1,1), ncol=3, byrow=TRUE))

writeRaster(mast_rast, 'data/process/mast_rast.tif', overwrite=T)
# ------------------------------------------------------------------------------
mast_rast <- rast('data/process/mast_rast.tif')
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Filter out fire-veg points that fall outside this more restrictive definition 
# of conifer-dominated forest. 
fire_veg_is_forest <- terra::extract(mast_rast, st_coordinates(fire_veg_sf))

fire_veg_sf <- fire_veg_sf %>% 
  mutate(is_forest = fire_veg_is_forest$frs_spatial) %>% 
  filter(!is.na(is_forest))

fire_veg_df <- fire_veg_sf %>% 
  st_drop_geometry() %>% 
  cbind(st_coordinates(fire_veg_sf))

# ------------------------------------------------------------------------------
# Climate 
# Read in MACA climate data and extract from unburned and burned forested points
# ------------------------------------------------------------------------------
# Read in AET and DEF TerraClim data for historical period
climate_stack <- list.files('data/terraclimate/', full.names = T) %>% 
  set_names(gsub("\\..*","",list.files('data/terraclimate/'))) %>% 
  stack() %>% 
  terra::rast()

# Burned points 
# Extract values from climate rasters
climate_fire_pts <- terra::extract(climate_stack, st_coordinates(fire_veg_sf))

# C-Bind to the fire-veg dataframe... 
pyroclimate_input_df <- cbind(fire_veg_df, climate_fire_pts) %>% 
  # Probably due to reprojection etc, 5 points from CA land outside the climate raster
  filter(aet_hist > 0)
saveRDS(pyroclimate_input_df, 'data/process/pyroclimate_input_df.Rdata')

# All forest points
# Get centroids
forest_pts_sp <- as.points(mast_rast)
# Extract climate data
forest_climate <- terra::extract(climate_stack, forest_pts_sp, xy = TRUE) 

# ------------------------------------------------------------------------------
# Calculate fire rotation period
# ------------------------------------------------------------------------------
# Done -------------------------------------------------------------------------
# Read in a rasterized version of MTBS data. Pre-processing in QGIS, based on this
# thread: https://gis.stackexchange.com/questions/297274/create-raster-out-of-overlapping-buffers-within-the-same-shapefile-in-qgis-3-0
# - Identify portions of polygons that overlap (aka reburns)
# - Rasterize to mast_rast grid
# - Mask out areas of non-forest (rasterize just matches grid)
reburn_rast <- rast('data/process/forest_reburns.tif')
# Change NAs to 0
reburn_rast <- classify(reburn_rast, cbind(NA, 0))
# Mask out non-forest (change those areas back to NA)
reburn_rast <- mask(reburn_rast, mast_rast)

# Record whether forest points burned or not
forest_burned <- terra::extract(reburn_rast, forest_pts_sp)

# Calculate the area of each cell (its a lat/long grid so area varies)
forest_cell_area <- cellSize(reburn_rast)
forest_area_df <- terra::extract(forest_cell_area, forest_pts_sp)

# Gut check on these values against Abatzoglous et al. observed data: https://datadryad.org/stash/landing/show?id=doi%3A10.6071%2FM3WQ1R
sum(forest_burned$forest_reburns*forest_area_df$area)*0.000001
sum(read_csv('data/forest/abatzoglou_comparison.csv')$area)
# Wow! That is crazy close! Good.

# Transform to dataframe and filter out non-forest cells, unrealistic AET vals
forest_clim_df <- forest_pts_sp %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  # Combine with climate data extracted for each point
  cbind(forest_climate,
        ., 
        'burned' = forest_burned$forest_reburns, 
        'area' = forest_area_df$area) %>% 
  # Leftover layer name from creating the mast_rast from Jen's FRS product
  dplyr::select(-frs_spatial) %>% 
  # Probably due to reprojection etc., a small number of points fall just off the climate grid, removing them
  filter(aet_hist > 0) 

saveRDS(forest_clim_df, 'data/process/forest_clim_df.Rdata')
# ------------------------------------------------------------------------------
forest_clim_df <- readRDS('data/process/forest_clim_df.Rdata')
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Some checks
# ------------------------------------------------------------------------------
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


