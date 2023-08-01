# ------------------------------------------------------------------------------
# Description: Prepare dataset for analysis from published data, save for use in 
# subsequent scripts. Burn severity, productivity, fire occurrence, and climate 
# data are harmonized.
#
# Input data that are publicly available elsewhere are not re-published in 
# this archive. The code in this script illustrates how the dataset we used for 
# analysis was created from publicly available data.
# A minimal dataset, which can be used to reproduce our analysis,
# is read in at the end of the script. 
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load packages
# ------------------------------------------------------------------------------
library(tidyverse)
library(raster)
library(terra)
library(sf)
select <- dplyr::select
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
# Create master grid 
# ------------------------------------------------------------------------------
# Stevens et al. FRS gridded dataset has the most constrained extent, so use
# this to create a master grid for masking out proceeding layers.
# frs_spatial is available here: https://www.sciencebase.gov/catalog/item/5e39f54ee4b0a79317e15e73

# This produces a coarser version of the Stevens layer 
frs_rast <- rast('data/composition/frs_spatial.tif')
frs_rast <- terra::project(frs_rast, y = "epsg:4326")
mast_rast <- terra::aggregate(frs_rast,
                              fact = c((1/96)/xres(frs_rast),(1/96)/yres(frs_rast)),
                              fun = 'modal', na.rm = T) %>%
  classify(., matrix(c(0,1,1), ncol=3, byrow=TRUE))

writeRaster(mast_rast, 'data/process/mast_rast.tif', overwrite=T)
mast_rast <- rast('data/process/mast_rast.tif')
# ------------------------------------------------------------------------------
# Filter out fire-veg points that fall outside this more restrictive definition 
# of conifer-dominated forest. 
# ------------------------------------------------------------------------------
fire_veg_is_forest <- terra::extract(mast_rast, st_coordinates(fire_veg_sf))

fire_veg_sf <- fire_veg_sf %>% 
  mutate(is_forest = fire_veg_is_forest$frs_spatial) %>% 
  filter(!is.na(is_forest))

# Save as a non-spatial dataframe, to add additional variables to
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

# Get centroids of all forest points
forest_pts_sp <- as.points(mast_rast)
# Extract climate data
forest_climate <- terra::extract(climate_stack, forest_pts_sp, xy = TRUE) 

# ------------------------------------------------------------------------------
# Calculate fire rotation period
# ------------------------------------------------------------------------------
# Read in a rasterized version of MTBS data. 
# Pre-processing in QGIS, based on this thread: https://gis.stackexchange.com/questions/297274/create-raster-out-of-overlapping-buffers-within-the-same-shapefile-in-qgis-3-0
# - Identify portions of polygons that overlap (aka reburns)
# - Rasterize to mast_rast grid
# - Mask out areas of non-forest (rasterize just matches grid)
reburn_rast <- rast('data/process/forest_reburns.tif')
# Change NAs to 0
reburn_rast <- classify(reburn_rast, cbind(NA, 0))
# Mask out non-forest (change those areas back to NA)
reburn_rast <- mask(reburn_rast, mast_rast)

# Aggregate to a 13.6 km grid (1/8th degree in both directions)
# This will take a "mean" of the number of times each small cell burned.
# This is equivalent to the proportion of the larger cell that burned.
reburn_frp_grid <- terra::aggregate(reburn_rast,
                                    fact = c((1/8)/xres(reburn_rast),
                                             (1/8)/yres(reburn_rast)),
                                    fun = mean, na.rm = T)

reburn_fire_pts <- terra::extract(reburn_frp_grid, st_coordinates(fire_veg_sf))

# Aggregate the mast_rast grid (indicating forested area) to 1/8th degree grid,
# indicating proportion of coarse grid that is forested.
forest_frp_grid <-  mast_rast %>%
  classify(., matrix(c(0,1,1), ncol=3, byrow=TRUE)) %>%
  classify(., cbind(NA, 0)) %>%
  terra::aggregate(.,
                   fact = c((1/8)/xres(.),
                            (1/8)/yres(.)),
                   fun = mean, na.rm = T)

forest_prop_fire_pts <- terra::extract(forest_frp_grid, st_coordinates(fire_veg_sf))

# Produce a dataframe with FRP values for all pyroclimate points
frp_df <- 
  tibble('prop_burn' = reburn_fire_pts$forest_reburns,
         'prop_forest' = forest_prop_fire_pts$frs_spatial) %>% 
  mutate(prop_for_burn = (prop_burn*prop_forest),
         prop_for_burn = ifelse(prop_for_burn < 0.036, runif(1,0.036,0.72), prop_for_burn),
         frp_pt = (36 / prop_for_burn),
         log_frp_pt = log(frp_pt))

# Gut check these values against Abatzoglous et al. observed data for same time period 
# I removed everything except 1984-2019!
# https://datadryad.org/stash/landing/show?id=doi%3A10.6071%2FM3WQ1R
forest_area_burned <- (reburn_frp_grid*forest_frp_grid)*(cellSize(reburn_frp_grid)*0.000001)
analog_estimate <- global(forest_area_burned, 'sum', na.rm = T)
abatzoglou_estimate <- read_csv('data/forest/abatzoglou_comparison.csv') %>% 
     dplyr::select(area) %>% 
     sum()

abatzoglou_estimate/36
analog_estimate/36
abatzoglou_estimate-analog_estimate  
# Difference of 777.2 ha - that's very close, 0.6% lower!

# Write a raster of estimated FRPs based on this method
prop_burn <- forest_frp_grid*reburn_frp_grid
m <- c(0, 0.036, 0.036)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
prop_burn <- classify(prop_burn, rclmat, include.lowest=TRUE)
frp_grid <- 36/prop_burn
frp_grid <- mask(frp_grid, mast_rast)

writeRaster(frp_grid, 'data/process/frp_grid_yrs_1-8.tiff', overwrite=TRUE)

# ------------------------------------------------------------------------------
# Combine things into one dataframe
# ------------------------------------------------------------------------------
pyroclimate_input_df <- cbind(fire_veg_df, climate_fire_pts, frp_df) %>% 
  # Probably due to reprojection, 5 points from CA land outside the climate raster
  filter(aet_hist > 0)

saveRDS(pyroclimate_input_df, 'data/process/pyroclimate_input_df.Rdata')

# Transform to dataframe and filter out non-forest cells, unrealistic AET vals
forest_clim_df <- forest_pts_sp %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  cbind(forest_climate,
        .) %>% 
  as_tibble() %>% 
  # Leftover layer name from creating the mast_rast from Jen's FRS product
  dplyr::select(-frs_spatial) %>% 
  # Probably due to reprojection etc., a small number of points fall just off the climate grid, removing them
  filter(aet_hist > 0) 

saveRDS(forest_clim_df, 'data/process/forest_clim_df.Rdata')
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Read in minimal dataset, to run proceeding steps, here
# ------------------------------------------------------------------------------
forest_clim_df <- readRDS('data/process/forest_clim_df.Rdata')
pyroclimate_input_df <-  readRDS('data/process/pyroclimate_input_df.Rdata')
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


