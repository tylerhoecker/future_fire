# ------------------------------------------------------------------------------
# This script compares TerraClimate data with ClimateNA using the burned points
# dataset provided by Sean Parks
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
library(tidyverse)
library(terra)
library(raster)
library(sf)

# ------------------------------------------------------------------------------
# Read in points from Sean, small subset of them
# ------------------------------------------------------------------------------
fire_veg_df_full <- list.files('../data/parks_data/', full.names = T) %>% 
  set_names(gsub("sev\\..*","",list.files('../data/parks_data/'))) %>% 
  map_df(., read_csv, .id = 'region') 

# Optionally, subset a proporiton of points
subprop <- .001

# Subset this 
fire_veg_df <- fire_veg_df_full %>% 
  group_by(region) %>% 
  slice_sample(., prop = subprop) %>% 
  dplyr::select(region, x, y, year = fire.year, fire.id, cbi, ndvi) %>% 
  filter(cbi > 0.1) %>% 
  ungroup()

# ------------------------------------------------------------------------------
# Prepare input for ClimateNA stand-alone program
# ------------------------------------------------------------------------------
# To get multiple points CSV set up 
# Read in a high-resolution DEM
dem_tif <- rast('../../../Work/GIS/DEM/dem90_hf.tif')

# In case of re-start, make sure the sample of points is the same
fire_veg_sf <- st_as_sf(fire_veg_df, coords = c('x', 'y'), crs = crs(dem_tif))

dem_pts <- terra::extract(dem_tif, st_coordinates(fire_veg_sf), xy = T)

burn_pts <- data.frame('ID1' = seq(1,length(dem_pts$dem90_hf)),
                       'ID2' = rep('.',length(dem_pts$dem90_hf)),
                       'lat' = dem_pts$y,
                       'long' = dem_pts$x,
                       'el' = dem_pts$dem90_hf) 
write_csv(burn_pts, 'burn_input_pts.csv')  

# Prepare a DEM in ASCII format
dem_tif <- raster('../../../Work/GIS/DEM/USGS_1_n49w114_20210607.tif')
dem_tif <- terra::aggregate(dem_tif,
                            fact = 3,
                            fun = 'modal', na.rm = T)

raster::NAvalue(dem_tif) <- -9999
writeRaster(dem_tif, '../../../Work/GIS/DEM/glac_tile.asc',
            overwrite=TRUE)

# In case of re-start, make sure the sample of points is the same
fire_veg_sf <- st_as_sf(burn_pts, coords = c('long', 'lat'), crs = crs(dem_tif))
# ------------------------------------------------------------------------------
# Extract TerraClimate data for sample of points
# ------------------------------------------------------------------------------
# TerraClimate data
climate_stack <- list.files('data/terraclimate/', full.names = T) %>% 
  set_names(gsub("\\..*","",list.files('data/terraclimate/'))) %>% 
  raster::stack() %>% 
  terra::rast()

# Extract values from TerraClimate rasters for burned points
terraclim <- terra::extract(climate_stack, st_coordinates(fire_veg_sf))


# ------------------------------------------------------------------------------
# Read in ClimateNA data for sample of points
# ------------------------------------------------------------------------------
climateNA <- read_csv('data/burn_input_pts_Normal_1961_1990Y.csv')


# ------------------------------------------------------------------------------
# Compare
# ------------------------------------------------------------------------------
# Combine into one dataframe

climate_hist_df <- data.frame(burn_pts$long,
                              burn_pts$lat,
                              def_terra = terraclim$def_hist,
                              cmd_climNA = climateNA$CMD,
                              cmi_climNA = climateNA$CMI)

ggplot(climate_hist_df, aes(x = def_terra, y = cmd_climNA)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = 'lm') +
  geom_abline(intercept = 0, slope = 1, color = 'red', size = 1) +
  theme_bw() +
  labs(x = 'Climatic water deficit - TerraClimate (mm)', 
       y = "Hargreaves moisture deficit - ClimateNA (mm)",
       title = 'Comparing moisture deficit metrics - 1961-1990 normals')

ggsave('TerraClimDef_vs_ClimateNAHargreaves.png', width = 5, height = 5, unit = 'in')

ggplot(climate_hist_df, aes(x = def_terra, y = cmd_climNA)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = 'lm') +
  geom_abline(intercept = 0, slope = 1, color = 'red', size = 1) +
  theme_bw() +
  labs(x = 'Climatic water deficit - TerraClimate', 
       y = "Hogg's climate (soil) moisture deficit - ClimateNA",
       title = 'Comparing moisture deficit metrics - 1961-1990 normals')





