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










# OLD


# REPLACED BY REBURN & CLIMATE BIN VERSION ------------------------------------

# ------------------------------------------------------------------------------
# Calculate fire return interval (fri) from burned area raster
# ------------------------------------------------------------------------------
# Read in raster of burned areas. This was created by rasterizing MTBS
# polygons to a raster grid
# Done -------------------------------------------------------------------------
burned_rast <- rast('data/forest/mtbs_rast.tif')

# Resample burned area raster to identical grid as mast_rast
burned_rast <- terra::resample(burned_rast, mast_rast, method = 'near')

# Now mask the burned area rast to the forest (it includes burned non-forest)
burned_forest_rast <- mask(burned_rast, mast_rast)
burned_forest_rast <- rast('data/process/burned_forest_rast.tif')

# Turn this grid into 0s (unburned) and 1s (burned), so that the mean of the aggregate cell is the 
# proportion of the larger cell that was burned
burned_forest_rast <- classify(burned_forest_rast, cbind(0, 1))
burned_forest_rast <- classify(burned_forest_rast, cbind(NA, 0))

# Write out for quickness or reference later
writeRaster(burned_forest_rast, 'data/process/burned_forest_rast.tif', overwrite = T)
# ------------------------------------------------------------------------------
burned_forest_rast <- rast('data/process/burned_forest_rast.tif')
# ------------------------------------------------------------------------------

# Aggregate to a 13.6 km grid (1/8th degree in both directions)
# Done -------------------------------------------------------------------------
burn_frp_grid <- terra::aggregate(burned_forest_rast,
                                  fact = c((1/8)/xres(burned_forest_rast),
                                           (1/8)/yres(burned_forest_rast)),
                                  fun = mean, na.rm = T)

# Aggregate the mast_rast grid (indicating forested area) to same grid
forest_frp_grid <-  mast_rast %>%
  classify(., matrix(c(0,1,1), ncol=3, byrow=TRUE)) %>%
  classify(., cbind(NA, 0)) %>%
  terra::aggregate(.,
                   fact = c((1/8)/xres(burned_forest_rast),
                            (1/8)/yres(burned_forest_rast)),
                   fun = mean, na.rm = T)

# Save these for reference
writeRaster(burn_frp_grid, 'data/process/frp_grid.tiff', overwrite = T)
writeRaster(forest_frp_grid, 'data/process/forest_grid.tiff', overwrite = T)
# ------------------------------------------------------------------------------
burn_frp_grid <- rast('data/process/frp_grid.tiff')
forest_frp_grid <- rast('data/process/forest_grid.tiff')
# ------------------------------------------------------------------------------

# Calculate the FRP, as a proportion (of forested area that burned)
frp_prop_rast <- (burn_frp_grid/forest_frp_grid)
# Change zero and very small proportions (0-0.05) to a small proportion (to avoid Inf FRPs)
frp_prop_rast <- classify(frp_prop_rast, matrix(c(0,0.04,0.04), ncol=3, byrow=TRUE),
                          include.lowest = TRUE, right = NA)

frp_prop_rast <- classify(frp_prop_rast, cbind(Inf,0.04))

# Convert this proporiton to years; FRPs will be constrained from 35-700
frp_yrs_rast <- 36/frp_prop_rast

# Extract for the study area points
fri_pts <- extract(frp_yrs_rast, st_coordinates(fire_veg_sf))[,'mtbs_rast']
# ------------------------------------------------------------------------------

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

saveRDS(pyroclimate_input_df, 'data/process/pyroclimate_input_df.Rdata')


