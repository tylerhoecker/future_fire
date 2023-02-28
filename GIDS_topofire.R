#######################################################################################
# CODE ADAPTED FROM:
# Rodman, K C., et al. 2020. 
# A Changing Climate is Snuffing Out Post-Fire Recovery in Montane Forests. 
# Global Ecology and Biogeography.

#######################################################################################
## This code performs a multi-step procedure that spatially downscales gridded climate
# data using GIDS (Gradient and Inverse Distance-Squared) of Nalder and Wein (1998)and
# Flint and Flint (2012)

###################################################################################
## Import necessary packages
###################################################################################
# Install missing packages and load them
using <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs,require,character.only=TRUE))
  need <- libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}

using('purrr','furrr','terra','spatialEco','geosphere')

# TESTING
#-------------------------------------------------------------------------------
# Focal testing area (Glacier)
bounds <- ext(-114, -113.5, 48.5, 48.8)

# Holden 250-m CWD normal raster
def_hist_fine <- rast('data/topofire/def_1981_2010.tif')
zoom(def_hist_fine, e = bounds)

# Raster that will be downscaled.
# Use this as  template for coarse raster to everything is aligned
def_future_coarse <- rast('data/terraclimate/def_2C.1961_1990.tif')
zoom(def_future_coarse, e = bounds)

# Create a coarse (4 km) version, which will align with future data, to regress against the fine data
# Resample, as a form of aggregation, the fine historical data to the desired coarse grid
def_hist_coarse <- resample(def_hist_fine, def_future_coarse, 'bilinear')
# Not sure why this is necessary... but it is
crs(def_hist_coarse) <- def_future_coarse
zoom(def_hist_coarse, e = bounds)
# Check
compareGeom(def_future_coarse, def_hist_coarse)

#-------------------------------------------------------------------------------

# DEM - elevation and associated metrics (HLI, TPI)
#-------------------------------------------------------------------------------
# Aggregation of 90-m DEM to 250-m (matching fine TopoFire data) complete, don't re-run:
# dem_tif <- rast('../../../hoecker/Work/GIS/DEM/dem90_hf.tif')
# dem_tif <- project(dem_tif, def_hist_fine)
# dem_fine <- resample(dem_tif, def_hist_fine)
# writeRaster(dem_fine, 'dem_250_topo.tif')

# Read in DEM that matches fine-scale climate data exactly
dem_fine <- rast('data/topofire/dem_250_topo.tif')
zoom(dem_fine, e = bounds)

# Coarse data
# Resample, as a form of aggregation
dem_coarse <- resample(dem_fine, def_future_coarse, 'bilinear')
# Not sure why this is necessary... but it is
crs(dem_coarse) <- def_future_coarse
zoom(dem_coarse, e = bounds)
# Yes - they align!
compareGeom(dem_coarse, def_future_coarse)

# DEM-based predictors
hli_coarse_conus <- hli(dem_coarse)
tpi_coarse_conus <- tpi(dem_coarse)

# Fine data
# This slow calculation is complete, don't re-run:
# hli_fine <- hli(dem_fine)
# writeRaster(hli_fine, 'data/topofire/hli_250_topo.tif')
# tpi_fine <- tpi(dem_fine)
# writeRaster(tpi_fine, 'data/topofire/tpi_250_topo.tif')
hli_fine <- rast('data/topofire/hli_250_topo.tif')
tpi_fine <- rast('data/topofire/tpi_250_topo.tif')


#-------------------------------------------------------------------------------

# Crop and stack coarse and fine predictor rasters
#-------------------------------------------------------------------------------
# Turn extent object into a polygon/vector
bound_object <- as.polygons(bounds, crs="+proj=longlat")

# Fine rasters
clim_fine <- crop(def_hist_fine, bounds)
elev_fine <- crop(dem_fine, bounds)
heat_fine <- crop(hli_fine, bounds)
topo_fine <- crop(tpi_fine, bounds)

# Coarse rasters
ds_coarse <- crop(def_future_coarse, bounds)
clim_coarse <- crop(def_hist_coarse, bounds) #crop(def_hist_coarse_fine, bounds)
hli_coarse <- crop(hli_coarse_conus, bounds)
tpi_coarse <- crop(tpi_coarse_conus, bounds)
elev_coarse <- crop(dem_coarse, bounds)

coarse_stack <- rast(list('ds' = ds_coarse,
                          'clim' = clim_coarse,
                          'elev' = elev_coarse,
                          'hli' = hli_coarse,
                          'tpi' = tpi_coarse))


fine_stack <- rast(list('clim' = clim_fine,
                        'elev' = elev_fine,
                        'hli' = heat_fine,
                        'tpi' = topo_fine))

# Convert raster to vector (points)
fine_pts <- as.points(clim_fine)

# Beginning from a small rectangular subset stack and centroids of each cell
#-------------------------------------------------------------------------------
downscale <- function(pt_id, fine_pts, fine_stack, coarse_stack){
  
  print(paste('Running point', pt_id, 'of', length(fine_pts), 'points'))
  
  # Pull out focal point
  focal_pt <- fine_pts[pt_id]
  # Plots for testing
  # plot(clim_fine)
  # plot(fine_pts[pt_id], add = T)

  # Create 1000-m buffer 
  focal_buff <- buffer(focal_pt, width = 15000)
  # Plot for testing
  # plot(focal_buff, add = T)
  
  # Extact predictor values within buffer
  fine_df <- terra::extract(fine_stack, focal_pt, cells = T, xy = T) 
  coarse_df <- terra::extract(coarse_stack, focal_buff, cells = T, xy = T)
  
  # Calculate distances between focal point and predictor points
  dists <- c(geosphere::distm(crds(focal_pt),coarse_df[,c('x','y')]))
  
  # Create index of points outside of nugget distance (the size of a coarse cell)
  nug_idx <- dists > 4000

  # Fine-grid information from fine-grid focal location
  X <- fine_df$x
  Y <- fine_df$y
  C <- fine_df$clim
  #E <- fine_df$elev
  #H <- fine_df$hli
  #P <- fine_df$tpi

  # Coarse-grid information from coarse grid points within buffer distance
  model_df <- data.frame(
    'Zi' = coarse_df$ds[nug_idx],
    'Xi' = coarse_df$x[nug_idx],
    'Yi' = coarse_df$y[nug_idx],
    'Ci' = coarse_df$clim[nug_idx],
    #'Ei' = coarse_df$elev[nug_idx],
    #'Hi' = coarse_df$hli[nug_idx],
    #'Pi' = coarse_df$tpi[nug_idx],
    # Distances between fine-grid focal location and coarse centroids
    'di' = dists[nug_idx]
    )
  
  # Remove incomplete cases... these are edges
  model_df <- na.omit(model_df)
  # Create a linear regression among coarse data
  # The lm() representation of the model: lm_mod <- lm(Zi ~ Xi + Yi + Ei + Hi + Pi)
  # Instead, use the base-level version to speed things up (thanks ChatGPT for translating!)
  # Fit the model using lm.fit() and model.matrix()
  x_lm <- model.matrix(Zi ~ Xi+Yi+Ci, data = model_df) # +Ei+Hi+Pi
  y_lm <- model_df$Zi
  lm_mod <- .lm.fit(x_lm, y_lm)
  
  # Extract coefficients
  Cx <- lm_mod$coefficients[2]
  Cy <- lm_mod$coefficients[3]
  Cc <- lm_mod$coefficients[4]
  Ce <- lm_mod$coefficients[5]
  Ch <- lm_mod$coefficients[6]
  Cp <- lm_mod$coefficients[7]
   
  # This forumula is provided on page 5 of Flint and Flint 2021
  # +(E-model_df$Ei)*Ce +(H-model_df$Hi)*Ch +(P-model_df$Pi)*Cp
  sum1 <- sum((model_df$Zi +(X-model_df$Xi)*Cx +(Y-model_df$Yi)*Cy +(C-model_df$Ci)*Cc) / model_df$di^2)
  sum2 <- sum(1/model_df$di^2)
  Z = sum1/sum2
  
  # Return as dataframe
  return(data.frame('value' = Z))
  
}

library(tictoc)
tic()
downscaled_df <- seq_along(1:length(fine_pts)) %>% 
  map_df(~ downscale(pt_id = .x, fine_pts, fine_stack, coarse_stack),.id = 'cell')
toc()
((expanse(def_hist_fine)/expanse(bound_object))*(204.57/60/60))/15

downscaled_rast <- cbind(crds(fine_pts), downscaled_df$value)
downscaled_rast <- rast(downscaled_rast, type = 'xyz', crs = crs(fine_pts))
plot(downscaled_rast)
