# ------------------------------------------------------------------------------
# Prepared 2/27/2023 by Tyler Hoecker: https://github.com/tylerhoecker
#
# Adapted from Rodman, K C., et al. 2020. 
# A Changing Climate is Snuffing Out Post-Fire Recovery in Montane Forests. 
# Global Ecology and Biogeography.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
## This code performs a multi-step procedure that spatially downscales gridded climate
# data using GIDS (Gradient and Inverse Distance-Squared) of Nalder and Wein (1998)
# as described in Flint and Flint (2012)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Import  packages
# ------------------------------------------------------------------------------
using <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs,require,character.only=TRUE))
  need <- libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}

using('purrr','furrr','terra','geosphere')
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------------------------
# Use correct paths for CyVerse...
cogs_path <- file.path('analyses','gids','cogs')
# For now local...
cogs_path <- file.path('data','cogs')

# Layer names
template_suffix <- '_1981_2010_cog.tif'
ds_suffix <- '_2C.1961_1990_cog.tif'

# Climate variables (should be one for each of template and ds)
vars <- list('def','aet')
# Climate variable of interest (may not generalize...)
vars |> 
  map(function(var){
    
    # Fine-scale template 
    template_fine <- rast(file.path(cogs_path,paste0(var,template_suffix)))
    
    # Raster that will be downscaled.
    ds_coarse <- rast(file.path(cogs_path,paste0(var,ds_suffix)))

    # Create a coarse (4 km) version, which will align with future data, to regress against the fine data
    # Resample, as a form of aggregation, the fine historical data to the desired coarse grid
    template_coarse <- resample(template_fine, ds_coarse, 'bilinear')
    # Not sure why this is necessary... but it is
    crs(template_coarse) <- ds_coarse
    # Check
    if(compareGeom(template_coarse, ds_coarse) == F){
      return(NULL)
    } 
    
    # The input data
    #fine <- rast('clim' = clim_fine)
    coarse <- rast(list('ds' = ds_coarse,'clim' = template_coarse))
    
    # Convert raster to vector (points)
    fine <- as.points(template_fine)
    
    # Beginning from a small rectangular subset stack and centroids of each cell
    #-------------------------------------------------------------------------------
    downscale <- function(pt_id){
      
      print(paste('Running point', pt_id, 'of', length(global_pts), 'points'))
      
      # Pull out focal point
      focal_pt <- global_pts[pt_id]
 
      # Create 1000-m buffer 
      focal_buff <- buffer(focal_pt, width = 15000)
      
      # Extract predictor values within buffer
      fine_df <- terra::extract(fine, focal_pt, cells = T, xy = T) 
      coarse_df <- terra::extract(coarse, focal_buff, cells = T, xy = T)
      
      # Calculate distances between focal point and predictor points
      dists <- c(distm(crds(focal_pt),coarse_df[,c('x','y')]))
      
      # Create index of points outside of nugget distance (the size of a coarse cell)
      nug_idx <- dists > 4000
      
      # Fine-grid information from fine-grid focal location
      X <- fine_df$x
      Y <- fine_df$y
      C <- fine_df$clim

      # Coarse-grid information from coarse grid points within buffer distance
      model_df <- data.frame(
        'Zi' = coarse_df$ds,
        'Xi' = coarse_df$x,
        'Yi' = coarse_df$y,
        'Ci' = coarse_df$clim,
        # Distances between fine-grid focal location and coarse centroids
        'di' = dists)[nug_idx] |> 
        # Remove incomplete cases... these are edges
        na.omit()
      
      # Create a linear regression among coarse data
      # Fit the model using lm.fit() and model.matrix(), faster than lm() 
      x_lm <- model.matrix(Zi ~ Xi+Yi+Ci, data = model_df) # +Ei+Hi+Pi
      y_lm <- model_df$Zi
      lm_mod <- .lm.fit(x_lm, y_lm)
      
      # Extract coefficients
      Cx <- lm_mod$coefficients[2]
      Cy <- lm_mod$coefficients[3]
      Cc <- lm_mod$coefficients[4]

      # This forumula is provided on page 5 of Flint and Flint 2021
      sum1 <- sum((model_df$Zi+(X-model_df$Xi)*Cx+(Y-model_df$Yi)*Cy+(C-model_df$Ci)*Cc) 
                  /model_df$di^2)
      sum2 <- sum(1/model_df$di^2)
      Z = sum1/sum2
      
      # Return as dataframe
      return(data.frame('value' = Z))
      
    }
    
    downscaled_df <- seq_along(1:length(fine_pts)) %>% 
      map_df(~ downscale(pt_id = .x, 
                         global_points = fine_pts, 
                         fine = fine_stack, 
                         coarse = coarse_stack))
    
    downscaled_rast <- cbind(crds(fine_pts), downscaled_df$value)
    downscaled_rast <- rast(downscaled_rast, type = 'xyz', crs = crs(fine_pts))
    
        
  })




# HOW TO SAVE JUST THE DOWNSCALED RAST?









