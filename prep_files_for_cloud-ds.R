# ------------------------------------------------------------------------------
# Prepared 2/27/2023 by Tyler Hoecker: https://github.com/tylerhoecker
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# This script prepares the geospatial layers needed to downscale future climate projections
#
# 250-m climate data from TopoFire (Holden et al.) will be used to downscale
# 4-km projected data, which was downscaled and bias-corrected from GCM output 
# by TerraClimate (Abatzoglou et al.)
# The downscaling process will be conducted in a cloud computing environment (CyVerse).

# The required files, all in cloud-optimized geotiff (COG) format, are:
# - Fine-scale "template" (in this case, 250-m climate data from TopoFire [Holden et al.])
# - Coarse-scale "template" that will be regressed against coarse-scale target layer
# - Coarse-scale "target" layer that will be downscaled
# - Multiply above by the number of climate variables / target layers to be downscaled.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load packages
using <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs,require,character.only=TRUE))
  need <- libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}
using('terra','gdalUtilities','purrr')
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Data paths (may not generalize)
template_path <- file.path('data','topofire')
ds_path <- file.path('data','terraclimate')
cogs_path <- file.path('data','cogs')

# Layer names
template_suffix <- '_1981_2010.tif'
ds_suffix <- '_2C.1961_1990.tif'

# Climate variables (should be one for each of template and ds)
vars <- list('def','aet')
# Climate variable of interest (may not generalize...)
vars |> 
  map(function(var){
    # Fine-scale template
    template_fine <- file.path(template_path,paste0(var,template_suffix))
    names(template_fine) <- names(rast(template_fine))
    # Raster that will be downscaled.
    ds_coarse <- file.path(ds_path,paste0(var,ds_suffix))
    # Crop it and save a temporary file
    temp_path <- file.path(cogs_path,paste0(var,'temp.tif'))
    writeRaster(crop(rast(ds_coarse), rast(template_fine)),temp_path)
    # Reset the ds path to this
    ds_coarse <- temp_path
    names(ds_coarse) <- names(rast(ds_coarse))
    
    list(template_fine, ds_coarse) |> 
      walk(~ gdal_translate(
        src_dataset = .x,
        dst_dataset = file.path(cogs_path,paste0(names(.x),'_cog.tif')),
        co = matrix(c("TILED=YES","COPY_SRC_OVERVIEWS=YES","COMPRESS=DEFLATE"),ncol = 1)
      ))
    unlink(temp_path)
  })



  




#----------------------
# Coarse template
# Resample, as a form of aggregation, the fine historical data to the desired coarse grid
template_coarse <- resample(template_fine, ds_coarse, 'bilinear')
# Not sure why this is necessary... but it is
crs(template_coarse) <- ds_coarse

# Check
compareGeom(template_coarse, ds_coarse)

