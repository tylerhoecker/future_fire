library(tidyverse) # Carpentry tools
library(furrr)     # Parralell computing


# ------------------------------------------------------------------------------
# The pyroclimate analog comparison without EMD 
# ------------------------------------------------------------------------------
rescale01 <- function(x){round(((x-min(x))/(max(x)-min(x))),3)}

calc_emds <- function(ref_df, fut_df, cell_pairs, index){
  
  print(paste0('Cell pair ',index, ': ', cell_pairs[1],',',cell_pairs[2]))
  
  # Sequence of conditions on output for scenarios with no change, no fire in reference bin, or no fire in future bin
  
  # Start with 0's in everything as default and update when conditions are met
  output_df <- data.frame(
    'cell_hist'= cell_pairs[['cell_hist']],
    'cell_2C'  = cell_pairs[['cell_2C']],
    "cbi_ref"  = 0,
    "ndvi_ref" = 0,
    "frp_ref"  = 0,
    "pca1_ref" = 0,
    "pca2_ref" = 0,
    "cbi_fut"  = 0,
    "ndvi_fut" = 0,
    "frp_fut"  = 0,
    "pca1_fut" = 0,
    "pca2_fut" = 0,
    "cbi_dir"  = 0,
    "ndvi_dir" = 0,
    "frp_dir"  = 0,
    "pca1_dir" = 0,
    "pca2_dir" = 0,
    'case'     = 0 # 0 = NO CHANGE
  )
  
  # Skip pairs where reference and future are the same cell, keep default output
  if(output_df$cell_hist == output_df$cell_2C){
    return(output_df)
    
    # No observed fire in reference pyroclimate
  } else if(nrow(ref_df) == 0 & nrow(fut_df) != 0) {
    output_df <- output_df %>%  
      mutate(across(c(cbi_ref, ndvi_ref, frp_ref, 
                      cbi_fut, ndvi_fut, frp_fut,
                      cbi_dir, ndvi_dir, frp_dir),
                    ~ case_when(. == 0 ~ NA)),
             case = 1) # 1 = NO FOREST FIRE REFERENCE, YES FIRE FUTURE
    
    # No observed fire in future pyroclimate
  } else if(nrow(ref_df) != 0 & nrow(fut_df) == 0) {
    output_df <- output_df %>% 
      mutate(across(c(cbi_ref, ndvi_ref, frp_ref, 
                      cbi_fut, ndvi_fut, frp_fut,
                      cbi_dir, ndvi_dir, frp_dir),
                    ~ case_when(. == 0 ~ NA)),
             case = 2) # 2 = YES FOREST FIRE REFERENCE, NO FIRE FUTURE
    
  } else if(nrow(ref_df) == 0 & nrow(fut_df) == 0) {
    output_df <- output_df %>% 
      mutate(across(c(cbi_ref, ndvi_ref, frp_ref, 
                      cbi_fut, ndvi_fut, frp_fut,
                      cbi_dir, ndvi_dir, frp_dir),
                    ~ case_when(. == 0 ~ NA)),
             case = 3) # 3 = NO FOREST FIRE REFERENCE OR FUTURE ("forest" but never was, never will be flammable)
    
    # Otherwise, there are data to compare... calculate EMD between ref and fut
  } else {
    
    # Calculate direction as well
    output_df <- output_df %>% 
      mutate(cbi_ref = mean(ref_df$cbi), 
             ndvi_ref = mean(ref_df$ndvi), 
             frp_ref = mean(exp(ref_df$log_frp_cell)), 
             cbi_fut = mean(fut_df$cbi), 
             ndvi_fut = mean(fut_df$ndvi), 
             frp_fut = mean(exp(fut_df$log_frp_cell)),
             cbi_dir  = (mean(fut_df$cbi) - mean(ref_df$cbi))/mean(ref_df$cbi),
             ndvi_dir = (mean(fut_df$ndvi) - mean(ref_df$ndvi))/mean(ref_df$ndvi),
             frp_dir  = mean(exp(fut_df$log_frp_cell)) - mean(exp(ref_df$log_frp_cell)),
             pca1_dir = mean(rescale01(fut_df$pca1)) - mean(rescale01(ref_df$pca1)),
             pca2_dir = mean(rescale01(fut_df$pca2)) - mean(rescale01(ref_df$pca2)),
             case = 4) # 4 = "TYPICAL" CASE - observations of fire in reference and future
  }
}

emd_chunk_df <- c(1:40) %>%  #40
  map_df(function(x){
    print(paste0('Running chunk ', x))
    ref_df_x <- readRDS(paste0('data/process/chunks/ref_df_list_',x,'.Rdata'))
    fut_df_x <- readRDS(paste0('data/process/chunks/fut_df_list_',x,'.Rdata'))
    cell_pairs_x <- readRDS(paste0('data/process/chunks/cell_pairs',x,'.Rdata'))
    
    spatial_df <- future_pmap_dfr(list(ref_df_x, fut_df_x, cell_pairs_x, seq_along(cell_pairs_x)), 
                              calc_emds)
    saveRDS(spatial_df, paste0('data/process/complete/spatial_df_',x,'.Rdata'))
    return(spatial_df)
  })

# Write to spatial objects
# Read in the master grid as raster!
mast_rast <- raster::raster('data/process/mast_rast.tif')

# Join to dataframe of reference and future climate bins for all forested cells
forest_pyrome_df <- readRDS('data/process/forest_pyrome_df.Rdata')
pyrome_emd_df <- full_join(forest_pyrome_df, emd_chunk_df) %>% 
  filter(!is.na(case))

# Convert to spatial object (sf)
pyrome_emd_sf <- sf::st_as_sf(pyrome_emd_df, coords = c('x', 'y'), crs = 4326) 

list("cbi_ref","ndvi_ref","frp_ref", "cbi_fut", "ndvi_fut","frp_fut") %>%  
  walk(function(x) {
    print(paste0('Rasterizing ',x,'...'))
    r <- terra::rasterize(pyrome_emd_sf, mast_rast, 
                          field = x, fun = mean, na.rm = T)
    print(paste0('Writing data/emd/',x,'.tiff'))
    terra::writeRaster(r, filename = paste0('data/emd/forest_',x,'.tiff'), overwrite = T)
    return(r)
  }) 

# ------------------------------------------------------------------------------
# CALCULATE BURNED FOR REFERENCE AND FUTURE PERIODS FOR COMPARISON 
# ------------------------------------------------------------------------------

# Read in FRP 
library(terra)
frp_ref <- rast('data/emd/forest_frp_ref.tiff')
plot(frp_ref)
frp_fut <- rast('data/emd/forest_frp_fut.tiff')
plot(frp_fut)

frp_diff <- frp_fut - frp_ref
plot(frp_diff)
hist(frp_diff)

# Calculate burned area reference and future from these
mast_rast <- rast('data/process/mast_rast.tif')
forest_pts_sp <- as.points(mast_rast)
cell_area_df <- terra::extract(cellSize(mast_rast), forest_pts_sp)

frp_ref_df <- terra::extract(frp_ref, forest_pts_sp) %>% 
  as_tibble() %>% 
  mutate(area = cell_area_df$area,
         prop_ref = 36/forest_frp_ref,
         area_burned = area*prop_ref) %>%  
  filter(forest_frp_ref > 0) 
(sum(frp_ref_df$area_burned)*0.000001)

abatzoglous_obs <- read_csv('data/forest/abatzoglou_comparison.csv') %>% 
  select(area) %>% 
  sum()

(sum(frp_ref_df$area_burned)*0.000001)
abatzoglous_obs/36

(abatzoglous_obs - sum(frp_ref_df$area_burned)*0.000001)/abatzoglous_obs

frp_fut_df <- terra::extract(frp_fut, forest_pts_sp)

# Future
frp_fut_df <- terra::extract(frp_fut, forest_pts_sp) %>% 
  as_tibble() %>% 
  mutate(area = cell_area_df$area,
         prop_fut = 35/forest_frp_fut,
         area_burned = area*prop_fut) %>%  
  filter(forest_frp_fut > 0) 
sum(frp_fut_df$area_burned)*0.000001




