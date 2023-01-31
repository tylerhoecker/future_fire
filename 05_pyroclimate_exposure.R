# This step performs the Earth Mover's Distance calculation
# EMD results are written out as rasters
# EMD results are pair with fire resistance trait data, then are
# summarized into firesheds and written out

library(tidyverse) # Carpentry tools
library(furrr)     # Parralell computing
library(emdist)    # Earth Movers Distance
library(tictoc)
library(infer)

# ------------------------------------------------------------------------------
# Measure 1-dimensional Earth Mover's Distance / Wasserstein Metric
# ------------------------------------------------------------------------------
# Testing: 
# x = 10 
# ref_df_x <- readRDS(paste0('data/process/chunks/ref_df_list_',x,'.Rdata'))
# fut_df_x <- readRDS(paste0('data/process/chunks/fut_df_list_',x,'.Rdata'))
# cell_pairs_x <- readRDS(paste0('data/process/chunks/cell_pairs',x,'.Rdata'))
# index = 1; ref_df = ref_df_x[[index]]; fut_df = fut_df_x[[index]]; cell_pairs = cell_pairs_x[[index]]
rescale01 <- function(x){round(((x-min(x))/(max(x)-min(x))),3)}

calc_emds <- function(ref_df, fut_df, cell_pairs, index){
  
  print(paste0('Cell pair ',index, ': ', cell_pairs[1],',',cell_pairs[2]))
  
  # Sequence of conditions on output for scenarios with no change, no fire in reference bin, or no fire in future bin
  
  # Start with 0's in everything as default and update when conditions are met
  output_df <- data.frame(
    'cell_hist'= cell_pairs[['cell_hist']],
    'cell_2C'  = cell_pairs[['cell_2C']],
    "multi_emd_50" = 0,
    "multi_emd_10" = 0,
    "multi_emd_90" = 0,
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
      mutate(across(c(multi_emd_50, multi_emd_10, multi_emd_90, 
                      cbi_dir, ndvi_dir, frp_dir),
                    ~ case_when(. == 0 ~ NA)),
             case = 1) # 1 = NO FOREST FIRE REFERENCE, YES FIRE FUTURE
    
    # No observed fire in future pyroclimate
  } else if(nrow(ref_df) != 0 & nrow(fut_df) == 0) {
    output_df <- output_df %>% 
      mutate(across(c(multi_emd_50, multi_emd_10, multi_emd_90, 
                      cbi_dir, ndvi_dir, frp_dir),
                    ~ case_when(. == 0 ~ NA)),
             case = 2) # 2 = YES FOREST FIRE REFERENCE, NO FIRE FUTURE
    
  } else if(nrow(ref_df) == 0 & nrow(fut_df) == 0) {
    output_df <- output_df %>% 
      mutate(across(c(multi_emd_50, multi_emd_10, multi_emd_90, 
                      cbi_dir, ndvi_dir, frp_dir),
                    ~ case_when(. == 0 ~ NA)),
             case = 3) # 3 = NO FOREST FIRE REFERENCE OR FUTURE ("forest" but never was, never will be flammable)
    
    # Otherwise, there are data to compare... calculate EMD between ref and fut
  } else {
    
    # # Univariate EMD -----------------------------------------------------------
    # # Function that calculates EMD distance after reworking data into matrices
    # emd_er <- function(var){
    # 
    #   dat_mats <- map(list(ref_df, fut_df), function(df){
    #     dat <- df %>%
    #       pull(var) %>%
    #       density(., bw = .05)
    #     dat <- cbind(dat$x, dat$y)
    #     return(dat)
    #   })
    # 
    #   emd <- emdist::emd(dat_mats[[1]], dat_mats[[2]], max.iter = 5000) %>%
    #     round(., 5)
    # 
    #   return(emd)
    # }
    # 
    # # Apply EMD calculation function to each variable
    # emds <- map_df(list('cbi_emd'='z_cbi',
    #                     'ndvi_emd'='z_ndvi',
    #                     'fri_emd'='z_fri',
    #                     'frs_emd'='z_frs'),
    #                emd_er)
    # # So BASEic
    # output_df[,colnames(emds)] <- emds
    # # END Univariate EMD -----------------------------------------------------------
    
    # Multi-dimensional EMD ----------------------------------------------------
    # How many replicates of 500 to subsample for multi EMD calc
    n_reps <- 10

    # Function to prepare dataframes for EMD function
    rep_prep <- function(rep_df){
      rep_dens <- Compositional::mkde(as.matrix(rep_df), h = 1, thumb = 'none')
      rep_mat <- cbind(rep_dens,as.matrix(rep_df))
      return(rep_mat)
    }

    # Create list of replicated samples, prepared for EMD function
    sample_dfs <- list(ref_df, fut_df) %>%
      set_names('ref','fut') %>%
      map(function(df){
        df %>%
          dplyr::select(pca1, pca2) %>% 
          infer::rep_slice_sample(n = 500, reps = n_reps) %>%
          group_map(~rep_prep(.x))
      })

    # Calculate multi EMD for each replicate
    rep_emds <- c(1:n_reps) %>%
      map_dbl(function(rep){
        emdist::emd(sample_dfs[['ref']][[rep]], sample_dfs[['fut']][[rep]], max.iter = 500)
      })
   # END multi-dimensional EMD -------------------------------------------------
    
    # Calculate direction as well
    output_df <- output_df %>% 
      mutate(multi_emd_50 = quantile(rep_emds, 0.50, na.rm = T),
             multi_emd_10 = quantile(rep_emds, 0.10, na.rm = T),
             multi_emd_90 = quantile(rep_emds, 0.90, na.rm = T),
             cbi_dir  = (mean(fut_df$cbi) - mean(ref_df$cbi))/mean(ref_df$cbi),
             ndvi_dir = (mean(fut_df$ndvi) - mean(ref_df$ndvi))/mean(ref_df$ndvi),
             frp_dir  = mean(exp(fut_df$log_frp_pt)) - mean(exp(ref_df$log_frp_pt)),
             pca1_dir = mean(rescale01(fut_df$pca1)) - mean(rescale01(ref_df$pca1)),
             pca2_dir = mean(rescale01(fut_df$pca2)) - mean(rescale01(ref_df$pca2)),
             case = 4) # 4 = "TYPICAL" CASE - observations of fire in reference and future
  }
}


# Run the function over all pairs

# Use parralell processing
plan(multisession, workers = 15)

# Start timer
#tic()
# There are 40 chunks... run the function over all pairs in all chunks
emd_chunk_df <- c(1:40) %>%  #40
  map_df(function(x){
    print(paste0('Running chunk ', x))
    ref_df_x <- readRDS(paste0('data/process/chunks/ref_df_list_',x,'.Rdata'))
    fut_df_x <- readRDS(paste0('data/process/chunks/fut_df_list_',x,'.Rdata'))
    cell_pairs_x <- readRDS(paste0('data/process/chunks/cell_pairs',x,'.Rdata'))
    
    emd_df <- future_pmap_dfr(list(ref_df_x, fut_df_x, cell_pairs_x, seq_along(cell_pairs_x)), 
                              calc_emds)
    saveRDS(emd_df, paste0('data/process/complete/emd_df_',x,'.Rdata'))
    return(emd_df)
  })
# Stop timer
#toc()

plan(sequential)

# Read/combine the saved chunks
# emd_chunk_df <- c(1:40) %>%
#   map_df(function(x){
#     print(paste0('Reading chunk ', x))
#     emd_df <- readRDS(paste0('data/process/complete/emd_df_',x,'.Rdata'))
#     return(emd_df)
#   })
# 
saveRDS(emd_chunk_df, 'data/emd/emd_2023_01_30.Rdata')
#emd_chunk_df <- readRDS('data/emd/emd_2023_01_25.Rdata')

# ------------------------------------------------------------------------------
# Temp. additions until can add above - update cases to ID dry vs wet
# ------------------------------------------------------------------------------
# emd_chunk_df <- emd_chunk_df %>% 
#   rowwise() %>% 
#   mutate(across(c(cbi_dir, ndvi_dir, fri_dir), ~ ifelse(case %in% c(1,2,3), NA, .x))) %>% 
#   ungroup()

# ------------------------------------------------------------------------------
# Rasterize output 
# ------------------------------------------------------------------------------
library(terra)
library(raster)
library(sf)

# Read in the master grid as raster!
mast_rast <- raster('data/process/mast_rast.tif')

# Join to dataframe of reference and future climate bins for all forested cells
forest_pyrome_df <- readRDS('data/process/forest_pyrome_df.Rdata')
pyrome_emd_df <- full_join(forest_pyrome_df, emd_chunk_df) %>% 
  filter(!is.na(case))

# Convert to spatial object (sf)
pyrome_emd_sf <- sf::st_as_sf(pyrome_emd_df, coords = c('x', 'y'), crs = 4326) 
# %>% 
#   as( "Spatial") %>% 
#   vect()

# Write out as shapefile
# sf::st_write(pyrome_emd_sf, 'C:/Users/hoecker/GitHub/future_fire/data/pyrome/emd_forest_pts.shp',
#              append = FALSE)

# Rasterize and write out the TIFFs
list('multi_emd_50','cbi_dir','ndvi_dir','frp_dir','pca1_dir','pca2_dir','case') %>%  
  walk(function(x) {
    print(paste0('Rasterizing ',x,'...'))
    r <- terra::rasterize(pyrome_emd_sf, mast_rast, 
                          field = x, fun = mean, na.rm = T)
    print(paste0('Writing data/emd/',x,'.tiff'))
    writeRaster(r, filename = paste0('data/emd/forest_',x,'.tiff'), overwrite = T)
    return(r)
  }) 

# ------------------------------------------------------------------------------
# Summarize into firesheds 
# ------------------------------------------------------------------------------
library(sf)
# Join to dataframe of reference and future climate bins for all forested cells
forest_pyrome_df <- readRDS('data/process/forest_pyrome_df.Rdata')
emd_df <- emd_chunk_df #readRDS('data/emd/emd_2023_01_25.Rdata')

length(distinct(emd_df, cell_hist, cell_2C)$cell_hist)

forest_emd_df <- full_join(forest_pyrome_df, emd_df) %>% 
  filter(!is.na(multi_emd_50))

forest_emd_sf <- st_as_sf(forest_emd_df, coords = c('x', 'y'), crs = 4326) %>% 
  st_transform(., crs = 5070)

firesheds_sf <- st_read('data/ecoregions/firesheds_5070.shp') %>% 
  st_crop(., forest_emd_sf) %>% 
  st_cast()

firesheds_emd_sf <- st_join(forest_emd_sf, firesheds_sf, join = st_within) %>% 
  filter(!is.na(multi_emd_50), !is.na(Fireshed_I)) %>% 
  group_by(Fireshed_I) %>% 
  summarize(adapt_dist = mean(multi_emd_50),
            units = n()) %>% 
  st_drop_geometry() %>% 
  filter(units > 500)

# This is a bit weird but preserves multipolygon geometry which gets mixed up otherwise
firesheds_emd_poly <- firesheds_sf %>% 
  filter(Fireshed_I %in% firesheds_emd_sf$Fireshed_I) %>% 
  left_join(., firesheds_emd_sf)

# Add frs --------------
# FRS rast
frs_rast <- terra::rast('data/composition/frs_spatial.tif')
# Crop for speed
frs_rast <- terra::project(frs_rast, y = "epsg:5070")
# Firesheds now as Terra vector
firesheds_vect <- terra::vect('data/ecoregions/firesheds_5070.shp')
# Filter out firesheds in forested study area
firesheds_within <- firesheds_vect[firesheds_vect$Fireshed_I %in% unique(firesheds_emd_sf$Fireshed_I)]
# Extract and summarise fire resistance scores
firesheds_frs <- terra::extract(frs_rast, firesheds_within, fun = mean, na.rm = T) 

# Summarize adaptation distance within these firesheds and add FRS values
firesheds_emd_frs_sf <- firesheds_emd_poly %>% 
  # Divide into groups for bivariate chloropleth categorization
  mutate(frs = firesheds_frs$frs_spatial,
         cat1 = case_when(adapt_dist <= quantile(firesheds_emd_poly$adapt_dist, 0.33) ~ 'A',
                          adapt_dist > quantile(firesheds_emd_poly$adapt_dist, 0.33) & 
                          adapt_dist <= quantile(firesheds_emd_poly$adapt_dist, 0.66) ~ 'B',
                          adapt_dist > quantile(firesheds_emd_poly$adapt_dist, 0.66) ~ 'C'),
         cat2 = case_when(frs <= quantile(firesheds_frs$frs_spatial, 0.33) ~ '3',
                          frs > quantile(firesheds_frs$frs_spatial, 0.33) & 
                          frs <= quantile(firesheds_frs$frs_spatial, 0.66) ~ '2',
                          frs > quantile(firesheds_frs$frs_spatial, 0.66) ~ '1')) %>% 
  unite('bicat', cat1, cat2) %>% 
  dplyr::select(Area_HA, Fireshed_I, Fireshed_N, adapt_dist, frs, bicat)
  
st_write(firesheds_emd_frs_sf, 'data/emd/bivariate_firesheds.shp', append= F)
















