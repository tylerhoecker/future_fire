# ------------------------------------------------------------------------------
# Description: The core of the analysis is done here, calculating dissimilarity
# between contemporary and projected future fire regime attributes using the 
# earth mover's distance.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load packages
# ------------------------------------------------------------------------------
library(tidyverse) # Carpentry tools
library(furrr)     # Parralell computing
library(emdist)    # Earth Movers Distance
library(infer)
library(Compositional)
library(sf)

calc_emds <- function(ref_df, fut_df, cell_pairs, index){
  
  print(paste0('Cell pair ',index, ': ', cell_pairs[1],',',cell_pairs[2]))
  
  # Sequence of conditions on output for scenarios with no change, no fire in reference bin, or no fire in future bin
  
  # Start with 0's in everything as default and update when conditions are met
  output_df <- data.frame(
    'cell_hist'= cell_pairs[['cell_hist']],
    'cell_2C'  = cell_pairs[['cell_2C']],
    emd_mean = 0,
    emd_05 = 0,
    emd_10 = 0,
    emd_90 = 0,
    emd_95 = 0,
    null_mean = 0,
    null_05 = 0,
    null_10 = 0,
    null_90 = 0,
    null_95 = 0,
    null_sum = 0,
    cbi_dir  = 0,
    cbi_pct  = 0,
    ndvi_dir = 0,
    ndvi_pct = 0,
    frp_dir  = 0,
    frp_pct  = 0,
    case = 0
  )
  
  # Skip pairs where reference and future are the same cell, keep default output
  if(output_df$cell_hist == output_df$cell_2C){
    return(output_df)
    
    # No observed fire in reference pyroclimate
  } else if(nrow(ref_df) == 0 & nrow(fut_df) != 0) {
    output_df <- output_df %>%  
      mutate(across(c(emd_mean,emd_05,emd_10,emd_90,emd_95,
                      null_mean,null_05,null_10,null_90,null_95,null_sum,
                      cbi_dir,cbi_pct,ndvi_dir,ndvi_pct, frp_dir, frp_pct),
                    ~ case_when(. == 0 ~ NA)),
             case = 1) # 1 = NO FOREST FIRE REFERENCE, YES FIRE FUTURE
    
    # No observed fire in future pyroclimate
  } else if(nrow(ref_df) != 0 & nrow(fut_df) == 0) {
    output_df <- output_df %>% 
      mutate(across(c(emd_mean,emd_05,emd_10,emd_90,emd_95,
                      null_mean,null_05,null_10,null_90,null_95,null_sum,
                      cbi_dir,cbi_pct,ndvi_dir,ndvi_pct, frp_dir, frp_pct),
                    ~ case_when(. == 0 ~ NA)),
             case = 2) # 2 = YES FOREST FIRE REFERENCE, NO FIRE FUTURE
    
  } else if(nrow(ref_df) == 0 & nrow(fut_df) == 0) {
    output_df <- output_df %>% 
      mutate(across(c(emd_mean,emd_05,emd_10,emd_90,emd_95,
                      null_mean,null_05,null_10,null_90,null_95,null_sum,
                      cbi_dir,cbi_pct,ndvi_dir,ndvi_pct, frp_dir, frp_pct),
                    ~ case_when(. == 0 ~ NA)),
             case = 3) # 3 = NO FOREST FIRE REFERENCE OR FUTURE ("forest" but never was, never will be flammable)
    
    # Otherwise, there are data to compare... calculate EMD between ref and fut
  } else {
    
    # Multi-dimensional EMD ----------------------------------------------------
    # How many replicates of 500 to subsample for multi EMD calc
    n_reps <- 100
    n_samp <- 500

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
          # Changed from PCA to Z-scores here on 4/19
          # dplyr::select(pca1, pca2) %>% 
          dplyr::select(z_cbi, z_ndvi, z_frp_pt) %>% 
          infer::rep_slice_sample(n = n_samp, reps = n_reps) %>%
          group_map(~rep_prep(.x))
      })
    
    null_dfs <- ref_df %>%
          # Changed from PCA to Z-scores here on 4/19
          # dplyr::select(pca1, pca2) %>% 
          dplyr::select(z_cbi, z_ndvi, z_frp_pt) %>% 
          infer::rep_slice_sample(n = n_samp, reps = n_reps) %>%
          group_map(~rep_prep(.x))

    # Visualize examples    
    # ggplot() +
    #   stat_ecdf(aes(x = sample_dfs[['ref']][[rep]][,c('z_cbi')]), size = 1, color = 'red') +
    #   stat_ecdf(aes(x = sample_dfs[['fut']][[rep]][,c('z_cbi')]), size = 1, color = 'blue') 
    

    # Calculate multi EMD for each replicate
    rep_emds <- c(1:n_reps) %>%
      map_dbl(function(rep){
        emdist::emd(sample_dfs[['ref']][[rep]][,c('z_cbi','z_ndvi','z_frp_pt')], 
                    sample_dfs[['fut']][[rep]][,c('z_cbi','z_ndvi','z_frp_pt')], 
                    max.iter = 500)
      })
    
    null_emds <- c(1:n_reps) %>%
      map_dbl(function(rep){
        emdist::emd(sample_dfs[['ref']][[rep]][,c('z_cbi','z_ndvi','z_frp_pt')], 
                    null_dfs[[rep]][,c('z_cbi','z_ndvi','z_frp_pt')], 
                    max.iter = 500)
      })
    
    
    # EMD multi-dimensional EMD ------------------------------------------------
    # Calculate direction as well
    output_df <- output_df %>% 
      mutate(emd_mean = mean(rep_emds),
             emd_05 = quantile(rep_emds, 0.05),
             emd_10 = quantile(rep_emds, 0.10),
             emd_90 = quantile(rep_emds, 0.90),
             emd_95 = quantile(rep_emds, 0.95),
             null_mean = mean(null_emds),
             null_05 = quantile(null_emds, 0.05),
             null_10 = quantile(null_emds, 0.10),
             null_90 = quantile(null_emds, 0.90),
             null_95 = quantile(null_emds, 0.95),
             null_sum = sum(null_emds > min(rep_emds)),
             cbi_dir  = (mean(fut_df$cbi) - mean(ref_df$cbi)),
             cbi_pct  = (mean(fut_df$cbi) - mean(ref_df$cbi))/mean(ref_df$cbi),
             ndvi_dir = (mean(fut_df$ndvi) - mean(ref_df$ndvi)),
             ndvi_pct = (mean(fut_df$ndvi) - mean(ref_df$ndvi))/mean(ref_df$ndvi),
             frp_dir  = mean(fut_df$frp_pt) - mean(ref_df$frp_pt),
             frp_pct  = mean(fut_df$frp_pt) - mean(ref_df$frp_pt)/mean(ref_df$frp_pt),
             case = 4) # 4 = "TYPICAL" CASE - observations of fire in reference and future
  }
}

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
saveRDS(emd_chunk_df, 'data/emd/emd_2023_06_13.Rdata')

plan(sequential)

# Read/combine the saved chunks
emd_chunk_df <- c(1:40) %>%
  map_df(function(x){
    print(paste0('Reading chunk ', x))
    emd_df <- readRDS(paste0('data/process/complete/emd_df_',x,'.Rdata'))
    return(emd_df)
  })

emd_chunk_df <- readRDS('data/emd/emd_2023_06_13.Rdata')

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

# Write out as shapefile
# sf::st_write(pyrome_emd_sf, 'C:/Users/hoecker/GitHub/future_fire/data/pyrome/emd_forest_pts.shp',
#              append = FALSE)

# Rasterize and write out the TIFFs
list('emd_mean','null_mean','null_sum','cbi_dir','ndvi_dir','frp_dir') %>%  #,'multi_emd_50','cbi_dir','ndvi_dir','frp_dir','pca1_dir','pca2_dir','case'
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
# Join to dataframe of reference and future climate bins for all forested cells
forest_pyrome_df <- readRDS('data/process/forest_pyrome_df.Rdata')
emd_df <- readRDS('data/emd/emd_2023_06_13.Rdata') #emd_chunk_df 

length(distinct(emd_df, cell_hist, cell_2C)$cell_hist)

forest_emd_df <- full_join(forest_pyrome_df, emd_df) %>% 
  filter(!is.na(emd_mean))

forest_emd_sf <- st_as_sf(forest_emd_df, coords = c('x','y'), crs = 4326) %>% 
  st_transform(., crs = 5070)

# Firesheds data are available here: https://www.fs.usda.gov/research/rmrs/projects/firesheds#overview
firesheds_sf <- st_read('data/ecoregions/firesheds_5070.shp') %>% 
  st_crop(., forest_emd_sf) %>% 
  st_cast()

firesheds_emd_sf <- st_join(forest_emd_sf, firesheds_sf, join = st_within) %>% 
  filter(!is.na(emd_mean), !is.na(Fireshed_I)) 

firesheds_info_sf <- firesheds_emd_sf |> 
  group_by(Fireshed_I) %>% 
  summarize(adapt_dist  = median(emd_mean),
            null_dist   = median(null_mean),
            null_pval   = median(null_sum),
            adapt_05    = median(emd_05),
            adapt_95    = median(emd_95),
            null_05     = median(null_05),
            null_95     = median(null_95),
            cbi_change  = median(cbi_dir),
            frp_change  = median(frp_dir),
            ndvi_change = median(ndvi_dir),
            cbi_pct  = median(cbi_pct),
            frp_pct  = median(frp_pct),
            ndvi_pct = median(ndvi_pct),
            units       = n()) %>% 
  st_drop_geometry() %>% 
  filter(units > 500)

# This is a bit weird but preserves multipolygon geometry which gets mixed up otherwise
firesheds_emd_poly <- firesheds_sf %>% 
  filter(Fireshed_I %in% firesheds_info_sf$Fireshed_I) %>% 
  left_join(., firesheds_info_sf)


# Add frs --------------
# FRS rast. These data are availabe here: https://www.sciencebase.gov/catalog/item/5e39f54ee4b0a79317e15e73
frs_rast <- terra::rast('data/composition/frs_spatial.tif')
# Crop for speed
frs_rast <- terra::project(frs_rast, y = "epsg:5070")
# Firesheds now as Terra vector
firesheds_vect <- terra::vect('data/ecoregions/firesheds_5070.shp')
# Filter out firesheds in forested study area
firesheds_within <- firesheds_vect[firesheds_vect$Fireshed_I %in% unique(firesheds_info_sf$Fireshed_I)]
# Extract and summarise fire resistance scores
firesheds_frs <- terra::extract(frs_rast, firesheds_within, fun = mean, na.rm = T) 

# Summarize adaptation distance within these firesheds and add FRS values
firesheds_emd_frs_sf <- firesheds_emd_poly %>% 
  # Divide into groups for bivariate chloropleth categorization
  mutate(frs = firesheds_frs$frs_spatial,
         cat1 = case_when(adapt_dist <= quantile(firesheds_emd_poly$adapt_dist, 0.33, na.rm = T) ~ 'A',
                          adapt_dist > quantile(firesheds_emd_poly$adapt_dist, 0.33, na.rm = T) & 
                          adapt_dist <= quantile(firesheds_emd_poly$adapt_dist, 0.66, na.rm = T) ~ 'B',
                          adapt_dist > quantile(firesheds_emd_poly$adapt_dist, 0.66, na.rm = T) ~ 'C'),
         cat2 = case_when(frs <= quantile(firesheds_frs$frs_spatial, 0.33, na.rm = T) ~ '3',
                          frs > quantile(firesheds_frs$frs_spatial, 0.33, na.rm = T) & 
                          frs <= quantile(firesheds_frs$frs_spatial, 0.66, na.rm = T) ~ '2',
                          frs > quantile(firesheds_frs$frs_spatial, 0.66, na.rm = T) ~ '1')) %>%
  unite('bicat', cat1, cat2) %>% 
  rowwise() |> 
  mutate(null_overlap = ifelse(null_05 <= adapt_95 && adapt_05 >= null_95, 1, 0)) |> 
  ungroup() |> 
  mutate('numcat' = as.numeric(factor(bicat))) |> 
  dplyr::select(Area_HA, Fireshed_I, Fireshed_N, 
                adapt_dist, adapt_05, adapt_95, 
                null_overlap,
                null_dist, null_pval, null_05, null_95, 
                cbi_change, cbi_pct, frp_change, frp_pct, ndvi_change, ndvi_pct,
                frs, bicat, numcat)
  
st_write(firesheds_emd_frs_sf, 'data/emd/bivariate_firesheds.shp', append= F)
