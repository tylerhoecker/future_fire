library(tidyverse) # Carpentry tools
library(tidymodels)
library(furrr)     # Parralell computing
library(emdist)    # Earth Movers Distance
library(tictoc)

# ------------------------------------------------------------------------------
# Measure 1-dimensional Earth Mover's Distance / Wasserstein Metric
# ------------------------------------------------------------------------------
# Testing: index = 10; ref_df = ref_df_x[[index]]; fut_df = fut_df_x[[index]]; cell_pairs = cell_pairs_x[[index]]

calc_emds <- function(ref_df, fut_df, cell_pairs, index){
  
  print(paste0('Cell pair ',index, ': ', cell_pairs[1],',',cell_pairs[2]))
  
  # Sequence of conditions on output for scenarios with no change, no fire in reference bin, or no fire in future bin
  # Start with 0's in everything and update when conditions are met
  output_df <- data.frame(
    'cell_hist'= cell_pairs[['cell_hist']],
    'cell_2C'  = cell_pairs[['cell_2C']],
    "multi_emd_50" = 0,
    "multi_emd_10" = 0,
    "multi_emd_90" = 0,
    'cbi_emd'  = 0,
    'ndvi_emd' = 0,
    'fri_emd'  = 0,
    'frs_emd'  = 0,
    'cbi_dir'  = 0,
    'ndvi_dir' = 0,
    'fri_dir'  = 0,
    'frs_dir'  = 0,
    'case'     = 0 # 0 = NO CHANGE
  )
  
  # Skip pairs where reference and future are the same cell
  if(output_df$cell_hist == output_df$cell_2C){
    return(output_df)
    
    # No observed fire in reference pyroclimate
  } else if(nrow(ref_df) == 0) {
    output_df <- output_df %>% 
      mutate(across(c(multi_emd_50,multi_emd_10,multi_emd_90,
                      cbi_emd,ndvi_emd,fri_emd,frs_emd,
                      cbi_dir,ndvi_dir,fri_dir,frs_dir),
                    ~ case_when(. == 0 ~ NA)),
             case = 1) # 1 = NO REFERENCE
    
    # No observed fire in future pyroclimate
  } else if(nrow(fut_df) == 0) {
    output_df <- output_df %>% 
      mutate(across(c(multi_emd_50,multi_emd_10,multi_emd_90,
                      cbi_emd,ndvi_emd,fri_emd,frs_emd,
                      cbi_dir,ndvi_dir,fri_dir,frs_dir),
                    ~ case_when(. == 0 ~ NA)),
             case = 2) # 2 = '0 future'
    
    # Otherwise, there are data to compare... calculate EMD between ref and fut
  } else {
    
    # Function that calculates EMD distance after reworking data into matrices
    emd_er <- function(var){
      
      #print(var)
      
      dat_mats <- map(list(ref_df, fut_df), function(df){
        dat <- df %>%
          pull(var) %>%
          density(., bw = .05)
        dat <- cbind(dat$x, dat$y)
        return(dat)
      })
      
      emd <- emdist::emd(dat_mats[[1]], dat_mats[[2]], max.iter = 5000) %>%
        round(., 5)
      
      return(emd)
    }
    
    # Apply EMD calculation function to each variable
    emds <- map_df(list('cbi_emd'='z_cbi',
                        'ndvi_emd'='z_ndvi',
                        'fri_emd'='z_fri',
                        'frs_emd'='z_frs'), 
                   emd_er)
    # So BASEic
    output_df[,colnames(emds)] <- emds
    
    
    # Multi-dimensional EMD ----------------------------------------------------
    # # How many replicates of 500 to subsample for multi EMD calc
    # n_reps <- 5
    # 
    # # Function to prepare dataframes for EMD function
    # rep_prep <- function(rep_df){
    #   rep_dens <- Compositional::mkde(as.matrix(rep_df), h = 1, thumb = 'none')
    #   rep_mat <- cbind(rep_dens,as.matrix(rep_df))
    #   return(rep_mat)
    # }
    # 
    # # Create list of replicated samples, prepared for EMD function 
    # sample_dfs <- list(ref_df, fut_df) %>%
    #   set_names('ref','fut') %>% 
    #   map(function(df){
    #     df %>% 
    #       dplyr::select(z_cbi, z_ndvi, z_fri, z_frs) %>% 
    #       rep_slice_sample(n = 500, reps = n_reps) %>% 
    #       group_map(~rep_prep(.x))
    #   })
    # 
    # # Calculate multi EMD for each replicate
    # rep_emds <- c(1:n_reps) %>% 
    #   map_dbl(function(rep){
    #     emdist::emd(sample_dfs[['ref']][[rep]], sample_dfs[['fut']][[rep]], max.iter = 500)
    #   })
   # --------------------------------------------------------------------------
    
    # Calculate direction as well. Could probably `map` this but leaving for now...
    output_df <- output_df %>% 
      mutate(#multi_emd_50 = quantile(rep_emds, 0.50, na.rm = T),
             #multi_emd_10 = quantile(rep_emds, 0.10, na.rm = T),
             #multi_emd_90 = quantile(rep_emds, 0.90, na.rm = T),
             cbi_dir  = round(mean(fut_df$cbi) - mean(ref_df$cbi), 2),
             ndvi_dir = round(mean(fut_df$ndvi) - mean(ref_df$ndvi), 2),
             fri_dir  = round(mean(fut_df$fri) - mean(ref_df$fri), 2),
             frs_dir  = round(mean(fut_df$frs) - mean(ref_df$frs), 2),
             case = 3)
  }
}


# Run the function over all pairs

# Use parralell processing
plan(multisession, workers = 15)

# Start timer
tic()
# There are 40 chunks... run the function over all pairs in all chunks
emd_chunk_df <- c(1:40) %>% 
  map_df(function(x){
    print(paste0('Running chunk ', x))
    ref_df_x <- readRDS(paste0('data/process/chunks/ref_df_list_',x,'.Rdata'))
    fut_df_x <- readRDS(paste0('data/process/chunks/fut_df_list_',x,'.Rdata'))
    cell_pairs_x <- readRDS(paste0('data/process/chunks/cell_pairs',x,'.Rdata'))
    
    emd_df <- future_pmap_dfr(list(ref_df_x, fut_df_x, cell_pairs_x, seq_along(cell_pairs_x)), 
                              calc_emds)
    saveRDS(emd_df, paste0('data/process/complete/uni_emd_df_',x,'.Rdata'))
    return(emd_df)
  })
toc()

plan(sequential)

# Read/combine the saved chunks
emd_chunk_df <- c(1:40) %>%
  map_df(function(x){
    print(paste0('Reading chunk ', x))
    emd_df <- readRDS(paste0('data/process/complete/uni_emd_df_',x,'.Rdata'))
    return(emd_df)
  })

saveRDS(emd_chunk_df, 'data/emd/uni_emd_2022_03_24.Rdata')

emd_chunk_df <- readRDS('data/emd/emd_2022_03_24.Rdata')

# Temp fix
emd_chunk_df <- emd_chunk_df %>% 
  rowwise() %>% 
  mutate(multi_emd_50 = ifelse(case %in% c(1,2), NA, multi_emd_50))

# ------------------------------------------------------------------------------
# Rasterize output 
# ------------------------------------------------------------------------------
library(terra)
library(raster)
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
'multi_emd_50','cbi_emd','ndvi_emd','fri_emd','frs_emd',
'cbi_dir','ndvi_dir','fri_dir','frs_dir'
list('case') %>%  
  walk(function(x) {
    print(paste0('Rasterizing ',x,'...'))
    r <- terra::rasterize(pyrome_emd_sf, mast_rast, 
                          field = 'case', fun = min, na.rm = T)
    print(paste0('Writing data/emd/',x,'.tiff'))
    writeRaster(r, filename = paste0('data/emd/forest_','case_min','.tiff'), overwrite = T)
    return(r)
  }) 


