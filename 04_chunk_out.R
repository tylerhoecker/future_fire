library(tidyverse)

# ------------------------------------------------------------------------------
# Rescale pyroclimate dimensions, 
# ------------------------------------------------------------------------------
pyrome_df <- readRDS('data/process/pyrome_df.Rdata')

# Use this function to rescale data from 0-1, so comparisons among dimensions make sense
rescale0100 <- function(x){round(((x-min(x))/(max(x)-min(x)))*100,3)}

# Create a list of unique reference-future climate cell pairs (its 1723) - based on all possible forested climate bins
forest_pyrome_df <- readRDS('data/process/forest_pyrome_df.Rdata')

cell_pairs <- forest_pyrome_df %>% 
  dplyr::select(cell_hist, cell_2C) %>% 
  distinct() %>% 
  split(seq(nrow(.)))

rescale_pyrome_df <- pyrome_df %>% 
  dplyr::select(cell_hist, cell_2C, cbi, ndvi, log_fri, frs, pca1, pca2) %>% 
  mutate(across(c(cbi, ndvi, log_fri, frs, pca1, pca2), ~ rescale0100(.x), .names = 'z_{.col}'))

# ------------------------------------------------------------------------------
# Create and save chunks of data to optimize parallel processing 
# ------------------------------------------------------------------------------
# Add columns for all variables that are rescaled from 0-100
ref_df_list <- cell_pairs %>% 
  map(function(pair){
    df <- filter(rescale_pyrome_df, cell_hist == pair[['cell_hist']]) 
    return(df)
  })

fut_df_list <- cell_pairs %>% 
  map(function(pair){
    df <- filter(rescale_pyrome_df, cell_hist == pair[['cell_2C']]) 
    return(df)
  })

chunks <- seq(from = 1, to = length(cell_pairs), by = 100) 
chunks[length(chunks) + 1] <- length(cell_pairs)

for(i in seq_along(chunks)){
  
  if(i == length(chunks)) break
  
  saveRDS(ref_df_list[chunks[i]:chunks[i+1]], 
          paste0('data/process/chunks/ref_df_list_',i,'.Rdata'))
  
  saveRDS(fut_df_list[chunks[i]:chunks[i+1]], 
          paste0('data/process/chunks/fut_df_list_',i,'.Rdata'))
  
  saveRDS(cell_pairs[chunks[i]:chunks[i+1]],
          paste0('data/process/chunks/cell_pairs',i,'.Rdata'))
}


rm(ref_df_list, fut_df_list)
gc()

