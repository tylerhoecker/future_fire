library(tidyverse)
select <- dplyr::select



library(sf)
library(terra)
library(raster)




# ------------------------------------------------------------------------------
# Read in adaptation distance, ref and future climate, fire severity, 
# ------------------------------------------------------------------------------
adapt_dist <- readRDS('data/emd/emd_2022_06_28.Rdata') %>% 
  select(cell_hist, cell_2C, multi_emd_50, case) %>% 
  filter(!is.na(multi_emd_50))

pyrome_df <- readRDS('data/process/pyrome_df.Rdata') %>% 
  select(cell_hist, cell_2C, 
         aet_hist, aet_2C, aet_diff, def_hist, def_2C, def_diff, 
         frs, log_fri, cbi, ndvi, fires)

adapt_model_df <- left_join(adapt_dist, pyrome_df)

# ------------------------------------------------------------------------------
# Develop additional predictors
# ------------------------------------------------------------------------------

# Ecoregion

# Management agency 

# Hazardous fuels reduction

# Timber harvest








