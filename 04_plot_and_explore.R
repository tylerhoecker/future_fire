library(tidyverse)
library(raster)
library(sf)
library(tmap)

# ------------------------------------------------------------------------------
# Plot and explore
# ------------------------------------------------------------------------------
emd_stack <- list.files('data/pyrome/', pattern = '.*.tiff$', full.names = T) %>% 
  stack()
  
  
boundary <- read_sf("C:\\Users\\hoecker\\Work\\GIS\\cb_2018_us_state_20m\\cb_2018_us_state_20m.shp") %>%
  st_transform(., crs = 4326)

tm_shape(boundary, bbox = emd_stack) +
  tm_borders() +
  tm_fill(col = 'grey10') +
  tm_shape(emd_stack) +
  tm_raster() +
  tm_facets(free.scales = TRUE) 

stars::write_stars(pyrome_emd_rast, layer =  'cbi_emd', dsn = 'pyrome_emd.tiff')

