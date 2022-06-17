library(tidyverse)
library(raster)
library(sf)
library(tmap)
library(ggridges)


# ------------------------------------------------------------------------------
# Read in EMD data
# ------------------------------------------------------------------------------
pyrome_df <- readRDS('data/process/pyrome_df.Rdata')
emd_chunk_df <- readRDS('data/emd/emd_2022_03_24.Rdata')

# emd_chunk_df %>% 
#   group_by(case) %>% 
#   summarize(climates = length(unique(cell_hist))) %>% 
#   as_tibble()

# There are a few thousands rows where reference and future data exist, but the ref-fut 
# did not occur in the observed fire dataset... I think... only explanation for NA in region I can think of
# Thus, left join.
burned_emd_df <- left_join(pyrome_df, emd_chunk_df)  


# Ecoregion names
ecoregs <- st_read('data/ecoregions/ecoregions_edc.shp') %>% 
  #dplyr::select(ECO_NUM, ECO_NAME) %>% 
  st_drop_geometry() %>% 
  mutate(region = paste0('r',as.character(Region_ID)))

burned_emd_df <- burned_emd_df %>% 
  full_join(., ecoregs)

# Sort regions by EMD
burned_emd_df <- burned_emd_df %>% 
  filter(ECO_NAME != 'Northern Great Plains Steppe') %>% 
  mutate(ECO_NAME = as.factor(ECO_NAME),
         ECO_NAME = fct_reorder(ECO_NAME, multi_emd_50, median),
         one_group = 1) 
# %>% 
#   filter(ECO_NAME %in% c('Southern Rocky Mountains', 'East Cascades - Modoc Plateau'))
  
# ------------------------------------------------------------------------------
# Plot and explore
# ------------------------------------------------------------------------------

# Ridge-style histograms
ggplot(burned_emd_df, aes(x = multi_emd_50, y = ECO_NAME, fill = stat(x))) +
  geom_density_ridges_gradient() +
  scale_fill_viridis_c(option = "C", guide = 'none') +
  xlim(0,50) +
  labs(y = '', x = 'Adaptation distance')

ggplot(burned_emd_df, aes(x = fri_dir, y = one_group, fill = stat(x))) +
  geom_density_ridges_gradient(bandwidth = 2) +
  scale_fill_gradient2(low = '#d7191c', mid = '#ffffbf', high = "#2c7bb6", guide = 'none') +
  #xlim(-200,200) +
  labs(y = '', x = 'Change in frequency')


ggplot(burned_emd_df, aes(x = case, y = region)) +
  geom_density_ridges() 


plot_df <- burned_emd_df %>% 
  group_by(region) %>% 
  slice_sample(n = 1000)


ggplot(plot_df) +
  geom_point(aes(x = ndvi_dir, y = cbi_dir, fill = cbi_emd)) +
  facet_wrap(~region) +
  #coord_cartesian(ylim = c(0,0.3)) +
  #scale_fill_brewer(palette = 1) +
  theme_bw()  





# ------------------------------------------------------------------------------
# Mapping
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

