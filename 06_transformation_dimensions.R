# Explores results in each pyroclimate dimension
library(tidyverse)
library(ggridges)
library(sf)
library(terra)
library(raster)
select <- dplyr::select

# ------------------------------------------------------------------------------
# EMD
# ------------------------------------------------------------------------------
emd_chunk_df <- readRDS('data/emd/emd_2023_01_30.Rdata')

# Join to dataframe of reference and future climate bins for all forested cells
forest_pyrome_df <- readRDS('data/process/forest_pyrome_df.Rdata')

# Just examine pyromes where exposure can be estimated (filter out shifts into or out of forest space)
pyrome_emd_df <- full_join(forest_pyrome_df, emd_chunk_df) %>%
  #filter(case == 0)
  filter(case %in% c(0,4))

# pyrome_emd_sf <- st_as_sf(pyrome_emd_df, coords = c('x', 'y'), crs = 4326)
# 
# ecoregs_sf <- st_read('data/ecoregions/ecoregions_edc.shp') %>%
#   st_transform(., crs = 4326) %>%
#   select(ECO_NAME)
# 
# pyrome_emd_df <- pyrome_emd_sf %>%
#   st_intersection(., ecoregs_sf) %>%
#   st_drop_geometry()
# 
# mean(pyrome_emd_df$multi_emd_50, na.rm = T)
# quant075 <- function(x) quantile(x, probs = 0.75, na.rm = T)

# Sort regions by EMD
pyrome_emd_df <- pyrome_emd_df %>% 
  mutate(#ECO_NAME = as.factor(ECO_NAME),
         #ECO_NAME = fct_reorder(ECO_NAME, multi_emd_50, quant075, .desc = T),
         one_group = 1) 

ggplot(pyrome_emd_df, aes(x = multi_emd_50, y = one_group, fill = stat(x))) +
  geom_density_ridges_gradient(bandwidth = .05) +
  #geom_vline(aes(xintercept = mean(multi_emd_50, na.rm = T)), linetype = 'dashed') +
  #scale_fill_viridis_c(option = "A", guide = 'none', limits = c(0,3), oob = scales::squish) +
  scale_fill_distiller(palette = 'Spectral', limits = c(0,1.25), oob = scales::squish,
                       guide = 'none') +
  coord_cartesian(xlim = c(-0.1,2.5)) +
  labs(x = 'Pyroclimate exposure')  +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = 'black'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(rect = element_rect(fill = "transparent"),
        panel.background = element_blank(),
        panel.grid = element_blank())

ggsave('emd_density.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs',
       height = 3, width = 6, dpi = 600)

# ------------------------------------------------------------------------------
# Severity change 
# ------------------------------------------------------------------------------
pyrome_emd_df %>% 
  filter(cbi_dir < -0.1) %>% 
  nrow()/nrow(pyrome_emd_df)

ggplot(pyrome_emd_df, aes(x = cbi_dir*100, y = one_group, fill = stat(x))) +
  geom_density_ridges_gradient(bandwidth = 2) +
  #geom_vline(aes(xintercept = 0), linetype = 'dashed') +
  #geom_vline(aes(xintercept = -0.3), linetype = 'dashed') +
  #geom_vline(aes(xintercept = 0.3), linetype = 'dashed') +
  scale_fill_distiller(palette = 'RdYlGn', limits = c(-20,20), oob = scales::squish,
                       guide = 'none') +
  #coord_cartesian(xlim = c(-0.6, 0.6)) +
  scale_x_continuous(breaks = c(-50,-25,0,25,50,75)) +
  labs(y = '', x = 'Percent change in burn severity') +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = 'black'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(rect = element_rect(fill = "transparent"),
        panel.background = element_blank(),
        panel.grid = element_blank())

ggsave('cbi_density.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs',
       height = 3, width = 6, dpi = 600)

# ------------------------------------------------------------------------------
# FRP change 
# ------------------------------------------------------------------------------
ggplot(pyrome_emd_df, aes(x = frp_dir, y = one_group, fill = stat(x))) +
  geom_density_ridges_gradient(bandwidth = 5) +
  #geom_vline(aes(xintercept = 0), linetype = 'dashed') +
  #geom_vline(aes(xintercept = -0.3), linetype = 'dashed') +
  #geom_vline(aes(xintercept = 0.3), linetype = 'dashed') +
  scale_fill_distiller(palette = 'RdYlBu', direction = 1, 
                       limits = c(-100,100), oob = scales::squish, guide = 'none') +
  coord_cartesian(xlim = c(-200, 200)) +
  #scale_x_continuous(breaks = seq(-300,300,100)) +
  labs(y = '', x = 'Change in FRP (yrs)') +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = 'black'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(rect = element_rect(fill = "transparent"),
        panel.background = element_blank(),
        panel.grid = element_blank())

ggsave('frp_density.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs',
       height = 3, width = 6, dpi = 600)

# ------------------------------------------------------------------------------
# NDVI change 
# ------------------------------------------------------------------------------
ggplot(pyrome_emd_df, aes(x = ndvi_dir, y = one_group, fill = stat(x))) +
  geom_density_ridges_gradient(bandwidth = 0.008) +
  #geom_vline(aes(xintercept = 0), linetype = 'dashed') +
  #geom_vline(aes(xintercept = -0.3), linetype = 'dashed') +
  #geom_vline(aes(xintercept = 0.3), linetype = 'dashed') +
  scale_fill_gradient2(low = '#a6611a', mid = '#ebf9b3', high = '#046235',
                       limits = c(-0.10,0.10), oob = scales::squish,
                       guide = 'none') +
  #coord_cartesian(xlim = c(-0.3, 0.3)) +
  scale_x_continuous(breaks = seq(-3,3,1)/10) +
  labs(y = '', x = 'Change in productivity (NDVI)') +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = 'black'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(rect = element_rect(fill = "transparent"),
        panel.background = element_blank(),
        panel.grid = element_blank())

 ggsave('ndvi_density.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs',
       height = 3, width = 6, dpi = 600)

# ------------------------------------------------------------------------------
# Various relationships
# ------------------------------------------------------------------------------
plot_df <- pyrome_emd_df %>% 
   group_by(cell_hist) %>% 
   slice_head(n = 1)
 
 ggplot(plot_df, aes(y = cbi_dir, x = frp_dir)) +
   geom_point(alpha = 0.8, aes(color = multi_emd_50)) +
   geom_smooth(span = 0.7, color = 'black', se = F, size = 1.5) +
   coord_cartesian(ylim = c(-.5,.5), xlim = c(-200,200)) +
   scale_color_distiller(palette = 'Spectral', limits = c(0,1.25), oob = scales::squish,
                         guide = 'none') +
   theme_bw()
   
 
  

# ------------------------------------------------------------------------------
# PCA1 change 
# ------------------------------------------------------------------------------
ggplot(pyrome_emd_df, aes(x = pca2_dir, y = one_group, fill = stat(x))) +
  geom_density_ridges_gradient(bandwidth = 0.1) +
  geom_vline(aes(xintercept = 0), linetype = 'dashed') +
  scale_fill_distiller(palette = 'BrBG', limits = c(-1.5,1.5), oob = scales::squish,
                       guide = 'none') +
  #coord_cartesian(xlim = c(-0.6, 0.6)) +
  #scale_x_continuous(breaks = c(-0.4,0,0.4)) +
  labs(y = '', x = 'Change in burn severity') +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = 'black'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(rect = element_rect(fill = "transparent"),
        panel.background = element_blank(),
        panel.grid = element_blank())


ggsave('cbi_delta_ridges.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs/',
       height = 7, width = 5, dpi = 600)

ggsave('cbi_density.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs',
       height = 3, width = 6, dpi = 600)



# ------------------------------------------------------------------------------
# Severity change - from burned points
# ------------------------------------------------------------------------------
change_df <- pyrome_df %>% 
  group_by(cell_hist) %>% 
  mutate(across(c(ndvi, log_fri, frs, cbi), mean, 
                .names = "ref_mean_{.col}")) %>% 
  group_by(cell_2C) %>% 
  mutate(across(c(ndvi, log_fri, frs, cbi), mean, 
                .names = "fut_mean_{.col}")) %>% 
  ungroup() %>% 
  mutate(delta_ndvi = fut_mean_ndvi - ref_mean_ndvi,
         delta_log_fri = fut_mean_log_fri - ref_mean_log_fri,
         delta_frs = fut_mean_frs - ref_mean_frs,
         delta_cbi = fut_mean_cbi - ref_mean_cbi) %>% 
  full_join(., ecoregs) 


plot_df <- change_df %>% 
  # filter(ECO_NAME %in% c('Southern Rocky Mountains',
  #                        'Utah-Wyoming Rocky Mountains',
  #                        'Klamath Mountains',
  #                        'East Cascades - Modoc Plateau',
  #                        'Middle Rockies - Blue Mountains',
  #                        'Canadian Rocky Mountains',
  #                        'Okanagan','Utah High Plateaus')) %>% 
  slice_sample(n = 1000000) %>% 
  mutate(ECO_NAME = as.factor(ECO_NAME),
         ECO_NAME = fct_reorder(ECO_NAME, delta_cbi, median, .desc = T)) 

ggplot(plot_df, aes(x = delta_cbi, y = ECO_NAME, fill = stat(x))) +
  geom_density_ridges_gradient(bandwidth = 0.02) +
  geom_vline(aes(xintercept = 0), linetype = 'dashed') +
  scale_fill_distiller(palette = 'RdYlGn', limits = c(-0.3,0.3), oob = scales::squish,
                       guide = 'none') +
  coord_cartesian(xlim = c(-0.6, 0.6)) +
  scale_x_continuous(breaks = c(-0.3,0,0.3)) +
  labs(y = '', x = 'Change in severity (CBI)') +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = 'black')) 


ggsave('cbi_delta_ridges.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
       height = 7, width = 7, dpi = 600)

# ------------------------------------------------------------------------------
# EMD Cases
# ------------------------------------------------------------------------------
mast_rast <- raster('data/process/mast_rast.tif')

# Join to dataframe of reference and future climate bins for all forested cells
forest_pyrome_df <- readRDS('data/process/forest_pyrome_df.Rdata')
pyrome_emd_df <- full_join(forest_pyrome_df, emd_chunk_df) %>% 
  filter(!is.na(case)) %>% 
  mutate(loss_case = case_when(case %in% c(1,2) & def_2C > 500 ~ 1))

emd_chunk_df <- emd_chunk_df %>%
  mutate(loss_case = case_when(case == 2 & def_2C > 500 ~ 1), # Was flammable forest, now fuel-limited (ie savanna)
         loss_case = case_when(case == 2 & def_2C < 500 ~ 2), # Was flammable forest, now climate-limited (ie PNW)
         loss_case = case_when(case == 3 & def_2C < 500 ~ 3), # Was fuel-limited, still is
         loss_case = case_when(case == 3 & def_2C > 500 ~ 4)) # Was climate-limited, still is



# Convert to spatial object (sf)
pyrome_emd_sf <- sf::st_as_sf(pyrome_emd_df, coords = c('x', 'y'), crs = 4326) 

r <- terra::rasterize(pyrome_emd_sf, mast_rast, 
                      field = 'loss_case', fun = modal, na.rm = T)
writeRaster(r, filename = paste0('data/emd/forest_','loss_case','.tiff'), overwrite = T)


 
