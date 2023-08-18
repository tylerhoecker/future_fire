# ------------------------------------------------------------------------------
# Description: Explore, summarize, and plot the fire-regime change exposure 
# results and the component dimensions.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load packages
# ------------------------------------------------------------------------------
library(tidyverse)
library(ggridges)
library(sf)
library(terra)
library(raster)
select <- dplyr::select

# ------------------------------------------------------------------------------
# EMD
# ------------------------------------------------------------------------------
emd_chunk_df <- readRDS('data/emd/emd_2023_06_13.Rdata')

# Join to dataframe of reference and future climate bins for all forested cells
forest_pyrome_df <- readRDS('data/process/forest_pyrome_df.Rdata')

# Just examine pyromes where exposure can be estimated (filter out shifts into or out of forest space)
pyrome_emd_df <- full_join(forest_pyrome_df, emd_chunk_df) %>%
  select(-null_sum)
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
median(pyrome_emd_df$emd_mean[pyrome_emd_df$emd_mean > 0], na.rm = T)
summary(pyrome_emd_df$frp_dir, na.rm = T)

# quant075 <- function(x) quantile(x, probs = 0.75, na.rm = T)

expanse(rast('data/process/mast_rast.tif'))

areas_pts <- extract(area(mast_rast), pyrome_emd_sf)

pyrome_emd_sf |> 
  st_drop_geometry() |> 
  select(ID, emd_mean, null_05, null_95, emd_05, emd_95) |> 
  mutate(area = areas_pts) |>
  rowwise() |> 
  mutate(null_overlap = ifelse(emd_mean  > 0 && null_05 <= emd_95 && emd_05 >= null_95, 1, 0),
         null_overlap = ifelse(emd_mean  == 0, 2, null_overlap)) |> 
  ungroup() |> 
  distinct(ID, .keep_all = T) |> 
  group_by(null_overlap) |> 
  summarize(n = n(),
            area = sum(area),
            pct = area/766934)

emd_chunk_df |> 
  select(emd_mean, null_05, null_95, emd_05, emd_95) |> 
  rowwise() |> 
  mutate(null_overlap = ifelse(emd_mean  > 0 && null_05 <= emd_95 && emd_05 >= null_95, 1, 0),
         null_overlap = ifelse(emd_mean  == 0, 2, null_overlap)) |> 
  group_by(null_overlap) |> 
  summarize(n = n(),
            pct = n()/3999)

summary(emd_chunk_df$emd_mean)

# Sort regions by EMD
pyrome_emd_df <- pyrome_emd_df %>% 
  mutate(#ECO_NAME = as.factor(ECO_NAME),
         #ECO_NAME = fct_reorder(ECO_NAME, multi_emd_50, quant075, .desc = T),
         one_group = 1) 

ggplot(pyrome_emd_df, aes(x = emd_mean, y = one_group, fill = stat(x))) +
  geom_density_ridges_gradient(bandwidth = 1) +
  #geom_vline(aes(xintercept = mean(multi_emd_50, na.rm = T)), linetype = 'dashed') +
  scale_fill_viridis_c(option = "B", guide = 'none', limits = c(0,16), oob = scales::squish) +
  # scale_fill_distiller(palette = 'Spectral', limits = c(0,25), oob = scales::squish,
  #                      guide = 'none') +
  coord_cartesian(xlim = c(-1,30)) +
  labs(x = 'Exposure')  +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(color = 'black'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        rect = element_rect(fill = "transparent"),
        panel.background = element_blank(),
        panel.border = element_blank()) 
# +
#   theme(,
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_blank())

ggsave('emd_density.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs',
       height = 1.3, width = 2, dpi = 600, type = "cairo-png")


# ------------------------------------------------------------------------------
# Significance null test
# ------------------------------------------------------------------------------
firesheds <- st_read('data/emd/bivariate_firesheds.shp') |> 
  st_drop_geometry() |> 
  mutate(name = fct_reorder(as.factor(Frshd_I), adpt_95)) |> 
  filter(adpt_ds > 0) |> 
  slice_sample(n = 100)
  

ggplot(firesheds, aes(y = name)) +
  geom_linerange(aes(xmin = adpt_05, xmax = adpt_95, color = 'Contemporary-Future')) +
  geom_linerange(aes(xmin = null_05, xmax = null_95, color = 'Contemporary-Contemporary \n(null model)')) +
  geom_point(aes(x = nll_dst, color ='Contemporary-Contemporary \n(null model)')) +
  geom_point(aes(x = adpt_ds, color =  'Contemporary-Future')) +
  scale_color_manual(name = 'Fire-regime comparison', values = c('red','black'))+
  theme_bw(base_size = 8) +
  labs(x = "Dissimilarity (Earth mover's distance)", y = "Fireshed ID")

ggsave('SI_null_model_firesheds.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Main_Figs',
       height = 10, width = 7, dpi = 600)

# Severity vs exposure vs productivity at fireshed level
firesheds <- st_read('data/emd/bivariate_firesheds.shp') |> 
  st_drop_geometry() 

ggplot(firesheds, aes(x = frp_chn, y = cb_chng)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_bw()

cor(firesheds$cb_chng, firesheds$frp_chn)
cor(firesheds$cb_chng, firesheds$ndv_chn)


ggplot(emd_chunk_df, aes(x = ndvi_dir, y = cbi_dir)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_bw()

post_mod <- lm(cbi_dir ~ ndvi_dir*frp_dir, data = pyrome_emd_df) #pyrome_emd_df
summary(post_mod)

cor(pyrome_emd_df$cbi_dir, pyrome_emd_df$ndvi_dir)
cor(pyrome_emd_df$cbi_dir, pyrome_emd_df$frp_dir)

# ------------------------------------------------------------------------------
# Severity change 
# ------------------------------------------------------------------------------
pyrome_emd_df %>% 
  filter(cbi_dir < -0.1) %>% 
  nrow()/nrow(pyrome_emd_df)

ggplot(pyrome_emd_df, aes(x = cbi_pct*100, y = one_group, fill = stat(x))) +
  geom_density_ridges_gradient(bandwidth = 2) +
  #geom_vline(aes(xintercept = 0), linetype = 'dashed') +
  #geom_vline(aes(xintercept = -0.3), linetype = 'dashed') +
  #geom_vline(aes(xintercept = 0.3), linetype = 'dashed') +
  scale_fill_distiller(palette = 'RdYlGn', limits = c(-20,20), oob = scales::squish,
                       guide = 'none') +
  #coord_cartesian(xlim = c(-0.6, 0.6)) +
  scale_x_continuous(breaks = c(-50,-25,0,25,50,75)) +
  labs(y = '', x = 'Percent change in burn severity') +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(color = 'black'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        rect = element_rect(fill = "transparent"),
        panel.background = element_blank(),
        panel.border = element_blank()) 

ggsave('cbi_density.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs',
       height = 1.3, width = 2, dpi = 500)

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
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(color = 'black'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        rect = element_rect(fill = "transparent"),
        panel.background = element_blank(),
        panel.border = element_blank()) 

ggsave('frp_density.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs',
       height = 1.3, width = 2, dpi = 500)

# ------------------------------------------------------------------------------
# NDVI change 
# ------------------------------------------------------------------------------
ggplot(pyrome_emd_df, aes(x = ndvi_pct*100, y = one_group, fill = stat(x))) +
  geom_density_ridges_gradient(bandwidth = 0.8) +
  #geom_vline(aes(xintercept = 0), linetype = 'dashed') +
  #geom_vline(aes(xintercept = -0.3), linetype = 'dashed') +
  #geom_vline(aes(xintercept = 0.3), linetype = 'dashed') +
  scale_fill_gradient2(low = '#a6611a', mid = '#ebf9b3', high = '#046235',
                       limits = c(-10,10), oob = scales::squish,
                       guide = 'none') +
  #coord_cartesian(xlim = c(-0.3, 0.3)) +
  scale_x_continuous(breaks = seq(-30,30,10)) +
  labs(y = '', x = 'Percent change in productivity') +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(color = 'black'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        rect = element_rect(fill = "transparent"),
        panel.background = element_blank(),
        panel.border = element_blank()) 

 ggsave('ndvi_density.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs',
        height = 1.3, width = 2, dpi = 500)

# ------------------------------------------------------------------------------
# Various relationships
# ------------------------------------------------------------------------------
plot_df <- pyrome_emd_df %>% 
   #group_by(cell_hist) %>% 
   slice_sample(n = 100000)
 
 ggplot(plot_df, aes(y = cbi_dir, x = frp_dir)) +
   geom_point(alpha = 0.8, aes(color = multi_emd_50)) +
   geom_smooth(method = 'lm', color = 'black', se = F, size = 1.5) +
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


 
