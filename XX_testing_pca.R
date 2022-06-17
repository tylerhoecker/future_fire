library(tidyverse)
library(factoextra)
library(sf)
library(tmap)

select <- dplyr::select

# PCA --------------------------------------------------------------------------
pyrome_df <- readRDS('data/process/pyrome_df.Rdata') %>% 
  filter(log_fri != '-Inf')

plot_df <- pyrome_df %>% 
  group_by(region) %>% 
  slice_sample(n = 1000) %>% 
  ungroup()

ggplot(plot_df, aes(x = cbi, y = log_fri)) +
  geom_point() +
  geom_smooth(method = 'lm')

pca_df <- pyrome_df %>% 
  select(cbi, ndvi, frs, log_fri)

pca_result <- prcomp(pca_df, scale = TRUE)

# Visualize results
fviz_eig(pca_result)

fviz_pca_var(pca_result,
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T)

# Mean of each dimension in original units for reference
pca_df %>% summarize(across(everything(), mean))

pca_ind <- get_pca_ind(pca_result)

n_fires <- pyrome_df %>% 
  select(cell_hist, fires)

pyrome_pca_df <- pyrome_df %>% 
  select(cell_hist, cell_2C, cbi, ndvi, log_fri, frs, def_hist, aet_hist) %>% 
  cbind(., pca_ind$coord[,c('Dim.1','Dim.2')], 'fires' = pyrome_df$fires)

# Pick a random subset of pyroclimates for plotting (there are too many to visualize clearly)
random_big_pyros <- pyrome_pca_df %>% 
  group_by(cell_hist) %>% 
  mutate(pyro_size = n()) %>% 
  filter(pyro_size > 1000) %>% 
  summarize() %>% 
  slice_sample(n = 100)


plot_pyros_pca <- pyrome_pca_df %>% 
  filter(cell_hist %in% random_big_pyros[['cell_hist']]) %>% 
  group_by(cell_hist) %>% 
  # slice_sample(n = 200) %>% 
  summarize(across(c(Dim.1,Dim.2), list('mean' = mean, 'sd' = sd)),
            fires = first(fires),
            aet_avg = mean(aet_hist),
            def_avg = mean(def_hist)) %>%
  mutate(Dim.1_lo = Dim.1_mean - (1.96*Dim.1_sd/fires),
         Dim.1_hi = Dim.1_mean + (1.96*Dim.1_sd/fires),
         Dim.2_lo = Dim.2_mean - (1.96*Dim.2_sd/fires),
         Dim.2_hi = Dim.2_mean + (1.96*Dim.2_sd/fires),
         cell_hist = as.factor(cell_hist))

# Extract loadings for plot
loadings <- data.frame(dimension = rownames(pca_result$rotation), pca_result$rotation)

# COLOR BY DEFICIT

ggplot(plot_pyros_pca, aes(x = Dim.1_mean, y = Dim.2_mean)) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = (PC1*2), yend = (PC2*2)), 
               arrow = arrow(length = unit(1/2, "picas")),
               #color = c('#e41a1c','#4daf4a','#984ea3','#a65628'), 
               size = 1) +
  geom_errorbar(aes(ymin = Dim.2_lo, ymax = Dim.2_hi), width = 0, size = 0.75, alpha = 0.6) +
  geom_errorbar(aes(xmin = Dim.1_lo, xmax = Dim.1_hi), width = 0, size = 0.75, alpha = 0.6) +
  geom_point(aes(fill = def_avg), color = 'black', shape = 21, size = 2) +
  
  geom_label(data = loadings, aes(x = (PC1*2), y = (PC2*2)+0.12), 
           label = c('Severity', 'Productivity','Resistance','FRI'),
           fill = c('grey'), 
           alpha = 0.8) +
  scale_fill_gradientn(colors = c('#01665e','#f6e8c3','#8c510a')) +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-1.5,1.5)) +
  theme_bw()


pyro_pca_means <- pyrome_pca_df %>% 
  filter(cell_hist %in% random_big_pyros[['cell_hist']]) %>% 
  group_by(cell_hist) %>% 
  summarize(pca1_mean = mean(Dim.1),
            pca2_mean = mean(Dim.2)) %>% 
  mutate(cell_hist = as.factor(cell_hist))

ggplot(plot_pyros_pca, aes(x = Dim.1, y = Dim.2, fill = cell_hist)) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_point(shape = 21, alpha = 0.7, color = 'black', size = 2) +
  geom_point(data = pyro_pca_means, 
             aes(x = pca1_mean, y = pca2_mean, fill = cell_hist), 
             shape = 21, color = 'black', size = 5) +
  scale_fill_brewer(palette = 'Set1') +
  theme_bw()



ref_pca_df <- pyrome_df %>% 
  group_by(cell_hist) %>% 
  summarise(ref_sev = mean(Dim.1),
            ref_freq = mean(Dim.2))

fut_pca_df <- pyrome_df %>% 
  group_by(cell_2C) %>% 
  summarise(fut_sev = mean(Dim.1),
            fut_freq = mean(Dim.2))

cell_pairs <- forest_pyrome_df %>% 
  group_by(cell_hist, cell_2C) %>% 
  tally() %>% 
  filter(n > 1000)

regime_change_df <- full_join(cell_pairs, ref_pca_df)

regime_change_df <- full_join(regime_change_df, fut_pca_df) 

regime_change_df <- regime_change_df %>% 
  filter(!is.na(ref_sev), !is.na(ref_freq), !is.na(fut_sev), !is.na(fut_freq))

ggplot(regime_change_df) +
  geom_segment(aes(x = ref_sev, xend = fut_sev, y = ref_freq, yend = fut_freq), 
               arrow = arrow(length = unit(0.2,"cm")),
               alpha = 0.5) +
  xlim(c(-1,1)) + ylim(c(-1,1))

pyrome_pca_df <- pyrome_pca_df %>% 
  mutate(regime = case_when(Dim.1 < 0 & Dim.2 < 0 ~ 1, #'Severe-Inrequent',
                            Dim.1 < 0 & Dim.2 > 0 ~ 2, #'Severe-Frequent',
                            Dim.1 > 0 & Dim.2 > 0 ~ 3, #'Insevere-Frequent',
                            Dim.1 > 0 & Dim.2 < 0 ~ 4)) #'Insevere-Infrequent'))

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

regimes <- pyrome_pca_df %>% 
  group_by(regime) %>% 
  summarise(across(c(cbi,ndvi,fri,frs), mean))

hist(regimes$mode_regime)


# WRITE OUT --------------------------------------------------------------------
forest_pyrome_df <- readRDS('data/process/forest_pyrome_df.Rdata')

forest_pca_df <- full_join(forest_pyrome_df, regimes)

forest_pca_sf <- sf::st_as_sf(forest_pca_df, coords = c('x', 'y'), crs = 4326) 

mast_rast <- raster('data/process/mast_rast.tif')

r <- terra::rasterize(forest_pca_sf, mast_rast, field = 'mode_regime', fun = modal, na.rm = T)

writeRaster(r, filename = paste0('data/emd/regime_pca','.tiff'), overwrite = T)




# ------------------------------ Experimenting
test_df <- as_tibble(pyrome_df) %>% 
  select(cell_hist, cell_2C, cbi, ndvi, fri, frs) 

# Euclidean distance
529 
562

ref_test <- test_df %>% 
  filter(cell_hist == 529) %>% 
  select(cbi, ndvi, fri, frs) %>% 
  slice_head(n = 1)

fut_test <- test_df %>% 
  filter(cell_hist == 529) %>% 
  summarise(across(everything(), mean)) %>% 
  select(cbi, ndvi, fri, frs)

dist(rbind(ref_test, fut_test))

# MAP ---------------------------------------------------------------------------
boundary <- read_sf("C:\\Users\\hoecker\\Work\\GIS\\cb_2018_us_state_20m\\cb_2018_us_state_20m.shp") %>%
  st_transform(., crs = 4326)

plot_df <- pyrome_df %>% 
  slice_sample(n = 2000)

plot_sf <- st_as_sf(plot_df, coords = c('x', 'y'), crs = 4326)

tm_shape(boundary, bbox = pyrome_sf) +
  tm_borders() +
  tm_fill(col = 'grey10') +
  tm_shape(plot_sf) +
  tm_dots(title = 'Fire regime',
          col = 'regime',
          size = 0.03,
          #palette = c('#8c510a','#f6e8c3','#01665e'),
          #style = 'cont',
          shape = 15) +
  tm_layout(legend.outside = TRUE)

