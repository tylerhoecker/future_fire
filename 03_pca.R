library(tidyverse)
library(factoextra)

select <- dplyr::select

# PCA --------------------------------------------------------------------------
pyrome_df <- readRDS('data/process/pyrome_df.Rdata') %>% 
  filter(log_fri != '-Inf')

pca_df <- pyrome_df %>% 
  select(cbi, ndvi, frs, log_fri)

pca_result <- prcomp(pca_df, scale = TRUE)

# Visualize results
fviz_eig(pca_result)

fviz_pca_var(pca_result,
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T)

pca_ind <- get_pca_ind(pca_result)

pyrome_pca_df <- pyrome_df %>%
  mutate(pca1 = pca_ind$coord[,c('Dim.1')],
         pca2 = pca_ind$coord[,c('Dim.2')],
         'fires' = pyrome_df$fires)

# Save this as pyrome_df for use in the adaptation script
saveRDS(pyrome_pca_df, 'data/process/pyrome_df.Rdata')

# Summarize PCA dimensions and climate for reference pyroclimates
ref_pyros_pca <- pyrome_pca_df %>% 
  group_by(cell_hist) %>% 
  summarize(ref_pca1_mean = mean(Dim.1),
            ref_pca2_mean = mean(Dim.2),
            ref_pca1_sd = sd(Dim.1),
            ref_pca2_sd = sd(Dim.2),
            ref_aet_mean = mean(aet_hist),
            ref_def_mean = mean(def_hist),
            ref_aet_sd = sd(aet_hist),
            ref_def_sd = sd(def_hist),
            fires = first(fires)) %>%
  mutate(ref_pca1_lo = ref_pca1_mean - (1.96*ref_pca1_sd/fires),
         ref_pca1_hi = ref_pca1_mean + (1.96*ref_pca1_sd/fires),
         ref_pca2_lo = ref_pca2_mean - (1.96*ref_pca2_sd/fires),
         ref_pca2_hi = ref_pca2_mean + (1.96*ref_pca2_sd/fires))

# Summarize PCA dimensions and climate for future pyroclimates
fut_pyros_pca <- pyrome_pca_df %>% 
  group_by(cell_2C) %>% 
  summarize(fut_pca1_mean = mean(Dim.1),
            fut_pca2_mean = mean(Dim.2),
            fut_pca1_sd = sd(Dim.1),
            fut_pca2_sd = sd(Dim.2),
            fut_aet_mean = mean(aet_2C),
            fut_def_mean = mean(def_2C),
            fut_aet_sd = sd(aet_2C),
            fut_def_sd = sd(def_2C),
            fires = first(fires)) %>%
  mutate(fut_pca1_lo = fut_pca1_mean - (1.96*fut_pca1_sd/fires),
         fut_pca1_hi = fut_pca1_mean + (1.96*fut_pca1_sd/fires),
         fut_pca2_lo = fut_pca2_mean - (1.96*fut_pca2_sd/fires),
         fut_pca2_hi = fut_pca2_mean + (1.96*fut_pca2_sd/fires))

# Pick a random subset of pyroclimates for plotting (there are too many to visualize clearly)
random_big_pyros <- pyrome_pca_df %>% 
  group_by(cell_hist) %>% 
  mutate(pyro_size = n()) %>% 
  filter(pyro_size > 1000) %>% 
  summarize() %>% 
  slice_sample(n = 100)


# Extract loadings for plot
loadings <- data.frame(dimension = rownames(pca_result$rotation), pca_result$rotation)

# Refernce period pyromes colored by deficit
plot_pyros_pca <- pyrome_pca_df %>%
  select(cell_hist, cell_2C) %>%
  distinct() %>%
  full_join(., ref_pyros_pca) %>%
  full_join(., fut_pyros_pca) %>% 
  filter(cell_hist %in% random_big_pyros[['cell_hist']]) %>% 
  mutate(cell_hist = as.factor(cell_hist),
         cell_2C = as.factor(cell_2C))

ggplot(plot_pyros_pca, aes(x = ref_pca1_mean, y = ref_pca2_mean)) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = (PC1), yend = (PC2)), 
               arrow = arrow(length = unit(1/2, "picas")),
               size = 1) +
  geom_errorbar(aes(ymin = ref_pca2_lo, ymax = ref_pca2_hi), width = 0, size = 0.5, alpha = 0.7) +
  geom_errorbar(aes(xmin = ref_pca1_lo, xmax = ref_pca1_hi), width = 0, size = 0.5, alpha = 0.7) +
  geom_point(aes(fill = ref_def_mean), color = 'black', shape = 21, size = 2) +
  geom_label_repel(data = loadings, aes(x = PC1, y = PC2), 
           label = c('Fire severity', 'Vegetation\nproductivity','Fire resistance','FRI'),
           fill = c('grey'), alpha = 0.8, min.segment.length = 1) +
  scale_fill_gradientn('Climatic\nwater\ndeficit', colors = c('#01665e','#f6e8c3','#8c510a')) +
  coord_cartesian(xlim = c(-1.25,1.25), ylim = c(-1.25,1.25)) +
  theme_bw() +
  theme(axis.text = element_text(color = 'black'),
        legend.position = c(0.9,0.2),
        legend.background = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)) +
  labs(x = 'Principle component 1 (31%)', y = 'Principle component 2 (26%)')

ggsave('pca_vectors.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire',
       height = 5, width = 5, dpi = 500)

# Refernce and future pyromes connected as vectors
ggplot(plot_pyros_pca) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_point(aes(x = ref_pca1_mean, y = ref_pca2_mean, fill = ref_def_mean), 
             color = 'black', shape = 21, size = 2) +
  # geom_point(aes(x = fut_pca1_mean, y = fut_pca2_mean, fill = fut_def_mean), 
  #            color = 'black', shape = 21, size = 2) +
  geom_segment(aes(x = ref_pca1_mean, y = ref_pca2_mean, 
                   xend = fut_pca1_mean, yend = fut_pca2_mean),
               arrow = arrow(length = unit(1, "picas")),
               size = 0.75) +
  geom_label(data = loadings, aes(x = (PC1*2), y = (PC2*2)+0.12), 
             label = c('Fire severity', 'Vegetation\nproductivity','Fire resistance','FRI'),
             fill = c('grey'), 
             alpha = 0.8) +
  scale_fill_gradientn('Climatic\nwater\ndeficit', colors = c('#01665e','#f6e8c3','#8c510a')) +
  #coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-1.5,1.5)) +
  coord_cartesian(xlim = c(-2,2), ylim = c(-2,2)) +
  
  theme_bw() +
  theme(axis.text = element_text(color = 'black'),
        legend.position = c(0.9,0.2),
        legend.background = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)) +
  labs(x = 'Principle component 1 (31%)', y = 'Principle component 2 (26%)')



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

