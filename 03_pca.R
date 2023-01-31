# Performs a PCA on pyroclimate dimensions to reduce from 3 to 2
# Creates plots
# Saves output for next step

library(tidyverse)
library(factoextra)
library(ggrepel)

select <- dplyr::select

# PCA --------------------------------------------------------------------------
pyrome_df <- readRDS('data/process/pyrome_df.Rdata') 

pca_df <- pyrome_df %>% 
  select(cbi, ndvi, log_frp_pt)

# Some stats which are useful for intepreting PCA which is rescaled
# c('mean' = mean(pca_df$cbi), '1 sd' = sd(pca_df$cbi))
# c('mean' = mean(pca_df$ndvi), '1 sd' = sd(pca_df$ndvi))
# c('mean' = exp(mean(pca_df$log_frp_cell)), '1 sd' = exp(sd(pca_df$log_frp_cell)))

pca_result <- prcomp(pca_df, scale = T) 

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

# Save this as pyrome_df for use in the exposure script
saveRDS(pyrome_pca_df, 'data/process/pyrome_df.Rdata')
pyrome_pca_df <- readRDS('data/process/pyrome_df.Rdata')
# Plotting ---------------------------------------------------------------------
# Summarize PCA dimensions and variables for reference pyroclimates
ref_pyros_pca <- pyrome_pca_df %>% 
  group_by(cell_hist) %>% 
  summarize(across(c(pca1, pca2, cbi, ndvi, frp_pt),
                   list(mean = mean, se = ~ sd(.x)/ sqrt(sum(fires)), 
                        low = ~ mean(.x) - sd(.x),
                        hi = ~ mean(.x) + sd(.x))))
# Extract loadings for plot
loadings <- data.frame(dimension = rownames(pca_result$rotation), pca_result$rotation)

# # Some stats which are useful for intepreting PCA which is rescaled
# c('mean' = mean(ref_pyros_pca$cbi_mean), '1 sd' = sd(ref_pyros_pca$cbi_mean))
# c('mean' = mean(ref_pyros_pca$ndvi_mean), '1 sd' = sd(ref_pyros_pca$ndvi_mean))
# c('mean' = mean(ref_pyros_pca$frp_cell_mean), '1 sd' = sd(ref_pyros_pca$frp_cell_mean))
# 
# 
# ref_to_clust <- ref_pyros_pca %>% 
#   select(pca1_mean, pca2_mean)
# 
# fviz_nbclust(ref_to_clust, FUNcluster=cluster::pam, k.max = 7)
# 
# pca_clust <- eclust(ref_to_clust, "kmeans", hc_metric="eucliden", k = 4)
# 
# ref_pyros_pca <- ref_pyros_pca %>% 
#   mutate(cluster = pca_clust$cluster)


# # Colored by cluster
# ggplot(ref_pyros_pca, aes(x = pca1_mean, y = pca2_mean)) +
#   geom_hline(aes(yintercept = 0)) +
#   geom_vline(aes(xintercept = 0)) +
#   geom_segment(data = loadings, 
#                aes(x = 0, y = 0, xend = (PC1), yend = (PC2)), 
#                arrow = arrow(length = unit(1/2, "picas")),
#                size = 1) +
#   # geom_errorbar(aes(xmin = pca1_low, xmax = pca1_hi), width = 0, size = 0.5, alpha = 0.7) +
#   # geom_errorbar(aes(ymin = pca2_low, ymax = pca2_hi), width = 0, size = 0.5, alpha = 0.7) +
#   geom_point(aes(fill = as.factor(cluster)), color = 'black', shape = 21, size = 2) +
#   geom_label_repel(data = loadings, aes(x = PC1, y = PC2),
#            label = c('Burn severity', 'Vegetation productivity','FRP'),
#            fill = c('grey'), alpha = 0.8, min.segment.length = 1) +
#   #scale_fill_gradientn('Climatic\nwater\ndeficit', colors = c('#01665e','#f6e8c3','#8c510a')) +
#   #coord_fixed(ratio = .94, xlim = c(-1.5,1.5), ylim = c(-1.5,1.5)) +
#   theme_bw(base_size = 14) +
#   theme(axis.text = element_text(color = 'black'),
#         legend.position = c(0.9,0.2),
#         legend.background = element_blank(),
#         legend.title = element_text(size = 9),
#         legend.text = element_text(size = 8)) +
#   labs(x = 'Principal component 1 (36.6%)', y = 'Principal component 2 (32.4%)')
# 
# ggsave('pca_vectors.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Main_Figs/',
#        height = 6, width = 6, dpi = 500)

# Colored by 3 main variables --------------------------------------------------
# FRP
frp_pca <- 
  ggplot(ref_pyros_pca, aes(x = pca1_mean, y = pca2_mean)) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  ggstar::geom_star(aes(fill = frp_pt_mean), color = 'black', starshape = "hexagon", size = 2) +
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = (PC1), yend = (PC2)), 
               arrow = arrow(length = unit(1/2, "picas")),
               size = 1) +
  geom_label_repel(data = loadings, aes(x = PC1, y = PC2),
                   label = c('Burn severity', 'Vegetation productivity','FRP'),
                   fill = c('grey'), alpha = 0.8, min.segment.length = 0.5, size = 2.5) +
  scale_fill_distiller('Fire rotation \nperiod (yrs)', palette = 'PuOr', direction =1,
                       limits = c(50,250), oob = scales::squish) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = c(0.92,0.2),
        legend.background = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)) +
  labs(x = 'Principal component 1 (36.1%)', y = 'Principal component 2 (33.7%)')
frp_pca

# CBI
cbi_pca <- 
  ggplot(ref_pyros_pca, aes(x = pca1_mean, y = pca2_mean)) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  ggstar::geom_star(aes(fill = cbi_mean), color = 'black', starshape = "hexagon", size = 2) +
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = (PC1), yend = (PC2)), 
               arrow = arrow(length = unit(1/2, "picas")),
               size = 1) +
  geom_label_repel(data = loadings, aes(x = PC1, y = PC2),
                   label = c('Burn severity', 'Vegetation productivity','FRP'),
                   fill = c('grey'), alpha = 0.8, min.segment.length = 0.5, size = 2.5) +
  scale_fill_distiller('Burn severity \n(CBI)', palette = 'RdYlGn',
                       limits = c(1.3,2), oob = scales::squish) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = c(0.92,0.2), #
        legend.background = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)) +
  labs(x = 'Principal component 1 (36.1%)', y = 'Principal component 2 (33.7%)')
# NDVI
ndvi_pca <- 
  ggplot(ref_pyros_pca, aes(x = pca1_mean, y = pca2_mean)) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  ggstar::geom_star(aes(fill = ndvi_mean), color = 'black', starshape = "hexagon", size = 2) +
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = (PC1), yend = (PC2)), 
               arrow = arrow(length = unit(1/2, "picas")),
               size = 1) +
  geom_label_repel(data = loadings, aes(x = PC1, y = PC2),
                   label = c('Burn severity', 'Vegetation productivity','FRP'),
                   fill = c('grey'), alpha = 0.8, min.segment.length = 0.5, size = 2.5) +
  scale_fill_distiller('Productivity \n(NDVI)', palette = 'BrBG', direction = 1,
                       limits = c(0.45,0.80), oob = scales::squish) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = c(0.92,0.2),
        legend.background = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)) +
  labs(x = 'Principal component 1 (36.1%)', y = 'Principal component 2 (33.7%)')

panel_hex <- cowplot::plot_grid(frp_pca, cbi_pca, ndvi_pca, nrow = 1)
panel_hex
ggsave('three_panel_pca.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs/',
       height = 5, width = 15, dpi = 500)


# Different version of similar for connected vectors
fut_pca <- pyrome_pca_df %>% 
  group_by(cell_2C) %>% 
  summarize(across(c(pca1, pca2), list(fut_mean = mean)))

ref_clusts <- ref_pyros_pca %>% 
  select(cell_hist, cluster)

ref_pca <- pyrome_pca_df %>% 
  group_by(cell_hist) %>% 
  summarize(across(c(pca1, pca2), list(ref_mean = mean))) %>% 
  left_join(., ref_clusts)

plot_pca <- pyrome_pca_df %>% 
  distinct(cell_hist, cell_2C) %>%
  left_join(., fut_pca) %>%
  left_join(., ref_pca) %>%
  filter(cell_hist %in% ref_pyros_pca$cell_hist) %>% 
  group_by(cluster) %>% 
  summarize(pca1_ref_regime = mean(pca1_ref_mean),
            pca2_ref_regime = mean(pca2_ref_mean),
            pca1_fut_regime = mean(pca1_fut_mean),
            pca2_fut_regime = mean(pca2_fut_mean))

# Refernce and future pyromes connected as vectors
ggplot(plot_pca) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = (PC1), yend = (PC2)), 
               arrow = arrow(length = unit(1/2, "picas")),
               size = 1) +
  geom_point(data = ref_pyros_pca, aes(x = pca1_mean, y = pca2_mean, fill = as.factor(cluster)),
             color = 'black', shape = 21, size = 2) +
  geom_point(aes(x = pca1_ref_regime, y = pca2_ref_regime), 
             color = 'black', fill = 'grey50', shape = 21, size = 4) +
  geom_segment(aes(x = pca1_ref_regime, y = pca2_ref_regime, 
                   xend = pca1_fut_regime, yend = pca2_fut_regime),
               arrow = arrow(length = unit(1, "picas")),
               size = 0.75) +
  geom_label_repel(data = loadings, aes(x = PC1, y = PC2),
                   label = c('Burn severity', 'Vegetation productivity','FRP'),
                   fill = c('grey'), alpha = 0.8, min.segment.length = 1) +
  coord_cartesian(xlim = c(-3,3), ylim = c(-3,3)) +
  theme_bw() +
  theme(axis.text = element_text(color = 'black'),
        legend.position = c(0.9,0.2),
        legend.background = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)) +
  labs(x = 'Principle component 1 (36%)', y = 'Principle component 2 (34%)')







