library(tidyverse)

# Abatzoglou data future projections: https://datadryad.org/stash/landing/show?id=doi%3A10.6071%2FM3WQ1R
# For unknown reason these do not align with the time series in the published paper. 
# Values here DO seem to align with the values in the appendix for the CWD model, but
# not the values shown in the body of the paper for the VPD model...

abatzoglou_strongfading <- 
  list.files('data/abatzoglou/', full.names = T, pattern = '.*_strong-fading.csv') %>% 
    map_df(., read_csv, .id = 'gcm') 

abatzoglou_future_yr <- abatzoglou_strongfading %>% 
  pivot_longer(cols = ens1:ens1000, names_to = 'ens', values_to = 'area') %>% 
  group_by(gcm, year) %>% 
  summarise(area = rollmean(area, k = 11, align = c("right"))) 
  

ggplot(abatzoglou_future_yr, aes(x = year, y = area)) +
  geom_line() +
  geom_smooth(span = 0.5) +
  theme_bw() +
  coord_cartesian(ylim = c(0,10000))

# ---
abatzoglou_static <- 
  list.files('data/abatzoglou/', full.names = T, pattern = '.*_static.csv') %>% 
  map_df(., read_csv, .id = 'gcm') 

abatzoglou_future_yr <- abatzoglou_static %>% 
  pivot_longer(cols = ens1:ens1000, names_to = 'ens', values_to = 'area') %>% 
  group_by(year) %>% 
  summarise(area = median(area)) 

ggplot(abatzoglou_future_yr, aes(x = year, y = area)) +
  geom_line() +
  geom_smooth(span = 0.5) +
  theme_bw() +
  coord_cartesian(ylim = c(0,10000))


# Analog model -----------------------------------------------------------------
forest_clim_df <- readRDS('data/process/forest_clim_df.Rdata') %>% 
  as_tibble()









sum(forest_clim_df$burned * forest_clim_df$area) * 0.000001


