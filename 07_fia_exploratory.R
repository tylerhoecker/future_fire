library(tidyverse)
library(sf)
library(terra)
#library(raster)
library(ggrepel)
select <- dplyr::select

# ------------------------------------------------------------------------------
# Fireshed-based method (to ensure best agreement with map)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Wrangle FIA plot data to determine dominant species 
# ------------------------------------------------------------------------------
# Read in a compiled table for CONUS of FIA data as of 2019: 
#      https://www.fs.usda.gov/rds/archive/Catalog/RDS-2019-0026
fia_tree <- read_csv('data/composition/Tree_table_CONUS.csv',
                     col_names = c('tl_id','CN','INVYR','STATECD',
                                   'State_Abbreviation','UNITCD','COUNTYCD',
                                   'PLOT','SUBP','TREE','STATUSCD','SPCD','DIA',
                                   'HT','ACTUALHT','CR','TPA_UNADJ','TreeVisitID',
                                   'TreeID')) %>% 
  # Reduce down to states in our study area
  filter(State_Abbreviation %in% c('WA','OR','ID','MT','CA','NV','UT','WY','CO','AZ','NM'))

# Link to species and forest type descriptions/names
spp_ref <- read_csv('data/composition/REF_SPECIES.csv')
# Make a simplified df with code and common name, species, genus
spp_names <- spp_ref %>% 
  dplyr::select(sp = SPCD, tree_group = W_SPGRPCD, name = COMMON_NAME) %>% #, name = COMMON_NAME, GENUS, SPECIES
  mutate(#latin_name = paste(GENUS, SPECIES),
    tree_group = factor(tree_group),
    tree_group = fct_recode(tree_group,
                            'Douglas-fir'='10',   
                            'Ponderosa pine'='11',
                            'True firs'='12', 
                            'Western hemlock'='13',  
                            'Sugar pine'='14',  
                            'Western white pine'='15', 
                            'Redwood'='16',   
                            'Sitka spruce'='17',  
                            'Engelmann spruce'='18',
                            'Western larch'='19',
                            'Incense-cedar'='20',
                            'Lodgepole pine'='21',  
                            'Western redcedar'='22' ,  
                            'Junipers and Pinyon pines'='23' ,
                            'Other western softwoods'='24',
                            'Aspen'='44',
                            'Alder'='45',
                            'Other west hardwoods'='46',
                            'Other west hardwoods'='47',
                            'Other west hardwoods'='48'))


# Identify the dominant species, by basal area, at each plot (CN)
fia_tree_basal <- fia_tree %>% 
  select(tl_id, fia_plot = CN, sp = SPCD, dbh = DIA) %>% 
  left_join(., spp_names) %>%
  mutate(basal = pi*((dbh/2)^2)) %>% 
  group_by(fia_plot) %>% 
  mutate(plot_count = n()) %>% 
  group_by(fia_plot,tree_group) %>%
  summarize(tree_basal = sum(basal, na.rm = T)) %>% 
  group_by(fia_plot) %>% 
  mutate(plot_basal = sum(tree_basal, na.rm = T),
         tree_prop = tree_basal/plot_basal) %>% 
  slice_max(n = 1, tree_prop) 

rm(fia_tree)

# ------------------------------------------------------------------------------
# Intersect firesheds with FIA plot data to determine dominant species 
# ------------------------------------------------------------------------------
# Read in the shapefile of firesheds with mean adaptation distance and mean fire 
# resistance score already calculated in previous steps (01-06)
firesheds <- read_sf('data/emd/bivariate_firesheds.shp')

firesheds %>% 
  st_drop_geometry() %>% 
  group_by(bicat) %>%
  summarize(n = n(),
            area = sum(Area_HA)) %>% 
  mutate(pct = round(n/sum(n)*100,1),
         pct_area = round(area/sum(area)*100,1))
            

# Load in TreeMap raster - this will link FIA plot data to space
treemap <- rast('data/composition/national_c2014_tree_list.tif')
# Load the key that connects TreeMap values to FIA plot IDs ('CN')
treemap_key <- read_csv('data/composition/TL_CN_Lookup.txt') %>% 
  select(fia_plot = CN, tl_id)



# Do this process again but summarize by tree group, to see if that changes anything
fireshed_treegroup_df <- split(vect(firesheds), 'Fireshed_I') %>% 
  map_dfr(function(fireshed){
    
    print(paste('Running fireshed', fireshed$Fireshed_I, fireshed$Fireshed_N))
    
    # Extract TreeMap values within fireshed polygon
    treemap_fs <- terra::crop(treemap, fireshed)
    tree_sheds <- terra::extract(treemap_fs, fireshed)
    tree_sheds$ID <- fireshed$Fireshed_I
    
    # Now summarize FIA data associated with this fireshed
    modal_dominant <- tree_sheds %>% 
      as_tibble() %>% 
      rename(tl_id = national_c2014_tree_list) %>% 
      filter(!is.na(tl_id)) %>% 
      left_join(., treemap_key) %>% 
      left_join(.,fia_tree_basal) %>% 
      group_by(ID) %>% 
      distinct() %>%
      mutate(fireshed_count = n()) %>% 
      group_by(ID, tree_group) %>% 
      summarise(tree_count_prop = n()/fireshed_count) %>% 
      ungroup() %>% 
      slice_max(n = 1, tree_count_prop, with_ties = F) %>% 
      rename(modal_dominant = tree_group)
    
    majority_basal <- tree_sheds %>% 
      as_tibble() %>% 
      rename(tl_id = national_c2014_tree_list) %>% 
      filter(!is.na(tl_id)) %>% 
      left_join(., treemap_key) %>% 
      left_join(.,fia_tree_basal) %>% 
      group_by(ID) %>% 
      distinct() %>% 
      mutate(shed_basal = sum(plot_basal, na.rm = T)) %>% 
      group_by(ID, tree_group) %>% 
      summarize(tree_basal_prop = sum(tree_basal, na.rm = T)/max(shed_basal)) %>% 
      ungroup() %>% 
      slice_max(n = 1, tree_basal_prop, with_ties = F) %>% 
      rename(majority_basal = tree_group)
    
    result <- full_join(modal_dominant, majority_basal)
    
    return(result)
  })


fireshed_treegroup_df 

fireshed_df <- firesheds %>% 
  st_drop_geometry() %>% 
  full_join(., fireshed_treegroup_df, by = c('Fireshed_I'='ID')) %>% 
  rename(tree_group = majority_basal) %>% 
  filter(tree_group != 'Aspen') %>% 
  mutate(tree_group = fct_relevel(as.factor(tree_group), 
                                  levels = c('Douglas-fir',
                                             'Ponderosa pine',
                                             'Junipers and Pinyon pines',
                                             'Lodgepole pine',
                                             'True firs',
                                             'Engelmann spruce',
                                             'Redwood'))) 

# fireshed_df %>% 
#   group_by(bicat, tree_group) %>% 
#   summarise(bicat_mode = n()) %>% 
#   slice_max(n = 3, bicat_mode) %>% 
#   View()

forest_type_egs <- fireshed_df %>% 
  filter(Fireshed_I %in% c(633, 1326, 1983, 21, 822, 1915, 260, 164, 102, 462,
                           524, 759, 1228, 1057, 1186, 1763, 2140, 1622))

# Colored by bivariate group
ggplot(fireshed_df, aes(x = adapt_dist, y = frs)) +
  geom_vline(aes(xintercept = 15)) +
  geom_vline(aes(xintercept = 20)) +
  geom_hline(aes(yintercept = 0.38)) +
  geom_hline(aes(yintercept = 0.54)) +
  # geom_errorbar(aes(xmin = lo_adapt, xmax = hi_adapt), width = 0, size = 0.5, alpha = 0.7) +
  # geom_errorbar(aes(ymin = lo_frs, ymax = hi_frs), width = 0, size = 0.5, alpha = 0.7) +
  geom_point(aes(fill = bicat), color = 'black', shape = 21, size = 2) +
  geom_label_repel(data = forest_type_egs,
                   aes(label = tree_group, fill = bicat),
                   alpha = 0.9, min.segment.length = 0.1, size = 3) +
  scale_fill_manual(values = c('A_1'='#e8e8e8','A_2'='#cbb8d7','A_3'='#9972af',
                               'B_1'='#e4d9ac','B_2'='#c8ada0','B_3'='#976b82',
                               'C_1'='#c8b35a','C_2'='#af8e53','C_3'='#804d36'),
                    guide = 'none') +
  scale_y_reverse() +
  #coord_cartesian(xlim = c(5,23), ylim = c(0.65,0.3)) +
  labs(x = 'Pyroclimate exposure', y = 'Fire resistance') +
  theme_bw(base_size = 12) +
  theme(rect = element_rect(fill = "transparent"))
 
 ggsave('fia_fireshed_bivar.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs',
        height = 5, width = 5, dpi = 600)
 

 ggplot(fireshed_df) +
   geom_histogram(aes(x = tree_group, fill = bicat), stat="count", color = 'black') +
   scale_fill_manual(values = c('A_1'='#e8e8e8','A_2'='#cbb8d7','A_3'='#9972af',
                                'B_1'='#e4d9ac','B_2'='#c8ada0','B_3'='#976b82',
                                'C_1'='#c8b35a','C_2'='#af8e53','C_3'='#804d36'),
                     guide = 'none') +
   theme_bw(base_size = 12) +
   labs(x = 'Forest type', y = 'Number of firesheds') +
   scale_y_continuous(expand = c(.01, 0)) +
   theme(axis.text.y = element_text(color = 'black'),
         axis.title.x = element_blank(),
         axis.text.x = element_text(angle = 10, margin = margin(t = 6.5)),
         axis.line = element_line(colour = "black"),
         panel.grid = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         rect = element_rect(fill = "transparent"))
 
 ggsave('bivar_histogram.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs',
        height = 7.5, width = 6.5, dpi = 600)
 
# Colored by forest type
forest_type_1 <- 
 forest_type_egs %>% 
   group_by(tree_group) %>% 
   slice_sample(n = 1) 
 
 ggplot(fireshed_df, aes(x = adapt_dist, y = frs)) +
   geom_vline(aes(xintercept = 15)) +
   geom_vline(aes(xintercept = 20)) +
   geom_hline(aes(yintercept = 0.38)) +
   geom_hline(aes(yintercept = 0.54)) +
   geom_point(aes(fill = tree_group), 
              alpha = 0.7, color = 'black', shape = 21, size = 2) +
   geom_label_repel(data = forest_type_1,
                    aes(label = tree_group, fill = tree_group),
                    alpha = 0.9, min.segment.length = 0.1, size = 3, show.legend = FALSE) +
   scale_fill_brewer('Forest type', palette = 'Dark2', guide = 'none') +
   scale_y_reverse() +
   labs(x = 'Pyroclimate exposure', y = 'Fire resistance') +
   theme_bw(base_size = 12) 

 ggsave('bivar_scatter_type.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs',
        height = 5, width = 5.5, dpi = 600)
 

# Previous methods:
 # # Extract and summarize FIA data by fireshed
 # fireshed_sp_df <- split(vect(firesheds), 'Fireshed_I') %>% 
 #   map_dfr(function(fireshed){
 #     
 #     print(paste('Running fireshed', fireshed$Fireshed_I, fireshed$Fireshed_N))
 #     
 #     # Extract TreeMap values within fireshed polygon
 #     treemap_fs <- terra::crop(treemap, fireshed)
 #     tree_sheds <- terra::extract(treemap_fs, fireshed)
 #     tree_sheds$ID <- fireshed$Fireshed_I
 #     
 #     # Now summarize FIA data associated with this fireshed
 #     modal_dominant <- tree_sheds %>% 
 #       as_tibble() %>% 
 #       rename(tl_id = national_c2014_tree_list) %>% 
 #       filter(!is.na(tl_id)) %>% 
 #       left_join(., treemap_key) %>% 
 #       left_join(.,fia_sp_basal) %>% 
 #       group_by(ID) %>% 
 #       distinct() %>%
 #       mutate(fireshed_count = n()) %>% 
 #       group_by(ID, sp) %>% 
 #       summarise(sp_count_prop = n()/fireshed_count) %>% 
 #       ungroup() %>% 
 #       slice_max(n = 1, sp_count_prop, with_ties = F) %>% 
 #       rename(modal_dominant = sp)
 #     
 #     majority_basal <- tree_sheds %>% 
 #       as_tibble() %>% 
 #       rename(tl_id = national_c2014_tree_list) %>% 
 #       filter(!is.na(tl_id)) %>% 
 #       left_join(., treemap_key) %>% 
 #       left_join(.,fia_sp_basal) %>% 
 #       group_by(ID) %>% 
 #       distinct() %>% 
 #       mutate(shed_basal = sum(plot_basal, na.rm = T)) %>% 
 #       group_by(ID, sp) %>% 
 #       summarize(sp_basal_prop = sum(sp_basal, na.rm = T)/max(shed_basal)) %>% 
 #       ungroup() %>% 
 #       slice_max(n = 1, sp_basal_prop, with_ties = F) %>% 
 #       rename(majority_basal = sp)
 #     
 #     result <- full_join(modal_dominant, majority_basal)
 #     
 #     return(result)
 #   })
 # 
 # fireshed_sp_df
 fireshed_sp_treegroup_df <- firesheds %>% 
   st_drop_geometry() %>% 
   full_join(., fireshed_sp_df, by = c('Fireshed_I'='ID')) %>% 
   left_join(., spp_names, by = c('majority_basal'='sp')) %>% 
   filter(tree_group != 'Aspen') %>% 
   mutate(tree_group = fct_relevel(as.factor(tree_group), 
                                   levels = c('Douglas-fir',
                                              'Ponderosa pine',
                                              'Lodgepole pine',
                                              'Junipers and Pinyon pines',
                                              'Engelmann spruce',
                                              'True firs',
                                              'Redwood'))) 
 
 fireshed_sp_treegroup_df %>% 
   group_by(bicat, name) %>% 
   summarise(bicat_mode = n()) %>% 
   slice_max(n = 3, bicat_mode) %>% 
   View()
 
 forest_type_egs <- fireshed_sp_treegroup_df %>% 
   filter(Fireshed_I %in% c(633, 1326, 1983, 21, 822, 1915, 260, 164, 102, 462,
                            524, 759, 1228, 1057, 1186, 1763, 2140, 1622))
 
 fia_sp_basal <- fia_tree %>% 
   select(tl_id, fia_plot = CN, sp = SPCD, dbh = DIA) %>% 
   mutate(basal = pi*((dbh/2)^2)) %>% 
   group_by(fia_plot) %>% 
   mutate(plot_count = n()) %>% 
   group_by(fia_plot,sp) %>%
   summarize(sp_basal = sum(basal, na.rm = T)) %>% 
   group_by(fia_plot) %>% 
   mutate(plot_basal = sum(sp_basal, na.rm = T),
          sp_prop = sp_basal/plot_basal) %>% 
   slice_max(n = 1, sp_prop) 
 
# ------------------------------------------------------------------------------
# EMD DATA for all forested points
# ------------------------------------------------------------------------------
emd_chunk_df <- readRDS('data/emd/emd_2022_08_02.Rdata') %>% 
  group_by(cell_hist, cell_2C) %>% 
  summarize(multi_emd_50 = mean(multi_emd_50)) %>% 
  filter(!is.na(multi_emd_50))

# Join to dataframe of reference and future climate bins for all forested cells
forest_pyrome_df <- readRDS('data/process/forest_pyrome_df.Rdata')

pyrome_emd_df <- left_join(forest_pyrome_df, emd_chunk_df) %>% 
  #filter(!is.na(multi_emd_50)) %>% 
  select(x, y, multi_emd_50)

pyrome_emd_sf <- st_as_sf(pyrome_emd_df, coords = c('x', 'y'), crs = 4326) 

# ------------------------------------------------------------------------------
# FIA Approach
# ------------------------------------------------------------------------------
# Read in a compiled table for CONUS of FIA data as of 2019: https://www.fs.usda.gov/rds/archive/Catalog/RDS-2019-0026
fia_tree <- read_csv('data/composition/Tree_table_CONUS.csv',
                     col_names = c('tl_id','CN','INVYR','STATECD',
                                   'State_Abbreviation','UNITCD','COUNTYCD',
                                   'PLOT','SUBP','TREE','STATUSCD','SPCD','DIA',
                                   'HT','ACTUALHT','CR','TPA_UNADJ','TreeVisitID',
                                   'TreeID')) %>% 
  # Reduce down to states in our study area
  filter(State_Abbreviation %in% c('WA','OR','ID','MT','CA','NV','UT','WY','CO','AZ','NM'))


# Identify the dominant species, by basal area, at each plot (CN)
fia_basal <- fia_tree %>% 
  select(tl_id, fia_plot = CN, sp = SPCD, dbh = DIA) %>% 
  mutate(basal = pi*((dbh/2)^2)) %>% 
  group_by(fia_plot,sp) %>%
  summarize(sp_basal = sum(basal)) %>% 
  group_by(fia_plot) %>% 
  arrange(desc(sp_basal)) %>%
  slice_head(n = 1)

# Link to species and forest type descriptions/names
spp_ref <- read_csv('data/composition/REF_SPECIES.csv')
# Make a simplified df with code and common name, species, genus
spp_names <- spp_ref %>% 
  dplyr::select(sp = SPCD, tree_group = W_SPGRPCD) %>% #, name = COMMON_NAME, GENUS, SPECIES
  mutate(#latin_name = paste(GENUS, SPECIES),
         tree_group = factor(tree_group),
         tree_group = fct_recode(tree_group,
                                 'Douglas-fir'='10',   
                                 'Ponderosa and Jeffrey pines'='11',
                                 'True firs'='12', 
                                 'Western hemlock'='13',  
                                 'Sugar pine'='14',  
                                 'Western white pine'='15', 
                                 'Redwood'='16',   
                                 'Sitka spruce'='17',  
                                 'Engelmann and other spruces'='18',
                                 'Western larch'='19',
                                 'Incense-cedar'='20',
                                 'Lodgepole pine'='21',  
                                 'Western redcedar'='22' ,  
                                 'Pinyon pines and Junipers'='23' ,
                                 'Other western softwoods'='24',
                                 'Aspen'='44',
                                 'Alder'='45',
                                 'Other west hardwoods'='46',
                                 'Other west hardwoods'='47',
                                 'Other west hardwoods'='48'))

# Check how many records of this species are present in the FIA databases from these states 
left_join(fia_basal, spp_names) %>% 
  group_by(tree_group) %>%  #sp, name, latin_name
  tally(sort = TRUE)

fia_groups <- left_join(fia_basal, spp_names) %>% 
  select(fia_plot, sp, tree_group)

# I transformed this to the WGS 84 projection in QGIS 
treemap <- terra::rast('data/composition/treemap.tif')
# Resample tree map to mast_rast

# <-  mast_rast <- raster('data/process/mast_rast.tif')
# <- tree_mast <- terra::resample(treemap, mast_rast, method = 'near')


# Read in TreeMap raster linking to FIA plot numbers
treemap_key <- read_csv('data/composition/TL_CN_Lookup.txt') %>% 
  select(fia_plot = CN, tl_id)

# Extract FRS
frs_rast <- terra::rast('data/process/frs_coarse.tif')
frs_pts <- terra::extract(frs_rast, vect(pyrome_emd_sf))$frs_spatial

# Extract Tree map code 
tlid_pts <- terra::extract(tree_mast, vect(pyrome_emd_sf))$treemap

pyrome_preds_df <- pyrome_emd_df %>% 
  as_tibble() %>% 
  mutate(frs = frs_pts,
         tl_id = tlid_pts) %>%  
  filter(!is.na(frs), !is.na(tl_id)) %>% 
  left_join(treemap_key) %>% 
  left_join(., fia_groups)

tree_emd_df <- pyrome_preds_df %>% 
  group_by(tree_group) %>% 
  #filter(n() > 5000) %>% 
  summarise(mean_adapt = mean(multi_emd_50),
            lo_adapt = mean(multi_emd_50)-(1.96*(sd(multi_emd_50)/sqrt(650))),
            hi_adapt = mean(multi_emd_50)+(1.96*(sd(multi_emd_50)/sqrt(650))),
            mean_frs = mean(frs),
            lo_frs = mean(frs)-(1.96*(sd(frs)/sqrt(650))),
            hi_frs = mean(frs)+(1.96*(sd(frs)/sqrt(650))),
            n = n()) %>% 
  mutate(cat1 = case_when(mean_adapt <= 12 ~ 'A',
                          mean_adapt > 12 & mean_adapt <= 20 ~ 'B',
                          mean_adapt > 20 ~ 'C'),
         cat2 = case_when(mean_frs <= 0.38 ~ '3',
                          mean_frs > 0.38 & mean_frs <= 0.54 ~ '2',
                          mean_frs > 0.54 ~ '1')) %>% 
  unite('bicat', cat1, cat2) %>% 
  filter(!tree_group %in% c('Other west hardwoods','Other western softwoods',
                            'Aspen','Alder'))

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
ggplot(tree_emd_df, aes(x = mean_adapt, y = mean_frs)) +
  geom_vline(aes(xintercept = 12)) +
  geom_vline(aes(xintercept = 20)) +
  geom_hline(aes(yintercept = 0.38)) +
  geom_hline(aes(yintercept = 0.54)) +
  geom_errorbar(aes(xmin = lo_adapt, xmax = hi_adapt), width = 0, size = 0.5, alpha = 0.7) +
  geom_errorbar(aes(ymin = lo_frs, ymax = hi_frs), width = 0, size = 0.5, alpha = 0.7) +
  geom_point(aes(fill = bicat), color = 'black', shape = 21, size = 2) +
  geom_label_repel(aes(label = tree_group, fill = bicat), 
                   alpha = 0.9, min.segment.length = 0.1, size = 3) +
  scale_fill_manual(values = c('A_1'='#e8e8e8','A_2'='#cbb8d7','A_3'='#9972af',
                               'B_1'='#e4d9ac','B_2'='#c8ada0','B_3'='#976b82',
                               'C_1'='#c8b35a','C_2'='#af8e53','C_3'='#804d36'),
                    guide = 'none') +
  scale_y_reverse() +
  #coord_cartesian(xlim = c(5,23), ylim = c(0.65,0.3)) +
  labs(x = 'Pyroclimate exposure', y = 'Fire resistance') +
  theme_bw(base_size = 12) 

ggsave('fia_bivar_test1`.png', path = 'C:/Users/hoecker/Work/Postdoc/Future_Fire/Progress_Figs',
       height = 5, width = 5, dpi = 600)




# ------------------------------------------------------------------------------
# Add FRS, and Landfire EVT
# ------------------------------------------------------------------------------
frs_rast <- terra::rast('data/composition/frs_spatial.tif')
frs_rast <- terra::project(frs_rast, y = "epsg:4326")
lf_evt <- rast('data/composition/landfire_evt_wgs.tif')
landfire_key <- read_csv('data/composition/landfire_key.csv')
modified_keys <- read_csv('data/composition/landfire_types.csv')

veg_type <- full_join(landfire_key,modified_keys) %>% 
  select(VALUE, Combined_Type) %>% 
  filter(!is.na(Combined_Type))

pyrome_preds_df <- pyrome_emd_df %>% 
  mutate(frs = terra::extract(frs_rast, vect(pyrome_emd_sf))$frs_spatial,
         evt = terra::extract(lf_evt, vect(pyrome_emd_sf))$landfire_evt_wgs) %>%  
  filter(!is.na(frs)) %>% 
  left_join(., veg_type, by = c('evt' = 'VALUE'))

evt_emd_df <- pyrome_preds_df %>% 
  group_by(Combined_Type) %>% 
  #filter(n() > 5000) %>% 
  summarise(mean_adapt = mean(multi_emd_50),
            lo_adapt = mean(multi_emd_50)-(1.96*(sd(multi_emd_50)/sqrt(650))),
            hi_adapt = mean(multi_emd_50)+(1.96*(sd(multi_emd_50)/sqrt(650))),
            mean_frs = mean(frs),
            lo_frs = mean(frs)-(1.96*(sd(frs)/sqrt(650))),
            hi_frs = mean(frs)+(1.96*(sd(frs)/sqrt(650))),
            n = n()) %>% 
  mutate(cat1 = case_when(mean_adapt <= 12 ~ 'A',
                          mean_adapt > 12 & mean_adapt <= 20 ~ 'B',
                          mean_adapt > 20 ~ 'C'),
         cat2 = case_when(mean_frs <= 0.38 ~ '3',
                          mean_frs > 0.38 & mean_frs <= 0.54 ~ '2',
                          mean_frs > 0.54 ~ '1')) %>% 
  unite('bicat', cat1, cat2) 



filter(bicat=='A_1' & tree_group=='Redwood'|
         bicat=='A_1' & tree_group=='Douglas-fir'|
         bicat=='A_2' & tree_group=='Ponderosa and Jeffrey pines'|
         #bicat=='A_2' & tree_group=='Lodgepole pine'|
         bicat=='A_3' & tree_group=='Lodgepole pine'|
         bicat=='A_3' & tree_group=='Engelmann and other spruces'|
         bicat=='B_1' & tree_group=='Ponderosa and Jeffrey pines'|
         bicat=='B_1' & tree_group=='Douglas-fir'|
         bicat=='B_2' & tree_group=='Lodgepole pine'|
         bicat=='B_2' & tree_group=='Douglas-fir'|
         bicat=='B_3' & tree_group=='Lodgepole pine'|
         bicat=='B_3' & tree_group=='Pinyon pines and Junipers'|
         bicat=='C_1' & tree_group=='Ponderosa and Jeffrey pines'|
         bicat=='C_1' & tree_group=='Ponderosa and Jeffrey pines'|
         bicat=='C_2' & tree_group=='True firs'|
         bicat=='C_2' & tree_group=='Douglas-fir'|
         bicat=='C_3' & tree_group=='Pinyon pines and Junipers'|
         bicat=='C_3' & tree_group=='Engelmann and other spruces') %>% 
  
