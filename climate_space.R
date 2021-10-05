library(tmap)
library(tidyverse)
library(sf)
library(ncdf4)
library(hexbin)

# Trying out the climate space binning...
region9_pointgrid <- read_sf("/Users/tylerhoecker/Box/Work/PhD/GIS/doi_12_unified_regions_20180801_shapefile/region9_pointgrid.shp")

point_sample <- region9_pointgrid %>% 
  st_transform(., crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  select(Longitude = X, Latitude = Y) %>% 
  sample_n(1000) %>% 
  split(.,1:nrow(.))

# This function based on code from TerraClimate website: http://www.climatologylab.org/uploads/2/2/1/3/22133936/read_terraclimate_point.r
terraclim_dl <- function(coords, variable, product){
  
  # enter in longitude, latitude here
  coords <- c(coords[1,"Longitude"], coords[1,"Latitude"])
  # enter in variable you want to download see: http://thredds.northwestknowledge.net:8080/thredds/terraclimate_aggregated.html
  # Updated this part to use climatologies with URL from here: http://thredds.northwestknowledge.net:8080/thredds/catalog/TERRACLIMATE_ALL/summaries/catalog.html
  # '#fillmismatch' to fix netCDF update that broke code
  base_url <- paste0('http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/summaries/')
  prod_url <- paste0(base_url, paste0(product,'_',variable,'.nc#fillmismatch'))
  nc <- nc_open(prod_url)
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  flat = match(abs(lat - coords[2]) < 1/48, 1)
  latindex = which(flat %in% 1)
  flon = match(abs(lon - coords[1]) < 1/48, 1)
  lonindex = which(flon %in% 1)
  start <- c(lonindex, latindex, 1)
  count <- c(1, 1, -1)
  
  # read in the full period of record using aggregated files
  data <- as.numeric(ncvar_get(nc, varid = variable, start = start, count))
  nc_close(nc)
  # Calculate mean of summer months only 
  summer_mean <- mean(data[5:7])
  
  # Put all info indo dataframe
  result <- data.frame('lon' = coords[1],
                       'lat' = coords[2],
                       'period' = product,
                       'variable' = variable,
                       'stat' = 'summer_mean',
                       'value' = summer_mean)
  return(result)
}

variables <- list('aet','def')
periods <- list('TerraClimate19812010','TerraClimate2C')

terra_clim_df <- 
  map_df(periods, function(product){
    map_df(variables, function(variable){
      map_df(point_sample, function(coords){
        terraclim_dl(coords,variable,product)
      })
    })
  })

saveRDS(terra_clim_df, 'terra_clim_n1000.R')
terra_clim_df <- readRDS('terra_clim_n1000.R')

# Conver to wide format to plot in bivariate climate space  
terra_wide <- terra_clim_df %>% 
  # Some spurious values...
  filter(value < 200) %>% 
  pivot_wider(names_from = variable, values_from = value) %>% 
  group_by(lon,lat) %>% 
  mutate(ptID = cur_group_id())

n_bins <- 5
xrange <- c(0,200)
yrange <- c(0,110)

# Create hexbins and save the bin IDs
h <- hexbin(x = terra_wide$def, y = terra_wide$aet, 
            xbins= n_bins, 
            shape = 1, 
            IDs = TRUE, 
            xbnds = xrange, ybnds = yrange)

hexdf <- data.frame(hcell2xy (h),  hexID = h@cell, counts = h@count)

terra_wide <- terra_wide %>% 
  mutate(cell = hexbin(x = def, y = aet, 
                       xbins= n_bins, 
                       shape = 1, 
                       IDs = TRUE, 
                       xbnds = xrange, ybnds = yrange)@cID) %>% 
  group_by(ptID) %>% 
  mutate(change = as.factor(ifelse(length(unique(cell)) > 1, 1, 0)))

terra_plot <- arrange(terra_wide, ptID)[300:700,]

# Plot in bivariate climate space 
ggplot(terra_plot, aes(x = def, y = aet)) +
  geom_hex(data = hexdf, aes(x = x, y = y, fill = hexID), color = 'black', alpha = 0.8, stat = 'identity') +
  geom_text(data = hexdf, aes(x = x, y = y, label = hexID)) +
  geom_point(aes(color = change), alpha = 0.5, size = 1) +
  geom_line(aes(group = ptID, color = change), alpha = 0.5) +
  scale_fill_gradient2('', low = '#8c510a', mid = '#f6e8c3', high = '#01665e', midpoint = 18) +
  scale_color_manual('+2C group shift?', values = c('white','black'), labels = c('No','Yes')) +
  labs(y = 'AET', x = 'DEF', title = 'Mean Summer Climate') +
  theme_bw() +
  coord_cartesian(xlim = xrange, ylim = yrange)

ggsave('climate_hex.png', height = 5, width = 7, units = 'in', dpi = 300)

terra_sf <- terra_wide %>% 
  filter(period == 'TerraClimate19812010') %>% 
  st_as_sf(coords = c('lon','lat'), crs = 4326)

# Make these climate data spatial
boundary <- read_sf("/Users/tylerhoecker/Box/Work/PhD/GIS/doi_12_unified_regions_20180801_shapefile/region9_boundary.shp") %>% 
  st_transform(., crs = 4326) 

# Plot on map
climate_map <- 
  tm_shape(boundary) +
  tm_borders() +
  tm_fill(col = 'grey90') +
  tm_shape(terra_sf) +
  tm_dots(title = 'Climate Bin',
          col = 'cell',
          size = 0.25, 
          palette = c('#8c510a','#f6e8c3','#01665e'),
          style = 'cont',
          shape = 15) +
  tm_text('cell', size = 0.5)


tmap_save(climate_map, 'climate_map.png', height = 5, width = 7, units = 'in', dpi = 300)


# Centroids of MTBS perimeters for an arbitrary rectangle around the NW
# These are not an even grid, the centroids for every fire perimeter, not spatially even
# mtbs_centroids <- read_sf("/Users/tylerhoecker/Box/Work/PhD/GIS/fire_data/mtbs_1984_2018/nw_mtbs_centroids.shp")
# mtbs_latlon <- st_coordinates(mtbs_centroids) %>% 
#   as.data.frame() %>% 
#   mutate(name = 'fire')
# 
# # Get a spatially even subset
# mtbs_4km <- 
#   spThin::thin(mtbs_latlon, lat.col = 'Y', long.col = 'X',
#                spec.col = 'name',
#                # KEY ARG: THIN RECORDS WITHIN 250 m (0.25 km)
#                thin.par = 4,
#                reps = 1,
#                write.files = F,
#                locs.thinned.list.return = T) 
# # It's a 1-long list, simplify
# mtbs_4km <- mtbs_4km[[1]]
# mtbs_4km_l <- split(mtbs_4km, 1:nrow(mtbs_4km))
# saveRDS(mtbs_4km_l, 'mtbs_4km_l.R')
mtbs_4km_l <- readRDS(file = 'mtbs_4km_l.R')
