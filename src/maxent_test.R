library(xtractomatic)
library(dismo)
library(sp)
library(dplyr)
library(ggplot2)

# Colony location and foraging range
wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
hi_aea_prj <- CRS('+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-157 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
kpc_col <- data.frame(x = -159.40, y = 22.23)
kpc_forage <- SpatialPoints(kpc_col,
                            proj4string = wgs84_prj) %>%
  spTransform(hi_aea_prj) %>%
  rgeos::gBuffer(width = 250e3) %>%
  spTransform(wgs84_prj)

# Environmental data
sst_id <- 'jplG1SST'
chla_id <- 'erdMH1pp1day'
bathy_id <- 'ETOPO180'
env_x <- bbox(kpc_forage)['x',]
env_y <- bbox(kpc_forage)['y',] 
env_t <- c('2016-06-01', '2016-06-01')
sst_data <- xtracto_3D(sst_id, env_x, env_y, env_t)
chl_data <- xtracto_3D(chla_id, env_x, env_y, env_t)
bathy_data <- xtracto_3D(bathy_id, env_x, env_y, env_t)

# Utility functions for ggplot compatibility
fortify_cube <- function(x, y, t, values) {
  if(!all.equal(dim(values), c(length(x), length(y), length(t)))) 
    stop('!all.equal(dim(values), c(length(x), length(y), length(t))')
  i = seq_along(x)
  j = seq_along(y)
  k = seq_along(t)
  
  if('POSIXlt' %in% class(t))
    t <- as.POSIXct(t)
  
  expand.grid(i = i, j = j, k = k) %>%
    mutate(x = x[i],
           y = y[j],
           t = t[k],
           value = mapply(function(i, j, k) values[i, j, k], i, j, k)) %>%
    select(-i, -j, -k)
}

fortify_xtracto <- function(xtracto_result) {
  with(xtracto_result, fortify_cube(longitude, latitude, time, data))
}

# Plot variables
fortify_xtracto(sst_data) %>%
  ggplot(aes(x = x, y = y, fill = value)) + 
  geom_raster()

fortify_xtracto(chl_data) %>%
  ggplot(aes(x = x, y = y, fill = value)) + 
  geom_raster()

fortify_xtracto(bathy_data) %>%
  ggplot(aes(x = x, y = y, fill = value)) + 
  geom_raster()
