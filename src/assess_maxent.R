library(sp)
library(dplyr)
library(ggplot2)
library(lubridate)
library(raster)
library(foreach)

# Spatial data
## Colony location and foraging range
wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
hi_aea_prj <- CRS('+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-163 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
leh_col <- data.frame(x = -160.1, y = 22.0)
leh_forage <- SpatialPoints(leh_col,
                            proj4string = wgs84_prj) %>%
  spTransform(hi_aea_prj) %>%
  rgeos::gBuffer(width = 250e3) %>%
  spTransform(wgs84_prj)
mcb_col <- data.frame(x = -157.7, y = 21.5)
mcb_forage <- SpatialPoints(mcb_col,
                            proj4string = wgs84_prj) %>%
  spTransform(hi_aea_prj) %>%
  rgeos::gBuffer(width = 250e3) %>%
  spTransform(wgs84_prj)

# Presence records
## Tracks
### Read from CSV
leh_mcb_tracks <- readr::read_csv('data/LEH_MCB_tracks.csv') %>%
  filter(Year > 2013)

### Sample 5 points per bird, per day, per deployment session
set.seed(2025)
leh_mcb_sample <- leh_mcb_tracks %>%
  group_by(Year, DepSess, DeployID, as.Date(TimestampUTC)) %>%
  filter(n() >= 20) %>%
  sample_n(5) %>%
  ungroup
leh_mcb_sample %>%
  group_by(SubColCode, Year, DepSess) %>%
  do(write = readr::write_csv(transmute(.,
                                        Species = 'RFBO',
                                        Longitude,
                                        Latitude), 
                              sprintf('data/out/Presences/%s/%i_%i.csv',
                                      first(.$SubColCode),
                                      first(.$Year),
                                      first(.$DepSess))))

# Aggregate environmental variables over deployment sessions
## SST & SST anomaly ncdf files
sst_nc_2014 <- RNetCDF::open.nc('data/sst/2014sst.nc')
sst_nc_2015 <- RNetCDF::open.nc('data/sst/2015sst.nc')
time_dim <- 'time'     # seconds since 1970-01-01 00:00:00.000 UTC
lat_dim <- 'latitude'  # degrees north
lon_dim <- 'longitude' # degrees east
sst_var <- 'sst'
anom_var <- 'anom'
time_sst_2014 <- RNetCDF::var.get.nc(sst_nc_2014, time_dim)
time_sst_2015 <- RNetCDF::var.get.nc(sst_nc_2015, time_dim)
t_max_2014 <- max(time_sst_2014)
t_min_2014 <- min(time_sst_2014)
t_max_2015 <- max(time_sst_2015)
t_min_2015 <- min(time_sst_2015)
lat_sst <- RNetCDF::var.get.nc(sst_nc_2014, lat_dim)
y_max <- max(lat_sst)
y_min <- min(lat_sst)
lon_sst <- RNetCDF::var.get.nc(sst_nc_2014, lon_dim)
x_max <- max(lon_sst)
x_min <- min(lon_sst)
sst_arr_2014 <- RNetCDF::var.get.nc(sst_nc_2014, sst_var)
anom_arr_2014 <- RNetCDF::var.get.nc(sst_nc_2014, anom_var)
sst_arr_2015 <- RNetCDF::var.get.nc(sst_nc_2015, sst_var)
anom_arr_2015 <- RNetCDF::var.get.nc(sst_nc_2015, anom_var)

## Deployment session dates
leh_depsess <- list(list(label = '2014_1',
                         dates = seq(ymd('2014-05-13', tz = 'UTC'), 
                                     ymd('2014-05-18', tz = 'UTC'), 
                                     by = '1 day')),
                    list(label = '2014_2',
                         dates = seq(ymd('2014-06-13', tz = 'UTC'), 
                                     ymd('2014-06-18', tz = 'UTC'), 
                                     by = '1 day')),
                    list(label = '2014_3',
                         dates = seq(ymd('2014-07-14', tz = 'UTC'), 
                                     ymd('2014-07-20', tz = 'UTC'), 
                                     by = '1 day')),
                    list(label = '2015_1',
                         dates = seq(ymd('2015-05-26', tz = 'UTC'), 
                                     ymd('2015-06-04', tz = 'UTC'), 
                                     by = '1 day')),
                    list(label = '2015_2',
                         dates = seq(ymd('2015-06-27', tz = 'UTC'), 
                                     ymd('2015-07-08', tz = 'UTC'), 
                                     by = '1 day')))
mcb_depsess <- list(list(label = '2014_1',
                         dates = seq(ymd('2014-06-01', tz = 'UTC'), 
                                     ymd('2014-06-07', tz = 'UTC'), 
                                     by = '1 day')),
                    list(label = '2015_1',
                         dates = seq(ymd('2015-06-17', tz = 'UTC'), 
                                     ymd('2015-07-06', tz = 'UTC'), 
                                     by = '1 day')),
                    list(label = '2015_2',
                         dates = seq(ymd('2015-06-29', tz = 'UTC'), 
                                     ymd('2015-07-09', tz = 'UTC'), 
                                     by = '1 day')))

## Clear out existing files
clear_leh <- foreach(depsess = leh_depsess) %do% {
  depsess_folder <- sprintf('data/out/Environment/LEH/%s', depsess$label)
  dir(depsess_folder, full.names = TRUE, pattern = 'asc') %>%
    sapply(file.remove)
}
clear_mcb <- foreach(depsess = leh_depsess) %do% {
  depsess_folder <- sprintf('data/out/Environment/MCB/%s', depsess$label)
  dir(depsess_folder, full.names = TRUE, pattern = 'asc') %>%
    sapply(file.remove)
}

## Create environmental aggregates
### Lehua
leh_template <- raster('data/out/EnergyLandscapes/all/LEH_20140513_in.tif',
                       crs = hi_aea_prj) %>%
  projectRaster(crs = wgs84_prj)
res(leh_template) = 0.25
leh_env <- foreach(depsess = leh_depsess) %do% {
  depsess_extent <- leh_forage
  template <- leh_template
  col_loc <- leh_col
  if(year(depsess$dates[1]) == 2014) {
    time_sst <- time_sst_2014
    t_max <- t_max_2014
    t_min <- t_min_2014
    sst_arr <- sst_arr_2014
    anom_arr <- anom_arr_2014
  } else {
    time_sst <- time_sst_2015
    t_max <- t_max_2015
    t_min <- t_min_2015
    sst_arr <- sst_arr_2015
    anom_arr <- anom_arr_2015
  }
  posix_depsess <- as.POSIXct(depsess$dates, 
                              origin = ymd('1970-01-01', tz = 'UTC'), 
                              tz = 'UTC')
  depsess_t_min <- as.numeric(first(posix_depsess))
  depsess_t_max <- as.numeric(last(posix_depsess))
  sst_t_min <- round((length(time_sst) - 1) * (depsess_t_min - t_min) / (t_max - t_min)) + 1
  sst_t_max <- round((length(time_sst) - 1) * (depsess_t_max - t_min) / (t_max - t_min)) + 1
  sst_r <- sst_arr[,,seq(sst_t_min, sst_t_max)] %>%
    apply(c(1,2), mean, na.rm = TRUE) %>%
    t %>% apply(2, rev) %>%
    raster(xmn = x_min - 360, xmx = x_max - 360, 
           ymn = y_min, ymx = y_max, 
           crs = wgs84_prj) %>%
    resample(template) %>%
    crop(template)
  anom_r <- anom_arr[,,seq(sst_t_min, sst_t_max)] %>%
    apply(c(1,2), mean, na.rm = TRUE) %>%
    t %>% apply(2, rev) %>%
    raster(xmn = x_min - 360, xmx = x_max - 360, 
           ymn = y_min, ymx = y_max, 
           crs = wgs84_prj) %>%
    resample(template) %>%
    crop(template)
  depsess_folder <- sprintf('data/out/Environment/LEH/%s', depsess$label)
  writeRaster(sst_r, 
              sprintf('%s/sst', depsess_folder),
              'ascii',
              overwrite = TRUE)
  writeRaster(anom_r, 
              sprintf('%s/anom', depsess_folder), 
              'ascii',
              overwrite = TRUE)
}
### MCB
mcb_template <- raster('data/out/EnergyLandscapes/all/MCB_20140601_in.tif',
                       crs = hi_aea_prj) %>%
  projectRaster(crs = wgs84_prj)
res(mcb_template) = 0.25
mcb_env <- foreach(depsess = mcb_depsess) %do% {
  depsess_extent <- mcb_forage
  col_loc <- mcb_col
  if(year(depsess$dates[1]) == 2014) {
    time_sst <- time_sst_2014
    t_max <- t_max_2014
    t_min <- t_min_2014
    sst_arr <- sst_arr_2014
    anom_arr <- anom_arr_2014
  } else {
    time_sst <- time_sst_2015
    t_max <- t_max_2015
    t_min <- t_min_2015
    sst_arr <- sst_arr_2015
    anom_arr <- anom_arr_2015
  }
  posix_depsess <- as.POSIXct(depsess$dates, 
                              origin = ymd('1970-01-01', tz = 'UTC'), 
                              tz = 'UTC')
  depsess_t_min <- as.numeric(first(posix_depsess))
  depsess_t_max <- as.numeric(last(posix_depsess))
  sst_t_min <- round((length(time_sst) - 1) * (depsess_t_min - t_min) / (t_max - t_min)) + 1
  sst_t_max <- round((length(time_sst) - 1) * (depsess_t_max - t_min) / (t_max - t_min)) + 1
  sst_r <- sst_arr[,,seq(sst_t_min, sst_t_max)] %>%
    apply(c(1,2), mean, na.rm = TRUE) %>%
    t %>% apply(2, rev) %>%
    raster(xmn = x_min - 360, xmx = x_max - 360, 
           ymn = y_min, ymx = y_max, 
           crs = wgs84_prj) %>%
    resample(template) %>%
    crop(template)
  anom_r <- anom_arr[,,seq(sst_t_min, sst_t_max)] %>%
    apply(c(1,2), mean, na.rm = TRUE) %>%
    t %>% apply(2, rev) %>%
    raster(xmn = x_min - 360, xmx = x_max - 360, 
           ymn = y_min, ymx = y_max, 
           crs = wgs84_prj) %>%
    resample(template) %>%
    crop(template)
  depsess_folder <- sprintf('data/out/Environment/MCB/%s', depsess$label)
  writeRaster(sst_r, 
              sprintf('%s/sst', depsess_folder), 
              'ascii',
              overwrite = TRUE)
  writeRaster(anom_r, 
              sprintf('%s/anom', depsess_folder), 
              'ascii',
              overwrite = TRUE)
}

# Aggregate Energy Landscapes
leh_el <- foreach(depsess = leh_depsess) %do% {
  template <- leh_template
  depsess_folder <- sprintf('data/out/Environment/LEH/%s', depsess$label)
  depsess_extent <- leh_forage
  col_loc <- leh_col
  posix_depsess <- as.POSIXct(depsess$dates, 
                              origin = ymd('1970-01-01', tz = 'UTC'), 
                              tz = 'UTC')
  
  ert_stack <- sprintf('data/out/EnergyLandscapes/all/LEH_%s_rt.tif', 
                       format(posix_depsess, '%Y%m%d')) %>% 
    stack
  projection(ert_stack) <- hi_aea_prj
  mean_ert <- mean(ert_stack)
  projectRaster(mean_ert, crs = wgs84_prj) %>% 
    resample(template) %>%
    crop(template) %>%
    writeRaster(sprintf('%s/Ert', depsess_folder), 
                'ascii',
                overwrite = TRUE)
}
mcb_el <- foreach(depsess = mcb_depsess) %do% {
  template <- mcb_template
  depsess_folder <- sprintf('data/out/Environment/MCB/%s', depsess$label)
  depsess_extent <- mcb_forage
  col_loc <- mcb_col
  posix_depsess <- as.POSIXct(depsess$dates, 
                              origin = ymd('1970-01-01', tz = 'UTC'), 
                              tz = 'UTC')
  
  ert_stack <- sprintf('data/out/EnergyLandscapes/all/MCB_%s_rt.tif', 
                       format(posix_depsess, '%Y%m%d')) %>% 
    stack
  projection(ert_stack) <- hi_aea_prj
  mean_ert <- mean(ert_stack)
  projectRaster(mean_ert, crs = wgs84_prj) %>% 
    resample(template) %>%
    crop(template) %>%
    writeRaster(sprintf('%s/Ert', depsess_folder),
                'ascii',
                overwrite = TRUE)
}

# D2C rasters
leh_d2c <- foreach(depsess = leh_depsess) %do% {
  template <- leh_template
  depsess_folder <- sprintf('data/out/Environment/LEH/%s', depsess$label)
  leh_col %>%
    SpatialPoints(proj4string = wgs84_prj) %>%
    distanceFromPoints(template, .) %>%
    resample(template) %>%
    crop(template) %>%
    writeRaster(sprintf('%s/D2C', depsess_folder),
                'ascii',
                overwrite = TRUE)
}
mcb_d2c <- foreach(depsess = mcb_depsess) %do% {
  template <- mcb_template
  depsess_folder <- sprintf('data/out/Environment/MCB/%s', depsess$label)
  mcb_col %>%
    SpatialPoints(proj4string = wgs84_prj) %>%
    distanceFromPoints(template, .) %>%
    resample(template) %>%
    crop(template) %>%
    writeRaster(sprintf('%s/D2C', depsess_folder),
                'ascii',
                overwrite = TRUE)
}

# Bathymetry rasters
leh_bathy <- foreach(depsess = leh_depsess) %do% {
  template <- leh_template
  depsess_folder <- sprintf('data/out/Environment/LEH/%s', depsess$label)
  raster('data/bathy/etopo1.tif') %>%
    resample(template) %>%
    crop(template) %>%
    writeRaster(sprintf('%s/bathy', depsess_folder),
                'ascii',
                overwrite = TRUE)
}
mcb_bathy <- foreach(depsess = mcb_depsess) %do% {
  template <- mcb_template
  depsess_folder <- sprintf('data/out/Environment/MCB/%s', depsess$label)
  raster('data/bathy/etopo1.tif') %>%
    resample(template) %>%
    crop(template) %>%
    writeRaster(sprintf('%s/bathy', depsess_folder),
                'ascii',
                overwrite = TRUE)
}

# CRW UD rasters
leh_crw <- foreach(depsess = leh_depsess) %do% {
  template <- leh_template
  depsess_folder <- sprintf('data/out/Environment/LEH/%s', depsess$label)
  raster('data/out/CyberBirds/LEH_CRW_UD.tif') %>%
    projectRaster(crs = wgs84_prj) %>%
    resample(template) %>%
    crop(template) %>%
    writeRaster(sprintf('%s/UD', depsess_folder),
                'ascii',
                overwrite = TRUE)
}
mcb_crw <- foreach(depsess = mcb_depsess) %do% {
  template <- mcb_template
  depsess_folder <- sprintf('data/out/Environment/MCB/%s', depsess$label)
  raster('data/out/CyberBirds/MCB_CRW_UD.tif') %>%
    projectRaster(crs = wgs84_prj) %>%
    resample(template) %>%
    crop(template) %>%
    writeRaster(sprintf('%s/UD', depsess_folder),
                'ascii',
                overwrite = TRUE)
}

# Align rasters
align_leh <- foreach(depsess = leh_depsess) %do% {
  template <- leh_template
  depsess_folder <- sprintf('data/out/Environment/LEH/%s', depsess$label)
  dir(depsess_folder, pattern = 'asc', full.names = TRUE) %>%
    lapply(function(r) {
      raster(r) %>%
        resample(template) %>%
        crop(template) %>%
        writeRaster(r, 'ascii', overwrite = TRUE)
    })
}
align_mcb <- foreach(depsess = mcb_depsess) %do% {
  template <- mcb_template
  depsess_folder <- sprintf('data/out/Environment/MCB/%s', depsess$label)
  dir(depsess_folder, pattern = 'asc', full.names = TRUE) %>%
    lapply(function(r) {
      raster(r) %>%
        resample(template) %>%
        crop(template) %>%
        writeRaster(r, 'ascii', overwrite = TRUE)
    })
}

