library(dplyr)
library(lubridate)
library(raster)
library(foreach)
source('src/energy_landscape.R')

# Simple test
## Test variables
origin <- c(0, 0)
radius <- 100
res <- 10
wind_u <- matrix(seq(from = -5, to = 10, length.out = 100), nrow = 10) %>%
  raster(xmn = origin[1] - radius, xmx = origin[1] + radius, 
         ymn = origin[2] - radius, ymx = origin[2] + radius)
wind_v <- matrix(seq(from = 5, to = -10, length.out = 100), nrow = 10) %>%
  raster(xmn = origin[1] - radius, xmx = origin[1] + radius, 
         ymn = origin[2] - radius, ymx = origin[2] + radius)

## Movement models
load('data/out/Models/EnergyModels.Rdata')
load('data/out/Models/SgTCModel.RData')
energy_mod <- function(a, m) {
  predict(oam, 
          newdata = data.frame(WindAngle = a,
                               WindSpd = m,
                               DeployID = 1145),
          exclude = 's(DeployID)',
          type = 'response')
}
dur_mod <- function(a, m, d) {
  t <- a * cos(a)
  spd <- predict(SgTCModel, 
                 re.form = NA,
                 newdata = data.frame(Tailwind = t))
  d / spd
}

## Barriers
poly <- Polygon(cbind(c(-80, -40, -40, -80), c(-80, -80, -40, -40)))
polys <- Polygons(srl = list(poly), ID = 1)
test_barr <- SpatialPolygons(list(polys))

## Run test
test_el <- energy_landscape(origin, radius, res, wind_u, wind_v,
                            energy_mod, dur_mod, test_barr)
print(plot_energy_landscape(test_el, origin, barriers = test_barr))

# Real data test
## Movement models
load('data/out/Models/EnergyModels.Rdata')
load('data/out/Models/SgTCModel.RData')
energy_mod <- function(a, m) {
  predict(oam, 
          newdata = data.frame(WindAngle = a,
                               WindSpd = m,
                               DeployID = 1145),
          exclude = 's(DeployID)',
          type = 'response')
}
dur_mod <- function(a, m, d) {
  t <- a * cos(a)
  spd <- predict(SgTCModel, 
                 re.form = NA,
                 newdata = data.frame(Tailwind = t))
  d / spd
}

wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
hi_aea_prj <- CRS('+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-163 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
project_col <- function(lat, lon) {
  data.frame(x = lon, y = lat) %>%
    SpatialPoints(wgs84_prj) %>%
    spTransform(hi_aea_prj) %>%
    as.data.frame %>%
    as.numeric
}
kpc_col <- project_col(22.2, -159.4)
leh_col <- project_col(22.0, -160.1)
mcb_col <- project_col(21.5, -157.7)

radius <- 250e3
## Resolution should be mean transit length (mean transit speed * 20 min)
res <- 11.4 * 20 * 60

## Wind data
## Each day's wind is the mean wind between the hours of 7am and 7pm HST
## This is when the birds are most active
daily_wind <- function(t, uv) {
  wind_path <- 'data/wind/WRF_HI/2016_WRF_Hawaii_Regional_Atmospheric_Model_best.ncd.nc'
  if(!(uv %in% c('u', 'v')))
    stop('uv must be "u" or "v"')
  
  wind_nc <- RNetCDF::open.nc(wind_path)
  time_dim <- 'time'   # hours since 2010-05-14 00:00:00.000 UTC
  lat_dim <- 'lat'     # degrees north
  lon_dim <- 'lon'     # degrees east
  u_var <- 'Uwind'
  v_var <- 'Vwind'
  time_wind_origin <- ymd('2010-05-14', tz = 'UTC')
  time_wind <- hours(RNetCDF::var.get.nc(wind_nc, time_dim)) + time_wind_origin
  lat_wind <- RNetCDF::var.get.nc(wind_nc, lat_dim)
  lon_wind <- RNetCDF::var.get.nc(wind_nc, lon_dim)
  u_wind <- RNetCDF::var.get.nc(wind_nc, u_var)
  v_wind <- RNetCDF::var.get.nc(wind_nc, v_var)
  wind_arr <- if(uv == 'u') u_wind else v_wind 
  
  wind_times <- sprintf('%s %i', as.Date(t, tz = 'US/Hawaii'), 7:19) %>% 
    ymd_h(tz = 'US/Hawaii')
  ks <- sapply(wind_times, function(t) which.min(abs(t - time_wind)))
  wgs84_wind <- apply(wind_arr[,,ks], c(1,2), mean, na.rm = TRUE) %>%
    t %>% apply(2, rev) %>%
    raster(xmn = min(lon_wind),
           xmx = max(lon_wind),
           ymn = min(lat_wind),
           ymx = max(lat_wind),
           crs = wgs84_prj)
  
  hi_aea_template <- projectExtent(wgs84_wind, hi_aea_prj)
  res(hi_aea_template) <- 6e3
  projectRaster(wgs84_wind, hi_aea_template)
}

hi_land <- rgdal::readOGR('data/coastline/hi_land', 'hi_land') %>%
  spTransform(hi_aea_prj)

## try one date
el_0528 <- energy_landscape(kpc_col, radius, res, 
                            daily_wind(ymd('2016-05-28', tz = 'US/Hawaii'), 'u'),
                            daily_wind(ymd('2016-05-28', tz = 'US/Hawaii'), 'v'),
                            energy_mod, dur_mod, hi_land);beepr::beep();
plot_energy_landscape(el_0528, kpc_col, 'rt_cost', hi_land)

## try lehua
el_leh_0615 <- energy_landscape(leh_col, radius, res, 
                                daily_wind(ymd('2016-06-15', tz = 'US/Hawaii'), 'u'),
                                daily_wind(ymd('2016-06-15', tz = 'US/Hawaii'), 'v'),
                                energy_mod, dur_mod, hi_land);beepr::beep();
plot_energy_landscape(el_leh_0615, leh_col, 'rt_cost')

## all dates, KPC, rasters + plots
foreach(el_date = ymd('2016-05-28', tz = 'US/Hawaii') + days(0:50)) %do% {
  el <- energy_landscape(location, radius, res, 
                         daily_wind(el_date, 'u'),
                         daily_wind(el_date, 'v'),
                         energy_mod, dur_mod, hi_land)
  writeRaster(el[['out_cost']], 
              sprintf('data/out/EnergyLandscapes/KPC/%s_out.tif',
                      format(el_date, '%Y%m%d')),
              'GTiff')
  writeRaster(el[['in_cost']], 
              sprintf('data/out/EnergyLandscapes/KPC/%s_in.tif',
                      format(el_date, '%Y%m%d')),
              'GTiff')
  writeRaster(el[['rt_cost']], 
              sprintf('data/out/EnergyLandscapes/KPC/%s_rt.tif',
                      format(el_date, '%Y%m%d')),
              'GTiff')
  p <- plot_energy_landscape(el, location, 'out_cost', hi_land)
  ggsave(sprintf('data/out/EnergyLandscapes/KPC/%s_out.png',
                 format(el_date, '%Y%m%d')),
         p, width = 5, height = 5, units = 'in')
  p <- plot_energy_landscape(el, location, 'in_cost', hi_land)
  ggsave(sprintf('data/out/EnergyLandscapes/KPC/%s_in.png',
                 format(el_date, '%Y%m%d')),
         p, width = 5, height = 5, units = 'in')
  p <- plot_energy_landscape(el, location, 'rt_cost', hi_land)
  ggsave(sprintf('data/out/EnergyLandscapes/KPC/%s_rt.png',
                 format(el_date, '%Y%m%d')),
         p, width = 5, height = 5, units = 'in')
}
## all dates, LEH, rasters + plots
daily_wind_2014 <- function(t, uv) {
  wind_path <- 'data/wind/WRF_HI/2014_WRF_Hawaii_Regional_Atmospheric_Model_best.ncd.nc'
  wind_nc <- RNetCDF::open.nc(wind_path)
  time_dim <- 'time'   # hours since 2010-05-14 00:00:00.000 UTC
  lat_dim <- 'lat'     # degrees north
  lon_dim <- 'lon'     # degrees east
  u_var <- 'Uwind'
  v_var <- 'Vwind'
  time_wind_origin <- ymd('2010-05-14', tz = 'UTC')
  time_wind <- hours(RNetCDF::var.get.nc(wind_nc, time_dim)) + time_wind_origin
  lat_wind <- RNetCDF::var.get.nc(wind_nc, lat_dim)
  lon_wind <- RNetCDF::var.get.nc(wind_nc, lon_dim)
  u_wind <- RNetCDF::var.get.nc(wind_nc, u_var)
  v_wind <- RNetCDF::var.get.nc(wind_nc, v_var)
  if(!any(between(t, min(time_wind), max(time_wind))))
    stop('t out of bounds')
  if(!(uv %in% c('u', 'v')))
    stop('uv must be "u" or "v"')
  wind_arr <- if(uv == 'u') u_wind else v_wind 
  wind_times <- sprintf('%s %i', as.Date(t, tz = 'US/Hawaii'), 7:19) %>% 
    ymd_h(tz = 'US/Hawaii')
  ks <- sapply(wind_times, function(t) which.min(abs(t - time_wind)))
  wgs84_wind <- apply(wind_arr[,,ks], c(1,2), mean, na.rm = TRUE) %>%
    t %>%
    raster(xmn = min(lon_wind),
           xmx = max(lon_wind),
           ymn = min(lat_wind),
           ymx = max(lat_wind),
           crs = wgs84_prj)
  
  hi_aea_template <- projectExtent(wgs84_wind, hi_aea_prj)
  res(hi_aea_template) <- 6e3
  projectRaster(wgs84_wind, hi_aea_template)
}
daily_wind_2015 <- function(t, uv) {
  wind_path <- 'data/wind/WRF_HI/2015_WRF_Hawaii_Regional_Atmospheric_Model_best.ncd.nc'
  wind_nc <- RNetCDF::open.nc(wind_path)
  time_dim <- 'time'   # hours since 2010-05-14 00:00:00.000 UTC
  lat_dim <- 'lat'     # degrees north
  lon_dim <- 'lon'     # degrees east
  u_var <- 'Uwind'
  v_var <- 'Vwind'
  time_wind_origin <- ymd('2010-05-14', tz = 'UTC')
  time_wind <- hours(RNetCDF::var.get.nc(wind_nc, time_dim)) + time_wind_origin
  lat_wind <- RNetCDF::var.get.nc(wind_nc, lat_dim)
  lon_wind <- RNetCDF::var.get.nc(wind_nc, lon_dim)
  u_wind <- RNetCDF::var.get.nc(wind_nc, u_var)
  v_wind <- RNetCDF::var.get.nc(wind_nc, v_var)
  if(!any(between(t, min(time_wind), max(time_wind))))
    stop('t out of bounds')
  if(!(uv %in% c('u', 'v')))
    stop('uv must be "u" or "v"')
  wind_arr <- if(uv == 'u') u_wind else v_wind 
  wind_times <- sprintf('%s %i', as.Date(t, tz = 'US/Hawaii'), 7:19) %>% 
    ymd_h(tz = 'US/Hawaii')
  ks <- sapply(wind_times, function(t) which.min(abs(t - time_wind)))
  wgs84_wind <- apply(wind_arr[,,ks], c(1,2), mean, na.rm = TRUE) %>%
    t %>%
    raster(xmn = min(lon_wind),
           xmx = max(lon_wind),
           ymn = min(lat_wind),
           ymx = max(lat_wind),
           crs = wgs84_prj)
  
  hi_aea_template <- projectExtent(wgs84_wind, hi_aea_prj)
  res(hi_aea_template) <- 6e3
  projectRaster(wgs84_wind, hi_aea_template)
}
leh_2014_depsess1 <- seq(ymd('2014-05-13', tz = 'US/Hawaii'), 
                         ymd('2014-05-18', tz = 'US/Hawaii'), 
                         by = '1 day')
leh_2014_depsess2 <- seq(ymd('2014-06-13', tz = 'US/Hawaii'), 
                         ymd('2014-06-18', tz = 'US/Hawaii'), 
                         by = '1 day')
leh_2014_depsess3 <- seq(ymd('2014-07-14', tz = 'US/Hawaii'), 
                         ymd('2014-07-20', tz = 'US/Hawaii'), 
                         by = '1 day')
leh_2015_depsess1 <- seq(ymd('2015-05-26', tz = 'US/Hawaii'), 
                         ymd('2015-06-04', tz = 'US/Hawaii'), 
                         by = '1 day')
leh_2015_depsess2 <- seq(ymd('2015-06-27', tz = 'US/Hawaii'), 
                         ymd('2015-07-05', tz = 'US/Hawaii'), 
                         by = '1 day')
leh_dates = c(leh_2014_depsess1, leh_2014_depsess2, leh_2014_depsess3,
              leh_2015_depsess1, leh_2015_depsess2)
foreach(el_date = leh_dates) %do% {
  wind_fun <- if(year(el_date) == 2014) daily_wind_2014 else daily_wind_2015
  el <- energy_landscape(leh_col, radius, res, 
                         wind_fun(el_date, 'u'),
                         wind_fun(el_date, 'v'),
                         energy_mod, dur_mod, hi_land)
  writeRaster(el[['out_cost']], 
              sprintf('data/out/EnergyLandscapes/LEH/%s_out.tif',
                      format(el_date, '%Y%m%d')),
              'GTiff')
  writeRaster(el[['in_cost']], 
              sprintf('data/out/EnergyLandscapes/LEH/%s_in.tif',
                      format(el_date, '%Y%m%d')),
              'GTiff')
  writeRaster(el[['rt_cost']], 
              sprintf('data/out/EnergyLandscapes/LEH/%s_rt.tif',
                      format(el_date, '%Y%m%d')),
              'GTiff')
  p <- plot_energy_landscape(el, location, 'out_cost', hi_land)
  ggsave(sprintf('data/out/EnergyLandscapes/LEH/%s_out.png',
                 format(el_date, '%Y%m%d')),
         p, width = 5, height = 5, units = 'in')
  p <- plot_energy_landscape(el, location, 'in_cost', hi_land)
  ggsave(sprintf('data/out/EnergyLandscapes/LEH/%s_in.png',
                 format(el_date, '%Y%m%d')),
         p, width = 5, height = 5, units = 'in')
  p <- plot_energy_landscape(el, location, 'rt_cost', hi_land)
  ggsave(sprintf('data/out/EnergyLandscapes/LEH/%s_rt.png',
                 format(el_date, '%Y%m%d')),
         p, width = 5, height = 5, units = 'in')
}
## all dates, MCB, rasters + plots
mcb_2014_depsess1 <- seq(ymd('2014-06-01', tz = 'US/Hawaii'), 
                         ymd('2014-06-07', tz = 'US/Hawaii'), 
                         by = '1 day')
mcb_2015_depsess1 <- seq(ymd('2015-06-17', tz = 'US/Hawaii'), 
                         ymd('2015-07-06', tz = 'US/Hawaii'), 
                         by = '1 day')
mcb_2015_depsess2 <- seq(ymd('2015-06-29', tz = 'US/Hawaii'), 
                         ymd('2015-07-08', tz = 'US/Hawaii'), 
                         by = '1 day')
mcb_dates = c(mcb_2014_depsess1, mcb_2015_depsess1, mcb_2015_depsess2)
foreach(el_date = mcb_dates) %do% {
  wind_fun <- if(year(el_date) == 2014) daily_wind_2014 else daily_wind_2015
  el <- energy_landscape(mcb_col, radius, res, 
                         wind_fun(el_date, 'u'),
                         wind_fun(el_date, 'v'),
                         energy_mod, dur_mod, hi_land)
  writeRaster(el[['out_cost']], 
              sprintf('data/out/EnergyLandscapes/MCB/%s_out.tif',
                      format(el_date, '%Y%m%d')),
              'GTiff')
  writeRaster(el[['in_cost']], 
              sprintf('data/out/EnergyLandscapes/MCB/%s_in.tif',
                      format(el_date, '%Y%m%d')),
              'GTiff')
  writeRaster(el[['rt_cost']], 
              sprintf('data/out/EnergyLandscapes/MCB/%s_rt.tif',
                      format(el_date, '%Y%m%d')),
              'GTiff')
  p <- plot_energy_landscape(el, location, 'out_cost', hi_land)
  ggsave(sprintf('data/out/EnergyLandscapes/MCB/%s_out.png',
                 format(el_date, '%Y%m%d')),
         p, width = 5, height = 5, units = 'in')
  p <- plot_energy_landscape(el, location, 'in_cost', hi_land)
  ggsave(sprintf('data/out/EnergyLandscapes/MCB/%s_in.png',
                 format(el_date, '%Y%m%d')),
         p, width = 5, height = 5, units = 'in')
  p <- plot_energy_landscape(el, location, 'rt_cost', hi_land)
  ggsave(sprintf('data/out/EnergyLandscapes/MCB/%s_rt.png',
                 format(el_date, '%Y%m%d')),
         p, width = 5, height = 5, units = 'in')
}
