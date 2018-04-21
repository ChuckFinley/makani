library(dplyr)
library(lubridate)
library(raster)
library(foreach)
source('src/energy_landscape.R')

# Movement models
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

# Spatial data
wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
hi_aea_prj <- CRS('+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-163 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
project_col <- function(latlon) {
  data.frame(x = latlon[2], y = latlon[1]) %>%
    SpatialPoints(wgs84_prj) %>%
    spTransform(hi_aea_prj) %>%
    as.data.frame %>%
    as.numeric
}
kpc_wgs84 <- c(22.2, -159.4)
leh_wgs84 <- c(22.0, -160.1)
mcb_wgs84 <- c(21.5, -157.7)
kpc_hi_aea <- project_col(kpc_wgs84)
leh_hi_aea <- project_col(leh_wgs84)
mcb_hi_aea <- project_col(mcb_wgs84)
hi_land <- rgdal::readOGR('data/coastline/hi_land', 'hi_land') %>%
  spTransform(hi_aea_prj)

# Energy landscape extent and resolution
el_radius <- 250e3
## Resolution should be mean transit length (mean transit speed * 20 min)
el_res <- 11.4 * 20 * 60

# Wind data
## Each day's wind is the mean wind between the hours of 7am and 7pm HST
## This is when the birds are most active
daily_wind <- function(t, uv) {
  wind_path <- sprintf('data/wind/WRF_HI/%i_WRF_Hawaii_Regional_Atmospheric_Model_best.ncd.nc',
                       year(t))
  
  if(!file.exists(wind_path))
    stop('wind file doesn\'t exist')
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

# Temporal range
## For KPC, valid dates are May 28 - July 17, 2016 (see ch. 1 fig. 2)
kpc_dates <- seq(ymd('2016-05-28'), ymd('2016-07-17'), by = '1 day')
## There were 3 LEH deployments in 2014 and 2 in 2015. Generate ELs for
## all these dates plus the KPC date range
leh_depsess <- c(seq(ymd('2014-05-13'), ymd('2014-05-18'), by = '1 day'),
                 seq(ymd('2014-06-13'), ymd('2014-06-18'), by = '1 day'),
                 seq(ymd('2014-07-14'), ymd('2014-07-20'), by = '1 day'),
                 seq(ymd('2015-05-26'), ymd('2015-06-04'), by = '1 day'),
                 seq(ymd('2015-06-27'), ymd('2015-07-08'), by = '1 day'))
leh_dates <- c(leh_depsess, kpc_dates)
## There was one MCB deployment in 2014 and two in 2015
mcb_depsess <- c(seq(ymd('2014-06-01'), ymd('2014-06-07'), by = '1 day'),
                 seq(ymd('2015-06-17'), ymd('2015-07-06'), by = '1 day'),
                 seq(ymd('2015-06-29'), ymd('2015-07-09'), by = '1 day'))
mcb_dates <- c(mcb_depsess, kpc_dates)

# Generate energy landscapes
## Foreach colony...
foreach(col_name = list('KPC', 'LEH', 'MCB'),
        col_loc = list(kpc_hi_aea, leh_hi_aea, mcb_hi_aea),
        col_dates = list(kpc_dates, leh_dates, mcb_dates)) %do% {
  ## Foreach date...        
  foreach(d = col_dates) %do% {
    el <- energy_landscape(col_loc, el_radius, el_res, 
                           daily_wind(d, 'u'),
                           daily_wind(d, 'v'),
                           energy_mod, dur_mod, hi_land)
    ## Foreach direction...
    foreach(dir = c('out', 'in', 'rt'),
            dir_name = c('Out', 'In', 'Roundtrip')) %do% {      
      file_name <- function(loc, d, dir) {
        paste(loc, format(d, '%Y%m%d'), dir, sep = '_')
      }
      raster_path1 <- file.path('data/out/EnergyLandscapes/',
                                sprintf('%s/Rasters/%s/',
                                        col_name, dir_name),
                                sprintf('%s.tif', 
                                        file_name(col_name, d, dir)))
      raster_path2 <- file.path('data/out/EnergyLandscapes/all',
                                sprintf('%s.tif', 
                                        file_name(col_name, d, dir)))
      figure_path <- file.path('data/out/EnergyLandscapes/',
                               sprintf('%s/Figures/%s/',
                                       col_name, dir_name),
                               sprintf('%s.png', 
                                       file_name(col_name, d, dir)))
      ## Save raster
      writeRaster(el[[paste(dir, 'cost', sep = '_')]], raster_path1, 'GTiff')
      writeRaster(el[[paste(dir, 'cost', sep = '_')]], raster_path2, 'GTiff')
      ## Save figure 
      p <- plot_energy_landscape(el, col_loc, sprintf('%s_cost', dir), hi_land)
      ggsave(figure_path, p, width = 5.5, height = 5, units = 'in')
    }
  }
}
