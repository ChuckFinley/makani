library(ggplot2)
library(dplyr)
library(raster)

wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
hi_aea_prj <- CRS('+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-163 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
hi_land <- rgdal::readOGR('data/coastline/hi_land', 'hi_land')
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

plot_wind <- function(d) {
  u <- daily_wind(ymd('2016-06-07'), 'u')
  v <- daily_wind(ymd('2016-06-07'), 'v')
  template <- u
  res(template) <- 40e3
  u2 <- resample(u, template)
  v2 <- resample(v, template)
  fortify_raster <- function(r) {
    data.frame(i = seq(ncell(r))) %>%
      mutate(x = xFromCell(r, i),
             y = yFromCell(r, i),
             val = getValues(r))
  }
  uv <- left_join(fortify_raster(u2), fortify_raster(v2), by = 'i') %>% 
    transmute(i, x = x.x, y = y.x, u = val.x, v = val.y)
  
  ggplot(uv, aes(x, y)) + 
    geom_segment(aes(xend = x + 5e3*u, yend = y + 5e3*v), 
                 arrow = arrow(length = unit(0.1,"cm"))) + 
    geom_polygon(aes(long, lat, group = group), 
                 fortify(spTransform(hi_land, hi_aea_prj)), inherit.aes = FALSE)
}


