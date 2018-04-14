library(dplyr)
library(lubridate)
library(raster)
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

## Movement model
load('data/out/Models/EnergyModels.Rdata')
load('data/out/Models/SgTCModel.RData')
mvmt_mod <- function(a, m, d) {
  t <- a * cos(a)
  spd <- predict(SgTCModel, 
                 re.form = NA,
                 newdata = data.frame(Tailwind = t))
  odba <- predict(oam, 
                  newdata = data.frame(WindAngle = a,
                                       WindSpd = m,
                                       DeployID = 1145),
                  exclude = 's(DeployID)',
                  type = 'response')
  dur <- d / spd
  dur * odba
}

## Barriers
poly <- Polygon(cbind(c(-80, -40, -40, -80), c(-80, -80, -40, -40)))
polys <- Polygons(srl = list(poly), ID = 1)
test_barr <- SpatialPolygons(list(polys))

## Run test
test_el <- energy_landscape(origin, radius, res, wind_u, wind_v, mvmt_mod, 'rt_cost', test_barr)
print(plot_energy_landscape(test_el, origin))

# Real data test
wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
hi_aea_prj <- CRS('+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-163 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
kpc_col <- data.frame(x = -159.40, y = 22.23) %>%
  SpatialPoints(wgs84_prj) %>%
  spTransform(hi_aea_prj) %>%
  as.data.frame %>%
  as.numeric
radius <- 250e3
## Resolution should be mean transit length (mean transit speed * 20 min)
res <- 11.4 * 20 * 60
## Wind data
wind_path <- 'data/wind/WRF_HI/WRF_Hawaii_Regional_Atmospheric_Model_best.ncd.nc'
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
## Wind accessor (index of closest value to x in v)
get_wind <- function(x, y, t, uv) {
  if(!any(between(x, min(lon_wind), max(lon_wind))))
    stop('x out of bounds')
  if(!any(between(y, min(lat_wind), max(lat_wind))))
    stop('y out of bounds')
  if(!any(between(t, min(time_wind), max(time_wind))))
    stop('t out of bounds')
  if(!(uv %in% c('u', 'v')))
    stop('uv must be "u" or "v"')
  i <- sapply(x, function(x) which.min(abs(x - lon_wind)))
  j <- sapply(y, function(y) which.min(abs(y - lat_wind)))
  k <- sapply(t, function(t) which.min(abs(t - time_wind)))
  if(uv == 'u')
    wind_arr <- u_wind
  else
    wind_arr <- v_wind
  mapply(function(i, j, k) wind_arr[i, j, k], i, j, k)
}
## Each day's wind is the mean wind between the hours of 7am and 7pm HST
## This is when the birds are most active
daily_wind <- function(t, uv) {
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

hi_land <- rgdal::readOGR('data/coastline/hi_land', 'hi_land') %>%
  spTransform(hi_aea_prj)

sapply(ymd('2016-05-28', tz = 'US/Hawaii') + days(0:50),
       function(d) {
         el <- energy_landscape(kpc_col, radius, res, 
                                daily_wind(d, 'u'),
                                daily_wind(d, 'v'),
                                mvmt_mod, 'rt_cost', hi_land)
         writeRaster(el, 
                     sprintf('data/out/EnergyLandscapes/KPC/%s.tif',
                             format(d, '%Y%m%d')),
                     'GTiff')
         p <- plot_energy_landscape(el, kpc_col) +
           geom_contour(aes(z = val))
         ggsave(sprintf('data/out/EnergyLandscapes/KPC/%s.png',
                        format(d, '%Y%m%d')),
                p, width = 5, height = 5, units = 'in')
       })
el_0528 <- energy_landscape(kpc_col, radius, res, 
                            daily_wind(ymd('2016-05-28', tz = 'US/Hawaii'), 'u'),
                            daily_wind(ymd('2016-05-28', tz = 'US/Hawaii'), 'v'),
                            mvmt_mod, 'rt_cost', hi_land)
plot_energy_landscape(el_0528, kpc_col)
beepr::beep()
