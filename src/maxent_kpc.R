library(sp)
library(dplyr)
library(ggplot2)
library(lubridate)
library(raster)

# Spatial data
## Colony location and foraging range
wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
hi_aea_prj <- CRS('+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-163 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
kpc_col <- data.frame(x = -159.40, y = 22.23)
kpc_forage <- SpatialPoints(kpc_col,
                            proj4string = wgs84_prj) %>%
  spTransform(hi_aea_prj) %>%
  rgeos::gBuffer(width = 250e3) %>%
  spTransform(wgs84_prj)

## Land
hi_land <- rgdal::readOGR('data/coastline/hi_land', 'hi_land')

# Presence records
## DB connections
mhi_db <-  src_sqlite('data/MHI_GPS.sqlite')
tidytracks_db <- tbl(mhi_db, 'TidyTracks')

## Tracks
### Get from database
t_rng <- as.POSIXct(seq(ymd('2016-05-29'), ymd('2016-07-16'), by = '1 day'))
t1 <- as.numeric(first(t_rng))
t2 <- as.numeric(last(t_rng))
rfbo_tracks <- tidytracks_db %>%
  filter(TimestampUTC >= t1,
         TimestampUTC <= t2,
         DistToCol <= 250e3) %>%
  collect

### Sample 5 points per bird per day
set.seed(1705)
dbdate_to_posix <- function(t)
  as.POSIXct(t, origin = ymd('1970-01-01', tz = 'UTC'), tz = 'UTC')
rfbo_sample <- rfbo_tracks %>%
  filter(PositionLag <= 180) %>%
  mutate(TimestampUTC = dbdate_to_posix(TimestampUTC)) %>%
  group_by(DeployID, as.Date(TimestampUTC)) %>%
  filter(n() >= 20) %>%
  sample_n(5) %>%
  ungroup

# Background records
set.seed(1554)
background_kpc <- spsample(kpc_forage, 1e4, type = 'random') %>%
  as.data.frame %>%
  transmute(Species = 'RFBO',
            Longitude = x,
            Latitude = y,
            LocDate = sample(t_rng, 1e4, replace = TRUE))

# Annotation
## Annotate presences with environmental data
dbdate_to_posixct <- function(dbdate) {
  posix_origin <- ymd('1970-01-01', tz = 'utc')
  as.POSIXct(dbdate, tz = 'utc', origin = posix_origin)
}
sst_nc <- RNetCDF::open.nc('data/sst/2016sst.nc')
time_dim <- 'time'     # seconds since 1970-01-01 00:00:00.000 UTC
lat_dim <- 'latitude'  # degrees north
lon_dim <- 'longitude' # degrees east
sst_var <- 'sst'
anom_var <- 'anom'
time_sst <- RNetCDF::var.get.nc(sst_nc, time_dim)
t_max <- max(time_sst)
t_min <- min(time_sst)
lat_sst <- RNetCDF::var.get.nc(sst_nc, lat_dim)
y_max <- max(lat_sst)
y_min <- min(lat_sst)
lon_sst <- RNetCDF::var.get.nc(sst_nc, lon_dim)
x_max <- max(lon_sst)
x_min <- min(lon_sst)
sst_arr <- RNetCDF::var.get.nc(sst_nc, sst_var)
anom_arr <- RNetCDF::var.get.nc(sst_nc, anom_var)
extract_sst <- function(x, y, t) {
  t <- as.numeric(t)
  x <- (x + 360) %% 360
  t_i <- round((length(time_sst) - 1) * (t - t_min) / (t_max - t_min)) + 1
  x_i <- round((length(lon_sst) - 1) * (x - x_min) / (x_max - x_min)) + 1
  y_i <- round((length(lat_sst) - 1) * (y - y_min) / (y_max - y_min)) + 1
  sst_arr[cbind(x_i, y_i, t_i)]
}
extract_anom <- function(x, y, t) {
  t <- as.numeric(t)
  x <- (x + 360) %% 360
  t_i <- round((length(time_sst) - 1) * (t - t_min) / (t_max - t_min)) + 1
  x_i <- round((length(lon_sst) - 1) * (x - x_min) / (x_max - x_min)) + 1
  y_i <- round((length(lat_sst) - 1) * (y - y_min) / (y_max - y_min)) + 1
  anom_arr[cbind(x_i, y_i, t_i)]
}
bathy_r <- raster('data/bathy/etopo1.tif')
extract_bathy <- function(x, y) {
  raster::extract(bathy_r, cbind(x, y))
}
rfbo_env <- rfbo_sample %>%
  mutate(sst = extract_sst(Longitude, Latitude, TimestampUTC),
         anom = extract_anom(Longitude, Latitude, TimestampUTC),
         bathy = extract_bathy(Longitude, Latitude)) %>%
  transmute(Species = 'RFBO',
            Longitude, 
            Latitude, 
            LocDate = as.Date(TimestampUTC),
            sst, 
            anom, 
            bathy) %>%
  na.omit

## Annotate presences with accessibility data
extract_el <- function(x, y, t, col) {
  if(!all(col %in% c('KPC', 'LEH', 'MCB'))) 
    stop('invalid colony')
  wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
  hi_aea_prj <- CRS('+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-163 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
  extract_el2 <- function(x, y, t, col) {
    sprintf('data/out/EnergyLandscapesFR/all/%s_%s_rt.tif',
                 first(col),
                 format(first(t), '%Y%m%d')) %>%
      raster(crs = hi_aea_prj) %>%
      projectRaster(crs = wgs84_prj) %>%
      raster::extract(cbind(x, y))
  }
  result <- data.frame(x = x, y = y, t = t, col = col) %>%
    mutate(orig_order = row_number()) %>%
    group_by(col, t) %>%
    do(data.frame(orig_order = .$orig_order,
                  el = extract_el2(.$x, .$y, .$t, .$col))) %>%
    ungroup %>%
    arrange(orig_order)
  result$el
}
extract_ud <- function(x, y, col) {
  if(!all(col %in% c('KPC', 'LEH', 'MCB'))) 
    stop('invalid colony')
  wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
  hi_aea_prj <- CRS('+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-163 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
  extract_ud2 <- function(x, y, col) {
    sprintf('data/out/CyberBirds/%s_CRW_UD.tif', first(col)) %>%
      raster %>%
      projectRaster(crs = wgs84_prj) %>%
      raster::extract(cbind(x, y))
  }
  result <- data.frame(x = x, y = y, col = col) %>%
    mutate(orig_order = row_number()) %>%
    group_by(col) %>%
    do(data.frame(orig_order = .$orig_order,
                  ud = extract_ud2(.$x, .$y, .$col))) %>%
    ungroup %>%
    arrange(orig_order)
  result$ud
}
rfbo_acc <- mutate(rfbo_env,
                   Ert = extract_el(Longitude, Latitude, LocDate, 'KPC'),
                   D2C = geosphere::distGeo(cbind(Longitude, Latitude),
                                            kpc_col),
                   UD = extract_ud(Longitude, Latitude, 'KPC')) %>%
  na.omit
readr::write_csv(rfbo_acc, 'data/out/Presences/rfbo_accessible.csv')

## Annotate background with environmental data  
background_env <- background_kpc %>%
  mutate(sst = extract_sst(Longitude, Latitude, LocDate),
         anom = extract_anom(Longitude, Latitude, LocDate),
         bathy = extract_bathy(Longitude, Latitude)) %>%
  transmute(Species = 'RFBO',
            Longitude, 
            Latitude, 
            LocDate,
            sst, 
            anom, 
            bathy) %>%
  na.omit 
background_acc <- mutate(background_env,
                         Ert = extract_el(Longitude, Latitude, LocDate, 'KPC'),
                         D2C = geosphere::distGeo(cbind(Longitude, Latitude),
                                                  kpc_col),
                         UD = extract_ud(Longitude, Latitude, 'KPC')) %>%
  na.omit
readr::write_csv(background_acc, 'data/out/Presences/kpc_background.csv')
