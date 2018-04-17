library(sp)
library(dplyr)
library(ggplot2)
library(lubridate)

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

## Xtractomatic data ids
sst_id <- 'jplG1SST'
chla_id <- 'mbchla14day'
bathy_id <- 'ETOPO180'

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

### Sample 5 points per trip
set.seed(1705)
rfbo_sample <- rfbo_tracks %>%
  filter(PositionLag <= 180) %>%
  group_by(TripID) %>%
  filter(n() >= 20) %>%
  sample_n(5) %>%
  ungroup

# Annotation
## Annotate presences with environmental data
dbdate_to_ymd <- function(dbdate) {
  posix_origin <- ymd('1970-01-01', tz = 'utc')
  as.POSIXct(dbdate, tz = 'utc', origin = posix_origin) %>%
    format('%Y-%m-%d')
}
extract_sst <- function(x, y, t) {
  xtractomatic::xtracto(sst_id, x, y, dbdate_to_ymd(t))$`mean SST`
}
extract_chl <- function(x, y, t) {
  xtractomatic::xtracto(chla_id, x, y, dbdate_to_ymd(t))$`mean chlorophyll`
}
extract_bathy <- function(x, y) {
  xtractomatic::xtracto(sst_id, x, y)$`mean altitude`
}
rfbo_env <- rfbo_sample %>%
  mutate(sst = extract_sst(Longitude, Latitude, TimestampUTC),
         chla = extract_chl(Longitude, Latitude, TimestampUTC),
         bathy = extract_bathy(Longitude, Latitude)) %>%
  transmute(Species = 'RFBO',
            Longitude, 
            Latitude, 
            sst, 
            chla, 
            bathy) %>%
  na.omit

## Annotate presences with accessibility data
extract_el <- function(x, y, t, col) {
  if(!all(col %in% c('KPC', 'LEH', 'MCB'))) 
    stop('invalid colony')
  t2 <- ymd(t)
  raster_path <- sprintf('data/out/EnergyLandscapes/%s/%s_rt.tif',
                         col,
                         format(t2, '%Y%m%d'))
  if(!any(file.exists(raster_path)))
    stop('no raster exists')
  wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
  hi_aea_prj <- CRS('+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-163 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
  rescale <- function(r) {
    r_min = cellStats(r, "min")
    r_max = cellStats(r, "max")
    (r - r_min) / (r_max - r_min)
  }
  mapply(FUN = function(x, y, path) {
    r <- raster(path, crs = hi_aea_prj) %>%
      projectRaster(crs = wgs84_prj) %>%
      rescale
    if(!between(x, extent(r)[1], extent(r)[2]) ||
       !between(y, extent(r)[3], extent(r)[4]))
      stop('point out of extent')
    raster::extract(r, cellFromXY(r, c(x, y)))
  },
  x, y, raster_path)
}
extract_ud <- function(x, y, col) {
  if(!all(col %in% c('KPC', 'LEH', 'MCB'))) 
    stop('invalid colony')
  raster_path <- sprintf('data/out/CyberBirds/%s_CRW_UD.tif', col)
  if(!any(file.exists(raster_path)))
    stop('no raster exists')
  wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
  hi_aea_prj <- CRS('+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-163 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
  mapply(FUN = function(x, y, path) {
    r <- raster(path) %>%
      projectRaster(crs = wgs84_prj)
    if(!between(x, extent(r)[1], extent(r)[2]) ||
       !between(y, extent(r)[3], extent(r)[4]))
      stop('point out of extent')
    raster::extract(r, cellFromXY(r, c(x, y)))
  },
  x, y, raster_path)
}
rfbo_acc <- mutate(rfbo_env,
                   Ert = extract_el(Longitude, Latitude, LocDate, 'KPC'),
                   D2C = geosphere::distGeo(cbind(Longitude, Latitude),
                                            kpc_col),
                   UD = sample_ud(Longitude, Latitude, 'KPC'))

                   