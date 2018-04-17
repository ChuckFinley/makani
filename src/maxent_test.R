library(xtractomatic)
library(dismo)
library(sp)
library(dplyr)
library(ggplot2)

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

# Environmental data
sst_id <- 'jplG1SST'
#chla_id <- 'erdMH1pp1day'
chla_id <- 'mbchla14day'
bathy_id <- 'ETOPO180'
env_x <- bbox(kpc_forage)['x',]
env_y <- bbox(kpc_forage)['y',] 
env_t <- c('2016-06-01', '2016-06-04')
sst_data <- xtracto_3D(sst_id, env_x, env_y, env_t)
chl_data <- xtracto_3D(chla_id, env_x, env_y, env_t)
bathy_data <- xtracto_3D(bathy_id, env_x, env_y, env_t)

# Plot environment
## Utility functions for ggplot compatibility
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

## Plot variables
fortify_xtracto(sst_data) %>%
  filter(t == first(t)) %>%
  ggplot(aes(x = x, y = y, fill = value)) + 
  geom_raster() +
  scale_fill_gradientn(colors = colorRamps::matlab.like(4))

fortify_xtracto(chl_data) %>%
  filter(t == first(t)) %>%
  ggplot(aes(x = x, y = y, fill = log(value))) + 
  geom_raster() +
  scale_fill_gradientn(colors = colorRamps::matlab.like(4))

fortify_xtracto(bathy_data) %>%
  ggplot(aes(x = x, y = y, fill = value)) + 
  geom_raster() +
  scale_fill_gradientn(colors = grDevices::rainbow(6))

ggplot() +
  geom_polygon(aes(long, lat, group = group), 
               data = fortify(hi_land)) +
  geom_polygon(aes(long, lat, group = group),
               data = fortify(kpc_forage),
               fill = NA,
               color = 'black',
               linetype = 'dashed') +
  geom_point(aes(x, y),
             data = kpc_col,
             color = 'red') +
  coord_fixed()

# Presence records
## DB connections
mhi_db <-  src_sqlite('data/MHI_GPS.sqlite')
tidytracks_db <- tbl(mhi_db, 'TidyTracks')

## Tracks
### Get from database
num_time <- as.numeric(lubridate::ymd(env_t, tz = 'UTC'))
rfbo_wk1 <- tidytracks_db %>%
  filter(TimestampUTC >= num_time[1],
         TimestampUTC <= num_time[2],
         DistToCol <= 250e3) %>%
  collect

ggplot(rfbo_wk1, aes(Longitude, 
                     Latitude, 
                     color = factor(TripID))) + 
  geom_path() +
  guides(color = FALSE)

### Sample 5 points per trip
set.seed(1705)
rfbo_sample <- rfbo_wk1 %>%
  filter(PositionLag <= 180) %>%
  group_by(TripID) %>%
  filter(n() >= 20) %>%
  sample_n(5) %>%
  ungroup

ggplot(rfbo_sample, 
       aes(Longitude, 
           Latitude, 
           color = factor(TripID))) + 
  geom_point() +
  guides(color = FALSE)

# Test MaxEnt
## Annotate presence data with SST, chl-a, and bathymetry
dbdate_to_ymd <- function(dbdate) {
  posix_origin <- lubridate::ymd('1970-01-01', tz = 'utc')
  as.POSIXct(dbdate, tz = 'utc', origin = posix_origin) %>%
    format('%Y-%m-%d')
}
rfbo_env <- rfbo_sample %>%
  mutate(sst = xtracto(sst_id, 
                       Longitude,
                       Latitude,
                       dbdate_to_ymd(TimestampUTC))$`mean SST`,
         chla = xtracto(chla_id, 
                        Longitude,
                        Latitude,
                        dbdate_to_ymd(TimestampUTC))$`mean chlorophyll`,
         bathy = xtracto(bathy_id, 
                         Longitude,
                         Latitude,
                         dbdate_to_ymd(TimestampUTC))$`mean altitude`) %>%
  transmute(TripID, 
            LocDate = dbdate_to_ymd(TimestampUTC), 
            Longitude, 
            Latitude, 
            sst, 
            chla, 
            bathy)
readr::write_csv(rfbo_env, 'data/out/Presences/rfbo_sample.csv')

## Annotate presence with energy landscape
sample_el <- function(x, y, t, col) {
  if(!all(col %in% c('KPC', 'LEH', 'MCB'))) 
    stop('invalid colony')
  raster_path <- sprintf('data/out/EnergyLandscapes/%s/%s_rt.tif',
                         col,
                         format(t, '%Y%m%d'))
  if(!file.exists(raster_path))
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
      stop(sprintf('point out of extent',
                   x, y, path, 
                   extent(r)[1], extent(r)[3], extent(r)[2], extent(r)[4]))
    extract(r, cellFromXY(r, c(x, y)))
  })
}

rfbo_env2 <- mutate(rfbo_env,
                    Ert = sample_el, Longitude, Latitude, LocDate, 'KPC',
                    D2C = geosphere::distGeo(cbind(Longitude, Latitude),
                                             kpc_col),
                    UD = )

