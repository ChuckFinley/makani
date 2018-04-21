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
hi_land <- rgdal::readOGR('data/coastline/hi_land', 'hi_land') %>%
  spTransform(wgs84_prj)

# Presence records
## DB connections
mhi_db <-  src_sqlite('data/MHI_GPS.sqlite')
tidytracks_db <- tbl(mhi_db, 'TidyTracks')

## Tracks
### Get from database
t_rng <- as.POSIXct(seq(ymd('2016-05-29'), ymd('2016-07-16'), by = '1 day'))
t1 <- as.numeric(first(t_rng))
t2 <- as.numeric(last(t_rng))
posix_origin <- ymd('1970-01-01', tz = 'UTC')
rfbo_tracks <- tidytracks_db %>%
  filter(TimestampUTC >= t1,
         TimestampUTC <= t2,
         PositionLag < 180) %>%
  collect %>%
  mutate(TimestampUTC = as.POSIXct(TimestampUTC, 
                                   origin = posix_origin, 
                                   tz = 'UTC'), 
         day = format(TimestampUTC, '%Y%m%d'))

plot_daily_tracks <- function(tracks) {
  el <- raster(sprintf('data/out/EnergyLandscapes/all/KPC_%s_rt.tif',
                       first(tracks$day)),
               crs = hi_aea_prj) %>%
    projectRaster(crs = wgs84_prj)
  el_df <- el%>%
    rasterToPoints() %>%
    data.frame
  colnames(el_df) <- c('Longitude', 'Latitude', 'EL')
  
  crop_extent <- union(extent(el),
                       extent(SpatialPoints(dplyr::select(tracks, 
                                                          Longitude, 
                                                          Latitude))))
  cropped_hi <- crop(hi_land, crop_extent)
  
  ggplot() +
    geom_raster(aes(Longitude, Latitude, fill = EL),
                el_df) +
    geom_polygon(aes(long, lat, group = group),
                 fortify(cropped_hi)) +
    geom_path(aes(long, lat),
              fortify(kpc_forage),
              linetype = 'dashed') +
    geom_path(aes(Longitude, Latitude),
              tracks) +
    theme_bw() +
    coord_fixed() +
    scale_fill_gradientn(name = expression(italic(E[rt])),
                         colors = colorRamps::matlab.like(4),
                         limits = 0:1)
}

overlay_plots <- rfbo_tracks %>%
  group_by(day) %>%
  do(p = plot_daily_tracks(.))
