library(raster)
library(tidyverse)
library(lubridate)

el <- raster('data/out/EnergyLandscapesFR/all/KPC_20160610_rt.tif')
wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
hi_aea_prj <- CRS('+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-163 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
hi_land <- rgdal::readOGR('data/coastline/hi_land', 'hi_land') %>%
  spTransform(hi_aea_prj)
project_col <- function(lat, lon) {
  data.frame(x = lon, y = lat) %>%
    SpatialPoints(wgs84_prj) %>%
    spTransform(hi_aea_prj) %>%
    as.data.frame %>%
    as.numeric
}
kpc_col <- project_col(22.2, -159.4)
kpc_extent <- rgeos::gBuffer(SpatialPoints(cbind(kpc_col[1], kpc_col[2]), proj4string = hi_aea_prj), width = 250e3)
origin <- kpc_col
start_time <- ymd_hms('20160610 00:00:00', tz = 'US/Hawaii')
end_time <- ymd_hms('20160611 00:00:00', tz = 'US/Hawaii')
tracks_df_wgs <- read_csv('data/KPC_tracks.csv') %>%
  filter(SubColCode == 'KPC',
         between(TimestampUTC, start_time, end_time))
tracks_sp <- SpatialPoints(tracks_df_wgs[,c('Longitude', 'Latitude')], 
                           proj4string = wgs84_prj) %>%
  spTransform(hi_aea_prj)
tracks_df_hi <- data.frame(tracks_sp)

dist_r <- distanceFromPoints(el, origin) %>%
  mask(el)
el_p <- rasterToPoints(el)
dist_p <- rasterToPoints(dist_r)
el_dist <- cbind(el_p, dist_p[,3]) %>% 
  data.frame
colnames(el_dist) <- c('x', 'y', 'el', 'dist')
el_lm <- lm(el ~ dist, el_dist)
colnames(tracks_df_hi) <- c('x', 'y')
el_df <- mutate(el_dist, el_resid = resid(el_lm))

sample_el_resid <- function(x, y, df) {
  dplyr::select(df, x, y, el_resid) %>%
    as.matrix %>%
    rasterFromXYZ(res = res(el),
                  crs = hi_aea_prj) %>%
    raster::extract(cbind(x, y))
}

presence <- tracks_df_hi %>%
  sample_frac(0.2) %>%
  mutate(el_resid = sample_el_resid(x, y, el_df),
         type = 'Presence')
background <- sample_frac(el_df, 0.2) %>%
  mutate(type = 'Background')
rbind(select(presence, el_resid, type),
      select(background, el_resid, type)) %>%
  ggplot(aes(el_resid, color = type)) +
  geom_density()

ggplot() +
  geom_polygon(aes(long, lat, group = group), fortify(crop(hi_land, kpc_extent))) +
  geom_raster(aes(x, y, fill = el_resid), el_df) +
  geom_path(aes(x, y), tracks_df_hi) +
  scale_fill_gradient2(low = 'red', high = 'blue') +
  labs(x = '',
       y = '',
       title = 'Energy Landscape and Tracks for 2016-06-10',
       fill = 'Residual') +
  coord_fixed() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave('analysis/figures/TracksResiduals.png', height = 5, width = 5.3, units = 'in')

rasterToPoints(el) %>% 
  data.frame %>%
  ggplot() +
  geom_raster(aes(x, y, fill = KPC_20160610_rt)) +
  scale_fill_gradientn(colors = colorRamps::matlab.like(4))
