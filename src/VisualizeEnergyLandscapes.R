library(tidyverse)
library(sp)
library(raster)

wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
hi_aea_prj <- CRS('+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-163 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
load_el <- function(raster_path) {
  raster(raster_path, crs = hi_aea_prj) %>%
    projectRaster(crs = wgs84_prj)
}

template <- raster('data/out/EnergyLandscapesFR/all/KPC_20160601_out.tif',
                   crs = hi_aea_prj)
hi_land <- rgdal::readOGR('data/coastline/hi_land', 'hi_land') %>%
  spTransform(hi_aea_prj) %>%
  crop(template) %>%
  spTransform(wgs84_prj)
kpc_col <- cbind(-159.3997, 22.22888) %>%
  SpatialPoints(proj4string = wgs84_prj) %>%
  spTransform(hi_aea_prj)
d2c <- distanceFromPoints(template, kpc_col) %>%
  mask(template) %>%
  projectRaster(crs = wgs84_prj)

el_resid <- function(el) {
  el_df <- rasterToPoints(el) %>%
    data.frame
  d2c_df <- rasterToPoints(d2c) %>%
    data.frame
  el_lm <- lm(el_df[,3] ~ d2c_df[,3])
  data.frame(Longitude = el_df[,1],
             Latitude = el_df[,2],
             Residual = resid(el_lm))
}

# OUTBOUND
el_out_r <- load_el('data/out/EnergyLandscapesFR/all/KPC_20160601_out.tif')
el_out_df <- rasterToPoints(el_out_r) %>%
  data.frame
colnames(el_out_df) <- c('Longitude', 'Latitude', 'EL_out')
# Eout
ggplot() +
  geom_raster(aes(Longitude, Latitude, fill = EL_out),
              el_out_df) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               color = '#BBBBBB') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      limits = c(0, 1),
                      name = expression(E[out])) +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/EOUT.png',
       height = 4,
       width = 5,
       units = 'in')

# Out residual
resid_out_df <- el_resid(el_out_r)
ggplot() +
  geom_raster(aes(Longitude, Latitude, fill = Residual),
              resid_out_df) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               color = '#BBBBBB') +
  scale_fill_gradient2(low = 'red',
                       high = 'blue',
                       name = expression(R[out])) +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/ROUT.png',
       height = 4,
       width = 5,
       units = 'in')

# INBOUND
el_in_r <- load_el('data/out/EnergyLandscapesFR/all/KPC_20160601_in.tif')
el_in_df <- rasterToPoints(el_in_r) %>%
  data.frame
colnames(el_in_df) <- c('Longitude', 'Latitude', 'EL_in')
# Eout
ggplot() +
  geom_raster(aes(Longitude, Latitude, fill = EL_in),
              el_in_df) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               color = '#BBBBBB') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      limits = c(0, 1),
                      name = expression(E[`in`])) +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/EIN.png',
       height = 4,
       width = 5,
       units = 'in')

# In residual
resid_in_df <- el_resid(el_in_r)
ggplot() +
  geom_raster(aes(Longitude, Latitude, fill = Residual),
              resid_in_df) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               color = '#BBBBBB') +
  scale_fill_gradient2(low = 'red',
                       high = 'blue',
                       name = expression(R[`in`])) +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/RIN.png',
       height = 4,
       width = 5,
       units = 'in')

# ROUNDTRIP
el_rt_r <- load_el('data/out/EnergyLandscapesFR/all/KPC_20160601_rt.tif')
el_rt_df <- rasterToPoints(el_rt_r) %>%
  data.frame
colnames(el_rt_df) <- c('Longitude', 'Latitude', 'EL_rt')
# Ert
ggplot() +
  geom_raster(aes(Longitude, Latitude, fill = EL_rt),
              el_rt_df) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               color = '#BBBBBB') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      limits = c(0, 1),
                      name = expression(E[rt])) +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/ERT.png',
       height = 4,
       width = 5,
       units = 'in')

# RT residual
resid_rt_df <- el_resid(el_rt_r)
ggplot() +
  geom_raster(aes(Longitude, Latitude, fill = Residual),
              resid_rt_df) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               color = '#BBBBBB') +
  scale_fill_gradient2(low = 'red',
                       high = 'blue',
                       name = expression(R[rt])) +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/RRT.png',
       height = 4,
       width = 5,
       units = 'in')

# With tracks
start_time <- ymd_hms('20160601 00:00:00', tz = 'US/Hawaii')
end_time <- ymd_hms('20160602 00:00:00', tz = 'US/Hawaii')
tracks_df <- read_csv('data/KPC_tracks.csv') %>%
  filter(SubColCode == 'KPC',
         between(TimestampUTC, start_time, end_time))
ggplot() +
  geom_raster(aes(Longitude, Latitude, fill = Residual),
              resid_rt_df) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               color = '#BBBBBB') +
  geom_path(aes(Longitude, Latitude), 
            tracks_df) +
  scale_fill_gradient2(low = 'red', 
                       high = 'blue',
                       name = expression(R[rt])) +
  labs(x = 'Longitude',
       y = 'Latitude') +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/RRTwTracks.png',
       height = 4,
       width = 5,
       units = 'in')

# CRW null model
crw_df <- raster('data/out/CyberBirds/KPC_CRW_UD.tif') %>%
  resample(template) %>%
  mask(template) %>%
  projectRaster(crs = wgs84_prj) %>%
  rasterToPoints %>%
  data.frame
colnames(crw_df) <- c('x', 'y', 'crw')
crw_df <- mutate(crw_df, crw = scales::rescale())
ggplot(crw_df, aes(x, y, fill = crw)) +
  geom_raster() +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               color = '#BBBBBB',
               inherit.aes = FALSE) +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      limits = c(0, 1),
                      name = expression(CRW)) +
  theme_bw() +
  labs(x = '', y = '')
ggsave('analysis/figures/CRW.png',
       height = 4,
       width = 5,
       units = 'in')

# D2C null model
d2c_df <- rasterToPoints(d2c) %>%
  data.frame
colnames(d2c_df) <- c('x', 'y', 'D2C')
d2c_df <- mutate(d2c_df, D2C = D2C / 1000)
ggplot(d2c_df, aes(x, y, fill = D2C)) +
  geom_raster() +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               color = '#BBBBBB',
               inherit.aes = FALSE) +
  scale_fill_gradient(high = '#fee0d2',
                      low = '#de2d26',
                      name = 'D2C (km)') +
  theme_bw() +
  labs(x = '', y = '')
ggsave('analysis/figures/D2C.png',
       height = 4,
       width = 5,
       units = 'in')
