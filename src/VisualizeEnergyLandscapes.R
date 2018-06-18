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

el_resid <- function(el, d2c) {
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
  geom_contour(aes(Longitude, Latitude, z = EL_out),
               el_out_df,
               color = '#666666',
               linetype = '88888888') +
  geom_point(aes(x, y),
             data.frame(x = -159.3997, 
                        y = 22.22888),
             shape = 16,
             size = 3,
             color = 'black') +
  geom_point(aes(x, y),
             data.frame(x = -159.3997, 
                        y = 22.22888),
             shape = 1,
             size = 3,
             color = 'white') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      limits = c(0, 1),
                      name = expression(E[out])) +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/EOUT.png',
       height = 3,
       width = 3.75,
       units = 'in')

# Out residual
resid_out_df <- el_resid(el_out_r, d2c)
ggplot() +
  geom_raster(aes(Longitude, Latitude, fill = Residual),
              resid_out_df) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               color = '#BBBBBB') +
  geom_point(aes(x, y),
             data.frame(x = -159.3997, 
                        y = 22.22888),
             shape = 16,
             size = 3,
             color = 'black') +
  geom_point(aes(x, y),
             data.frame(x = -159.3997, 
                        y = 22.22888),
             shape = 1,
             size = 3,
             color = 'white') +
  scale_fill_gradient2(low = 'red',
                       high = 'blue',
                       name = expression(R[out])) +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/ROUT.png',
       height = 3,
       width = 3.75,
       units = 'in')

# INBOUND
el_in_r <- load_el('data/out/EnergyLandscapesFR/all/KPC_20160601_in.tif')
el_in_df <- rasterToPoints(el_in_r) %>%
  data.frame
colnames(el_in_df) <- c('Longitude', 'Latitude', 'EL_in')
# Ein
ggplot() +
  geom_raster(aes(Longitude, Latitude, fill = EL_in),
              el_in_df) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               color = '#BBBBBB') +
  geom_contour(aes(Longitude, Latitude, z = EL_in),
               el_in_df,
               color = '#666666',
               linetype = '88888888') +
  geom_point(aes(x, y),
             data.frame(x = -159.3997, 
                        y = 22.22888),
             shape = 16,
             size = 3,
             color = 'black') +
  geom_point(aes(x, y),
             data.frame(x = -159.3997, 
                        y = 22.22888),
             shape = 1,
             size = 3,
             color = 'white') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      limits = c(0, 1),
                      name = expression(E[`in`])) +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/EIN.png',
       height = 3,
       width = 3.75,
       units = 'in')

# In residual
resid_in_df <- el_resid(el_in_r, d2c)
ggplot() +
  geom_raster(aes(Longitude, Latitude, fill = Residual),
              resid_in_df) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               color = '#BBBBBB') +
  geom_point(aes(x, y),
             data.frame(x = -159.3997, 
                        y = 22.22888),
             shape = 16,
             size = 3,
             color = 'black') +
  geom_point(aes(x, y),
             data.frame(x = -159.3997, 
                        y = 22.22888),
             shape = 1,
             size = 3,
             color = 'white') +
  scale_fill_gradient2(low = 'red',
                       high = 'blue',
                       name = expression(R[`in`])) +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/RIN.png',
       height = 3,
       width = 3.75,
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
  geom_contour(aes(Longitude, Latitude, z = EL_rt),
               el_rt_df,
               color = '#666666',
               linetype = '88888888') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      limits = c(0, 1),
                      name = expression(E[rt])) +
  geom_point(aes(x, y),
             data.frame(x = -159.3997, 
                        y = 22.22888),
             shape = 16,
             size = 3,
             color = 'black') +
  geom_point(aes(x, y),
             data.frame(x = -159.3997, 
                        y = 22.22888),
             shape = 1,
             size = 3,
             color = 'white') +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/ERT.png',
       height = 3,
       width = 3.75,
       units = 'in')

# RT residual
resid_rt_df <- el_resid(el_rt_r, d2c)
ggplot() +
  geom_raster(aes(Longitude, Latitude, fill = Residual),
              resid_rt_df) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               color = '#BBBBBB') +
  scale_fill_gradient2(low = 'red',
                       high = 'blue',
                       name = expression(R[rt])) +
  geom_point(aes(x, y),
             data.frame(x = -159.3997, 
                        y = 22.22888),
             shape = 16,
             size = 3,
             color = 'black') +
  geom_point(aes(x, y),
             data.frame(x = -159.3997, 
                        y = 22.22888),
             shape = 1,
             size = 3,
             color = 'white') +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/RRT.png',
       height = 3,
       width = 3.75,
       units = 'in')

# RT + residual for LEH, MCB
leh_template <- raster('data/out/EnergyLandscapesFR/all/LEH_20160601_out.tif',
                   crs = hi_aea_prj)
leh_land <- rgdal::readOGR('data/coastline/hi_land', 'hi_land') %>%
  spTransform(hi_aea_prj) %>%
  crop(leh_template) %>%
  spTransform(wgs84_prj)
leh_col <- cbind(-160.1, 22.0) %>%
  SpatialPoints(proj4string = wgs84_prj) %>%
  spTransform(hi_aea_prj)
leh_d2c <- distanceFromPoints(leh_template, leh_col) %>%
  mask(leh_template) %>%
  projectRaster(crs = wgs84_prj)

leh_rt_r <- load_el('data/out/EnergyLandscapesFR/all/LEH_20160601_rt.tif')
leh_rt_df <- rasterToPoints(leh_rt_r) %>%
  data.frame
colnames(leh_rt_df) <- c('Longitude', 'Latitude', 'EL_rt')
# LEH Ert
ggplot() +
  geom_raster(aes(Longitude, Latitude, fill = EL_rt),
              leh_rt_df) +
  geom_polygon(aes(long, lat, group = group),
               fortify(leh_land),
               color = '#BBBBBB') +
  geom_contour(aes(Longitude, Latitude, z = EL_rt),
               leh_rt_df,
               color = '#666666',
               linetype = '88888888') +
  geom_point(aes(x, y),
             data.frame(x = -160.1, 
                        y = 22.0),
             shape = 16,
             size = 3,
             color = 'black') +
  geom_point(aes(x, y),
             data.frame(x = -160.1, 
                        y = 22.0),
             shape = 1,
             size = 3,
             color = 'white') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      limits = c(0, 1),
                      name = expression(E[rt])) +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/ERT_LEH.png',
       height = 3,
       width = 3.75,
       units = 'in')

# RT residual
leh_resid_df <- el_resid(leh_rt_r, leh_d2c)
ggplot() +
  geom_raster(aes(Longitude, Latitude, fill = Residual),
              leh_resid_df) +
  geom_polygon(aes(long, lat, group = group),
               fortify(leh_land),
               color = '#BBBBBB') +
  geom_point(aes(x, y),
             data.frame(x = -160.1, 
                        y = 22.0),
             shape = 16,
             size = 3,
             color = 'black') +
  geom_point(aes(x, y),
             data.frame(x = -160.1, 
                        y = 22.0),
             shape = 1,
             size = 3,
             color = 'white') +
  scale_fill_gradient2(low = 'red',
                       high = 'blue',
                       name = expression(R[rt])) +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/RRT_LEH.png',
       height = 3,
       width = 3.75,
       units = 'in')

mcb_template <- raster('data/out/EnergyLandscapesFR/all/MCB_20160601_out.tif',
                       crs = hi_aea_prj)
mcb_land <- rgdal::readOGR('data/coastline/hi_land', 'hi_land') %>%
  spTransform(hi_aea_prj) %>%
  crop(mcb_template) %>%
  spTransform(wgs84_prj)
mcb_col <- cbind(-157.7, 21.5) %>%
  SpatialPoints(proj4string = wgs84_prj) %>%
  spTransform(hi_aea_prj)
mcb_d2c <- distanceFromPoints(mcb_template, mcb_col) %>%
  mask(mcb_template) %>%
  projectRaster(crs = wgs84_prj)

mcb_rt_r <- load_el('data/out/EnergyLandscapesFR/all/MCB_20160601_rt.tif')
mcb_rt_df <- rasterToPoints(mcb_rt_r) %>%
  data.frame
colnames(mcb_rt_df) <- c('Longitude', 'Latitude', 'EL_rt')
# MCB Ert
ggplot() +
  geom_raster(aes(Longitude, Latitude, fill = EL_rt),
              mcb_rt_df) +
  geom_polygon(aes(long, lat, group = group),
               fortify(mcb_land),
               color = '#BBBBBB') +
  geom_contour(aes(Longitude, Latitude, z = EL_rt),
               mcb_rt_df,
               color = '#666666',
               linetype = '88888888') +
  geom_point(aes(x, y),
             data.frame(x = -157.7, 
                        y = 21.5),
             shape = 16,
             size = 3,
             color = 'black') +
  geom_point(aes(x, y),
             data.frame(x = -157.7, 
                        y = 21.5),
             shape = 1,
             size = 3,
             color = 'white') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      limits = c(0, 1),
                      name = expression(E[rt])) +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/ERT_MCB.png',
       height = 3,
       width = 3.75,
       units = 'in')

# RT residual
mcb_resid_df <- el_resid(mcb_rt_r, mcb_d2c)
ggplot() +
  geom_raster(aes(Longitude, Latitude, fill = Residual),
              mcb_resid_df) +
  geom_polygon(aes(long, lat, group = group),
               fortify(mcb_land),
               color = '#BBBBBB') +
  geom_point(aes(x, y),
             data.frame(x = -157.7, 
                        y = 21.5),
             shape = 16,
             size = 3,
             color = 'black') +
  geom_point(aes(x, y),
             data.frame(x = -157.7, 
                        y = 21.5),
             shape = 1,
             size = 3,
             color = 'white') +
  scale_fill_gradient2(low = 'red',
                       high = 'blue',
                       name = expression(R[rt])) +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/RRT_MCB.png',
       height = 3,
       width = 3.75,
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
                       name = expression(R[rt]),
                       lim = c(-0.17, 0.17)) +
  labs(x = 'Longitude',
       y = 'Latitude') +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/RRTwTracks0601.png',
       height = 3,
       width = 3.75,
       units = 'in')
start_time <- ymd_hms('20160615 00:00:00', tz = 'US/Hawaii')
end_time <- ymd_hms('20160616 00:00:00', tz = 'US/Hawaii')
tracks_df <- read_csv('data/KPC_tracks.csv') %>%
  filter(SubColCode == 'KPC',
         between(TimestampUTC, start_time, end_time))
resid_rt_df2 <- el_resid(load_el('data/out/EnergyLandscapesFR/all/KPC_20160615_rt.tif'), 
                         d2c)
ggplot() +
  geom_raster(aes(Longitude, Latitude, fill = Residual),
              resid_rt_df2) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               color = '#BBBBBB') +
  geom_path(aes(Longitude, Latitude), 
            tracks_df) +
  scale_fill_gradient2(low = 'red', 
                       high = 'blue',
                       name = expression(R[rt]),
                       lim = c(-0.17, 0.17)) +
  labs(x = 'Longitude',
       y = 'Latitude') +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/RRTwTracks0615.png',
       height = 3,
       width = 3.75,
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

# Ert variability
fortify_raster <- function(raster) {
  result <- raster %>%
    as('SpatialPixelsDataFrame') %>%
    as.data.frame
  colnames(result) <- c('value', 'x', 'y')
  result
}

## KPC
kpc_stack <- dir('data/out/EnergyLandscapesFR/all/', 
    pattern = 'KPC_2016[0-9]{4}_rt.tif',
    full.names = TRUE) %>%
  stack
crs(kpc_stack) <- hi_aea_prj
kpc_stack <- projectRaster(kpc_stack, crs = wgs84_prj)
kpc_mean <- calc(kpc_stack, fun = mean)
kpc_sd <- calc(kpc_stack, fun = sd)
ggplot() +
  geom_raster(aes(x, y, fill = value),
              fortify_raster(kpc_mean))  +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               color = '#BBBBBB') +
  geom_contour(aes(x, y, z = value),
               fortify_raster(kpc_mean),
               color = '#666666',
               linetype = '88888888') +
  geom_point(aes(x, y),
             data.frame(x = -159.3997, 
                        y = 22.22888),
             shape = 16,
             size = 3,
             color = 'black') +
  geom_point(aes(x, y),
             data.frame(x = -159.3997, 
                        y = 22.22888),
             shape = 1,
             size = 3,
             color = 'white') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      limits = c(0, 1),
                      name = expression(mean(E[rt]))) +
  theme_bw() +
  theme(legend.title = element_text(size = 9)) +
  coord_fixed() +
  labs(x = 'Longitude',
       y = 'Latitude')
ggsave('analysis/figures/KPC_mean.png',
       height = 3,
       width = 3.75,
       units = 'in')
ggplot() +
  geom_raster(aes(x, y, fill = value),
              fortify_raster(kpc_sd))  +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               color = '#BBBBBB') +
  geom_contour(aes(x, y, z = value),
               fortify_raster(kpc_sd),
               color = '#666666',
               linetype = '88888888') +
  geom_point(aes(x, y),
             data.frame(x = -159.3997, 
                        y = 22.22888),
             shape = 16,
             size = 3,
             color = 'black') +
  geom_point(aes(x, y),
             data.frame(x = -159.3997, 
                        y = 22.22888),
             shape = 1,
             size = 3,
             color = 'white') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      lim = c(0, 0.06),
                      breaks = seq(0, 0.06, by = 0.02),
                      name = expression(sd(E[rt]))) +
  theme_bw() +
  coord_fixed() +
  labs(x = 'Longitude',
       y = 'Latitude')
ggsave('analysis/figures/KPC_sd.png',
       height = 3,
       width = 3.75,
       units = 'in')
## MCB
mcb_stack <- dir('data/out/EnergyLandscapesFR/all/', 
                 pattern = 'MCB_2016[0-9]{4}_rt.tif',
                 full.names = TRUE) %>%
  stack
crs(mcb_stack) <- hi_aea_prj
mcb_stack <- projectRaster(mcb_stack, crs = wgs84_prj)
mcb_mean <- calc(mcb_stack, fun = mean)
mcb_sd <- calc(mcb_stack, fun = sd)
ggplot() +
  geom_raster(aes(x, y, fill = value),
              fortify_raster(mcb_mean))  +
  geom_polygon(aes(long, lat, group = group),
               fortify(mcb_land),
               color = '#BBBBBB') +
  geom_contour(aes(x, y, z = value),
               fortify_raster(mcb_mean),
               color = '#666666',
               linetype = '88888888') +
  geom_point(aes(x, y),
             data.frame(x = -157.7, 
                        y = 21.5),
             shape = 16,
             size = 3,
             color = 'black') +
  geom_point(aes(x, y),
             data.frame(x = -157.7, 
                        y = 21.5),
             shape = 1,
             size = 3,
             color = 'white') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      limits = c(0, 1),
                      name = expression(mean(E[rt]))) +
  theme_bw() +
  theme(legend.title = element_text(size = 9)) +
  coord_fixed() +
  labs(x = 'Longitude',
       y = 'Latitude')
ggsave('analysis/figures/MCB_mean.png',
       height = 3,
       width = 3.75,
       units = 'in')
ggplot() +
  geom_raster(aes(x, y, fill = value),
              fortify_raster(mcb_sd))  +
  geom_polygon(aes(long, lat, group = group),
               fortify(mcb_land),
               color = '#BBBBBB') +
  geom_contour(aes(x, y, z = value),
               fortify_raster(mcb_sd),
               color = '#666666',
               linetype = '88888888') +
  geom_point(aes(x, y),
             data.frame(x = -157.7, 
                        y = 21.5),
             shape = 16,
             size = 3,
             color = 'black') +
  geom_point(aes(x, y),
             data.frame(x = -157.7, 
                        y = 21.5),
             shape = 1,
             size = 3,
             color = 'white') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      lim = c(0, 0.06),
                      breaks = seq(0, 0.06, by = 0.02),
                      name = expression(sd(E[rt]))) +
  theme_bw() +
  coord_fixed() +
  labs(x = 'Longitude',
       y = 'Latitude')
ggsave('analysis/figures/MCB_sd.png',
       height = 3,
       width = 3.75,
       units = 'in')

## LEH
leh_stack <- dir('data/out/EnergyLandscapesFR/all/', 
                 pattern = 'LEH_2016[0-9]{4}_rt.tif',
                 full.names = TRUE) %>%
  stack
crs(leh_stack) <- hi_aea_prj
leh_stack <- projectRaster(leh_stack, crs = wgs84_prj)
leh_mean <- calc(leh_stack, fun = mean)
leh_sd <- calc(leh_stack, fun = sd)
ggplot() +
  geom_raster(aes(x, y, fill = value),
              fortify_raster(leh_mean))  +
  geom_polygon(aes(long, lat, group = group),
               fortify(leh_land),
               color = '#BBBBBB') +
  geom_contour(aes(x, y, z = value),
               fortify_raster(leh_mean),
               color = '#666666',
               linetype = '88888888') +
  geom_point(aes(x, y),
             data.frame(x = -160.1, 
                        y = 22.0),
             shape = 16,
             size = 3,
             color = 'black') +
  geom_point(aes(x, y),
             data.frame(x = -160.1, 
                        y = 22.0),
             shape = 1,
             size = 3,
             color = 'white') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      limits = c(0, 1),
                      name = expression(mean(E[rt]))) +
  theme_bw() +
  theme(legend.title = element_text(size = 9)) +
  coord_fixed() +
  labs(x = 'Longitude',
       y = 'Latitude')
ggsave('analysis/figures/LEH_mean.png',
       height = 3,
       width = 3.75,
       units = 'in')
ggplot() +
  geom_raster(aes(x, y, fill = value),
              fortify_raster(leh_sd))  +
  geom_polygon(aes(long, lat, group = group),
               fortify(leh_land),
               color = '#BBBBBB') +
  geom_contour(aes(x, y, z = value),
               fortify_raster(leh_sd),
               color = '#666666',
               linetype = '88888888') +
  geom_point(aes(x, y),
             data.frame(x = -160.1, 
                        y = 22.0),
             shape = 16,
             size = 3,
             color = 'black') +
  geom_point(aes(x, y),
             data.frame(x = -160.1, 
                        y = 22.0),
             shape = 1,
             size = 3,
             color = 'white') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      lim = c(0, 0.06),
                      breaks = seq(0, 0.06, by = 0.02),
                      name = expression(sd(E[rt]))) +
  theme_bw() +
  coord_fixed() +
  labs(x = 'Longitude',
       y = 'Latitude')
ggsave('analysis/figures/LEH_sd.png',
       height = 3,
       width = 3.75,
       units = 'in')

## Legends
ggplot() +
  geom_raster(aes(x, y, fill = value),
              fortify_raster(leh_mean))  +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      limits = c(0, 1),
                      name = expression(mean(E[rt]))) +
  theme_bw() +
  theme(legend.position = 'bottom')
ggsave('analysis/figures/legend_mean.png',
       height = 3,
       width = 3.75,
       units = 'in')
ggplot() +
  geom_raster(aes(x, y, fill = value),
              fortify_raster(leh_sd)) +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      lim = c(0, 0.06),
                      breaks = seq(0, 0.06, by = 0.02),
                      name = expression(sd(E[rt]))) +
  theme_bw() +
  theme(legend.position = 'bottom')
ggsave('analysis/figures/legend_sd.png',
       height = 3,
       width = 3.75,
       units = 'in')

## Tracks and residuals
extract_resid <- function(x, y, t, col) {
  date <- as.Date(t, tz = 'US/Hawaii')
  result <- data.frame(x = x, y = y, d = date) %>%
    group_by(d) %>%
    do({
      ert_meters <- raster(sprintf('data/out/EnergyLandscapesFR/all/KPC_%s_rt.tif',
                                   format(.$d[1], '%Y%m%d')),
                           crs = hi_aea_prj)
      col_meters <- spTransform(col, hi_aea_prj)
      d2c <- distanceFromPoints(ert_meters, col_meters) %>%
        mask(ert_meters) %>%
        projectRaster(crs = wgs84_prj)
      ert <- projectRaster(ert_meters, crs = wgs84_prj)
      rrt_lm <- data.frame(ert = getValues(ert), 
                           d2c = getValues(d2c)) %>%
        lm(ert ~ d2c, data = .)
      ert_pts <- raster::extract(ert, cbind(.$x, .$y))
      d2c_pts <- raster::extract(d2c, cbind(.$x, .$y))
      mutate(., rrt = ert_pts - predict(rrt_lm, newdata = data.frame(d2c = d2c_pts)))
    }) %>%
    ungroup
  result$rrt
}
kpc_col <- cbind(-159.3997, 22.22888) %>%
  SpatialPoints(proj4string = wgs84_prj)
tracks_rrt <- read_csv('data/KPC_tracks.csv') %>%
  filter(TimestampUTC < ymd('20160718')) %>%
  mutate(Rrt = extract_resid(Longitude, Latitude, TimestampUTC, kpc_col),
         D2C = geosphere::distGeo(cbind(Longitude, Latitude),
                                  cbind(-159.3997, 22.22888))) %>%
  group_by(TripID) %>%
  mutate(TimeElapsed = as.numeric(TimestampUTC - min(TimestampUTC), unit = 'secs')) %>%
  ungroup
farthest_point <- tracks_rrt %>%
  group_by(TripID) %>%
  do(arrange(., desc(D2C)) %>%
       slice(1)) %>%
  ungroup

mcb_col <- cbind(-157.7, 21.5) %>%
  SpatialPoints(proj4string = wgs84_prj)
leh_col <- cbind(-160.1, 22.0) %>%
  SpatialPoints(proj4string = wgs84_prj)
extract_resid2 <- function(x, y, t, col) {
  date <- as.Date(t, tz = 'US/Hawaii')
  result <- data.frame(x = x, y = y, d = date) %>%
    group_by(d) %>%
    do({
      ert_meters <- raster(sprintf('data/out/EnergyLandscapesFR/all/%s_%s_rt.tif',
                                   col,
                                   format(.$d[1], '%Y%m%d')),
                           crs = hi_aea_prj)
      if(col == 'MCB'){
        col_meters <- spTransform(mcb_col, hi_aea_prj)
      } else if(col == 'LEH') {
        col_meters <- spTransform(leh_col, hi_aea_prj)
      } else {
        stop('colony not found')
      }
      d2c <- distanceFromPoints(ert_meters, col_meters) %>%
        mask(ert_meters) %>%
        projectRaster(crs = wgs84_prj)
      ert <- projectRaster(ert_meters, crs = wgs84_prj)
      rrt_lm <- data.frame(ert = getValues(ert), 
                           d2c = getValues(d2c)) %>%
        lm(ert ~ d2c, data = .)
      ert_pts <- raster::extract(ert, cbind(.$x, .$y))
      d2c_pts <- raster::extract(d2c, cbind(.$x, .$y))
      mutate(., rrt = ert_pts - predict(rrt_lm, newdata = data.frame(d2c = d2c_pts)))
    }) %>%
    ungroup
  result$rrt
}
tracks_rrt_mcb <- read_csv('data/LEH_MCB_tracks.csv') %>%
  filter(SubColCode == 'MCB') %>%
  mutate(Rrt = extract_resid2(Longitude, Latitude, TimestampUTC, 'MCB'),
         D2C = geosphere::distGeo(cbind(Longitude, Latitude),
                                  cbind(-157.7, 21.5))) %>%
  group_by(TripID) %>%
  mutate(TimeElapsed = as.numeric(TimestampUTC - min(TimestampUTC), unit = 'secs')) %>%
  ungroup
farthest_point_mcb <- tracks_rrt_mcb %>%
  group_by(TripID) %>%
  do(arrange(., desc(D2C)) %>%
       slice(1)) %>%
  ungroup
tracks_rrt_leh <- read_csv('data/LEH_MCB_tracks.csv') %>%
  filter(SubColCode == 'LEH') %>%
  mutate(Rrt = extract_resid2(Longitude, Latitude, TimestampUTC, 'LEH'),
         D2C = geosphere::distGeo(cbind(Longitude, Latitude),
                                  cbind(-160.1, 22.0))) %>%
  group_by(TripID) %>%
  mutate(TimeElapsed = as.numeric(TimestampUTC - min(TimestampUTC), unit = 'secs')) %>%
  ungroup
farthest_point_leh <- tracks_rrt_leh %>%
  group_by(TripID) %>%
  do(arrange(., desc(D2C)) %>%
       slice(1)) %>%
  ungroup
# ggsave('analysis/figures/TripRRT_LEH.png',
#        height = 3,
#        width = 3,
#        units = 'in')

rbind(farthest_point,
      farthest_point_mcb,
      farthest_point_leh) %>% 
  mutate(SubColCode = factor(SubColCode, levels = c('KPC', 'MCB', 'LEH'))) %>%
  ggplot(aes(Rrt)) +
  geom_histogram(binwidth = 0.02,
                 boundary = 0) +
  geom_vline(xintercept = 0,
             color = 'red') +
  theme_bw() +
  labs(x = expression(R[RT]),
       y = 'Count of Trips') +
  facet_wrap(~ SubColCode,
             scales = 'free_y')
ggsave('analysis/figures/farthestpoint.png',
       height = 2.5,
       width = 5,
       units = 'in')

rbind(farthest_point,
      farthest_point_mcb,
      farthest_point_leh) %>% 
  group_by(SubColCode) %>%
  summarize(meanRrt = mean(Rrt, na.rm = TRUE),
            sdRrt = sd(Rrt, na.rm = TRUE))

rbind(farthest_point,
      farthest_point_mcb,
      farthest_point_leh) %>% 
  mutate(SubColCode = factor(SubColCode, levels = c('KPC', 'MCB', 'LEH'))) %>%
  ggplot(aes(D2C, Rrt)) +
  geom_point() +
  theme_bw() +
  labs(x = 'Dist. to Colony (m)',
       y = expression(R[RT])) +
  facet_wrap(~ SubColCode)

# KPC vs MCB
rbind(farthest_point,
      farthest_point_mcb) %>% 
  mutate(SubColCode = factor(SubColCode)) %>%
  t.test(Rrt ~ SubColCode, data = .)
# KPC vs LEH
rbind(farthest_point,
      farthest_point_leh) %>% 
  mutate(SubColCode = factor(SubColCode)) %>%
  t.test(Rrt ~ SubColCode, data = .)
# MCB vs LEH
rbind(farthest_point_mcb,
      farthest_point_leh) %>% 
  mutate(SubColCode = factor(SubColCode)) %>%
  t.test(Rrt ~ SubColCode, data = .)
# Greater than 0?
t.test(farthest_point$Rrt, mu = 0, alternative = 'greater')
t.test(farthest_point_mcb$Rrt, mu = 0, alternative = 'greater')
t.test(farthest_point_leh$Rrt, mu = 0, alternative = 'less')

## Presences
### KPC first half (May 28 - June 22)
kpc_tracks1 <- read_csv('data/out/Presences/rfbo_accessible.csv') %>%
  filter(LocDate < ymd('20160622'))
kpc_tracks1_sp <- SpatialPoints(dplyr::select(kpc_tracks1, Longitude, Latitude),
                                proj4string = wgs84_prj) %>%
  spTransform(hi_aea_prj)
template <- raster('data/out/EnergyLandscapesFR/all/KPC_20160528_rt.tif',
                   crs = hi_aea_prj)
kpc_ud1 <- adehabitatHR::kernelUD(kpc_tracks1_sp) %>%
  raster %>%
  resample(template) %>%
  mask(template) %>%
  projectRaster(crs = wgs84_prj)
hi_land <- rgdal::readOGR('data/coastline/hi_land', 'hi_land') %>%
  spTransform(hi_aea_prj) %>%
  crop(template) %>%
  spTransform(wgs84_prj)
ggplot(fortify_raster(kpc_ud1), aes(x, y, fill = value)) +
  geom_raster() +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               inherit.aes = FALSE) +
  labs(x = 'Longitude',
       y = 'Latitude',
       fill = 'UD') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26') +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/KPC_UD_1.png',
       height = 3,
       width = 3.75,
       units = 'in')

### KPC second half (June 23 - July 18)
kpc_tracks2 <- read_csv('data/out/Presences/rfbo_accessible.csv') %>%
  filter(LocDate >= ymd('20160622'))
kpc_tracks2_sp <- SpatialPoints(dplyr::select(kpc_tracks2, Longitude, Latitude),
                                proj4string = wgs84_prj) %>%
  spTransform(hi_aea_prj)
kpc_ud2 <- adehabitatHR::kernelUD(kpc_tracks2_sp) %>%
  raster %>%
  resample(template) %>%
  mask(template) %>%
  projectRaster(crs = wgs84_prj)
ggplot(fortify_raster(kpc_ud2), aes(x, y, fill = value)) +
  geom_raster() +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               inherit.aes = FALSE) +
  labs(x = 'Longitude',
       y = 'Latitude',
       fill = 'UD') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26') +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/KPC_UD_2.png',
       height = 3,
       width = 3.75,
       units = 'in')

## Mean energy landscape during those periods
kpc_el1 <- stack(dir('data/out/EnergyLandscapesFR/KPC/Rasters/Roundtrip',
                     full.names = TRUE)[2:27]) %>%
  calc(fun = mean)
crs(kpc_el1) <- hi_aea_prj
kpc_el1 <- projectRaster(kpc_el1, crs = wgs84_prj)
ggplot(fortify_raster(kpc_el1), aes(x, y, fill = value)) +
  geom_raster() +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               inherit.aes = FALSE) +
  geom_contour(aes(x, y, z = value),
               inherit.aes = FALSE,
               color = '#666666',
               linetype = '88888888') +
  labs(x = 'Longitude',
       y = 'Latitude') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      lim = c(0, 1),
                      name = expression(mean(E[RT]))) +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/KPC_EL_1.png',
       height = 3,
       width = 3.75,
       units = 'in')
kpc_el2 <- stack(dir('data/out/EnergyLandscapesFR/KPC/Rasters/Roundtrip',
                     full.names = TRUE)[28:52]) %>%
  calc(fun = mean)
crs(kpc_el2) <- hi_aea_prj
kpc_el2 <- projectRaster(kpc_el2, crs = wgs84_prj)
ggplot(fortify_raster(kpc_el2), aes(x, y, fill = value)) +
  geom_raster() +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               inherit.aes = FALSE) +
  geom_contour(aes(x, y, z = value),
               inherit.aes = FALSE,
               color = '#666666',
               linetype = '88888888') +
  labs(x = 'Longitude',
       y = 'Latitude') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      lim = c(0, 1),
                      name = expression(mean(E[RT]))) +
  theme_bw() +
  coord_fixed()
ggsave('analysis/figures/KPC_EL_2.png',
       height = 3,
       width = 3.75,
       units = 'in')
### MCB sessions
foreach(presences = dir('data/out/Presences/MCB', full.names = TRUE)) %do% {
  tracks <- read_csv(presences)
  tracks_sp <- SpatialPoints(dplyr::select(tracks, Longitude, Latitude),
                             proj4string = wgs84_prj) %>%
    spTransform(hi_aea_prj)
  ud <- adehabitatHR::kernelUD(tracks_sp) %>%
    raster %>%
    resample(template) %>%
    mask(template) %>%
    projectRaster(crs = wgs84_prj)
  template <- raster('data/out/EnergyLandscapesFR/all/MCB_20140601_rt.tif',
                     crs = hi_aea_prj)
  hi_land <- rgdal::readOGR('data/coastline/hi_land', 'hi_land') %>%
    spTransform(hi_aea_prj) %>%
    crop(template) %>%
    spTransform(wgs84_prj)
  p <- ggplot(fortify_raster(ud), aes(x, y, fill = value)) +
    geom_raster() +
    geom_polygon(aes(long, lat, group = group),
                 fortify(hi_land),
                 inherit.aes = FALSE) +
    labs(x = 'Longitude',
         y = 'Latitude',
         fill = 'UD') +
    scale_fill_gradient(low = '#fee0d2',
                        high = '#de2d26') +
    theme_bw() +
    coord_fixed()
  sess_name <- basename(presences) %>% substr(1, nchar(.) - 4)
  ggsave(sprintf('analysis/figures/UD/MCB_%s.png', sess_name),
         p,
         height = 3,
         width = 3.75,
         units = 'in')
}

### LEH sessions
foreach(presences = dir('data/out/Presences/LEH', full.names = TRUE)) %do% {
  tracks <- read_csv(presences)
  tracks_sp <- SpatialPoints(dplyr::select(tracks, Longitude, Latitude),
                             proj4string = wgs84_prj) %>%
    spTransform(hi_aea_prj)
  template <- predict_envonly
  template <- raster('data/out/EnergyLandscapesFR/all/LEH_20140513_rt.tif',
                     crs = hi_aea_prj)
  # ud <- adehabitatHR::kernelUD(tracks_sp) %>%
  #   raster %>%
  #   resample(template) %>%
  #   mask(template) %>%
  #   projectRaster(crs = wgs84_prj)
  ud <- adehabitatHR::kernelUD(tracks_sp) %>%
    raster %>%
    resample(template) %>%
    mask(template) %>%
    projectRaster(crs = wgs84_prj)
  hi_land <- rgdal::readOGR('data/coastline/hi_land', 'hi_land') %>%
    spTransform(hi_aea_prj) %>%
    crop(template) %>%
    spTransform(wgs84_prj)
  p <- ggplot(fortify_raster(ud), aes(x, y, fill = value)) +
    geom_raster() +
    geom_polygon(aes(long, lat, group = group),
                 fortify(hi_land),
                 inherit.aes = FALSE) +
    labs(x = 'Longitude',
         y = 'Latitude',
         fill = 'UD') +
    scale_fill_gradient(low = '#fee0d2',
                        high = '#de2d26') +
    theme_bw() +
    coord_fixed()
  sess_name <- basename(presences) %>% substr(1, nchar(.) - 4)
  ggsave(sprintf('analysis/figures/UD/LEH_%s.png', sess_name),
         p,
         height = 3,
         width = 3.75,
         units = 'in')
}

