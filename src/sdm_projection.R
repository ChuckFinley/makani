library(tidyverse)
library(sp)
library(raster)
library(lubridate)

# Colony data
wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
hi_aea_prj <- CRS('+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-163 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
project_col <- function(latlon) {
  data.frame(x = latlon[2], y = latlon[1]) %>%
    SpatialPoints(wgs84_prj) %>%
    spTransform(hi_aea_prj) %>%
    as.data.frame
}
kpc_wgs84 <- c(22.2, -159.4)
leh_wgs84 <- c(22.0, -160.1)
mcb_wgs84 <- c(21.5, -157.7)
kpc_hi_aea <- project_col(kpc_wgs84)
leh_hi_aea <- project_col(leh_wgs84)
mcb_hi_aea <- project_col(mcb_wgs84)
hi_land <- rgdal::readOGR('data/coastline/hi_land', 'hi_land')

# Around KPC for 2016-06-01
start_time <- ymd_hm('20160601 00:00', tz = 'US/Hawaii')
end_time <- ymd_hm('20160602 00:00', tz = 'US/Hawaii')
t <- as.numeric(start_time)
sst_nc <- RNetCDF::open.nc('data/sst/2016sst.nc')
time_sst <- RNetCDF::var.get.nc(sst_nc, 'time')
t_min <- min(time_sst)
t_max <- max(time_sst)
t_i <- round((length(time_sst) - 1) * (t - t_min) / (t_max - t_min)) + 1
sst_stack <- stack('data/sst/2016sst.nc', varname = 'sst')
anom_stack <- stack('data/sst/2016sst.nc', varname = 'anom')
writeRaster(sst_stack[[t_i]], 'data/SDM_env/kpc_20160601/sst.asc', 
            overwrite = TRUE)
writeRaster(anom_stack[[t_i]], 'data/SDM_env/kpc_20160601/anom.asc', 
            overwrite = TRUE)
bathy_r <- raster('data/bathy/etopo1.tif')
extent(bathy_r)[1:2] <- extent(bathy_r)[1:2] + 360
bathy_r %>%
  resample(sst_stack[[1]]) %>%
  crop(sst_stack[[1]]) %>% 
  writeRaster('data/SDM_env/kpc_20160601/bathy.asc', overwrite = TRUE)
ert_r <- raster('data/out/EnergyLandscapesFR/all/KPC_20160601_rt.tif',
       crs = hi_aea_prj) %>%
  projectRaster(crs = wgs84_prj)
extent(ert_r)[1:2] <- extent(ert_r)[1:2] + 360
ert_r %>%
  resample(sst_stack[[1]]) %>%
  crop(sst_stack[[1]]) %>% 
  writeRaster('data/SDM_env/kpc_20160601/ert.asc', overwrite = TRUE)
d2c_template <- raster('data/out/EnergyLandscapesFR/all/KPC_20160601_rt.tif',
                       crs = hi_aea_prj)
d2c_r <- distanceFromPoints(d2c_template, kpc_hi_aea) %>%
  projectRaster(crs = wgs84_prj)
hi_land2 <- crop(hi_land, d2c_r)
extent(d2c_r)[1:2] <- extent(d2c_r)[1:2] + 360
d2c_r %>%
  resample(sst_stack[[1]]) %>%
  crop(sst_stack[[1]]) %>% 
  writeRaster('data/SDM_env/kpc_20160601/d2c.asc', overwrite = TRUE)

system('java -cp tools/maxent/maxent.jar density.Project analysis/maxent/KPC/envonly/RFBO.lambdas data/SDM_env/kpc_20160601 analysis/maxent/KPC/projections/KPC_envonly.asc')
system('java -cp tools/maxent/maxent.jar density.Project analysis/maxent/KPC/ert/RFBO.lambdas data/SDM_env/kpc_20160601 analysis/maxent/KPC/projections/KPC_ert.asc')
system('java -cp tools/maxent/maxent.jar density.Project analysis/maxent/KPC/d2c/RFBO.lambdas data/SDM_env/kpc_20160601 analysis/maxent/KPC/projections/KPC_d2c.asc')

predict_envonly <- raster('analysis/maxent/KPC/projections/KPC_envonly.asc',
                          crs = wgs84_prj)
predict_envonly <- predict_envonly %>%
  resample(ert_r) %>%
  mask(ert_r) %>%
  rasterToPoints %>%
  data.frame %>%
  mutate(x = x - 360)
colnames(predict_envonly)[3] <- 'ROR'

predict_d2c <- raster('analysis/maxent/KPC/projections/KPC_d2c.asc',
                      crs = wgs84_prj)
predict_d2c <- predict_d2c %>%
  resample(ert_r) %>%
  mask(ert_r) %>%
  rasterToPoints %>%
  data.frame %>%
  mutate(x = x - 360)
colnames(predict_d2c)[3] <- 'ROR'

predict_ert <- raster('analysis/maxent/KPC/projections/KPC_ert.asc',
                      crs = wgs84_prj)
predict_ert <- predict_ert %>%
  resample(ert_r) %>%
  mask(ert_r) %>%
  rasterToPoints %>%
  data.frame %>%
  mutate(x = x - 360)
colnames(predict_ert)[3] <- 'ROR'

# Tracks
tracks_df <- read_csv('data/KPC_tracks.csv') %>%
  filter(SubColCode == 'KPC',
         between(TimestampUTC, start_time, end_time))

## Suitability only
ggplot() +
  geom_raster(aes(x, y, fill = ROR), predict_envonly) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land2),
               color = '#BBBBBB') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      name = 'ROR',
                      lim = c(0, 1)) +
  #geom_path(aes(Longitude, Latitude), tracks_df, alpha = 0.5) +
  coord_fixed() +
  theme_bw() +
  labs(x = '',
       y = '')
ggsave('analysis/figures/KPC_envonly_20160601.png',
       height = 3,
       width = 3.75)

## D2C
ggplot() +
  geom_raster(aes(x, y, fill = ROR), predict_d2c) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land2),
               color = '#BBBBBB') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      name = 'ROR',
                      lim = c(0, 1)) +
  #geom_path(aes(Longitude, Latitude), tracks_df, alpha = 0.5) +
  coord_fixed() +
  theme_bw() +
  labs(x = '',
       y = '')
ggsave('analysis/figures/KPC_d2c_20160601.png',
       height = 3,
       width = 3.75)
## Ert
ggplot() +
  geom_raster(aes(x, y, fill = ROR), predict_ert) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land2),
               color = '#BBBBBB') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      name = 'ROR',
                      lim = c(0, 1)) +
  #geom_path(aes(Longitude, Latitude), tracks_df, alpha = 0.5) +
  coord_fixed() +
  theme_bw() +
  labs(x = '',
       y = '')
ggsave('analysis/figures/KPC_ert_20160601.png',
       height = 3,
       width = 3.75)

# correlation
cor(predict_ert$ROR, predict_d2c$ROR)
cor(predict_ert$ROR, predict_envonly$ROR)
cor(predict_envonly$ROR, predict_d2c$ROR)

# Around LEH for 2015-05-28
start_time <- ymd_hm('20150528 00:00', tz = 'US/Hawaii')
end_time <- ymd_hm('20150529 00:00', tz = 'US/Hawaii')
t <- as.numeric(start_time)
sst_nc <- RNetCDF::open.nc('data/sst/2015sst.nc')
time_sst <- RNetCDF::var.get.nc(sst_nc, 'time')
t_min <- min(time_sst)
t_max <- max(time_sst)
t_i <- round((length(time_sst) - 1) * (t - t_min) / (t_max - t_min)) + 1
sst_stack <- stack('data/sst/2015sst.nc', varname = 'sst')
anom_stack <- stack('data/sst/2015sst.nc', varname = 'anom')
writeRaster(sst_stack[[t_i]], 'data/SDM_env/leh_20150528/sst.asc', overwrite = TRUE)
writeRaster(anom_stack[[t_i]], 'data/SDM_env/leh_20150528/anom.asc', overwrite = TRUE)
bathy_r <- raster('data/bathy/etopo1.tif')
extent(bathy_r)[1:2] <- extent(bathy_r)[1:2] + 360
bathy_r %>%
  resample(sst_stack[[1]]) %>%
  crop(sst_stack[[1]]) %>% 
  writeRaster('data/SDM_env/leh_20150528/bathy.asc', overwrite = TRUE)
ert_r <- raster('data/out/EnergyLandscapesFR/all/LEH_20150528_rt.tif',
                crs = hi_aea_prj) %>%
  projectRaster(crs = wgs84_prj)
extent(ert_r)[1:2] <- extent(ert_r)[1:2] + 360
ert_r %>%
  resample(sst_stack[[1]]) %>%
  crop(sst_stack[[1]]) %>% 
  writeRaster('data/SDM_env/leh_20150528/ert.asc', overwrite = TRUE)
d2c_template <- raster('data/out/EnergyLandscapesFR/all/LEH_20150528_rt.tif',
                       crs = hi_aea_prj)
d2c_r <- distanceFromPoints(d2c_template, leh_hi_aea) %>%
  projectRaster(crs = wgs84_prj)
hi_land2 <- crop(hi_land, d2c_r)
extent(d2c_r)[1:2] <- extent(d2c_r)[1:2] + 360
d2c_r %>%
  resample(sst_stack[[1]]) %>%
  crop(sst_stack[[1]]) %>% 
  writeRaster('data/SDM_env/leh_20150528/d2c.asc', overwrite = TRUE)

system('java -cp tools/maxent/maxent.jar density.Project analysis/maxent/KPC/envonly/RFBO.lambdas data/SDM_env/leh_20150528 analysis/maxent/LEH/projections/LEH_envonly.asc')
system('java -cp tools/maxent/maxent.jar density.Project analysis/maxent/KPC/ert/RFBO.lambdas data/SDM_env/leh_20150528 analysis/maxent/LEH/projections/LEH_ert.asc')
system('java -cp tools/maxent/maxent.jar density.Project analysis/maxent/KPC/d2c/RFBO.lambdas data/SDM_env/leh_20150528 analysis/maxent/LEH/projections/LEH_d2c.asc')

predict_envonly <- raster('analysis/maxent/LEH/projections/LEH_envonly.asc',
                          crs = wgs84_prj) %>%
  resample(ert_r) %>%
  mask(ert_r) %>%
  rasterToPoints %>%
  data.frame
colnames(predict_envonly)[3] <- 'ROR'
predict_envonly[, 'x'] <- predict_envonly[, 'x'] - 360
predict_d2c <- raster('analysis/maxent/LEH/projections/LEH_d2c.asc',
                      crs = wgs84_prj) %>%
  resample(ert_r) %>%
  mask(ert_r) %>%
  rasterToPoints %>%
  data.frame
colnames(predict_d2c)[3] <- 'ROR'
predict_d2c[, 'x'] <- predict_d2c[, 'x'] - 360
predict_ert <- raster('analysis/maxent/LEH/projections/LEH_ert.asc',
                      crs = wgs84_prj) %>%
  resample(ert_r) %>%
  mask(ert_r) %>%
  rasterToPoints %>%
  data.frame
colnames(predict_ert)[3] <- 'ROR'
predict_ert[, 'x'] <- predict_ert[, 'x'] - 360

# Tracks
tracks_df <- read_csv('data/LEH_MCB_tracks.csv') %>%
  filter(SubColCode == 'LEH',
         between(TimestampUTC, start_time, end_time))

## Suitability only
ggplot() +
  geom_raster(aes(x, y, fill = ROR), predict_envonly) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land2),
               color = '#BBBBBB') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      name = 'ROR',
                      lim = c(0, 1)) +
  #geom_path(aes(Longitude, Latitude), tracks_df, alpha = 0.5) +
  coord_fixed() +
  theme_bw() +
  labs(x = '',
       y = '')
ggsave('analysis/figures/LEH_envonly_20150528.png',
       height = 3,
       width = 3.75)

## D2C
ggplot() +
  geom_raster(aes(x, y, fill = ROR), predict_d2c) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land2),
               color = '#BBBBBB') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      name = 'ROR',
                      lim = c(0, 1)) +
  geom_path(aes(Longitude, Latitude), tracks_df, alpha = 0.5) +
  coord_fixed() +
  theme_bw() +
  labs(x = '',
       y = '')
ggsave('analysis/figures/LEH_d2c_20150528.png',
       height = 3,
       width = 3.75)
## Ert
ggplot() +
  geom_raster(aes(x, y, fill = ROR), predict_ert) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land2),
               color = '#BBBBBB') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      name = 'ROR',
                      lim = c(0, 1)) +
  geom_path(aes(Longitude, Latitude), tracks_df, alpha = 0.5) +
  coord_fixed() +
  theme_bw() +
  labs(x = '',
       y = '')
ggsave('analysis/figures/LEH_ert_20150528.png',
       height = 3,
       width = 3.75)

# Around MCB for 2014-06-05
start_time <- ymd_hm('20140605 00:00', tz = 'US/Hawaii')
end_time <- ymd_hm('20140606 00:00', tz = 'US/Hawaii')
t <- as.numeric(start_time)
sst_nc <- RNetCDF::open.nc('data/sst/2014sst.nc')
time_sst <- RNetCDF::var.get.nc(sst_nc, 'time')
t_min <- min(time_sst)
t_max <- max(time_sst)
t_i <- round((length(time_sst) - 1) * (t - t_min) / (t_max - t_min)) + 1
sst_stack <- stack('data/sst/2014sst.nc', varname = 'sst')
anom_stack <- stack('data/sst/2014sst.nc', varname = 'anom')
writeRaster(sst_stack[[t_i]], 'data/SDM_env/mcb_20140605/sst.asc')
writeRaster(anom_stack[[t_i]], 'data/SDM_env/mcb_20140605/anom.asc')
bathy_r <- raster('data/bathy/etopo1.tif')
extent(bathy_r)[1:2] <- extent(bathy_r)[1:2] + 360
bathy_r %>%
  resample(sst_stack[[1]]) %>%
  crop(sst_stack[[1]]) %>% 
  writeRaster('data/SDM_env/mcb_20140605/bathy.asc', overwrite = TRUE)
ert_r <- raster('data/out/EnergyLandscapesFR/all/MCB_20140605_rt.tif',
                crs = hi_aea_prj) %>%
  projectRaster(crs = wgs84_prj)
extent(ert_r)[1:2] <- extent(ert_r)[1:2] + 360
ert_r %>%
  resample(sst_stack[[1]]) %>%
  crop(sst_stack[[1]]) %>% 
  writeRaster('data/SDM_env/mcb_20140605/ert.asc', overwrite = TRUE)
d2c_template <- raster('data/out/EnergyLandscapesFR/all/MCB_20140605_rt.tif',
                       crs = hi_aea_prj)
d2c_r <- distanceFromPoints(d2c_template, mcb_hi_aea) %>%
  projectRaster(crs = wgs84_prj)
hi_land2 <- crop(hi_land, d2c_r)
extent(d2c_r)[1:2] <- extent(d2c_r)[1:2] + 360
d2c_r %>%
  resample(sst_stack[[1]]) %>%
  crop(sst_stack[[1]]) %>% 
  writeRaster('data/SDM_env/mcb_20140605/d2c.asc', overwrite = TRUE)

system('java -cp tools/maxent/maxent.jar density.Project analysis/maxent/KPC/envonly/RFBO.lambdas data/SDM_env/mcb_20140605 analysis/maxent/MCB/projections/MCB_envonly.asc')
system('java -cp tools/maxent/maxent.jar density.Project analysis/maxent/KPC/ert/RFBO.lambdas data/SDM_env/mcb_20140605 analysis/maxent/MCB/projections/MCB_ert.asc')
system('java -cp tools/maxent/maxent.jar density.Project analysis/maxent/KPC/d2c/RFBO.lambdas data/SDM_env/mcb_20140605 analysis/maxent/MCB/projections/MCB_d2c.asc')

predict_envonly <- raster('analysis/maxent/MCB/projections/MCB_envonly.asc',
                          crs = wgs84_prj) %>%
  resample(ert_r) %>%
  mask(ert_r) %>%
  rasterToPoints %>%
  data.frame
colnames(predict_envonly)[3] <- 'ROR'
predict_envonly[, 'x'] <- predict_envonly[, 'x'] - 360
predict_d2c <- raster('analysis/maxent/MCB/projections/MCB_d2c.asc',
                      crs = wgs84_prj) %>%
  resample(ert_r) %>%
  mask(ert_r) %>%
  rasterToPoints %>%
  data.frame
colnames(predict_d2c)[3] <- 'ROR'
predict_d2c[, 'x'] <- predict_d2c[, 'x'] - 360
predict_ert <- raster('analysis/maxent/MCB/projections/MCB_ert.asc',
                      crs = wgs84_prj) %>%
  resample(ert_r) %>%
  mask(ert_r) %>%
  rasterToPoints %>%
  data.frame
colnames(predict_ert)[3] <- 'ROR'
predict_ert[, 'x'] <- predict_ert[, 'x'] - 360

# Tracks
tracks_df <- read_csv('data/LEH_MCB_tracks.csv') %>% 
  filter(SubColCode == 'MCB',
         between(TimestampUTC, start_time, end_time))

## Suitability only
ggplot() +
  geom_raster(aes(x, y, fill = ROR), predict_envonly) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land2),
               color = '#BBBBBB') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      name = 'ROR',
                      lim = c(0,1)) +
  geom_path(aes(Longitude, Latitude), tracks_df, alpha = 0.5) +
  coord_fixed() +
  theme_bw() +
  labs(x = '',
       y = '')

ggsave('analysis/figures/MCB_envonly_20140605.png',
       height = 3,
       width = 3.75)

## D2C
ggplot() +
  geom_raster(aes(x, y, fill = ROR), predict_d2c) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land2),
               color = '#BBBBBB') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      name = 'ROR',
                      lim = c(0,1)) +
  geom_path(aes(Longitude, Latitude), tracks_df, alpha = 0.5) +
  coord_fixed() +
  theme_bw() +
  labs(x = '',
       y = '')
ggsave('analysis/figures/MCB_d2c_20140605.png',
       height = 3,
       width = 3.75)
## Ert
ggplot() +
  geom_raster(aes(x, y, fill = ROR), predict_ert) +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land2),
               color = '#BBBBBB') +
  scale_fill_gradient(low = '#fee0d2',
                      high = '#de2d26',
                      name = 'ROR',
                      lim = c(0,1)) +
  geom_path(aes(Longitude, Latitude), tracks_df, alpha = 0.5) +
  coord_fixed() +
  theme_bw() +
  labs(x = '',
       y = '')
ggsave('analysis/figures/MCB_ert_20140605.png',
       height = 3,
       width = 3.75)

