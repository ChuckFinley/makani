library(rgdal)
library(raster)
library(tidyverse)
library(moveHMM)
library(geosphere)
library(CircStats)
select <- dplyr::select

# Load tracking data
raw_tracks <- read_csv('data/tracking/RFBOtracksLEH.csv') %>%
  mutate_at(vars(DeployID, TripID, TripNumber, Behavior, Site, SubColCode),
            funs(factor)) %>%
  rename(ID = TripID,
         x = Longitude,
         y = Latitude) %>%
  arrange(ID, TimestampUTC) %>%
  as.data.frame

summary(raw_tracks)

# Prepare data and estimate initial parameters
# tracks_prepped <- prepData(raw_tracks, type = 'LL')
# Just the first ten trips from Lehua in 2015
tracks_prepped <- raw_tracks %>%
  filter(SubColCode == 'LEH',
         Year == 2015) %>%
  group_by(ID) %>%
  summarize() %>%
  slice(1:10) %>%
  semi_join(raw_tracks, .) %>% 
  select(x, y, ID, TimestampUTC) %>%
  prepData(type = 'LL')

summary(tracks_prepped)

# Choose parameter initial settings
params0 <- left_join(transmute(raw_tracks, as.character(ID), TimestampUTC, Behavior), 
                     transmute(tracks_prepped, as.character(ID), TimestampUTC, step, angle)) %>% 
  na.omit %>% 
  group_by(Behavior) %>%
  summarize(meanStep = mean(step),
            sdStep = sd(step),
            zeroMass = sum(step == 0, na.rm = TRUE) / n(),
            meanAngle = circ.mean(angle),
            concAngle = est.kappa(angle)) %>%
  as.data.frame

# NOTE: if some step lengths are 0, include zero mass
step_par <- c(params0$meanStep, params0$sdStep)
if(any(tracks_prepped$step == 0, na.rm = TRUE)) {
  step_par <- c(step_par, params0$zeroMass)
} 
angle_par <- c(params0$meanAngle, params0$concAngle)

# Fit model
tracksHMM <- fitHMM(data = tracks_prepped,
                    nbStates = 3,
                    stepPar0 = step_par,
                    anglePar0 = angle_par)

# Check model
plot(tracksHMM)
plotStates(tracksHMM, '95200001')
plotPR(tracksHMM)

# Simulate tracks
sim_tracks <- simData(nbAnimals = 100, 
                      model = tracksHMM, 
                      obsPerAnimal = c(340, 1058)) %>%
  mutate(ID = factor(ID))

colony_crw_ud <- function(crw, col_coords, radius, resolution) {
  # Convert to UD
  wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
  hi_aea_prj <- CRS('+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-163 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
  
  crw_sp <- crw %>%
    mutate(bearing = atan2(y, x) * 180 / pi,
           dist = sqrt(x^2 + y^2),
           lon = geosphere::destPoint(col_coords, bearing, dist*1000)[,1],
           lat = geosphere::destPoint(col_coords, bearing, dist*1000)[,2]) %>%
    # Select longitude, latitude columns
    select(lon, lat) %>%
    # Create SpatialPoints object
    SpatialPoints(proj4string = wgs84_prj) %>%
    # Project to Hawaii Albers Equal Area Conic
    spTransform(hi_aea_prj)
  crw_ud <- adehabitatHR::kernelUD(crw_sp)
  # Crop, re-sample, rescale
  col_extent <- SpatialPoints(matrix(col_coords, ncol = 2),
                              proj4string = wgs84_prj) %>%
    spTransform(hi_aea_prj) %>%
    rgeos::gBuffer(width = radius)
  col_template <- raster(col_extent)
  res(col_template) <- resolution
  rescale <- function(r) {
    r_min = cellStats(r, "min")
    r_max = cellStats(r, "max")
    (r - r_min) / (r_max - r_min)
  }
  raster(crw_ud) %>%
    resample(col_template) %>%
    crop(col_extent) %>%
    rescale
}

crw_res <- res(raster('data/out/EnergyLandscapes/KPC/20160528_in.tif'))
kpc_ud <- colony_crw_ud(sim_tracks, c(-159.40, 22.23), 250e3, crw_res)
ggplot(fortify_raster(kpc_crw_ud), aes(x, y, fill = val)) + 
  geom_raster() +
  scale_fill_gradientn(colors = colorRamps::matlab.like(4))
writeRaster(kpc_ud, 'data/out/CyberBirds/KPC_CRW_UD.tif', 'GTiff', overwrite = TRUE)

leh_ud <- colony_crw_ud(sim_tracks, c(-160.1, 22.0), 250e3, crw_res)
writeRaster(leh_ud, 'data/out/CyberBirds/LEH_CRW_UD.tif', 'GTiff', overwrite = TRUE)
mcb_ud <- colony_crw_ud(sim_tracks, c(-157.7, 21.5), 250e3, crw_res)
writeRaster(mcb_ud, 'data/out/CyberBirds/MCB_CRW_UD.tif', 'GTiff', overwrite = TRUE)
  