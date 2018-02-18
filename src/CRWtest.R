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
  prepData(type = 'LL')

summary(tracks_prepped)

# Choose parameter initial settings
params0 <- left_join(select(raw_tracks, ID, TimestampUTC, Behavior), 
                     select(tracks_prepped, ID, TimestampUTC, step, angle)) %>% 
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
sim_tracks <- simData(nbAnimals = 10, 
                          model = tracksHMM, 
                          obsPerAnimal = c(340, 1058)) %>%
  mutate(ID = factor(ID))
ggplot(sim_tracks, aes(x, y, color = ID)) +
  geom_path() +
  coord_fixed() +
  labs(x = 'X (km)',
       y = 'Y (km)',
       title = 'Simulated Correlated Random Walks') 
