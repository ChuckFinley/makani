library(tidyverse)

wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
hi_land <- rgdal::readOGR('data/coastline/hi_land', 'hi_land') %>%
  sp::spTransform(wgs84_prj) %>%
  crop(extent(-160, -158, 22, 23.5))

load('data/out/TransitSegments/SegSample.RData')

tripid <- 115300027
tracks <- read_csv('data/KPC_tracks.csv') %>% 
  filter(TripID == tripid)
segments <- filter(segSample, TripID == tripid)
out_tracks <- filter(tracks, 
                     between(TimestampUTC, 
                             segments$Start[1], 
                             segments$End[1]))
in_tracks <- filter(tracks, 
                    between(TimestampUTC, 
                            segments$Start[2], 
                            segments$End[2]))

# Map
ggplot() +
  geom_polygon(aes(long, lat, group = group),
               fortify(hi_land),
               fill = 'light grey') +
  geom_path(aes(Longitude, Latitude),
            tracks) +
  geom_path(aes(Longitude, Latitude),
            out_tracks,
            color = 'orange',
            size = 2) +
  geom_path(aes(Longitude, Latitude),
            in_tracks,
            color = 'green',
            size = 2) +
  xlim(-160, -158) +
  ylim(22, 23.5) +
  coord_fixed() +
  theme_bw() +
  labs(x = '',
       y = '')
ggsave('analysis/figures/ExampleTripSegments.png',
       height = 3.8,
       width = 5,
       units = 'in')

# Projection
flight_x <- segments$meanSpd[1] * cos(segments$bearing[1] * pi / 180)
flight_y <- segments$meanSpd[1] * sin(segments$bearing[1] * pi / 180)
projection_df <- data.frame(x0 = c(0,0),
                            y0 = c(0, 0),
                            x1 = c(flight_x, segments$uWind[1]),
                            y1 = c(flight_y, segments$vWind[1]),
                            label = c('Flight Vector', 'Wind Vector'))
ggplot() + 
  geom_segment(aes(x0, y0, xend = x1, yend = y1, color = label),
               projection_df,
               arrow = arrow(length = unit(0.05, "npc")),
               size = 1.2) +
  scale_color_manual(values = c('green', 'purple')) +
  theme_bw() +
  coord_fixed() +
  labs(x = 'x (m/s)',
       y = 'y (m/s)',
       color = '') +
  theme(legend.position = 'bottom')
ggsave('analysis/figures/ExampleProjection.png',
       height = 5,
       width = 3.8,
       units = 'in')

# Acceleration
get.acc <- function(segID) {
  POSIX.origin <- ymd_hm('19700101 00:00', tz = 'UTC')
  mhi_db <- DBI::dbConnect(RSQLite::SQLite(), dbname = 'data/MHI_GPS.sqlite')
  this.acc_db <- tbl(mhi_db, 'SegmentACC')
  result <- this.acc_db %>%
    filter(SegmentID == segID) %>%
    collect(n = Inf) %>%
    transmute(DeployID,
              SegmentID,
              BurstID,
              ACCTimestampUTC = as.POSIXct(ACCTimestampUTC,
                                           tz = 'UTC',
                                           origin = POSIX.origin),
              X, Y, Z,
              DynX, DynY, DynZ,
              ODBA)
  DBI::dbDisconnect(mhi_db)
  result
}

acc_data_in <- get.acc(6131) %>%
  group_by(BurstID) %>%
  mutate(Time = as.numeric(ACCTimestampUTC - min(ACCTimestampUTC), unit = 'secs')) %>%
  ungroup

# glider
acc_data_in %>% 
  filter(BurstID == 39363) %>%
  gather(axis, acc, DynX, DynY, DynZ) %>%
  ggplot(aes(Time, acc)) +
  geom_line() +
  facet_wrap(~ axis) +
  theme_bw() +
  labs(x = 'Time (s)',
       y = 'Acceleration (m/s^2)')
ggsave('analysis/figures/Glider.png',
       height = 4,
       width = 6,
       units = 'in')
# flapper
acc_data_in %>% 
  filter(BurstID == 39365) %>%
  gather(axis, acc, DynX, DynY, DynZ) %>%
  ggplot(aes(Time, acc)) +
  geom_line() +
  facet_wrap(~ axis) +
  theme_bw() +
  labs(x = 'Time (s)',
       y = 'Acceleration (m/s^2)')
ggsave('analysis/figures/Flapper.png',
       height = 4,
       width = 6,
       units = 'in')

acc_data_out <- get.acc(6122) %>%
  group_by(BurstID) %>%
  mutate(Time = as.numeric(ACCTimestampUTC - min(ACCTimestampUTC), unit = 'secs')) %>%
  ungroup
acc_data_in %>%
  mutate(Leg = 'In',
         x = row_number()) %>%
  rbind(acc_data_out %>%
          mutate(Leg = 'Out',
                 x = row_number())) %>%
  mutate(Leg = factor(Leg, levels = c('Out', 'In'))) %>%
  filter(between(x, 500, 900)) %>%
  ggplot(aes(x, DynY, color = Leg)) +
  geom_line() +
  facet_wrap(~ Leg) +
  scale_color_manual(values = c('orange', 'green')) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = '',
       y = 'Y-axis Acceleration (m/s^2)')
ggsave('analysis/figures/ACCtwosegs.png',
       height = 4,
       width = 6,
       units = 'in')
