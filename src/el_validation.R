library(adehabitatHR)
library(tidyverse)
library(RSQLite)
library(rgeos)
library(rgdal)

# OBSERVED DISTRIBUTION
## Load KPC tracks
kpc_tracks <- src_sqlite('data/MHI_GPS.sqlite') %>%
  tbl('TidyTracks') %>%
  filter(DistToCol <= 250e3) %>%
  collect
## Generate utilization distribution
### Define projections (WGS84 and Hawaii Albers Equal Area Conic)
wgs84_prj <- CRS('+proj=longlat +datum=WGS84')
hi_aea_prj <- CRS('+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-157 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
### Define extent
kpc_col <- data.frame(x = -159.3997, y = 22.22888)
kpc_extent <- SpatialPoints(kpc_col,
                            proj4string = wgs84_prj) %>%
  spTransform(hi_aea_prj) %>%
  gBuffer(width = 250e3)
### Create UD from tracks
kpc_sp <- kpc_tracks %>%
  # Select longitude, latitude columns
  select(x = Longitude,
         y = Latitude) %>%
  # Create SpatialPoints object
  SpatialPoints(proj4string = wgs84_prj) %>%
  # Project to Hawaii Albers Equal Area Conic
  spTransform(hi_aea_prj)
kpc_ud <- kernelUD(kpc_sp, 
                   grid = 100)

# Visualize distribution and foraging range
coasts <- readOGR('data/coastline/ne_10m_coastline', 'ne_10m_coastline') %>%
  spTransform(hi_aea_prj)
hi_coast <- gIntersection(coasts, kpc_extent) %>%
  SpatialLinesDataFrame(data = coasts@data)
kpc_ud_df <- kpc_ud %>%
  raster %>%
  rasterToPoints(fun = function(x) { x > 0 }, spatial = TRUE) %>%
  as.data.frame
ggplot() +
  geom_raster(aes(x, y, fill = ud), data = kpc_ud_df) +
  geom_path(data = fortify(hi_coast),
            aes(long, lat, group = group)) +
  geom_path(data = fortify(kpc_extent), 
            aes(long, lat),
            linetype = 'dashed') +
  coord_fixed() +
  scale_fill_distiller(type = 'div', palette = 5) +
  labs(x = 'Easting',
       y = 'Northing')
