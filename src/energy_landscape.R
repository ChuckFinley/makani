library(tidyverse)
library(sp)
library(igraph)

select <- dplyr::select

# Create hexagonal grid 
create_grid <- function(origin, radius, spacing) {
  if(length(origin) != 2 || class(origin) != 'numeric')
    stop('origin must be x,y coordinates (numeric vector of length 2)')
  if(length(radius) != 1 || class(radius) != 'numeric')
    stop('radius must be a single number')
  if(length(spacing) != 1 || class(spacing) != 'numeric')
    stop('spacing must be a single number')
  
  o <- origin
  l <- spacing
  L <- 2 * radius
  h <- l * sin(pi / 3)
  H <- 2 * radius
  
  i <- seq(floor(L / l) + 1) - 1
  j <- seq(floor(H / h) + 1) - 1
  
  coords <- expand.grid(i = i, j = j) %>%
    mutate(x = i * l + (j %% 2) * l / 2,
           y = j * h,
           d2o = sqrt((x - o[1])^2 + (y - o[2])^2))
  
  o_error <- filter(coords, d2o == min(d2o)) %>%
    mutate(dx = x - o[1],
           dy = y - o[2])
  
  mutate(coords, 
         x = x - o_error$dx,
         y = y - o_error$dy,
         d2o = sqrt((x - o[1])^2 + (y - o[2])^2)) %>%
    filter(d2o <= radius) %>%
    mutate(id = row_number()) %>%
    select(id, i, j, x, y)
}

origin <- c(30, 30)
radius <- 30
spacing <- 15
grid_coords <- create_grid(origin, radius, spacing)
circ_coords <- data.frame(i = seq(0, 2*pi, length.out = 100)) %>%
  mutate(x = radius * cos(i) + origin[1],
         y = radius * sin(i) + origin[2])
ggplot(grid_coords, aes(x, y)) + 
  geom_text(aes(label = sprintf('(%i, %i)', i, j))) + 
  annotate('point', origin[1], origin[2], color = 'red') + 
  geom_path(data = circ_coords, linetype = 'dashed') + 
  coord_fixed()

# Connect neighboring grid points
connect_neighbors <- function(grid_coords) {
  neighbors <- function(i, j) {
    # left-to-right, top-to-bottom
    x_8 <- c(-1, 0, 1, -1, 1, -1, 0 , 1)
    y_8 <- c(1, 1, 1, 0, 0, -1, -1, -1)
    if((j %% 2) == 0) {
      x_6 <- x_8[c(-3, -8)]
      y_6 <- y_8[c(-3, -8)]
    } else {
      x_6 <- x_8[c(-1, -6)]
      y_6 <- y_8[c(-1, -6)]
    }
    data.frame(i = i,
               j = j,
               i2 = i + x_6,
               j2 = j + y_6)
  }
  
  grid_coords2 <- rename(grid_coords, id2 = id, x2 = x, y2 = y)
  
  grid_coords %>%
    rowwise %>%
    do(neighbors(.$i, .$j)) %>%
    right_join(grid_coords, by = c('i' = 'i', 'j' = 'j')) %>%
    right_join(grid_coords2, by = c('i2' = 'i', 'j2' = 'j')) %>%
    mutate(azimuth_out = atan2(y2 - y, x2 - x),
           azimuth_in = pi - azimuth_out,
           distance = sqrt((x2 - x)^2 + (y2 - y)^2))
}

connections <- connect_neighbors(grid_coords)
ggplot(grid_coords, aes(x, y)) + 
  geom_text(aes(label = sprintf('(%i, %i)', i, j))) + 
  annotate('point', 30, 30, color = 'red') + 
  geom_path(data = circ_coords, linetype = 'dashed') +
  geom_segment(aes(x, y, xend = x2, yend = y2), 
               data = connections,
               arrow = arrow(length = unit(0.03, 'npc'))) +
  coord_fixed()

# dummy wind
wind_u <- matrix(seq(from = -5, to = 10, length.out = 100), nrow = 10) %>%
  raster(xmn = 0, xmx = 60, ymn = 0, ymx = 60)
wind_v <- matrix(seq(from = 5, to = -10, length.out = 100), nrow = 10) %>%
  raster(xmn = 0, xmx = 60, ymn = 0, ymx = 60)

load('data/out/Models/EnergyModels.Rdata')
load('data/out/Models/SgTCModel.RData')

pred_spd <- function(t) {
  predict(SgTCModel, 
          re.form = NA,
          newdata = data.frame(Tailwind = t))
}
pred_energy <- function(a, m) {
  predict(fam, 
          newdata = data.frame(WindAngle = a,
                               WindSpd = m,
                               DeployID = 1145),
          exclude = 's(DeployID)',
          type = 'response')
}

# Annotate connections with wind conditions and trip predictions
annotate_wind <- function(connections, u, v, pred_spd, pred_energy) {
  sample_raster <- function(x1, y1, x2, y2, r) {
    len <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
    npoints <- len / xres(r) + 1
    points <- approx(c(x1, x2), c(y1, y2), n = npoints)
    cbind(points$x, points$y) %>%
      cellFromXY(r, .) %>%
      raster::extract(r, .) %>%
      mean
  }
  
  vector_angle <- function(u, v) {
    dot <- function(u, v) sum(u * v)
    norm <- function(u) sqrt(sum(u^2))
    acos(dot(u, v) / (norm(u) * norm(v))) %>%
      round(digits = 3)
  }
  
  mutate(connections,
         mean_u = mapply(FUN = sample_raster, 
                         x, y, x2, y2, 
                         MoreArgs = list(r = wind_u)),
         mean_v = mapply(FUN = sample_raster, 
                         x, y, x2, y2, 
                         MoreArgs = list(r = wind_v)),
         wind_angle_out = vector_angle(c(x2 - x, y2 - y), c(mean_u, mean_v)),
         wind_angle_in = vector_angle(c(x - x2, y - y2), c(mean_u, mean_v)),
         wind_speed = sqrt(mean_u^2 + mean_v^2),
         tailwind_out = wind_speed * cos(wind_angle_out),
         tailwind_in = wind_speed * cos(wind_angle_in),
         flight_speed_out = pred_spd(tailwind_out),
         duration_out = flight_speed_out / distance,
         flap_ratio_out = pred_energy(wind_angle_out, wind_speed),
         flap_out = duration_out * flap_ratio_out / 60,
         flight_speed_in = pred_spd(tailwind_in),
         duration_in = flight_speed_in / distance,
         flap_ratio_in = pred_energy(wind_angle_in, wind_speed),
         flap_in = duration_in * flap_ratio_in / 60)
}

wind_conns <- annotate_wind(connections, wind_u, wind_v, pred_spd, pred_energy)

# Calculate trip costs
## NOTE: uses duration as weight. This is hard-coded and needs future editing.
trip_costs <- function(grid, conns, origin) {
  origin_id <- filter(grid, x == origin[1], y == origin[2])$id
  
  out_mat <- transmute(conns, id, id2,
                       weight = duration_out) %>%
    igraph::graph_from_data_frame(directed = TRUE,
                                  vertices = grid) %>%
    igraph::distances(graph = .,
                      mode = 'out', 
                      weights = NULL, 
                      algorithm = 'dijkstra')
  out_cost <- out_mat[origin_id, ]
  
  in_mat <- transmute(conns, id, id2,
                      weight = duration_in) %>%
    igraph::graph_from_data_frame(directed = TRUE,
                                  vertices = grid_coords) %>%
    igraph::distances(graph = .,
                      mode = 'in', 
                      weights = NULL, 
                      algorithm = 'dijkstra')
  in_cost <- in_mat[origin_id, ]
  
  mutate(grid,
         out_cost = out_cost,
         in_cost = in_cost,
         rt_cost = out_cost + in_cost)
}

grid_costs <- trip_costs(grid_coords, wind_conns, origin)

landscape_mask <- data.frame(x = origin[1], y = origin[2]) %>%
  SpatialPoints %>%
  rgeos::gBuffer(width = radius)
landscape_template <- raster(landscape_mask, res = res(wind_u))
idw <- interpolate(landscape_template, 
                   gstat::gstat(formula = rt_cost ~ 1,
                                locations = ~ x + y,
                                data = grid_costs))

ggplot(grid_costs, aes(x, y, color = rt_cost)) +
  geom_point(size = 4) +
  scale_color_gradientn(colors = colorRamps::matlab.like(4))

fortify_raster <- function(r) {
  data.frame(i = seq(ncell(r))) %>%
    mutate(x = xFromCell(r, i),
           y = yFromCell(r, i),
           val = getValues(r))
}

ggplot(grid_costs) +
  geom_raster(aes(x, y, fill = val), data = fortify_raster(idw)) +
  scale_fill_gradientn(colors = colorRamps::matlab.like(4)) +
  geom_point(aes(x, y, color = rt_cost), size = 4) +
  scale_color_gradientn(colors = colorRamps::matlab.like(4))
  