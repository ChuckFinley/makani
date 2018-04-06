library(tidyverse)
library(sp)

energy_landscape <- function(mvmt_mod, wind_uv, col_loc, for_rng, grd_sz) {
  
}

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
    select(i, j, x, y)
}

grid_coords <- create_grid(c(30, 30), 30, 15)
ggplot(grid_coords, aes(x, y)) + 
  geom_text(aes(label = sprintf('(%i, %i)', i, j))) + 
  annotate('point', 30, 30, color = 'red') + 
  geom_path(data = circ_coords, linetype = 'dashed') + 
  coord_fixed()

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
  
  grid_coords2 <- rename(grid_coords, x2 = x, y2 = y)
  
  grid_coords %>%
    rowwise %>%
    do(neighbors(.$i, .$j)) %>%
    right_join(grid_coords, by = c('i' = 'i', 'j' = 'j')) %>%
    right_join(grid_coords2, by = c('i2' = 'i', 'j2' = 'j'))
}

neighbors <- connect_neighbors(grid_coords)
ggplot(grid_coords, aes(x, y)) + 
  geom_text(aes(label = sprintf('(%i, %i)', i, j))) + 
  annotate('point', 30, 30, color = 'red') + 
  geom_path(data = circ_coords, linetype = 'dashed') +
  geom_segment(aes(x, y, xend = x2, yend = y2), 
               data = neighbors,
               arrow = arrow(length = unit(0.03, 'npc'))) +
  coord_fixed()

