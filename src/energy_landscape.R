library(tidyverse)
library(sp)
library(igraph)
library(dplyr)
library(foreach)

select <- dplyr::select

# Create hexagonal grid 
create_grid <- function(origin, radius, res) {
  if(length(origin) != 2 || class(origin) != 'numeric')
    stop('origin must be x,y coordinates (numeric vector of length 2)')
  if(length(radius) != 1 || class(radius) != 'numeric')
    stop('radius must be a single number')
  if(length(res) != 1 || class(res) != 'numeric')
    stop('res must be a single number')
  
  o <- origin
  l <- res
  L <- 2 * radius
  h <- l * sin(pi / 3)
  H <- 2 * radius
  x0 <- origin[1] - radius
  y0 <- origin[2] - radius
  
  i <- seq(floor(L / l) + 1) - 1
  j <- seq(floor(H / h) + 1) - 1
  
  coords <- expand.grid(i = i, j = j) %>%
    mutate(x = x0 + i * l + (j %% 2) * l / 2,
           y = y0 + j * h,
           d2o = sqrt((x - o[1])^2 + (y - o[2])^2))
  
  o_error <- filter(coords, d2o == min(d2o)) %>%
    mutate(dx = x - o[1],
           dy = y - o[2])
  
  mutate(coords, 
         x = x - o_error$dx,
         y = y - o_error$dy,
         d2o = sqrt((x - o[1])^2 + (y - o[2])^2)) %>%
    filter(d2o <= radius) %>%
    mutate(id = as.character(row_number())) %>%
    select(id, i, j, x, y)
}

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
           dist = sqrt((x2 - x)^2 + (y2 - y)^2))
}

# Annotate connections with wind conditions and trip predictions
annotate_wind <- function(connections, u, v, mvmt_mod, barriers) {
  sample_raster <- function(x1, y1, x2, y2, r) {
    len <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
    npoints <- len / xres(r) + 1
    points <- approx(c(x1, x2), c(y1, y2), n = npoints)
    cbind(points$x, points$y) %>%
      cellFromXY(r, .) %>%
      raster::extract(r, .) %>%
      mean(na.rm = TRUE)
  }
  
  vector_angle <- function(u, v) {
    dot <- function(u, v) sum(u * v)
    norm <- function(u) sqrt(sum(u^2))
    result <- acos(dot(u, v) / (norm(u) * norm(v))) %>%
      round(digits = 3)
    result[is.nan(result)] <- 0
    result
  }
  
  if(!is.null(barriers)) {
    foreach(segment = iterators::iter(connections, by = 'row')) %do% {
      line <- with(segment, Line(rbind(c(x, y), c(x2, y2))))
      lines <- with(segment, Lines(list(line), ID = sprintf('%s_%s', id, id2)))
    } %>%
      SpatialLines() -> sp_segments
    
    seg_barr <- intersect(sp_segments, barriers)
    blocked_segs <- data.frame(ids = row.names(seg_barr),
                               length = SpatialLinesLengths(seg_barr)) %>%
      separate(ids, c('sp_segments_id', 'barriers_id'), sep = ' ') %>%
      filter(length >= 0.5 * connections$dist[1]) %>% 
      mutate(blocked = TRUE,
             id_id2 = row.names(sp_segments)[as.numeric(sp_segments_id)]) %>%
      separate(id_id2, c('id', 'id2'), sep = '_') %>%
      select(id, id2, blocked)
  } else {
    blocked_segs <- data.frame(id = character(0), id2 = character(0), 
                               blocked = logical(0), stringsAsFactors = FALSE)
  }
  
  connections %>%
    left_join(blocked_segs, by = c('id', 'id2')) %>%
    replace_na(list(blocked = FALSE)) %>%
    mutate(mean_u = mapply(FUN = sample_raster, 
                           x, y, x2, y2, 
                           MoreArgs = list(r = u)),
           mean_v = mapply(FUN = sample_raster, 
                           x, y, x2, y2, 
                           MoreArgs = list(r = v)),
           a_out = vector_angle(c(x2 - x, y2 - y), c(mean_u, mean_v)),
           a_in = vector_angle(c(x - x2, y - y2), c(mean_u, mean_v)),
           m = sqrt(mean_u^2 + mean_v^2),
           e_out = ifelse(blocked, 0, mvmt_mod(a_out, m, dist)),
           e_in = ifelse(blocked, 0, mvmt_mod(a_in, m, dist)))
}

# Calculate trip costs based on movement model
trip_costs <- function(grid, conns, origin) {
  origin_id <- filter(grid, x == origin[1], y == origin[2])$id
  
  out_mat <- transmute(conns, id, id2,
                       weight = e_out) %>%
    igraph::graph_from_data_frame(directed = TRUE,
                                  vertices = grid) %>%
    igraph::distances(graph = .,
                      mode = 'out', 
                      weights = NULL, 
                      algorithm = 'dijkstra')
  out_cost <- out_mat[origin_id, ]
  
  in_mat <- transmute(conns, id, id2,
                      weight = e_in) %>%
    igraph::graph_from_data_frame(directed = TRUE,
                                  vertices = grid) %>%
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

# Interpolate grid to raster
interp_grid <- function(grid, origin, res, dir = 'rt_cost') {
  if(!(dir %in% c('out_cost', 'in_cost', 'rt_cost'))) {
    stop('"dir" must be one of "out_cost", "in_cost", or "rt_cost"')
  }
  landscape_mask <- data.frame(x = origin[1], y = origin[2]) %>%
    SpatialPoints %>%
    rgeos::gBuffer(width = radius)
  landscape_template <- raster(landscape_mask, res = res)
  idw_mod <- gstat::gstat(formula = rt_cost ~ 1,
                          locations = ~ x + y,
                          data = filter(grid, !is.na(rt_cost)))
  interpolate(landscape_template, idw_mod) %>%
    mask(landscape_mask)
}

# Estimate energy landscape
energy_landscape <- function(origin, radius, res, u, v, mvmt_mod, dir, barriers = NULL) {
  grid_coords <- create_grid(origin, radius, res)
  connections <- connect_neighbors(grid_coords)
  wind_conns <- annotate_wind(connections, u, v, mvmt_mod, barriers)
  grid_costs <- trip_costs(grid_coords, wind_conns, origin)
  el_grid <- interp_grid(grid_costs, origin, res, dir)
  if(is.null(barriers))
    el_grid
  else
    mask(el_grid, barriers, inverse = TRUE)
}

plot_energy_landscape <- function(energy_landscape, origin) {
  # Utility function for plotting a raster in ggplot
  fortify_raster <- function(r) {
    data.frame(i = seq(ncell(r))) %>%
      mutate(x = xFromCell(r, i),
             y = yFromCell(r, i),
             val = getValues(r))
  }
  
  ggplot(fortify_raster(energy_landscape),
         aes(x, y, fill = val)) +
    geom_raster() +
    scale_fill_gradientn(colors = colorRamps::matlab.like(4)) +
    annotate(geom = 'point', origin[1], origin[2], color = 'black', size = 4) +
    coord_fixed()
}
