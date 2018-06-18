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
connect_neighbors <- function(grid_coords, res) {
  grid_dist <- function(id1, id2) {
    id1 <- as.numeric(id1)
    id2 <- as.numeric(id2)
    with(grid_coords, sqrt((x[id1] - x[id2])^2 + (y[id1] - y[id2])^2))
  }
  neighbor_dist <- sqrt(2 * res^2 - 2 * res^2 * cos(2 * pi / 3))
  grid_coords2 <- rename(grid_coords, id2 = id, x2 = x, y2 = y, i2 = i, j2 = j)
  tol <- 1.1
  expand.grid(id = grid_coords$id,
              id2 = grid_coords$id,
              stringsAsFactors = FALSE) %>%
    mutate(dist = grid_dist(id, id2)) %>% 
    filter(dist <= neighbor_dist * tol,
           dist > 0) %>% 
    left_join(grid_coords) %>% 
    left_join(grid_coords2)
}

# Annotate connections with wind conditions and trip predictions
annotate_wind <- function(connections, u, v, e_mod, barriers) {
  extract_line <- function(x1, y1, x2, y2, r) {
    len <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
    npoints <- ceiling(len / xres(r) + 1)
    interp_seg <- function(x1, y1, x2, y2, n) {
      list(x = seq(x1, x2, length.out = n),
           y = seq(y1, y2, length.out = n))
    }
    cbind(x1, y1, x2, y2, npoints) %>% 
      apply(1, function(row) {
        tryCatch({
          ptlist <- interp_seg(row[1], row[2], row[3], row[4], row[5])
          raster::extract(r, cbind(ptlist$x, ptlist$y)) %>%
            mean(na.rm = TRUE)
        }, error = function(e) browser())
      })
  }
  
  vector_angle <- function(ux, uy, vx, vy) {
    dot <- function(u, v) sum(u * v)
    norm <- function(u) sqrt(sum(u^2))
    angle <- function(u, v) acos(dot(u, v) / (norm(u) * norm(v))) %>%
      round(digits = 3)
    result <- apply(cbind(ux, uy, vx, vy), 1,
                    function(row) angle(row[1:2], row[3:4]))
    result[is.nan(result)] <- 0
    result
  }

  if(!is.null(barriers)) {
    sp_segments <- foreach(segment = iterators::iter(connections, by = 'row')) %do% {
      line <- with(segment, Line(rbind(c(x, y), c(x2, y2))))
      lines <- with(segment, Lines(list(line), ID = sprintf('%s_%s', id, id2)))
    } %>%
      SpatialLines()
    
    sp_segments_df <- SpatialLinesDataFrame(sp_segments, 
                                            select(connections, id, id2),
                                            match.ID = FALSE)
    sp_segments_df@data <- as.data.frame(sp_segments_df@data)
    
    projection(sp_segments_df) <- projection(barriers)
    seg_barr <- intersect(sp_segments_df, barriers)
    blocked_segs <- as.data.frame(seg_barr) %>%
      select(id, id2) %>%
      mutate(length = SpatialLinesLengths(seg_barr)) %>%
      filter(length >= 0.5 * connections$dist[1]) %>% 
      mutate(blocked = TRUE) %>%
      select(id, id2, blocked)
  } else {
    blocked_segs <- data.frame(id = character(0), id2 = character(0), 
                               blocked = logical(0), stringsAsFactors = FALSE)
  }
  
  connections %>%
    anti_join(blocked_segs, by = c('id', 'id2')) %>%
    mutate(mean_u = extract_line(x, y, x2, y2, u),
           mean_v = extract_line(x, y, x2, y2, v),
           angle = vector_angle(x2 - x, y2 - y, mean_u, mean_v),
           mag = sqrt(mean_u^2 + mean_v^2),
           energy = e_mod(angle, mag, dist)) %>% 
    filter(!is.nan(mean_u),
           !is.nan(mean_v))
}

# Calculate trip costs based on movement model
trip_costs <- function(grid, conns, origin) {
  origin_id <- filter(grid, x == origin[1], y == origin[2])$id
  
  energy_graph <- transmute(conns, id, id2,
                            weight = energy) %>%
    igraph::graph_from_data_frame(directed = TRUE,
                                  vertices = grid)
  
  mutate(grid,
         out_cost = igraph::distances(energy_graph,
                                      origin_id,
                                      mode = 'out')[1,],
         in_cost = igraph::distances(energy_graph,
                                     origin_id,
                                     mode = 'in')[1,],
         rt_cost = out_cost + in_cost)
}

# Interpolate grid to raster
interp_grid <- function(grid, origin, radius, res, barriers) {
  landscape_mask <- data.frame(x = origin[1], y = origin[2]) %>%
    SpatialPoints %>%
    rgeos::gBuffer(width = radius)
  landscape_template <- raster(landscape_mask, res = res)
  
  
  gstat_data <- select(grid, x, y, out_cost, in_cost, rt_cost) %>%
    filter(is.finite(out_cost), 
           is.finite(in_cost),
           is.finite(rt_cost)) %>%
    na.omit
  rasters <- foreach(dir = c('out_cost', 'in_cost', 'rt_cost')) %do% {
    gstat_form <- as.formula(sprintf('%s ~ 1', dir))
    idw_mod <- gstat::gstat(formula = gstat_form,
                            locations = ~ x + y,
                            data = gstat_data)
    result <- interpolate(landscape_template, idw_mod) %>%
      mask(landscape_mask)
    if(!is.null(barriers))
      result <- mask(result, barriers, inverse = TRUE)
    result
  }
  names(rasters) <- c('out_cost', 'in_cost', 'rt_cost')
  rasters
}

# Calculate residual of energy landscape ~ distance to colony
d2c_resid <- function(el, origin) {
  lapply(el, function(el) {
    d2c <- distanceFromPoints(el, cbind(origin[1], origin[2])) %>%
      mask(el)
    el_d2c_lm <- lm(getValues(el) ~ getValues(d2c), na.action = na.exclude)
    result <- el
    result[] <- resid(el_d2c_lm) / predict(el_d2c_lm)
    result
  })
}

# Estimate energy landscape
energy_landscape <- function(origin, radius, res, u, v, 
                             energy_mod, barriers = NULL) {
  grid_coords <- create_grid(origin, radius, res)
  connections <- connect_neighbors(grid_coords, res)
  wind_conns <- annotate_wind(connections, u, v, energy_mod, barriers)
  grid_costs <- trip_costs(grid_coords, wind_conns, origin)
  el_raster <- interp_grid(grid_costs, origin, radius, res, barriers)
  el_resid <- d2c_resid(el_raster, origin)
}

plot_energy_landscape <- function(energy_landscape, origin, 
                                  dir = 'rt_cost', barriers = NULL) {
  # Utility function for plotting a raster in ggplot
  fortify_raster <- function(r) {
    data.frame(i = seq(ncell(r))) %>%
      mutate(x = xFromCell(r, i),
             y = yFromCell(r, i),
             val = getValues(r))
  }
  
  plot_raster <- energy_landscape[[dir]]
  barriers2 <- crop(barriers, plot_raster)
  
  ggplot(fortify_raster(plot_raster),
         aes(x, y, fill = val)) +
    geom_raster() +
    scale_fill_gradientn(colors = colorRamps::matlab.like(4),
                         na.value = '#00000000') +
    annotate(geom = 'point', origin[1], origin[2], color = 'black', size = 4) +
    if(!is.null(barriers2)) {
      geom_polygon(aes(long, lat, group = group),
                   fortify(barriers2),
                   inherit.aes = FALSE)
    } else {
      NULL
    } +
    coord_fixed()
}
