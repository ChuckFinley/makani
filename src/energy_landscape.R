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
    ungroup %>%
    right_join(grid_coords, by = c('i' = 'i', 'j' = 'j')) %>%
    right_join(grid_coords2, by = c('i2' = 'i', 'j2' = 'j')) %>%
    mutate(azimuth_out = atan2(y2 - y, x2 - x),
           azimuth_in = pi - azimuth_out,
           dist = sqrt((x2 - x)^2 + (y2 - y)^2))
}

# Annotate connections with wind conditions and trip predictions
annotate_wind <- function(connections, u, v, e_mod, t_mod, barriers) {
  extract_line <- function(x1, y1, x2, y2, r) {
    len <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
    npoints <- len / xres(r) + 1
    cbind(x1, x2, y1, y2, npoints) %>% 
      apply(1, function(row) {
        ptlist <- approx(row[1:2], row[3:4], n = row[5])
        raster::extract(r, cbind(ptlist$x, ptlist$y)) %>%
          mean(na.rm = TRUE)
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
           a_out = vector_angle(x2 - x, y2 - y, mean_u, mean_v),
           a_in = vector_angle(x - x2, y - y2, mean_u, mean_v),
           m = sqrt(mean_u^2 + mean_v^2),
           e_out = scales::rescale(e_mod(a_out, m)),
           e_in = scales::rescale(e_mod(a_in, m)),
           t_out = scales::rescale(t_mod(a_out, m, dist)),
           t_in = scales::rescale(t_mod(a_out, m, dist)))
}

# Calculate trip costs based on movement model
trip_costs <- function(grid, conns, origin) {
  origin_id <- filter(grid, x == origin[1], y == origin[2])$id
  
  e_out_graph <- transmute(conns, id, id2,
                           weight = e_out * t_out) %>%
    filter(weight > 0) %>%
    igraph::graph_from_data_frame(directed = TRUE,
                                  vertices = grid)
  e_out_costs <- igraph::distances(e_out_graph,
                                   origin_id,
                                   mode = 'out') %>%
    as.numeric
  e_out_costs[is.infinite(e_out_costs)] <- NA
  
  e_in_graph <- transmute(conns, id, id2,
                           weight = e_in * t_in) %>%
    filter(weight > 0) %>%
    igraph::graph_from_data_frame(directed = TRUE,
                                  vertices = grid)
  e_in_costs <- igraph::distances(e_in_graph,
                                  origin_id,
                                  mode = 'in') %>%
    as.numeric
  e_in_costs[is.infinite(e_in_costs)] <- NA
  
  mutate(grid,
         out_cost = e_out_costs,
         in_cost = e_in_costs,
         rt_cost = out_cost + in_cost)
}

# Interpolate grid to raster
interp_grid <- function(grid, origin, radius, res, barriers) {
  landscape_mask <- data.frame(x = origin[1], y = origin[2]) %>%
    SpatialPoints %>%
    rgeos::gBuffer(width = radius)
  landscape_template <- raster(landscape_mask, res = res)
  
  rasters <- foreach(dir = c('out_cost', 'in_cost', 'rt_cost')) %do% {
    gstat_form <- as.formula(sprintf('%s ~ 1', dir))
    gstat_data <- select(grid, x, y, out_cost, in_cost, rt_cost) %>%
      na.omit
    idw_mod <- gstat::gstat(formula = gstat_form,
                            locations = ~ x + y,
                            data = gstat_data)
    result <- interpolate(landscape_template, idw_mod) %>%
      mask(landscape_mask)
    if(is.null(barriers))
      result
    else
      mask(result, barriers, inverse = TRUE)
  }
  names(rasters) <- c('out_cost', 'in_cost', 'rt_cost')
  rasters
}

# Rescale rasters to [0, 1]
rescale <- function(rasters) {
  lapply(rasters,
         function(r) {
           r_min = cellStats(r, "min")
           r_max = cellStats(r, "max")
           (r - r_min) / (r_max - r_min)
         })
}

# Estimate energy landscape
energy_landscape <- function(origin, radius, res, u, v, 
                             energy_mod, dur_mod, barriers = NULL) {
  grid_coords <- create_grid(origin, radius, res)
  connections <- connect_neighbors(grid_coords)
  wind_conns <- annotate_wind(connections, u, v, energy_mod, dur_mod, barriers)
  grid_costs <- trip_costs(grid_coords, wind_conns, origin)
  el_raster <- interp_grid(grid_costs, origin, radius, res, barriers)
  norm_raster <- rescale(el_raster)
}

plot_energy_landscape <- function(energy_landscape, origin, 
                                  dir = 'rt_cost', barriers = NULL) {
  rasterPoints <- rasterToPoints(energy_landscape[[dir]]) %>%
    as.data.frame
  colnames(rasterPoints) <- c('x', 'y', 'val')
  
  p <- ggplot(rasterPoints, aes(x, y, fill = val)) +
    geom_raster() +
    scale_fill_gradientn(colors = colorRamps::matlab.like(4),
                         na.value = '#00000000') +
    annotate(geom = 'point', origin[1], origin[2], color = 'black', size = 4) +
    coord_fixed()
  
  if(!is.null(barriers)) 
    p <- p + geom_polygon(aes(long, lat, group = group),
                          fortify(crop(barriers, energy_landscape[[dir]])),
                          inherit.aes = FALSE)
  
  p
}
