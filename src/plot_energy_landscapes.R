rescaled_kpc <- dir('data/out/EnergyLandscapes/KPC', pattern = '.tif', full.names = TRUE) %>%
  lapply(raster, crs = hi_aea_prj) %>%
  lapply(function(r) projectRaster(r, crs = wgs84_prj)) %>%
  lapply(function(r) {
    r_min = cellStats(r, "min")
    r_max = cellStats(r, "max")
    (r - r_min) / (r_max - r_min)
  })

raster_names_kpc <- dir('data/out/EnergyLandscapes/KPC', pattern = '.tif') %>%
  sub('.tif', '', ., fixed = TRUE) %>%
  paste('KPC', ., sep = '_')

rescaled_leh <- dir('data/out/EnergyLandscapes/KPC', pattern = '.tif', full.names = TRUE) %>%
  lapply(raster, crs = hi_aea_prj) %>%
  lapply(function(r) projectRaster(r, crs = wgs84_prj)) %>%
  lapply(function(r) {
    r_min = cellStats(r, "min")
    r_max = cellStats(r, "max")
    (r - r_min) / (r_max - r_min)
  })

raster_names_leh <- dir('data/out/EnergyLandscapes/LEH', pattern = '.tif') %>%
  sub('.tif', '', ., fixed = TRUE) %>%
  paste('LEH', ., sep = '_')

rescaled_mcb <- dir('data/out/EnergyLandscapes/MCB', pattern = '.tif', full.names = TRUE) %>%
  lapply(raster, crs = hi_aea_prj) %>%
  lapply(function(r) projectRaster(r, crs = wgs84_prj)) %>%
  lapply(function(r) {
    r_min = cellStats(r, "min")
    r_max = cellStats(r, "max")
    (r - r_min) / (r_max - r_min)
  })

raster_names_mcb <- dir('data/out/EnergyLandscapes/MCB', pattern = '.tif') %>%
  sub('.tif', '', ., fixed = TRUE) %>%
  paste('MCB', ., sep = '_')

colonies <- data.frame(label = c('KPC', 'LEH', 'MCB'),
                       lat = c(22.2, 22.0, 21.5),
                       lon = c(-159.4, -160.1, -157.7)) 
row.names(colonies) <- colonies$label

plot_r <- function(r, origin, dir) {
  # Utility function for plotting a raster in ggplot
  fortify_raster <- function(r) {
    data.frame(i = seq(ncell(r))) %>%
      mutate(x = xFromCell(r, i),
             y = yFromCell(r, i),
             val = getValues(r))
  }
  
  degree_labels <- function(breaks) {
    breaks_char <- as.character(breaks)
    do.call(expression, lapply(breaks_char, function(b) bquote(.(b)*degree)))
  }
  
  hi_land2 <- hi_land %>%
    spTransform(wgs84_prj) %>%
    crop(r)
  
  fill_lab <- bquote(italic(E[.(dir)]))
  
  ggplot(fortify_raster(r),
         aes(x, y, fill = val)) +
    geom_raster() +
    scale_fill_gradientn(colors = colorRamps::matlab.like(4),
                         na.value = '#00000000') +
    annotate(geom = 'point', origin[1], origin[2], color = 'black', size = 4) +
    geom_polygon(aes(long, lat, group = group),
                 fortify(hi_land2),
                 inherit.aes = FALSE) +
    coord_fixed() +
    scale_x_continuous(labels = degree_labels,
                       breaks = function(lim) seq(ceiling(lim[1]), floor(lim[2]), by = 1)) +
    scale_y_continuous(labels = degree_labels) +
    labs(x = NULL, y = NULL, fill = fill_lab) +
    theme_bw()
}

plot_r(rescaled[[1]], as.numeric(colonies[1,3:2]), 'in')
plot_r(rescaled[[2]], as.numeric(colonies[1,3:2]), 'out')
plot_r(rescaled[[3]], as.numeric(colonies[1,3:2]), 'rt')

foreach(r = rescaled[55:57], n = raster_names[55:57]) %do% {
  colony <- substr(n, 1, 3)
  dir <- sub('.*_.*_(.*)', '\\1', n)
  origin <- as.numeric(colonies[colony, 3:2])
  p <- plot_r(r, origin, dir)
  ggsave(sprintf('analysis/figures/EnergyLandscapes/%s.png', n),
         height = 3,
         width = 3.2,
         units = 'in')
}

foreach(r = rescaled_leh[55:57], n = raster_names_leh[55:57]) %do% {
  colony <- substr(n, 1, 3)
  dir <- sub('.*_.*_(.*)', '\\1', n)
  origin <- as.numeric(colonies[colony, 3:2])
  p <- plot_r(r, origin, dir)
  ggsave(sprintf('analysis/figures/EnergyLandscapes/%s.png', n),
         height = 3,
         width = 3.2,
         units = 'in')
}

foreach(r = rescaled_mcb[55:57], n = raster_names_mcb[55:57]) %do% {
  colony <- substr(n, 1, 3)
  dir <- sub('.*_.*_(.*)', '\\1', n)
  origin <- as.numeric(colonies[colony, 3:2])
  p <- plot_r(r, origin, dir)
  ggsave(sprintf('analysis/figures/EnergyLandscapes/%s.png', n),
         height = 3,
         width = 3.2,
         units = 'in')
}
