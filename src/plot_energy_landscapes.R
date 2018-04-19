colonies <- data.frame(label = c('KPC', 'LEH', 'MCB'),
                       lat = c(22.2, 22.0, 21.5),
                       lon = c(-159.4, -160.1, -157.7)) 
row.names(colonies) <- colonies$label

plot_r <- function(r, origin, dir) {
  # Re-project R to WGS84
  projection(r) <- hi_aea_prj
  r <- projectRaster(r, crs = wgs84_prj)
  
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
    scale_fill_gradientn(colors = rev(colorRamps::matlab.like(4)),
                         na.value = '#00000000',
                         limits = 0:1) +
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

# Plot each dir, each location for 2016-06-15
foreach(col_loc = iterators::iter(colonies, by = 'row')) %do%{
  foreach(dir = c('out', 'in', 'rt')) %do% {
    p <- sprintf('data/out/EnergyLandscapes2/all/%s_20160615_%s.tif', 
                 col_loc$label, dir) %>%
      raster %>%
      plot_r(as.numeric(col_loc[1,3:2]), dir)
    ggsave(sprintf('analysis/figures/EnergyLandscapes/%s_20160615_%s.png', 
                   col_loc$label, dir),
           height = 3,
           width = 3.2,
           units = 'in')
  }
}
