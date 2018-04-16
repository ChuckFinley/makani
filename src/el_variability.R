# Energy landscape variability across colonies
kpc_rescaled <- dir('data/out/EnergyLandscapes/KPC', pattern = '_rt.tif', full.names = TRUE) %>%
  lapply(raster, crs = hi_aea_prj) %>%
  lapply(function(r) projectRaster(r, crs = wgs84_prj)) %>%
  lapply(function(r) {
    r_min = cellStats(r, "min")
    r_max = cellStats(r, "max")
    (r - r_min) / (r_max - r_min)
  })
leh_rescaled<- dir('data/out/EnergyLandscapes/LEH', pattern = '_rt.tif', full.names = TRUE) %>%
  lapply(raster, crs = hi_aea_prj) %>%
  lapply(function(r) projectRaster(r, crs = wgs84_prj)) %>%
  lapply(function(r) {
    r_min = cellStats(r, "min")
    r_max = cellStats(r, "max")
    (r - r_min) / (r_max - r_min)
  })
mcb_rescaled<- dir('data/out/EnergyLandscapes/MCB', pattern = '_rt.tif', full.names = TRUE) %>%
  lapply(raster, crs = hi_aea_prj) %>%
  lapply(function(r) projectRaster(r, crs = wgs84_prj)) %>%
  lapply(function(r) {
    r_min = cellStats(r, "min")
    r_max = cellStats(r, "max")
    (r - r_min) / (r_max - r_min)
  })

kpc_stack <- stack(kpc_rescaled)
leh_stack <- stack(leh_rescaled)
mcb_stack <- stack(mcb_rescaled)
kpc_mean <- calc(kpc_stack, fun = mean)
leh_mean <- calc(leh_stack, fun = mean)
mcb_mean <- calc(mcb_stack, fun = mean)
range(cellStats(kpc_mean, range), 
      cellStats(leh_mean, range), 
      cellStats(mcb_mean, range))
kpc_sd <- calc(kpc_stack, fun = sd)
leh_sd <- calc(leh_stack, fun = sd)
mcb_sd <- calc(mcb_stack, fun = sd)
range(cellStats(kpc_sd, range), 
      cellStats(leh_sd, range), 
      cellStats(mcb_sd, range))

plot_r <- function(r, origin, stat) {
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
  
  fill_lab <- bquote(.(stat)(italic(E[rt])))
  fill_lim <- if(stat == 'sd') {
    c(0.07, 0.28)
  } else {
    c(0.1, 1)
  }
  
  hi_land2 <- hi_land %>%
    spTransform(wgs84_prj) %>%
    crop(r)
  
  ggplot(fortify_raster(r),
         aes(x, y, fill = val)) +
    geom_raster() +
    annotate(geom = 'point', origin[1], origin[2], color = 'black', size = 4) +
    geom_polygon(aes(long, lat, group = group),
                 fortify(hi_land2),
                 inherit.aes = FALSE) +
    coord_fixed() +
    scale_x_continuous(labels = degree_labels,
                       breaks = function(lim) seq(ceiling(lim[1]), floor(lim[2]), by = 1)) +
    scale_y_continuous(labels = degree_labels,
                       breaks = function(lim) seq(ceiling(lim[1]), floor(lim[2]), by = 1)) +
    labs(x = NULL, y = NULL, fill = fill_lab) +
    theme_bw() +
    scale_fill_gradientn(colors = colorRamps::matlab.like(4),
                         na.value = '#00000000',
                         limits = fill_lim)
}

p_kpc_mean <- plot_r(kpc_mean, as.numeric(colonies['KPC', 3:2]), 'mean')
p_leh_mean <- plot_r(leh_mean, as.numeric(colonies['LEH', 3:2]), 'mean')
p_mcb_mean <- plot_r(mcb_mean, as.numeric(colonies['MCB', 3:2]), 'mean')

p_kpc_sd <- plot_r(kpc_sd, as.numeric(colonies['KPC', 3:2]), 'sd')
p_leh_sd <- plot_r(leh_sd, as.numeric(colonies['LEH', 3:2]), 'sd')
p_mcb_sd <- plot_r(mcb_sd, as.numeric(colonies['MCB', 3:2]), 'sd')

ggsave('analysis/figures/KPC_mean.png', p_kpc_mean,
       width = 5, height = 3, units ='in')
ggsave('analysis/figures/LEH_mean.png', p_leh_mean,
       width = 5, height = 3, units ='in')
ggsave('analysis/figures/MCB_mean.png', p_mcb_mean,
       width = 5, height = 3, units ='in')

ggsave('analysis/figures/KPC_sd.png', p_kpc_sd,
       width = 5, height = 3, units ='in')
ggsave('analysis/figures/LEH_sd.png', p_leh_sd,
       width = 5, height = 3, units ='in')
ggsave('analysis/figures/MCB_sd.png', p_mcb_sd,
       width = 5, height = 3, units ='in')
