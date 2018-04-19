# Energy landscape variability across colonies
kpc_rt <- dir('data/out/EnergyLandscapes2/all/', 
              pattern = 'KPC_.*_rt.tif', 
              full.names = TRUE) %>%
  lapply(raster, crs = hi_aea_prj) %>%
  lapply(projectRaster, crs = wgs84_prj)
dates <- dir('data/out/EnergyLandscapes2/all/', 
             pattern = 'KPC_.*_rt.tif', 
             full.names = TRUE) %>% 
  sub('.*KPC_(.*)_rt.tif', '\\1', .)
leh_rt <- sprintf('data/out/EnergyLandscapes2/all/LEH_%s_rt.tif', 
                  dates) %>%
  lapply(raster, crs = hi_aea_prj) %>%
  lapply(projectRaster, crs = wgs84_prj)
mcb_rt <- sprintf('data/out/EnergyLandscapes2/all/MCB_%s_rt.tif', 
                  dates) %>%
  lapply(raster, crs = hi_aea_prj) %>%
  lapply(projectRaster, crs = wgs84_prj)

kpc_stack <- stack(kpc_rt)
leh_stack <- stack(leh_rt)
mcb_stack <- stack(mcb_rt)
kpc_mean <- calc(kpc_stack, fun = mean)
leh_mean <- calc(leh_stack, fun = mean)
mcb_mean <- calc(mcb_stack, fun = mean)
range(cellStats(kpc_mean, range), 
      cellStats(leh_mean, range), 
      cellStats(mcb_mean, range))
# [1] 0.04217044 0.89535627
kpc_sd <- calc(kpc_stack, fun = sd)
leh_sd <- calc(leh_stack, fun = sd)
mcb_sd <- calc(mcb_stack, fun = sd)
range(cellStats(kpc_sd, range), 
      cellStats(leh_sd, range), 
      cellStats(mcb_sd, range))
# [1] 0.07811175 0.27136536

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
  
  hi_land2 <- hi_land %>%
    spTransform(wgs84_prj) %>%
    crop(r)
  
  fill_lab <- bquote(.(stat)(italic(E[rt])))
  fill_lim <- if(stat == 'mean') c(0, 1) else c(0.07, 0.28)
  
  ggplot(fortify_raster(r),
         aes(x, y, fill = val)) +
    geom_raster() +
    scale_fill_gradientn(colors = rev(colorRamps::matlab.like(4)),
                         na.value = '#00000000',
                         limits = fill_lim) +
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

colonies <- data.frame(label = c('KPC', 'LEH', 'MCB'),
                       lat = c(22.2, 22.0, 21.5),
                       lon = c(-159.4, -160.1, -157.7)) 
row.names(colonies) <- colonies$label

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
