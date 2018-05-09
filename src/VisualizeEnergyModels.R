library(tidyverse)
library(rgl)

load('data/out/Models/EnergyModels.Rdata')

pred_f <- function(a, s) {
  predict(fam, newdata = data.frame(WindAngle = a,
                                    WindSpd = s,
                                    DeployID = 1145), 
          exclude = 's(DeployID)', type = 'response')
}
pred_o <- function(a, s) {
  predict(oam, newdata = data.frame(WindAngle = a,
                                    WindSpd = s,
                                    DeployID = 1145), 
          exclude = 's(DeployID)', type = 'response')
}

grain <- 15
angle <- seq(0, pi, length.out = grain)
spd <- seq(5, 11, length.out = grain)
odba <- outer(angle, spd, pred_o)
flap <- outer(angle, spd, pred_f)
color_fun <- colorRamp(c('blue', 'red'))
odba_col <- scales::rescale(odba) %>%
  color_fun %>%
  apply(1, function(row) rgb(row[1], row[2], row[3], maxColorValue = 255)) %>%
  matrix(ncol = grain)
flap_col <- scales::rescale(flap) %>%
  color_fun %>%
  apply(1, function(row) rgb(row[1], row[2], row[3], maxColorValue = 255)) %>%
  matrix(ncol = grain)

open3d()
persp3d(x = angle * 180/pi,
        y = spd,
        z = odba,
        color = odba_col,
        xlab = '',
        ylab = '',
        zlab = '',
        lit = FALSE)
surface3d(x = angle * 180/pi,
          y = spd,
          z = odba, 
          front = "lines")

open3d()
persp3d(x = angle * 180/pi,
        y = spd,
        z = flap,
        color = flap_col,
        xlab = '',
        ylab = '',
        zlab = '',
        lit = FALSE)
surface3d(x = angle * 180/pi,
          y = spd,
          z = flap, 
          front = "lines")
