
center_fire <- read.csv("../../data/Vertical velocity/xz_all.csv", row.names = 1)
center_fire <- as.matrix(center_fire)
center_fire <- as.matrix(center_fire) - 1.904

load("./emulation_fire_1.RData")
fire_approx <- fire_approx - 1.904
load("./emulation_PCA.RData")
fire_PCA <- fire_PCA - 1.904
stations <- expand.grid(x=1:37, z=1:600)

library(ggplot2)
library(ggh4x)
pal <- rev(RColorBrewer::brewer.pal(11,"Spectral"))
for(tt in c(10, 500, 900)){
  tmp_dat <- data.frame(x=stations$x, z=stations$z, fire = center_fire[,tt], model="Vertical velocity")
  tmp_dat <- rbind(tmp_dat, data.frame(x=stations$x, z=stations$z, fire = fire_approx[,tt], model="xVAE emulation (K = 100)"))
  tmp_dat <- rbind(tmp_dat, data.frame(x=stations$x, z=stations$z, fire = fire_PCA[,tt], model="POD emulation (K = 200)"))
  tmp_dat$model <- factor(
    tmp_dat$model,
    levels = c(
      "Vertical velocity",
      "xVAE emulation (K = 100)",
      "POD emulation (K = 200)"
    )
  )
  
  ggplot(tmp_dat) + geom_raster(aes(x = x, y = z, fill = fire)) +
    scale_fill_gradientn(colours = pal, name=paste("Time", tt), limits = c(-1.904, 15.13)) +
    scale_x_continuous(expand = c(0, 0), breaks=seq(0, 37,length.out = 5),labels=seq(0, 4000,length.out = 5)) +
    scale_y_continuous(expand = c(0, 0), breaks=seq(0, 600,length.out = 6),labels=seq(0, 7000,length.out = 6)) + labs(x='x (m)', y='z (m)') +
    force_panelsizes(rows = unit(2.4, "in"),
                     cols = unit(1.37, "in")) +
    facet_grid(cols = vars(model)) +
    theme(strip.text.x = element_text(size = 11), axis.title = element_text(size = 10))
  
  ggsave(paste0("./Figures/fire",tt+100,"_all_three_methods_add_data.png"), width = 7.2, height = 4)
}
