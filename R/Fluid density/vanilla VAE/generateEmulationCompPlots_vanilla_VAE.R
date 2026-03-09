
load("../../../data/Fluid density/center_fire_input.RData")
ncomp <- 50
load(paste0("./emulation_vanilla_VAE_50_ncomp_",ncomp,".RData"))
fire_vanilla_VAE_50 <- fire_vanilla_VAE
load(paste0("./emulation_vanilla_VAE__ncomp_",ncomp,".RData"))
stations <- expand.grid(x=1:198, z=1:500)

library(ggplot2)
library(ggh4x)
pal <- RColorBrewer::brewer.pal(11,"Spectral")
for(tt in c(3, 40, 90)){
  tmp_dat <- data.frame(x=stations$x, z=stations$z, fire = center_fire[,tt], model="Fluid density")
  tmp_dat <- rbind(tmp_dat, data.frame(x=stations$x, z=stations$z, fire = fire_vanilla_VAE_50[,tt], model="VAE, β = 0.5"))
  tmp_dat <- rbind(tmp_dat, data.frame(x=stations$x, z=stations$z, fire = fire_vanilla_VAE[,tt], model="VAE, β = 1"))
  
  if(abs(min(center_fire[,tt]))< 1e-3) {
    brks <- c(-2e-04, -1e-04); labels = c('-2e-04', '-1e-04')}else{
      brks <- seq(-4e-03, 0, length.out = 5); labels = c('-4e-03', '-3e-03', '-2e-03', '-1e-03',"0")
    }
  ggplot(tmp_dat) + geom_raster(aes(x = x, y = z, fill = fire)) +
    scale_fill_gradientn(colours = pal, name=paste("Time", tt), breaks = brks, labels = labels, limits = c(min(c(center_fire[,tt], fire_approx[,tt])),-min(center_fire[,tt])/40)) +
    scale_x_continuous(expand = c(0, 0), breaks=seq(0, 200,length.out = 5),labels=seq(0, 8000,length.out = 5)) +
    scale_y_continuous(expand = c(0, 0), breaks=seq(0, 500,length.out = 6),labels=seq(0, 5000,length.out = 6)) + labs(x='x (m)', y='z (m)') +
    force_panelsizes(rows = unit(2.4, "in"),
                     cols = unit(2.44, "in")) +
    facet_grid(cols = vars(model)) +
    theme(strip.text.x = element_text(size = 13), axis.title = element_text(size = 12), axis.text = element_text(size = 12))
  
  ggsave(paste0("./Figures/fire",tt+100,"_vanilla_VAE.png"), width = 12, height = 4)
}
