# setwd("~/Desktop/Turbulence/Another_analysis_for_review")
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Turbulence/Another_analysis_for_review")
load("./xz_all.RData")

stations <- expand.grid(x=1:37, z=1:600)
# xz_allt: 37 x 600 x 900
xz_mat <- matrix(xz_all, nrow = 37*600, ncol = 900)   # 22200 x 900
rm(xz_all)

library(ggplot2)
library(ggh4x)
pal <- rev(RColorBrewer::brewer.pal(11,"Spectral"))
for(tt in 1:ncol(xz_mat)){
  ggplot(stations) + geom_raster(aes(x = x, y = z, fill = xz_mat[,tt])) +
    scale_fill_gradientn(colours = pal, name=paste("Time", tt)) +
    scale_x_continuous(breaks=seq(0, 37,length.out = 3),labels=seq(0, 148,length.out = 3)) +
    scale_y_continuous(breaks=seq(0, 600,length.out = 6),labels=seq(0, 2400,length.out = 6)) + labs(x='x (m)', y='z (m)') +
    force_panelsizes(rows = unit(4.5, "in"),
                     cols = unit(1.5, "in"))
  ggsave(paste0("./Figures/fire",tt+100,".png"), width = 6, height = 4)
}

library(magick)
library(magrittr)
list.files(path='./Figures/', pattern = '*.png', full.names = TRUE) %>%
  image_read() %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=10) %>% # animates, can opt for number of loops
  image_write("./fire.gif") # write to current dir



xz_mat <- xz_mat + 1.904
write.csv(xz_mat, file="./xz_all.csv")
