library(torch)

setwd("~/Desktop/Turbulence/")
source("../GEV-GP_VAE/extCVAE/utils.R")

###### ---------------------------------------------------------------------- ######
###### ----------------------------- Load in Data  -------------------------- ######
###### ---------------------------------------------------------------------- ######
stations <- expand.grid(x=1:198, z=1:500)

load("./center_fire_input.RData")


ggplot(stations) + geom_raster(aes(x=x, y=z, fill=center_fire[,1])) +
  geom_hline(yintercept = 10) +
  geom_hline(yintercept = 250) +
  geom_hline(yintercept = 490) +
  scale_fill_gradientn(colours = topo.colors(100), name = paste0('Original replicate #', 1), na.value = "transparent")

H <- 1:10
z = 10
str2 <- rep(NA, length(H))
str3 <- rep(NA, length(H))
for(h.iter in 1:length(H)){
  h = H[h.iter]
  start <- seq(1, 198, by = 2*h)
  start <- start[start + h <= 198]
  wh <- which(stations$z == z)
  str2[h.iter] <- mean((center_fire[wh[start+h],] - center_fire[wh[start],])^2)
  str3[h.iter] <- mean((abs(center_fire[wh[start+h],] - center_fire[wh[start],]))^3)
}


z = 250
str2_250 <- rep(NA, length(H))
str3_250 <- rep(NA, length(H))
for(h.iter in 1:length(H)){
  h = H[h.iter]
  start <- seq(1, 198, by = 2*h)
  start <- start[start + h <= 198]
  wh <- which(stations$z == z)
  str2_250[h.iter] <- mean((center_fire[wh[start+h],] - center_fire[wh[start],])^2)
  str3_250[h.iter] <- mean((abs(center_fire[wh[start+h],] - center_fire[wh[start],]))^3)
}



z = 490
str2_490 <- rep(NA, length(H))
str3_490 <- rep(NA, length(H))
for(h.iter in 1:length(H)){
  h = H[h.iter]
  start <- seq(1, 198, by = 2*h)
  start <- start[start + h <= 198]
  wh <- which(stations$z == z)
  str2_490[h.iter] <- mean((center_fire[wh[start+h],] - center_fire[wh[start],])^2)
  str3_490[h.iter] <- mean((abs(center_fire[wh[start+h],] - center_fire[wh[start],]))^3)
}

pdf(file="./Figures/struc_fun.pdf", width = 8.9, height=8.8)
par(mfrow = c(3,3))
plot(H, log(str2), type='b', ylab=expression(log(delta[2])), ylim=c(-55,-15), main="Original", mgp=c(2.4,1,0), xlab='h')
points(H, log(str2_250), type='b', col='blue')
points(H, log(str2_490), type='b', col='red')
legend('left', col=c("black",'blue', 'red'), lty=1,  legend=c("z=10", "z=250", "z=490"))

plot(H, log(str3), type='b', ylab = expression(log(delta[3])), ylim=c(-80,-20), main="Original", mgp=c(2.4,1,0), xlab='h')
points(H, log(str3_250), type='b', col='blue')
points(H, log(str3_490), type='b', col='red')
# legend('bottomleft', col=c("black",'blue', 'red'), lty=1,  legend=c("z=10", "z=250", "z=490"))

plot(H, str3/str2^{3/2}, type='b', ylab=expression(delta[3]/delta[2]^{3/2}),  ylim=c(1,9), main="Original", mgp=c(2.4,1,0), xlab='h')
points(H, str3_250/str2_250^{3/2}, type='b', col='blue')
points(H, str3_490/str2_490^{3/2}, type='b', col='red')
# legend('bottomleft', col=c("black",'blue', 'red'), lty=1,  legend=c("z=10", "z=250", "z=490"))





load("./emulation_PCA_ncomp_7.RData")


H <- 1:10
z = 10
str_POD2 <- rep(NA, length(H))
str_POD3 <- rep(NA, length(H))
for(h.iter in 1:length(H)){
  h = H[h.iter]
  start <- seq(1, 198, by = 2*h)
  start <- start[start + h <= 198]
  wh <- which(stations$z == z)
  str_POD2[h.iter] <- mean((fire_PCA[wh[start+h],] - fire_PCA[wh[start],])^2)
  str_POD3[h.iter] <- mean((abs(fire_PCA[wh[start+h],] - fire_PCA[wh[start],]))^3)
}


z = 250
str_POD2_250 <- rep(NA, length(H))
str_POD3_250 <- rep(NA, length(H))
for(h.iter in 1:length(H)){
  h = H[h.iter]
  start <- seq(1, 198, by = 2*h)
  start <- start[start + h <= 198]
  wh <- which(stations$z == z)
  str_POD2_250[h.iter] <- mean((fire_PCA[wh[start+h],] - fire_PCA[wh[start],])^2)
  str_POD3_250[h.iter] <- mean((abs(fire_PCA[wh[start+h],] - fire_PCA[wh[start],]))^3)
}



z = 490
str_POD2_490 <- rep(NA, length(H))
str_POD3_490 <- rep(NA, length(H))
for(h.iter in 1:length(H)){
  h = H[h.iter]
  start <- seq(1, 198, by = 2*h)
  start <- start[start + h <= 198]
  wh <- which(stations$z == z)
  str_POD2_490[h.iter] <- mean((fire_PCA[wh[start+h],] - fire_PCA[wh[start],])^2)
  str_POD3_490[h.iter] <- mean((abs(fire_PCA[wh[start+h],] - fire_PCA[wh[start],]))^3)
}

plot(H, log(str_POD2), type='b', ylab=expression(log(delta[2])), ylim=c(-55,-15), main="POD", mgp=c(2.4,1,0), xlab='h')
points(H, log(str_POD2_250), type='b', col='blue')
points(H, log(str_POD2_490), type='b', col='red')
legend('left', col=c("black",'blue', 'red'), lty=1,  legend=c("z=10", "z=250", "z=490"))

plot(H, log(str_POD3), type='b', ylab = expression(log(delta[3])), ylim=c(-80,-20),main="POD", mgp=c(2.4,1,0), xlab='h')
points(H, log(str_POD3_250), type='b', col='blue')
points(H, log(str_POD3_490), type='b', col='red')
# legend('bottomleft', col=c("black",'blue', 'red'), lty=1,  legend=c("z=10", "z=250", "z=490"))

plot(H, str_POD3/str_POD2^{3/2}, type='b', ylab=expression(delta[3]/delta[2]^{3/2}),  ylim=c(1,9),main="POD", mgp=c(2.4,1,0), xlab='h')
points(H, str_POD3_250/str_POD2_250^{3/2}, type='b', col='blue')
points(H, str_POD3_490/str_POD2_490^{3/2}, type='b', col='red')
# legend('bottomleft', col=c("black",'blue', 'red'), lty=1,  legend=c("z=10", "z=250", "z=490"))




load("./emulation_PCA.RData")


H <- 1:10
z = 10
str_XVAE2 <- rep(NA, length(H))
str_XVAE3 <- rep(NA, length(H))
for(h.iter in 1:length(H)){
  h = H[h.iter]
  start <- seq(1, 198, by = 2*h)
  start <- start[start + h <= 198]
  wh <- which(stations$z == z)
  str_XVAE2[h.iter] <- mean((fire_PCA[wh[start+h],] - fire_PCA[wh[start],])^2)
  str_XVAE3[h.iter] <- mean((abs(fire_PCA[wh[start+h],] - fire_PCA[wh[start],]))^3)
}


z = 250
str_XVAE2_250 <- rep(NA, length(H))
str_XVAE3_250 <- rep(NA, length(H))
for(h.iter in 1:length(H)){
  h = H[h.iter]
  start <- seq(1, 198, by = 2*h)
  start <- start[start + h <= 198]
  wh <- which(stations$z == z)
  str_XVAE2_250[h.iter] <- mean((fire_PCA[wh[start+h],] - fire_PCA[wh[start],])^2)
  str_XVAE3_250[h.iter] <- mean((abs(fire_PCA[wh[start+h],] - fire_PCA[wh[start],]))^3)
}



z = 490
str_XVAE2_490 <- rep(NA, length(H))
str_XVAE3_490 <- rep(NA, length(H))
for(h.iter in 1:length(H)){
  h = H[h.iter]
  start <- seq(1, 198, by = 2*h)
  start <- start[start + h <= 198]
  wh <- which(stations$z == z)
  str_XVAE2_490[h.iter] <- mean((fire_PCA[wh[start+h],] - fire_PCA[wh[start],])^2)
  str_XVAE3_490[h.iter] <- mean((abs(fire_PCA[wh[start+h],] - fire_PCA[wh[start],]))^3)
}

plot(H, log(str_XVAE2), type='b',ylim=c(-55,-15), ylab=expression(log(delta[2])), main="xVAE", mgp=c(2.4,1,0), xlab='h')
points(H, log(str_XVAE2_250), type='b', col='blue')
points(H, log(str_XVAE2_490), type='b', col='red')
legend('left', col=c("black",'blue', 'red'), lty=1,  legend=c("z=10", "z=250", "z=490"))

plot(H, log(str_XVAE3), type='b', ylab = expression(log(delta[3])), ylim=c(-80,-20),  main="xVAE", mgp=c(2.4,1,0), xlab='h')
points(H, log(str_XVAE3_250), type='b', col='blue')
points(H, log(str_XVAE3_490), type='b', col='red')
# legend('bottomleft', col=c("black",'blue', 'red'), lty=1,  legend=c("z=10", "z=250", "z=490"))

plot(H, str_XVAE3/str_XVAE2^{3/2}, type='b', ylab=expression(delta[3]/delta[2]^{3/2}),  ylim=c(1,9), main="xVAE", mgp=c(2.4,1,0), xlab='h')
points(H, str_XVAE3_250/str_XVAE2_250^{3/2}, type='b', col='blue')
points(H, str_XVAE3_490/str_XVAE2_490^{3/2}, type='b', col='red')
# legend('bottomleft', col=c("black",'blue', 'red'), lty=1,  legend=c("z=10", "z=250", "z=490"))
dev.off()

