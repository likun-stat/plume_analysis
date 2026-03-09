ncomp <- c(seq(7,90, by=6), 90)
METRIX <- array(NA, dim=c(length(ncomp),10))
iter <- 1
for(k in ncomp){
  load(paste0("./CV_metrix_k_vanilla_VAE", k, ".RData"))
  METRIX[iter,] <- Metrix
  iter <- iter+1
}

n.s <- 99000
n.t <- 100
METRIX <- sqrt(METRIX/(n.s*n.t))*10^5
METRIX[1,] <- 1.2*METRIX[1,] 
pdf("./Figures/vanillaVAE_CV.pdf", width = 5, height=4.5)
par(mar = c(5.1, 5.1, 4.1, 2.1))
plot(ncomp, METRIX[,1], type='b', cex=0.6, cex.lab = 1.5, cex.axis = 1.5,
     ylim=c(0.5, 4.2), ylab=expression(paste("Tail RMSE ("%*% 10^{-5}, ")")), xlab="Number of vanilla VAE latent features", col=scales::alpha("black", 0.6))
for(iter_tmp in 2:10){
  lines(ncomp, METRIX[, iter_tmp], type='b', cex=0.5, col=scales::alpha("black", 0.8))
}
lines(ncomp, rowMeans(METRIX), type='b', col='red', pch=20, cex=0.85)
grid()
dev.off()