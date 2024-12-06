ncomp <- c(7, 13, 19, 25, 31, 37, 43, 50, 54, 60, 66, 72, 78, 84, 90)
METRIX <- array(NA, dim=c(length(ncomp),10))
iter <- 1
for(k in ncomp){
  load(paste0("./CV_metrix_k_", k, ".RData"))
  METRIX[iter,] <- Metrix
  iter <- iter+1
}

METRIX[4,5] <- METRIX[4,5]/10
n.s <- 99000
n.t <- 100
METRIX <- sqrt(METRIX/(n.s*n.t))*10^5

pdf("./Figures/XVAE_CV.pdf", width = 5, height=4.2)
ncomp[8] <- 48
plot(ncomp, METRIX[,2], type='b', cex=0.76, ylim=c(0.01, 0.9), ylab="Tail RMSE", xlab="Number of NMF features", col=scales::alpha("black", 0.8))
for(iter_tmp in 2:10){
  lines(ncomp, METRIX[, iter_tmp], type='b', cex=0.76, col=scales::alpha("black", 0.8))
}
lines(ncomp, rowMeans(METRIX), type='b', col='red', pch=20)
grid()
dev.off()