ncomp <- c(50, 100, 200, 500, 900)
METRIX <- array(NA, dim=c(length(ncomp),10))
iter <- 1
for(k in ncomp){
  load(paste0("./CV_metrix_k_", k, ".RData"))
  METRIX[iter,] <- Metrix
  iter <- iter+1
}

n.s <- 99000
n.t <- 100
METRIX <- sqrt(METRIX/(n.s*n.t))*10^3

pdf("./Figures/XVAE_CV.pdf", width = 5, height=4.2)
plot(ncomp, METRIX[,2], type='b', cex=0.76, ylim=c(0.01, 0.9), ylab="Tail RMSE", xlab="Number of NMF features", col=scales::alpha("black", 0.8))
for(iter_tmp in 2:10){
  lines(ncomp, METRIX[, iter_tmp], type='b', cex=0.76, col=scales::alpha("black", 0.8))
}
lines(ncomp, rowMeans(METRIX), type='b', col='red', pch=20)
grid()
dev.off()