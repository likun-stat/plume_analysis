
setwd("~/Desktop/Turbulence/")



###### ---------------------------------------------------------------------- ######
###### ----------------------------- Load in Data  -------------------------- ######
###### ---------------------------------------------------------------------- ######
stations <- expand.grid(x=1:198, z=1:500)
# 
# center_fire <- matrix(NA, nrow = nrow(stations), ncol=100)
# for(tt in 1:ncol(center_fire)){
#   tmp_file <- paste0("~/Documents/MATLAB/tmp",tt,".mat")
#   dat_mat <- readMat(tmp_file)
#   center_fire[,tt] <- as.numeric(dat_mat$aa)
# }
# 
# save(stations, center_fire, file="./center_fire_input.RData")

load("./center_fire_input.RData")
threshold <- quantile(center_fire,0.05)
# all(apply(center_fire,2, function(x) any(x<threshold))) # all times have threshold exceedances
center_fire_ori <- center_fire
center_fire <- log(-center_fire+1e-9)
pc_center <- colMeans(center_fire)
Center <- matrix(rep(pc_center, each=nrow(center_fire)), ncol=100)
center_fire_centered <- apply(center_fire, 2, function(x) x-mean(x))
pc <- prcomp(center_fire_centered, center = FALSE, scale. = FALSE)

eigenvalues <- pc$sdev^2
png(filename = "./Figures/percent_variation.png", width=400, height=320)
plot(eigenvalues/sum(eigenvalues), type = "b",
     xlab = "Principal Component",
     ylab = "Percentage of Variance Explained")
grid()
dev.off()



time <- 2
X_approx <- pc$x[,1:10] %*% pc$rotation[time, 1:10] + Center[,time]
# sum((center_fire[,time]- X_approx)^2)
# tmp_where <- center_fire[,time] < threshold
# sum((center_fire[tmp_where,time]- X_approx[tmp_where])^2)
Limits <- range(c(center_fire[,2], X_approx))
ggplot(stations) + geom_raster(aes(x=x, y=z, fill=X_approx)) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11,"Spectral"),
                       name = paste0('Approx time ', time), na.value = "transparent", limits=Limits)




## For the chi plot
fire_PCA <- pc$x[,1:50] %*% t(pc$rotation[, 1:50]) + Center
save(fire_PCA, file="./emulation_PCA.RData")



## For tail RMSE
threshold <- quantile(center_fire_ori,0.05)
metric_CV <- function(emulation, X, threshold){
  tmp_where <- (X < threshold)
  return(sum((emulation[tmp_where]-X[tmp_where])^2))
}

metric_CV(-exp(fire_PCA)+1e-9, center_fire_ori, threshold)
n.s <- 99000
n.t <- 100
sqrt(0.0009281847/(n.s*n.t))


###### ---------------------------------------------------------------------- ######
###### --------------------------- Cross validation ------------------------- ######
###### ---------------------------------------------------------------------- ######
metric_CV <- function(emulation, X, threshold){
  tmp_where <- (X < threshold)
  return(sum((emulation[tmp_where]-X[tmp_where])^2))
}

## CV to choose the optimum number of PCs
num_of_PCs <- round(seq(1,90, length.out=16))
num_of_PCs <- num_of_PCs[-1]
Metrix <- array(NA, dim=c(10,length(num_of_PCs))) 
for(fold in 1:10){
  cat("fold = ", fold, '\n')
  fold_ind <- seq(fold*10-9,fold*10,by=1)
  pc_tmp <- prcomp(center_fire_centered[,-fold_ind], center = FALSE, scale. = FALSE)
  for(iter in 1:length(num_of_PCs)){
    num <- num_of_PCs[iter]
    cat(" --- num_of_PC = ", num, '\n')
    bases <- pc_tmp$x[,1:num]
    coef <-  solve(t(bases)%*%bases)%*%t(bases)%*%center_fire_centered[,fold_ind]
    emulation_tmp <- bases %*% coef + Center[,fold_ind]
    emulation_t <- -exp(emulation_tmp) + 1e-9
    # fire_PCA <- emulation_tmp
    # save(fire_PCA, file=paste0("./emulation_PCA_FOLD",fold, "k_",num, ".RData"))
    Metrix[fold, iter] <- metric_CV(emulation_t, center_fire_ori[,fold_ind], threshold)
  }
}

# save(num_of_PCs, Metrix, file="./POD_results.RData")

load("./POD_results.RData")
n.s <- 99000
n.t <- 100
Metrix <- sqrt(Metrix/(n.s*n.t))*10^5
pdf("./Figures/POD_CV.pdf", width = 5, height=4.2)
plot(num_of_PCs, Metrix[1,], type='b', ylim=c(0.7, 3.45), mgp=c(2.8,1,0), xlab="Number of POD modes", ylab=expression(paste("Tail RMSE ("%*%10^{-5}, ")")), cex=0.76)
for(iter_tmp in 2:10){
  lines(num_of_PCs, Metrix[iter_tmp,], type='b', cex=0.76)
}
lines(num_of_PCs, colMeans(Metrix), type='b', col='red', pch=20)
grid()
dev.off()

par(mar = c(3.5,4,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
extRemes::qqplot(original_data[,ind], mon_max_emulation_smooth, 
                 ylab='Emulation', main='', regress = FALSE, xlab="")
title(xlab="Observations", line=1.9, cex.lab=1.2)
dev.off()


### ------------------- qq plots -------------------
### FOLD 1
fold=1
fold_ind <- seq(fold*10-9,fold*10,by=1)
pc_tmp <- prcomp(center_fire_centered[,-fold_ind], center = FALSE, scale. = FALSE)


num <- 7
bases <- pc_tmp$x[,1:num]
coef <-  solve(t(bases)%*%bases)%*%t(bases)%*%center_fire_centered[,fold_ind]
emulation_tmp <- bases %*% coef + Center[,fold_ind]
png(paste0("./Figures/qq_mon_POD_fold_", fold, "_num_of_PCS_", num, ".png"),width = 360, height = 380)
extRemes::qqplot(emulation_tmp[,1], center_fire[,fold_ind][,1], regress = FALSE,  
                 xlab=expression(paste('Observed ', X[1])), 
                 ylab=expression(paste('Emulated ', X[1])),
                 main = paste0("Fold # = ", fold, ", num_of_PCs = ", num))
dev.off()


num <- 48
bases <- pc_tmp$x[,1:num]
coef <-  solve(t(bases)%*%bases)%*%t(bases)%*%center_fire_centered[,fold_ind]
emulation_tmp <- bases %*% coef + Center[,fold_ind]
png(paste0("./Figures/qq_mon_POD_fold_", fold, "_num_of_PCS_", num, ".png"),width = 360, height = 380)
extRemes::qqplot(emulation_tmp[,1], center_fire[,fold_ind][,1], regress = FALSE,  
                 xlab=expression(paste('Observed ', X[1])), 
                 ylab=expression(paste('Emulated ', X[1])),
                 main = paste0("Fold # = ", fold, ", num_of_PCs = ", num))
dev.off()

num <- 90
bases <- pc_tmp$x[,1:num]
coef <-  solve(t(bases)%*%bases)%*%t(bases)%*%center_fire_centered[,fold_ind]
emulation_tmp <- bases %*% coef + Center[,fold_ind]
png(paste0("./Figures/qq_mon_POD_fold_", fold, "_num_of_PCS_", num, ".png"),width = 360, height = 380)
extRemes::qqplot(emulation_tmp[,1], center_fire[,fold_ind][,1], regress = FALSE,  
                 xlab=expression(paste('Observed ', X[1])), 
                 ylab=expression(paste('Emulated ', X[1])),
                 main = paste0("Fold # = ", fold, ", num_of_PCs = ", num))
dev.off()


### FOLD 10
fold=10
fold_ind <- seq(fold*10-9,fold*10,by=1)
pc_tmp <- prcomp(center_fire_centered[,-fold_ind], center = FALSE, scale. = FALSE)


num <- 7
bases <- pc_tmp$x[,1:num]
coef <-  solve(t(bases)%*%bases)%*%t(bases)%*%center_fire_centered[,fold_ind]
emulation_tmp <- bases %*% coef + Center[,fold_ind]
png(paste0("./Figures/qq_mon_POD_fold_", fold, "_num_of_PCS_", num, ".png"),width = 360, height = 380)
extRemes::qqplot(emulation_tmp[,1], center_fire[,fold_ind][,1], regress = FALSE,  
                 xlab=expression(paste('Observed ', X[1])), 
                 ylab=expression(paste('Emulated ', X[1])),
                 main = paste0("Fold # = ", fold, ", num_of_PCs = ", num))
dev.off()


num <- 48
bases <- pc_tmp$x[,1:num]
coef <-  solve(t(bases)%*%bases)%*%t(bases)%*%center_fire_centered[,fold_ind]
emulation_tmp <- bases %*% coef + Center[,fold_ind]
png(paste0("./Figures/qq_mon_POD_fold_", fold, "_num_of_PCS_", num, ".png"),width = 360, height = 380)
extRemes::qqplot(emulation_tmp[,1], center_fire[,fold_ind][,1], regress = FALSE,  
                 xlab=expression(paste('Observed ', X[1])), 
                 ylab=expression(paste('Emulated ', X[1])),
                 main = paste0("Fold # = ", fold, ", num_of_PCs = ", num))
dev.off()

num <- 90
bases <- pc_tmp$x[,1:num]
coef <-  solve(t(bases)%*%bases)%*%t(bases)%*%center_fire_centered[,fold_ind]
emulation_tmp <- bases %*% coef + Center[,fold_ind]
png(paste0("./Figures/qq_mon_POD_fold_", fold, "_num_of_PCS_", num, ".png"),width = 360, height = 380)
extRemes::qqplot(emulation_tmp[,1], center_fire[,fold_ind][,1], regress = FALSE,  
                 xlab=expression(paste('Observed ', X[1])), 
                 ylab=expression(paste('Emulated ', X[1])),
                 main = paste0("Fold # = ", fold, ", num_of_PCs = ", num))
dev.off()