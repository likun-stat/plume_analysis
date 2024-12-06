###############################################################################
###### -----------------------------  XVAE ---------------------------#########
###############################################################################
setwd('~/Desktop/Turbulence/')
load("./center_fire_input.RData")
load("./CV_metrix_k_7fold1.RData")

dim(center_fire)

pdf(file = "./Figures/qqplot_fold1_K_7.pdf", width = 4.5, height = 4.8)
lims <- range(center_fire[,1:10], fire_approx)
extRemes::qqplot(center_fire[,1:10], fire_approx, regress = FALSE,  
                 ylim=lims, xlim=lims,
                 ylab=expression(paste('Emulated ', bolditalic(X)[1:10])), 
                 xlab=expression(paste('Observed ', bolditalic(X)[1:10])), 
                 main = expression(italic(K)==7))
dev.off()


load("./CV_metrix_k_50fold1.RData")
pdf(file = "./Figures/qqplot_fold1_K_50.pdf", width = 4.5, height = 4.8)
lims <- range(center_fire[,1:10], fire_approx[,1:10])
extRemes::qqplot(center_fire[,1:10], fire_approx[,1:10], regress = FALSE,  
                 ylim=lims, xlim=lims,
                 ylab=expression(paste('Emulated ', bolditalic(X)[1:10])), 
                 xlab=expression(paste('Observed ', bolditalic(X)[1:10])), 
                 main = expression(italic(K)==50))
dev.off()


load("./CV_metrix_k_90fold1.RData")
pdf(file = "./Figures/qqplot_fold1_K_90.pdf", width = 4.5, height = 4.8)
lims <- range(center_fire[,1:10], fire_approx[,1:10])
extRemes::qqplot(center_fire[,1:10], fire_approx[,1:10], regress = FALSE,  
                 ylim=lims, xlim=lims,
                 ylab=expression(paste('Emulated ', bolditalic(X)[1:10])), 
                 xlab=expression(paste('Observed ', bolditalic(X)[1:10])), 
                 main = expression(italic(K)==90))
dev.off()








load("./CV_metrix_k_7fold10.RData")
pdf(file = "./Figures/qqplot_fold10_K_7.pdf", width = 4.5, height = 4.8)
lims <- range(center_fire[,91:100], fire_approx)
extRemes::qqplot(center_fire[,91:100], fire_approx, regress = FALSE,  
                 ylim=lims, xlim=lims,
                 ylab=expression(paste('Emulated ', bolditalic(X)[91:100])), 
                 xlab=expression(paste('Observed ', bolditalic(X)[91:100])), 
                 main = expression(italic(K)==7))
dev.off()


load("./CV_metrix_k_50fold10.RData")
pdf(file = "./Figures/qqplot_fold10_K_50.pdf", width = 4.5, height = 4.8)
lims <- range(center_fire[,91:100], fire_approx)
extRemes::qqplot(center_fire[,91:100], fire_approx, regress = FALSE,  
                 ylim=lims, xlim=lims,
                 ylab=expression(paste('Emulated ', bolditalic(X)[91:100])), 
                 xlab=expression(paste('Observed ', bolditalic(X)[91:100])), 
                 main = expression(italic(K)==50))
dev.off()


load("./CV_metrix_k_90fold10.RData")
pdf(file = "./Figures/qqplot_fold10_K_90.pdf", width = 4.5, height = 4.8)
lims <- range(center_fire[,91:100], fire_approx)
extRemes::qqplot(center_fire[,91:100], fire_approx, regress = FALSE,  
                 ylim=lims, xlim=lims,
                 ylab=expression(paste('Emulated ', bolditalic(X)[91:100])), 
                 xlab=expression(paste('Observed ', bolditalic(X)[91:100])), 
                 main = expression(italic(K)==90))
dev.off()




###############################################################################
###### -----------------------------  POD ----------------------------#########
###############################################################################
setwd('~/Desktop/Turbulence/')
load("./center_fire_input.RData")
load("./emulation_PCA_FOLD1k_7.RData")

dim(center_fire)
fire_PCA <- -exp(fire_PCA)+1e-9
pdf(file = "./Figures/qqplot_fold1_K_7_POD.pdf", width = 4.5, height = 4.8)
lims <- range(center_fire[,1:10], fire_PCA)
extRemes::qqplot(center_fire[,1:10], fire_PCA, regress = FALSE,  
                 ylim=lims, xlim=lims,
                 ylab=expression(paste('Emulated ', bolditalic(X)[1:10])), 
                 xlab=expression(paste('Observed ', bolditalic(X)[1:10])), 
                 main = expression(italic(K)==7))
dev.off()


load("./emulation_PCA_FOLD1k_50.RData")
fire_PCA <- -exp(fire_PCA)+1e-9
pdf(file = "./Figures/qqplot_fold1_K_50_POD.pdf", width = 4.5, height = 4.8)
lims <- range(center_fire[,1:10], fire_PCA)
extRemes::qqplot(center_fire[,1:10], fire_PCA, regress = FALSE,  
                 ylim=lims, xlim=lims,
                 ylab=expression(paste('Emulated ', bolditalic(X)[1:10])), 
                 xlab=expression(paste('Observed ', bolditalic(X)[1:10])), 
                 main = expression(italic(K)==50))
dev.off()


load("./emulation_PCA_FOLD1k_90.RData")
fire_PCA <- -exp(fire_PCA)+1e-9
pdf(file = "./Figures/qqplot_fold1_K_90_POD.pdf", width = 4.5, height = 4.8)
lims <- range(center_fire[,1:10], fire_PCA)
extRemes::qqplot(center_fire[,1:10], fire_PCA, regress = FALSE,  
                 ylim=lims, xlim=lims,
                 ylab=expression(paste('Emulated ', bolditalic(X)[1:10])), 
                 xlab=expression(paste('Observed ', bolditalic(X)[1:10])), 
                 main = expression(italic(K)==90))
dev.off()




load("./center_fire_input.RData")
load("./emulation_PCA_FOLD10k_7.RData")

dim(center_fire)
fire_PCA <- -exp(fire_PCA)+1e-9
pdf(file = "./Figures/qqplot_fold10_K_7_POD.pdf", width = 4.5, height = 4.8)
lims <- c(-5e-4,0)
extRemes::qqplot(center_fire[,91:100], fire_PCA, regress = FALSE,  
                 ylim=lims, xlim=lims,
                 ylab=expression(paste('Emulated ', bolditalic(X)[91:100])), 
                 xlab=expression(paste('Observed ', bolditalic(X)[91:100])), 
                 main = expression(italic(K)==7))
dev.off()


load("./emulation_PCA_FOLD10k_50.RData")
fire_PCA <- -exp(fire_PCA)+1e-9
pdf(file = "./Figures/qqplot_fold10_K_50_POD.pdf", width = 4.5, height = 4.8)
lims <- c(-5e-4,0)
extRemes::qqplot(center_fire[,91:100], fire_PCA, regress = FALSE,  
                 ylim=lims, xlim=lims,
                 ylab=expression(paste('Emulated ', bolditalic(X)[91:100])), 
                 xlab=expression(paste('Observed ', bolditalic(X)[91:100])), 
                 main = expression(italic(K)==50))
dev.off()


load("./emulation_PCA_FOLD10k_90.RData")
fire_PCA <- -exp(fire_PCA)+1e-9
pdf(file = "./Figures/qqplot_fold10_K_90_POD.pdf", width = 4.5, height = 4.8)
lims <- c(-5e-4,0)
extRemes::qqplot(center_fire[,91:100], fire_PCA, regress = FALSE,  
                 ylim=lims, xlim=lims,
                 ylab=expression(paste('Emulated ', bolditalic(X)[91:100])), 
                 xlab=expression(paste('Observed ', bolditalic(X)[91:100])), 
                 main = expression(italic(K)==90))
dev.off()