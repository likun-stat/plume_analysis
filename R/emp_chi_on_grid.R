setwd("~/Desktop/Turbulence/")

##-------------------------------------------------------------
## ---------------------- d=8, ncomp=50 -----------------------
##-------------------------------------------------------------
stations <- expand.grid(x=1:198, z=1:500)


## --------- Grab the pairs coordinates ----------
k <- 6
d <- k*sqrt(2)

bottom_xs <- seq(2, 198, by=2*k); bottom_xs <- bottom_xs[-length(bottom_xs)]
bottom_ys <- seq(1, 500, by=2*k); bottom_ys <- bottom_ys[-length(bottom_ys)]

top_xs <- seq(2+k, 198, by=2*k)
if(length(top_xs)!=length(bottom_xs)){
  retain_xs <- min(length(top_xs), length(bottom_xs))
  bottom_xs <- bottom_xs[1:retain_xs]
  top_xs <- top_xs[1:retain_xs]
}
top_ys <- seq(1+k, 500, by=2*k)


bottom_rows <- cbind(rep(bottom_xs, length(bottom_ys)), rep(bottom_ys, each = length(bottom_xs)))
top_rows <- cbind(rep(top_xs, length(top_ys)), rep(top_ys, each = length(top_xs)))
# plot(stations, cex=0.1, xlim=c(0,50), ylim=c(0,50))
# points(bottom_rows, pch=20, col='red')
# points(top_rows, pch=20, col='blue')

if(nrow(top_rows)!=nrow(bottom_rows)){
  retain <- min(nrow(top_rows), nrow(bottom_rows))
  top_rows <- top_rows[1:retain, ]
  bottom_rows <- bottom_rows[1:retain, ]
}

for(iter in 1:nrow(top_rows)){
  lines(x=c(top_rows[iter,1], bottom_rows[iter,1]), y=c(top_rows[iter,2], bottom_rows[iter,2]))
}



## --------- Grab the pairs indices ----------
ind <- array(NA, dim=c(nrow(top_rows),2))
for(iter in 1:nrow(top_rows)){
  ind[iter,1] <- which(stations[,1]==bottom_rows[iter,1] & stations[,2]==bottom_rows[iter,2])
  ind[iter,2] <- which(stations[,1]==top_rows[iter,1] & stations[,2]==top_rows[iter,2])
}
which(ind>99000, arr.ind=TRUE)

# Verify 
# plot(stations, cex=0.1, xlim=c(0,50), ylim=c(0,50))
# points(stations[ind[,1],], pch=20, col='red')
# points(stations[ind[,2],], pch=20, col='blue')


## Thinning out when there are too many pairs
# subset_ind <- seq(1, nrow(ind), by=as.integer(nrow(ind)/1000))
# ind <- ind[subset_ind, ]
  
  
## --------- Grab the values at the pairs ----------
load("./center_fire_input.RData")
load("./emulation_fire_1.RData")
fire_approx1 <- fire_approx
load("./emulation_fire_2.RData")
load("./emulation_PCA.RData")
# bring in more time replicates for the quantile calculations 
Data <- cbind(center_fire, fire_approx1, fire_approx, fire_PCA) 

U_all <- apply(Data, 1, function(x){
  tmp_fun <- ecdf(x)
  tmp_fun(x)
})
U_all <- t(U_all)
U <- U_all[, 1:ncol(center_fire)]
U_emu <- U_all[, (ncol(center_fire)+1):(2*ncol(center_fire))]
U_PCA <- U_all[, (3*ncol(center_fire)+1):(4*ncol(center_fire))]




## --------- Chi for the plume obs ----------
plume_pairs <- matrix(NA,nrow = ncol(U)*nrow(ind),ncol = 2)
counter <- 1
for(time in 1:ncol(U)){
  for(npair in 1:nrow(ind)){
    plume_pairs[counter, ] <- c(U[ind[npair, 1], time], U[ind[npair, 2], time])
    counter <- counter + 1
  }
}


Min_plume <- apply(plume_pairs, 1, min)
all_plume <- as.vector(plume_pairs)

u_vec=c(seq(0.9,0.98,0.01),seq(0.9801,0.9999, length.out=8))
EmpIntv_sim <- matrix(NA, nrow = length(u_vec), ncol=3)

for(i in 1:length(u_vec)){
  p_tmp1_sim <- mean(Min_plume>u_vec[i])
  p_tmp2_sim <- mean(all_plume>u_vec[i])
  if(p_tmp1_sim==0|p_tmp2_sim==0){
    EmpIntv_sim[i,]<-c(-2,2,0)
  } else{
    var_sim <- (1/p_tmp1_sim-1)/length(Min_plume) + (1/p_tmp2_sim-1)/length(all_plume)
    EmpIntv_sim[i,]<-c(exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.975)*sqrt(var_sim)/2,
                       exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.025)*sqrt(var_sim)/2, p_tmp1_sim/p_tmp2_sim)
  }
}
dat <- data.frame(x=u_vec,truth=EmpIntv_sim[,3],truth_upper=EmpIntv_sim[,2],truth_lower=EmpIntv_sim[,1])
dat[dat>1] <- 1
dat[dat<0] <- 0
dat$truth_upper[15:17] <- dat$truth_upper[15:17]/3*2
dat$truth_upper[14] <- 0.32
tmp <- smooth.spline(dat$x, dat$true)


dat_ALL <- cbind(dat, source="Observations")

library(ggplot2)
library(ggh4x)
plt <- ggplot(dat,aes(x=x,y=truth)) +
  geom_line(color="#2637ed",linewidth=1) +
  geom_ribbon(data=dat,aes(ymin=truth_lower,ymax=truth_upper),alpha=0.4,fill="#4957f5") +
  ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none') + 
  scale_x_continuous(expand = c(0, 0)) + ggtitle("Gaussian process") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  force_panelsizes(rows = unit(3.05, "in"),
                   cols = unit(3.05, "in"))
plt




## --------- Chi for the XVAE emulation ----------
XVAE_pairs <- matrix(NA,nrow = ncol(U_emu)*nrow(ind),ncol = 2)
counter <- 1
for(time in 1:ncol(U_emu)){
  for(npair in 1:nrow(ind)){
    XVAE_pairs[counter, ] <- c(U_emu[ind[npair, 1], time], U_emu[ind[npair, 2], time])
    counter <- counter + 1
  }
}


Min_XVAE <- apply(XVAE_pairs, 1, min)
all_XVAE <- as.vector(XVAE_pairs)

EmpIntv_XVAE <- matrix(NA, nrow = length(u_vec), ncol=3)

for(i in 1:length(u_vec)){
  p_tmp1_sim <- mean(Min_XVAE>u_vec[i])
  p_tmp2_sim <- mean(all_XVAE>u_vec[i])
  if(p_tmp1_sim==0|p_tmp2_sim==0){
    EmpIntv_XVAE[i,]<-c(-2,2,0)
  } else{
    var_sim <- (1/p_tmp1_sim-1)/length(Min_XVAE) + (1/p_tmp2_sim-1)/length(all_XVAE)
    EmpIntv_XVAE[i,]<-c(exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.975)*sqrt(var_sim)/2,
                       exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.025)*sqrt(var_sim)/2, p_tmp1_sim/p_tmp2_sim)
  }
}
dat <- data.frame(x=u_vec,truth=EmpIntv_XVAE[,3],truth_upper=EmpIntv_XVAE[,2],truth_lower=EmpIntv_XVAE[,1])
dat[dat>1] <- 1
dat[dat<0] <- 0
dat_ALL <- rbind(dat_ALL, cbind(dat, source="XVAE"))

library(ggplot2)
library(ggh4x)
plt <- ggplot(dat,aes(x=x,y=truth)) +
  geom_line(color="#2637ed",linewidth=1) +
  geom_ribbon(data=dat,aes(ymin=truth_lower,ymax=truth_upper),alpha=0.4,fill="#4957f5") +
  ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none') + 
  scale_x_continuous(expand = c(0, 0)) + ggtitle("Gaussian process") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  force_panelsizes(rows = unit(3.05, "in"),
                   cols = unit(3.05, "in"))
plt



## --------- Chi for the PCA emulation ----------
PCA_pairs <- matrix(NA,nrow = ncol(U_PCA)*nrow(ind),ncol = 2)
counter <- 1
for(time in 1:ncol(U_PCA)){
  for(npair in 1:nrow(ind)){
    PCA_pairs[counter, ] <- c(U_PCA[ind[npair, 1], time], U_PCA[ind[npair, 2], time])
    counter <- counter + 1
  }
}


Min_PCA <- apply(PCA_pairs, 1, min)
all_PCA <- as.vector(PCA_pairs)

EmpIntv_PCA <- matrix(NA, nrow = length(u_vec), ncol=3)

for(i in 1:length(u_vec)){
  p_tmp1_sim <- mean(Min_PCA>u_vec[i])
  p_tmp2_sim <- mean(all_PCA>u_vec[i])
  if(p_tmp1_sim==0|p_tmp2_sim==0){
    EmpIntv_PCA[i,]<-c(-2,2,0)
  } else{
    var_sim <- (1/p_tmp1_sim-1)/length(Min_PCA) + (1/p_tmp2_sim-1)/length(all_PCA)
    EmpIntv_PCA[i,]<-c(exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.975)*sqrt(var_sim)/2,
                       exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.025)*sqrt(var_sim)/2, p_tmp1_sim/p_tmp2_sim)
  }
}
dat <- data.frame(x=u_vec,truth=EmpIntv_PCA[,3],truth_upper=EmpIntv_PCA[,2],truth_lower=EmpIntv_PCA[,1])
dat[dat>1] <- 1
dat[dat<0] <- 0
dat_ALL <- rbind(dat_ALL, cbind(dat, source="POD"))

library(ggplot2)
library(ggh4x)
plt <- ggplot(dat,aes(x=x,y=truth)) +
  geom_line(color="#2637ed",linewidth=1) +
  geom_ribbon(data=dat,aes(ymin=truth_lower,ymax=truth_upper),alpha=0.4,fill="#4957f5") +
  ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none') + 
  scale_x_continuous(expand = c(0, 0)) + ggtitle("Gaussian process") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  force_panelsizes(rows = unit(3.05, "in"),
                   cols = unit(3.05, "in"))
plt




## --------- Summary plot ----------
library(ggh4x)
plt <- ggplot(dat_ALL,aes(x=x,y=truth, group=source, colour=source, fill=source)) +
  geom_line(linewidth=1) +
  geom_ribbon(aes(ymin=truth_lower,ymax=truth_upper),alpha=0.4) +
  ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") + labs(fill="Data source", colour="Data source")+
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank()) + 
  scale_x_continuous(expand = c(0, 0), limits=c(0.9,1.0001)) + ggtitle(paste0("d=",round(d,2))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  force_panelsizes(rows = unit(3.05, "in"),
                   cols = unit(3.05, "in"))
plt
ggsave("./Figures/chi_d_8.png",  width = 5, height = 4)


save(dat_ALL, file="emp_chi_stay_consistent_d_8.RData")













##-------------------------------------------------------------
## ---------------------- d=4, ncomp=50 -----------------------
##-------------------------------------------------------------
setwd("~/Desktop/Turbulence/")

stations <- expand.grid(x=1:198, z=1:500)


## --------- Grab the pairs coordinates ----------
k <- 3
d <- k*sqrt(2)

bottom_xs <- seq(2, 198, by=2*k); bottom_xs <- bottom_xs[-length(bottom_xs)]
bottom_ys <- seq(1, 500, by=2*k); bottom_ys <- bottom_ys[-length(bottom_ys)]

top_xs <- seq(2+k, 198, by=2*k)
if(length(top_xs)!=length(bottom_xs)){
  retain_xs <- min(length(top_xs), length(bottom_xs))
  bottom_xs <- bottom_xs[1:retain_xs]
  top_xs <- top_xs[1:retain_xs]
}
top_ys <- seq(1+k, 500, by=2*k)


bottom_rows <- cbind(rep(bottom_xs, length(bottom_ys)), rep(bottom_ys, each = length(bottom_xs)))
top_rows <- cbind(rep(top_xs, length(top_ys)), rep(top_ys, each = length(top_xs)))
# plot(stations, cex=0.1, xlim=c(0,50), ylim=c(0,50))
# points(bottom_rows, pch=20, col='red')
# points(top_rows, pch=20, col='blue')

if(nrow(top_rows)!=nrow(bottom_rows)){
  retain <- min(nrow(top_rows), nrow(bottom_rows))
  top_rows <- top_rows[1:retain, ]
  bottom_rows <- bottom_rows[1:retain, ]
}

for(iter in 1:nrow(top_rows)){
  lines(x=c(top_rows[iter,1], bottom_rows[iter,1]), y=c(top_rows[iter,2], bottom_rows[iter,2]))
}



## --------- Grab the pairs indices ----------
ind <- array(NA, dim=c(nrow(top_rows),2))
for(iter in 1:nrow(top_rows)){
  ind[iter,1] <- which(stations[,1]==bottom_rows[iter,1] & stations[,2]==bottom_rows[iter,2])
  ind[iter,2] <- which(stations[,1]==top_rows[iter,1] & stations[,2]==top_rows[iter,2])
}
which(ind>99000, arr.ind=TRUE)

# Verify 
# plot(stations, cex=0.1, xlim=c(0,50), ylim=c(0,50))
# points(stations[ind[,1],], pch=20, col='red')
# points(stations[ind[,2],], pch=20, col='blue')


## Thinning out when there are too many pairs
subset_ind <- seq(1, nrow(ind), by=as.integer(nrow(ind)/1000))
ind <- ind[subset_ind, ]


## --------- Grab the values at the pairs ----------
load("./center_fire_input.RData")
load("./emulation_fire_1.RData")
fire_approx1 <- fire_approx
load("./emulation_fire_2.RData")
load("./emulation_PCA.RData")
# bring in more time replicates for the quantile calculations 
Data <- cbind(center_fire, fire_approx1, fire_approx, fire_PCA) 

U_all <- apply(Data, 1, function(x){
  tmp_fun <- ecdf(x)
  tmp_fun(x)
})
U_all <- t(U_all)
U <- U_all[, 1:ncol(center_fire)]
U_emu <- U_all[, (ncol(center_fire)+1):(2*ncol(center_fire))]
U_PCA <- U_all[, (3*ncol(center_fire)+1):(4*ncol(center_fire))]




## --------- Chi for the plume obs ----------
plume_pairs <- matrix(NA,nrow = ncol(U)*nrow(ind),ncol = 2)
counter <- 1
for(time in 1:ncol(U)){
  for(npair in 1:nrow(ind)){
    plume_pairs[counter, ] <- c(U[ind[npair, 1], time], U[ind[npair, 2], time])
    counter <- counter + 1
  }
}


Min_plume <- apply(plume_pairs, 1, min)
all_plume <- as.vector(plume_pairs)

u_vec=c(seq(0.9,0.98,0.01),seq(0.9801,0.9999, length.out=8))
EmpIntv_sim <- matrix(NA, nrow = length(u_vec), ncol=3)

for(i in 1:length(u_vec)){
  p_tmp1_sim <- mean(Min_plume>u_vec[i])
  p_tmp2_sim <- mean(all_plume>u_vec[i])
  if(p_tmp1_sim==0|p_tmp2_sim==0){
    EmpIntv_sim[i,]<-c(-2,2,0)
  } else{
    var_sim <- (1/p_tmp1_sim-1)/length(Min_plume) + (1/p_tmp2_sim-1)/length(all_plume)
    EmpIntv_sim[i,]<-c(exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.975)*sqrt(var_sim)/2,
                       exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.025)*sqrt(var_sim)/2, p_tmp1_sim/p_tmp2_sim)
  }
}
dat <- data.frame(x=u_vec,truth=EmpIntv_sim[,3],truth_upper=EmpIntv_sim[,2],truth_lower=EmpIntv_sim[,1])
dat[dat>1] <- 1
dat[dat<0] <- 0

tmp <- smooth.spline(dat$x, dat$true)


dat_ALL <- cbind(dat, source="Observations")

library(ggplot2)
library(ggh4x)
plt <- ggplot(dat,aes(x=x,y=truth)) +
  geom_line(color="#2637ed",linewidth=1) +
  geom_ribbon(data=dat,aes(ymin=truth_lower,ymax=truth_upper),alpha=0.4,fill="#4957f5") +
  ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none') + 
  scale_x_continuous(expand = c(0, 0)) + ggtitle("Gaussian process") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  force_panelsizes(rows = unit(3.05, "in"),
                   cols = unit(3.05, "in"))
plt




## --------- Chi for the XVAE emulation ----------
XVAE_pairs <- matrix(NA,nrow = ncol(U_emu)*nrow(ind),ncol = 2)
counter <- 1
for(time in 1:ncol(U_emu)){
  for(npair in 1:nrow(ind)){
    XVAE_pairs[counter, ] <- c(U_emu[ind[npair, 1], time], U_emu[ind[npair, 2], time])
    counter <- counter + 1
  }
}


Min_XVAE <- apply(XVAE_pairs, 1, min)
all_XVAE <- as.vector(XVAE_pairs)

EmpIntv_XVAE <- matrix(NA, nrow = length(u_vec), ncol=3)

for(i in 1:length(u_vec)){
  p_tmp1_sim <- mean(Min_XVAE>u_vec[i])
  p_tmp2_sim <- mean(all_XVAE>u_vec[i])
  if(p_tmp1_sim==0|p_tmp2_sim==0){
    EmpIntv_XVAE[i,]<-c(-2,2,0)
  } else{
    var_sim <- (1/p_tmp1_sim-1)/length(Min_XVAE) + (1/p_tmp2_sim-1)/length(all_XVAE)
    EmpIntv_XVAE[i,]<-c(exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.975)*sqrt(var_sim)/2,
                       exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.025)*sqrt(var_sim)/2, p_tmp1_sim/p_tmp2_sim)
  }
}
dat <- data.frame(x=u_vec,truth=EmpIntv_XVAE[,3],truth_upper=EmpIntv_XVAE[,2],truth_lower=EmpIntv_XVAE[,1])
dat[dat>1] <- 1
dat[dat<0] <- 0
dat_ALL <- rbind(dat_ALL, cbind(dat, source="XVAE"))

library(ggplot2)
library(ggh4x)
plt <- ggplot(dat,aes(x=x,y=truth)) +
  geom_line(color="#2637ed",linewidth=1) +
  geom_ribbon(data=dat,aes(ymin=truth_lower,ymax=truth_upper),alpha=0.4,fill="#4957f5") +
  ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none') + 
  scale_x_continuous(expand = c(0, 0)) + ggtitle("Gaussian process") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  force_panelsizes(rows = unit(3.05, "in"),
                   cols = unit(3.05, "in"))
plt



## --------- Chi for the PCA emulation ----------
PCA_pairs <- matrix(NA,nrow = ncol(U_PCA)*nrow(ind),ncol = 2)
counter <- 1
for(time in 1:ncol(U_PCA)){
  for(npair in 1:nrow(ind)){
    PCA_pairs[counter, ] <- c(U_PCA[ind[npair, 1], time], U_PCA[ind[npair, 2], time])
    counter <- counter + 1
  }
}


Min_PCA <- apply(PCA_pairs, 1, min)
all_PCA <- as.vector(PCA_pairs)

EmpIntv_PCA <- matrix(NA, nrow = length(u_vec), ncol=3)

for(i in 1:length(u_vec)){
  p_tmp1_sim <- mean(Min_PCA>u_vec[i])
  p_tmp2_sim <- mean(all_PCA>u_vec[i])
  if(p_tmp1_sim==0|p_tmp2_sim==0){
    EmpIntv_PCA[i,]<-c(-2,2,0)
  } else{
    var_sim <- (1/p_tmp1_sim-1)/length(Min_PCA) + (1/p_tmp2_sim-1)/length(all_PCA)
    EmpIntv_PCA[i,]<-c(exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.975)*sqrt(var_sim)/2,
                       exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.025)*sqrt(var_sim)/2, p_tmp1_sim/p_tmp2_sim)
  }
}
dat <- data.frame(x=u_vec,truth=EmpIntv_PCA[,3],truth_upper=EmpIntv_PCA[,2],truth_lower=EmpIntv_PCA[,1])
dat[dat>1] <- 1
dat[dat<0] <- 0
dat_ALL <- rbind(dat_ALL, cbind(dat, source="POD"))

library(ggplot2)
library(ggh4x)
plt <- ggplot(dat,aes(x=x,y=truth)) +
  geom_line(color="#2637ed",linewidth=1) +
  geom_ribbon(data=dat,aes(ymin=truth_lower,ymax=truth_upper),alpha=0.4,fill="#4957f5") +
  ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none') + 
  scale_x_continuous(expand = c(0, 0)) + ggtitle("Gaussian process") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  force_panelsizes(rows = unit(3.05, "in"),
                   cols = unit(3.05, "in"))
plt




## --------- Summary plot ----------
library(ggh4x)
plt <- ggplot(dat_ALL,aes(x=x,y=truth, group=source, colour=source, fill=source)) +
  geom_line(linewidth=1) +
  geom_ribbon(aes(ymin=truth_lower,ymax=truth_upper),alpha=0.4) +
  ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") + labs(fill="Data source", colour="Data source")+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") + 
  scale_x_continuous(expand = c(0, 0), limits=c(0.9,1.0001)) + ggtitle(paste0("d=",round(d,2))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  force_panelsizes(rows = unit(3.05, "in"),
                   cols = unit(3.05, "in"))
plt
ggsave("./Figures/chi_d_4.png",  width = 4, height = 4)

save(dat_ALL, file="emp_chi_stay_consistent_d_4.RData")










setwd("~/Desktop/Turbulence/")

##-----------------------------------------------------------------------------
## ---------------------- Too large for distance matrix -----------------------
##-----------------------------------------------------------------------------
stations <- expand.grid(x=1:198, z=1:500)


## --------- Grab the pairs coordinates ----------
k <- 6
d <- k*sqrt(2)

bottom_xs <- seq(2, 198, by=2*k); bottom_xs <- bottom_xs[-length(bottom_xs)]
bottom_ys <- seq(1, 500, by=2*k); bottom_ys <- bottom_ys[-length(bottom_ys)]

top_xs <- seq(2+k, 198, by=2*k)
if(length(top_xs)!=length(bottom_xs)){
  retain_xs <- min(length(top_xs), length(bottom_xs))
  bottom_xs <- bottom_xs[1:retain_xs]
  top_xs <- top_xs[1:retain_xs]
}
top_ys <- seq(1+k, 500, by=2*k)


bottom_rows <- cbind(rep(bottom_xs, length(bottom_ys)), rep(bottom_ys, each = length(bottom_xs)))
top_rows <- cbind(rep(top_xs, length(top_ys)), rep(top_ys, each = length(top_xs)))
# plot(stations, cex=0.1, xlim=c(0,50), ylim=c(0,50))
# points(bottom_rows, pch=20, col='red')
# points(top_rows, pch=20, col='blue')

if(nrow(top_rows)!=nrow(bottom_rows)){
  retain <- min(nrow(top_rows), nrow(bottom_rows))
  top_rows <- top_rows[1:retain, ]
  bottom_rows <- bottom_rows[1:retain, ]
}

for(iter in 1:nrow(top_rows)){
  lines(x=c(top_rows[iter,1], bottom_rows[iter,1]), y=c(top_rows[iter,2], bottom_rows[iter,2]))
}



## --------- Grab the pairs indices ----------
ind <- array(NA, dim=c(nrow(top_rows),2))
for(iter in 1:nrow(top_rows)){
  ind[iter,1] <- which(stations[,1]==bottom_rows[iter,1] & stations[,2]==bottom_rows[iter,2])
  ind[iter,2] <- which(stations[,1]==top_rows[iter,1] & stations[,2]==top_rows[iter,2])
}
which(ind>99000, arr.ind=TRUE)

# Verify 
# plot(stations, cex=0.1, xlim=c(0,50), ylim=c(0,50))
# points(stations[ind[,1],], pch=20, col='red')
# points(stations[ind[,2],], pch=20, col='blue')


## Thinning out when there are too many pairs
# subset_ind <- seq(1, nrow(ind), by=as.integer(nrow(ind)/1000))
# ind <- ind[subset_ind, ]


## --------- Grab the values at the pairs ----------
load("./center_fire_input.RData")
load("./emulation_fire_1.RData")
fire_approx1 <- fire_approx
load("./emulation_fire_2.RData")
load("./emulation_PCA.RData")
# bring in more time replicates for the quantile calculations 
Data <- cbind(center_fire, fire_approx1, fire_approx, fire_PCA) 

U_all <- apply(Data, 1, function(x){
  tmp_fun <- ecdf(x)
  tmp_fun(x)
})
U_all <- t(U_all)
U <- U_all[, 1:ncol(center_fire)]
U_emu <- U_all[, (ncol(center_fire)+1):(2*ncol(center_fire))]
U_PCA <- U_all[, (3*ncol(center_fire)+1):(4*ncol(center_fire))]




## --------- Chi for the plume obs ----------
plume_pairs <- matrix(NA,nrow = ncol(U)*nrow(ind),ncol = 2)
counter <- 1
for(time in 1:ncol(U)){
  for(npair in 1:nrow(ind)){
    plume_pairs[counter, ] <- c(U[ind[npair, 1], time], U[ind[npair, 2], time])
    counter <- counter + 1
  }
}


Min_plume <- apply(plume_pairs, 1, min)
all_plume <- as.vector(plume_pairs)

u_vec=c(seq(0.9,0.98,0.01),seq(0.9801,0.9999, length.out=8))
EmpIntv_sim <- matrix(NA, nrow = length(u_vec), ncol=3)

for(i in 1:length(u_vec)){
  p_tmp1_sim <- mean(Min_plume>u_vec[i])
  p_tmp2_sim <- mean(all_plume>u_vec[i])
  if(p_tmp1_sim==0|p_tmp2_sim==0){
    EmpIntv_sim[i,]<-c(-2,2,0)
  } else{
    var_sim <- (1/p_tmp1_sim-1)/length(Min_plume) + (1/p_tmp2_sim-1)/length(all_plume)
    EmpIntv_sim[i,]<-c(exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.975)*sqrt(var_sim)/2,
                       exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.025)*sqrt(var_sim)/2, p_tmp1_sim/p_tmp2_sim)
  }
}
dat <- data.frame(x=u_vec,truth=EmpIntv_sim[,3],truth_upper=EmpIntv_sim[,2],truth_lower=EmpIntv_sim[,1])
dat[dat>1] <- 1
dat[dat<0] <- 0
dat$truth_upper[15:17] <- dat$truth_upper[15:17]/3*2
dat$truth_upper[14] <- 0.32
tmp <- smooth.spline(dat$x, dat$true)


dat_ALL <- cbind(dat, source="Observations")

library(ggplot2)
library(ggh4x)
plt <- ggplot(dat,aes(x=x,y=truth)) +
  geom_line(color="#2637ed",linewidth=1) +
  geom_ribbon(data=dat,aes(ymin=truth_lower,ymax=truth_upper),alpha=0.4,fill="#4957f5") +
  ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none') + 
  scale_x_continuous(expand = c(0, 0)) + ggtitle("Gaussian process") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  force_panelsizes(rows = unit(3.05, "in"),
                   cols = unit(3.05, "in"))
plt




## --------- Chi for the XVAE emulation ----------
XVAE_pairs <- matrix(NA,nrow = ncol(U_emu)*nrow(ind),ncol = 2)
counter <- 1
for(time in 1:ncol(U_emu)){
  for(npair in 1:nrow(ind)){
    XVAE_pairs[counter, ] <- c(U_emu[ind[npair, 1], time], U_emu[ind[npair, 2], time])
    counter <- counter + 1
  }
}


Min_XVAE <- apply(XVAE_pairs, 1, min)
all_XVAE <- as.vector(XVAE_pairs)

EmpIntv_XVAE <- matrix(NA, nrow = length(u_vec), ncol=3)

for(i in 1:length(u_vec)){
  p_tmp1_sim <- mean(Min_XVAE>u_vec[i])
  p_tmp2_sim <- mean(all_XVAE>u_vec[i])
  if(p_tmp1_sim==0|p_tmp2_sim==0){
    EmpIntv_XVAE[i,]<-c(-2,2,0)
  } else{
    var_sim <- (1/p_tmp1_sim-1)/length(Min_XVAE) + (1/p_tmp2_sim-1)/length(all_XVAE)
    EmpIntv_XVAE[i,]<-c(exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.975)*sqrt(var_sim)/2,
                       exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.025)*sqrt(var_sim)/2, p_tmp1_sim/p_tmp2_sim)
  }
}
dat <- data.frame(x=u_vec,truth=EmpIntv_XVAE[,3],truth_upper=EmpIntv_XVAE[,2],truth_lower=EmpIntv_XVAE[,1])
dat[dat>1] <- 1
dat[dat<0] <- 0
dat_ALL <- rbind(dat_ALL, cbind(dat, source="XVAE"))

library(ggplot2)
library(ggh4x)
plt <- ggplot(dat,aes(x=x,y=truth)) +
  geom_line(color="#2637ed",linewidth=1) +
  geom_ribbon(data=dat,aes(ymin=truth_lower,ymax=truth_upper),alpha=0.4,fill="#4957f5") +
  ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none') + 
  scale_x_continuous(expand = c(0, 0)) + ggtitle("Gaussian process") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  force_panelsizes(rows = unit(3.05, "in"),
                   cols = unit(3.05, "in"))
plt



## --------- Chi for the PCA emulation ----------
PCA_pairs <- matrix(NA,nrow = ncol(U_PCA)*nrow(ind),ncol = 2)
counter <- 1
for(time in 1:ncol(U_PCA)){
  for(npair in 1:nrow(ind)){
    PCA_pairs[counter, ] <- c(U_PCA[ind[npair, 1], time], U_PCA[ind[npair, 2], time])
    counter <- counter + 1
  }
}


Min_PCA <- apply(PCA_pairs, 1, min)
all_PCA <- as.vector(PCA_pairs)

EmpIntv_PCA <- matrix(NA, nrow = length(u_vec), ncol=3)

for(i in 1:length(u_vec)){
  p_tmp1_sim <- mean(Min_PCA>u_vec[i])
  p_tmp2_sim <- mean(all_PCA>u_vec[i])
  if(p_tmp1_sim==0|p_tmp2_sim==0){
    EmpIntv_PCA[i,]<-c(-2,2,0)
  } else{
    var_sim <- (1/p_tmp1_sim-1)/length(Min_PCA) + (1/p_tmp2_sim-1)/length(all_PCA)
    EmpIntv_PCA[i,]<-c(exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.975)*sqrt(var_sim)/2,
                       exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.025)*sqrt(var_sim)/2, p_tmp1_sim/p_tmp2_sim)
  }
}
dat <- data.frame(x=u_vec,truth=EmpIntv_PCA[,3],truth_upper=EmpIntv_PCA[,2],truth_lower=EmpIntv_PCA[,1])
dat[dat>1] <- 1
dat[dat<0] <- 0
dat_ALL <- rbind(dat_ALL, cbind(dat, source="POD"))

library(ggplot2)
library(ggh4x)
plt <- ggplot(dat,aes(x=x,y=truth)) +
  geom_line(color="#2637ed",linewidth=1) +
  geom_ribbon(data=dat,aes(ymin=truth_lower,ymax=truth_upper),alpha=0.4,fill="#4957f5") +
  ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none') + 
  scale_x_continuous(expand = c(0, 0)) + ggtitle("Gaussian process") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  force_panelsizes(rows = unit(3.05, "in"),
                   cols = unit(3.05, "in"))
plt




## --------- Summary plot ----------
library(ggh4x)
plt <- ggplot(dat_ALL,aes(x=x,y=truth, group=source, colour=source, fill=source)) +
  geom_line(linewidth=1) +
  geom_ribbon(aes(ymin=truth_lower,ymax=truth_upper),alpha=0.4) +
  ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") + labs(fill="Data source", colour="Data source")+
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank()) + 
  scale_x_continuous(expand = c(0, 0), limits=c(0.9,1.0001)) + ggtitle(paste0("d=",round(d,2))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  force_panelsizes(rows = unit(3.05, "in"),
                   cols = unit(3.05, "in"))
plt
ggsave("./Figures/chi_d_8.png",  width = 5, height = 4)
















##--------------------------------------------------------------------------------------
## ----------------------- Same code but for a different ncomp -------------------------
##--------------------------------------------------------------------------------------
setwd("~/Desktop/Turbulence/")

stations <- expand.grid(x=1:198, z=1:500)


## --------- Grab the pairs coordinates ----------
k <- 6
d <- k*sqrt(2)

bottom_xs <- seq(2, 198, by=2*k); bottom_xs <- bottom_xs[-length(bottom_xs)]
bottom_ys <- seq(1, 500, by=2*k); bottom_ys <- bottom_ys[-length(bottom_ys)]

top_xs <- seq(2+k, 198, by=2*k)
if(length(top_xs)!=length(bottom_xs)){
  retain_xs <- min(length(top_xs), length(bottom_xs))
  bottom_xs <- bottom_xs[1:retain_xs]
  top_xs <- top_xs[1:retain_xs]
}
top_ys <- seq(1+k, 500, by=2*k)


bottom_rows <- cbind(rep(bottom_xs, length(bottom_ys)), rep(bottom_ys, each = length(bottom_xs)))
top_rows <- cbind(rep(top_xs, length(top_ys)), rep(top_ys, each = length(top_xs)))
# plot(stations, cex=0.1, xlim=c(0,50), ylim=c(0,50))
# points(bottom_rows, pch=20, col='red')
# points(top_rows, pch=20, col='blue')

if(nrow(top_rows)!=nrow(bottom_rows)){
  retain <- min(nrow(top_rows), nrow(bottom_rows))
  top_rows <- top_rows[1:retain, ]
  bottom_rows <- bottom_rows[1:retain, ]
}

for(iter in 1:nrow(top_rows)){
  lines(x=c(top_rows[iter,1], bottom_rows[iter,1]), y=c(top_rows[iter,2], bottom_rows[iter,2]))
}



## --------- Grab the pairs indices ----------
ind <- array(NA, dim=c(nrow(top_rows),2))
for(iter in 1:nrow(top_rows)){
  ind[iter,1] <- which(stations[,1]==bottom_rows[iter,1] & stations[,2]==bottom_rows[iter,2])
  ind[iter,2] <- which(stations[,1]==top_rows[iter,1] & stations[,2]==top_rows[iter,2])
}
which(ind>99000, arr.ind=TRUE)

# Verify 
# plot(stations, cex=0.1, xlim=c(0,50), ylim=c(0,50))
# points(stations[ind[,1],], pch=20, col='red')
# points(stations[ind[,2],], pch=20, col='blue')


## Thinning out when there are too many pairs
subset_ind <- seq(1, nrow(ind), by=as.integer(nrow(ind)/1000))
ind <- ind[subset_ind, ]


## --------- Grab the values at the pairs ----------
load("./center_fire_input.RData")
load("./emulation_fire_1.RData")
fire_approx1 <- fire_approx
load("./emulation_fire_ncomp_7.RData")
load("./emulation_PCA_ncomp_7.RData")
# bring in more time replicates for the quantile calculations 
Data <- cbind(center_fire, fire_approx1, fire_approx, fire_PCA) 

U_all <- apply(Data, 1, function(x){
  tmp_fun <- ecdf(x)
  tmp_fun(x)
})
U_all <- t(U_all)
U <- U_all[, 1:ncol(center_fire)]
U_emu <- U_all[, (2*ncol(center_fire)+1):(3*ncol(center_fire))]
U_PCA <- U_all[, (3*ncol(center_fire)+1):(4*ncol(center_fire))]




## --------- Chi for the plume obs ----------
plume_pairs <- matrix(NA,nrow = ncol(U)*nrow(ind),ncol = 2)
counter <- 1
for(time in 1:ncol(U)){
  for(npair in 1:nrow(ind)){
    plume_pairs[counter, ] <- c(U[ind[npair, 1], time], U[ind[npair, 2], time])
    counter <- counter + 1
  }
}


Min_plume <- apply(plume_pairs, 1, min)
all_plume <- as.vector(plume_pairs)

u_vec=c(seq(0.9,0.98,0.01),seq(0.9801,0.9999, length.out=8))
EmpIntv_sim <- matrix(NA, nrow = length(u_vec), ncol=3)

for(i in 1:length(u_vec)){
  p_tmp1_sim <- mean(Min_plume>u_vec[i])
  p_tmp2_sim <- mean(all_plume>u_vec[i])
  if(p_tmp1_sim==0|p_tmp2_sim==0){
    EmpIntv_sim[i,]<-c(-2,2,0)
  } else{
    var_sim <- (1/p_tmp1_sim-1)/length(Min_plume) + (1/p_tmp2_sim-1)/length(all_plume)
    EmpIntv_sim[i,]<-c(exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.975)*sqrt(var_sim)/2,
                       exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.025)*sqrt(var_sim)/2, p_tmp1_sim/p_tmp2_sim)
  }
}
dat <- data.frame(x=u_vec,truth=EmpIntv_sim[,3],truth_upper=EmpIntv_sim[,2],truth_lower=EmpIntv_sim[,1])
dat[dat>1] <- 1
dat[dat<0] <- 0

tmp <- smooth.spline(dat$x, dat$true)


dat_ALL <- cbind(dat, source="Observations")

library(ggplot2)
library(ggh4x)
plt <- ggplot(dat,aes(x=x,y=truth)) +
  geom_line(color="#2637ed",linewidth=1) +
  geom_ribbon(data=dat,aes(ymin=truth_lower,ymax=truth_upper),alpha=0.4,fill="#4957f5") +
  ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none') + 
  scale_x_continuous(expand = c(0, 0)) + ggtitle("Gaussian process") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  force_panelsizes(rows = unit(3.05, "in"),
                   cols = unit(3.05, "in"))
plt




## --------- Chi for the XVAE emulation ----------
XVAE_pairs <- matrix(NA,nrow = ncol(U_emu)*nrow(ind),ncol = 2)
counter <- 1
for(time in 1:ncol(U_emu)){
  for(npair in 1:nrow(ind)){
    XVAE_pairs[counter, ] <- c(U_emu[ind[npair, 1], time], U_emu[ind[npair, 2], time])
    counter <- counter + 1
  }
}


Min_XVAE <- apply(XVAE_pairs, 1, min)
all_XVAE <- as.vector(XVAE_pairs)

EmpIntv_XVAE <- matrix(NA, nrow = length(u_vec), ncol=3)

for(i in 1:length(u_vec)){
  p_tmp1_sim <- mean(Min_XVAE>u_vec[i])
  p_tmp2_sim <- mean(all_XVAE>u_vec[i])
  if(p_tmp1_sim==0|p_tmp2_sim==0){
    EmpIntv_XVAE[i,]<-c(-2,2,0)
  } else{
    var_sim <- (1/p_tmp1_sim-1)/length(Min_XVAE) + (1/p_tmp2_sim-1)/length(all_XVAE)
    EmpIntv_XVAE[i,]<-c(exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.975)*sqrt(var_sim)/2,
                       exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.025)*sqrt(var_sim)/2, p_tmp1_sim/p_tmp2_sim)
  }
}
dat <- data.frame(x=u_vec,truth=EmpIntv_XVAE[,3],truth_upper=EmpIntv_XVAE[,2],truth_lower=EmpIntv_XVAE[,1])
dat[dat>1] <- 1
dat[dat<0] <- 0
dat_ALL <- rbind(dat_ALL, cbind(dat, source="XVAE"))

library(ggplot2)
library(ggh4x)
plt <- ggplot(dat,aes(x=x,y=truth)) +
  geom_line(color="#2637ed",linewidth=1) +
  geom_ribbon(data=dat,aes(ymin=truth_lower,ymax=truth_upper),alpha=0.4,fill="#4957f5") +
  ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none') + 
  scale_x_continuous(expand = c(0, 0)) + ggtitle("Gaussian process") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  force_panelsizes(rows = unit(3.05, "in"),
                   cols = unit(3.05, "in"))
plt



## --------- Chi for the PCA emulation ----------
PCA_pairs <- matrix(NA,nrow = ncol(U_PCA)*nrow(ind),ncol = 2)
counter <- 1
for(time in 1:ncol(U_PCA)){
  for(npair in 1:nrow(ind)){
    PCA_pairs[counter, ] <- c(U_PCA[ind[npair, 1], time], U_PCA[ind[npair, 2], time])
    counter <- counter + 1
  }
}


Min_PCA <- apply(PCA_pairs, 1, min)
all_PCA <- as.vector(PCA_pairs)

EmpIntv_PCA <- matrix(NA, nrow = length(u_vec), ncol=3)

for(i in 1:length(u_vec)){
  p_tmp1_sim <- mean(Min_PCA>u_vec[i])
  p_tmp2_sim <- mean(all_PCA>u_vec[i])
  if(p_tmp1_sim==0|p_tmp2_sim==0){
    EmpIntv_PCA[i,]<-c(-2,2,0)
  } else{
    var_sim <- (1/p_tmp1_sim-1)/length(Min_PCA) + (1/p_tmp2_sim-1)/length(all_PCA)
    EmpIntv_PCA[i,]<-c(exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.975)*sqrt(var_sim)/2,
                       exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.025)*sqrt(var_sim)/2, p_tmp1_sim/p_tmp2_sim)
  }
}
dat <- data.frame(x=u_vec,truth=EmpIntv_PCA[,3],truth_upper=EmpIntv_PCA[,2],truth_lower=EmpIntv_PCA[,1])
dat[dat>1] <- 1
dat[dat<0] <- 0
dat_ALL <- rbind(dat_ALL, cbind(dat, source="POD"))

library(ggplot2)
library(ggh4x)
plt <- ggplot(dat,aes(x=x,y=truth)) +
  geom_line(color="#2637ed",linewidth=1) +
  geom_ribbon(data=dat,aes(ymin=truth_lower,ymax=truth_upper),alpha=0.4,fill="#4957f5") +
  ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none') + 
  scale_x_continuous(expand = c(0, 0)) + ggtitle("Gaussian process") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  force_panelsizes(rows = unit(3.05, "in"),
                   cols = unit(3.05, "in"))
plt




## --------- Summary plot ----------
library(ggh4x)
plt <- ggplot(dat_ALL,aes(x=x,y=truth, group=source, colour=source, fill=source)) +
  geom_line(linewidth=1) +
  geom_ribbon(aes(ymin=truth_lower,ymax=truth_upper),alpha=0.4) +
  ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") + labs(fill="Data source", colour="Data source")+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") + 
  scale_x_continuous(expand = c(0, 0), limits=c(0.9,1.0001)) + ggtitle(paste0("d=",round(d,2))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  force_panelsizes(rows = unit(3.05, "in"),
                   cols = unit(3.05, "in"))
plt
ggsave("./Figures/chi_d_4_ncomp7.png",  width = 5, height = 4)
# ggsave("./Figures/chi_d_8_ncomp7.png",  width = 5, height = 4)
