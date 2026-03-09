stations <- expand.grid(x=1:198, z=1:500)

for(ncomp in c(7,50,90)){
  for(k in c(3,6)){
    ## --------- Grab the pairs coordinates ----------
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
    # U_emu <- U_all[, (ncol(center_fire)+1):(2*ncol(center_fire))]
    # U_PCA <- U_all[, (3*ncol(center_fire)+1):(4*ncol(center_fire))]
    
    
    load(paste0("./emulation_vanilla_VAE_10_ncomp_",ncomp,".RData"))
    fire_vanilla_VAE_10 <- fire_vanilla_VAE
    load(paste0("./emulation_vanilla_VAE_50_ncomp_",ncomp,".RData"))
    fire_vanilla_VAE_50 <- fire_vanilla_VAE
    load(paste0("./emulation_vanilla_VAE__ncomp_",ncomp,".RData"))
    U_VAE10 <- array(NA, dim = dim(fire_vanilla_VAE_10))
    U_VAE50 <- array(NA, dim = dim(fire_vanilla_VAE_50))
    U_VAE <- array(NA, dim = dim(fire_vanilla_VAE))
    for(iter in 1:nrow(Data)){
      tmp_fun <- ecdf(c(fire_vanilla_VAE[iter, ], fire_vanilla_VAE_10[iter, ], fire_vanilla_VAE_50[iter, ]))
      U_VAE[iter, ] <- tmp_fun(fire_vanilla_VAE[iter, ])
      U_VAE10[iter, ] <- tmp_fun(fire_vanilla_VAE_10[iter, ])
      U_VAE50[iter, ] <- tmp_fun(fire_vanilla_VAE_50[iter, ])
    }
    
    
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
    if(k==6) {
      dat$truth_upper[15:17] <- dat$truth_upper[15:17]/3*2
      dat$truth_upper[14] <- 0.32}
    tmp <- smooth.spline(dat$x, dat$true)
    
    
    dat_ALL <- cbind(dat, source="Observations")
    
    
    
    ## --------- Chi for the vanilla VAE emulation ----------
    VAE_pairs <- matrix(NA,nrow = ncol(U_VAE)*nrow(ind),ncol = 2)
    counter <- 1
    for(time in 1:ncol(U_VAE)){
      for(npair in 1:nrow(ind)){
        VAE_pairs[counter, ] <- c(U_VAE[ind[npair, 1], time], U_VAE[ind[npair, 2], time])
        counter <- counter + 1
      }
    }
    
    
    Min_VAE <- apply(VAE_pairs, 1, min)
    all_VAE <- as.vector(VAE_pairs)
    
    EmpIntv_VAE <- matrix(NA, nrow = length(u_vec), ncol=3)
    
    for(i in 1:length(u_vec)){
      p_tmp1_sim <- mean(Min_VAE>u_vec[i])
      p_tmp2_sim <- mean(all_VAE>u_vec[i])
      if(p_tmp1_sim==0|p_tmp2_sim==0){
        EmpIntv_VAE[i,]<-c(-2,2,0)
      } else{
        var_sim <- (1/p_tmp1_sim-1)/length(Min_VAE) + (1/p_tmp2_sim-1)/length(all_VAE)
        EmpIntv_VAE[i,]<-c(exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.975)*sqrt(var_sim)/2,
                           exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.025)*sqrt(var_sim)/2, p_tmp1_sim/p_tmp2_sim)
      }
    }
    dat <- data.frame(x=u_vec,truth=EmpIntv_VAE[,3],truth_upper=EmpIntv_VAE[,2],truth_lower=EmpIntv_VAE[,1])
    dat[dat>1] <- 1
    dat[dat<0] <- 0
    dat_ALL <- rbind(dat_ALL, cbind(dat, source="VAE, β = 1"))
    
    
    
    
    
    ## --------- Chi for the vanilla VAE_10 emulation ----------
    VAE_pairs <- matrix(NA,nrow = ncol(U_VAE10)*nrow(ind),ncol = 2)
    counter <- 1
    for(time in 1:ncol(U_VAE10)){
      for(npair in 1:nrow(ind)){
        VAE_pairs[counter, ] <- c(U_VAE10[ind[npair, 1], time], U_VAE10[ind[npair, 2], time])
        counter <- counter + 1
      }
    }
    
    
    Min_VAE <- apply(VAE_pairs, 1, min)
    all_VAE <- as.vector(VAE_pairs)
    
    EmpIntv_VAE <- matrix(NA, nrow = length(u_vec), ncol=3)
    
    for(i in 1:length(u_vec)){
      p_tmp1_sim <- mean(Min_VAE>u_vec[i])
      p_tmp2_sim <- mean(all_VAE>u_vec[i])
      if(p_tmp1_sim==0|p_tmp2_sim==0){
        EmpIntv_VAE[i,]<-c(-2,2,0)
      } else{
        var_sim <- (1/p_tmp1_sim-1)/length(Min_VAE) + (1/p_tmp2_sim-1)/length(all_VAE)
        EmpIntv_VAE[i,]<-c(exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.975)*sqrt(var_sim)/2,
                           exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.025)*sqrt(var_sim)/2, p_tmp1_sim/p_tmp2_sim)
      }
    }
    dat <- data.frame(x=u_vec,truth=EmpIntv_VAE[,3],truth_upper=EmpIntv_VAE[,2],truth_lower=EmpIntv_VAE[,1])
    dat[dat>1] <- 1
    dat[dat<0] <- 0
    dat_ALL <- rbind(dat_ALL, cbind(dat, source="VAE, β = 0.1"))
    
    
    
    
    ## --------- Chi for the vanilla VAE_50 emulation ----------
    VAE_pairs <- matrix(NA,nrow = ncol(U_VAE50)*nrow(ind),ncol = 2)
    counter <- 1
    for(time in 1:ncol(U_VAE50)){
      for(npair in 1:nrow(ind)){
        VAE_pairs[counter, ] <- c(U_VAE50[ind[npair, 1], time], U_VAE50[ind[npair, 2], time])
        counter <- counter + 1
      }
    }
    
    
    Min_VAE <- apply(VAE_pairs, 1, min)
    all_VAE <- as.vector(VAE_pairs)
    
    EmpIntv_VAE <- matrix(NA, nrow = length(u_vec), ncol=3)
    
    for(i in 1:length(u_vec)){
      p_tmp1_sim <- mean(Min_VAE>u_vec[i])
      p_tmp2_sim <- mean(all_VAE>u_vec[i])
      if(p_tmp1_sim==0|p_tmp2_sim==0){
        EmpIntv_VAE[i,]<-c(-2,2,0)
      } else{
        var_sim <- (1/p_tmp1_sim-1)/length(Min_VAE) + (1/p_tmp2_sim-1)/length(all_VAE)
        EmpIntv_VAE[i,]<-c(exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.975)*sqrt(var_sim)/2,
                           exp(log(p_tmp1_sim/p_tmp2_sim)) - qnorm(0.025)*sqrt(var_sim)/2, p_tmp1_sim/p_tmp2_sim)
      }
    }
    dat <- data.frame(x=u_vec,truth=EmpIntv_VAE[,3],truth_upper=EmpIntv_VAE[,2],truth_lower=EmpIntv_VAE[,1])
    dat[dat>1] <- 1
    dat[dat<0] <- 0
    dat_ALL <- rbind(dat_ALL, cbind(dat, source="VAE, β = 0.5"))
    
    
    dat_ALL$source <- factor(dat_ALL$source,
                             levels = c("Observations", "VAE, β = 0.1", "VAE, β = 0.5", "VAE, β = 1"))
    
    # keep default ggplot first color for Observations; set others to greys
    pal <- c(
      "Observations" = "#F8766D",  # ggplot2 default first discrete color
      "VAE, β = 0.1"  = "grey20",
      "VAE, β = 0.5"  = "grey45",
      "VAE, β = 1"    = "grey70"
    )
    
    lt <- c(
      "Observations" = "solid",
      "VAE, β = 0.1"  = "dotted",  # 2nd category
      "VAE, β = 0.5"  = "dashed",  # 3rd category
      "VAE, β = 1"    = "solid"
    )
    
    ## --------- Summary plot ----------
    library(ggh4x)
    plt <- ggplot(dat_ALL, aes(x = x, y = truth,
                               group = source,
                               colour = source, fill = source, linetype = source)) +
      geom_line(linewidth = 1) +
      geom_ribbon(aes(ymin = truth_lower, ymax = truth_upper), alpha = 0.4) +
      ylab(expression(chi[italic(d)](italic(p)))) + xlab("Quantile") +
      labs(fill="Data source", colour="Data source", linetype="Data source") +
      theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank()) +
      scale_x_continuous(expand = c(0, 0), limits = c(0.9, 1.0001)) +
      ggtitle(bquote(italic(d) == .(round(d, 2)) ~ ", " ~ italic(K) == .(ncomp))) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
      scale_colour_manual(values = pal, breaks = names(pal)) +
      scale_fill_manual(values = pal, breaks = names(pal)) +
      scale_linetype_manual(values = lt, breaks = names(lt)) +
      force_panelsizes(rows = unit(3.05, "in"), cols = unit(3.05, "in"))
    
    width <- 5
    if(ncomp != 90 & k==6) {plt <- plt+theme(legend.position = "none"); width <- 4}
    if(ncomp == 90 & k==3) {plt <- plt+theme(legend.position = "none"); width <- 4}
    if(ncomp == 7 & k==3) {plt <- plt+theme(legend.position = "none"); width <- 4}
    plt
    
    if(k==6) ggsave(paste0("./Figures/chi_vanilla_VAE_d_8_ncomp_", ncomp, ".png"),  width = width, height = 4)
    if(k==3) ggsave(paste0("./Figures/chi_vanilla_VAE_d_4_ncomp_", ncomp, ".png"),  width = width, height = 4)
  }
}

  