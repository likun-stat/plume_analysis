library(torch)
library(VGAM)
source("./utils.R")

###### ---------------------------------------------------------------------- ######
###### ----------------------------- Load in Data  -------------------------- ######
###### ---------------------------------------------------------------------- ######
stations <- expand.grid(x=1:37, z=1:600)

center_fire <- read.csv("../../data/Vertical velocity/xz_all.csv", row.names = 1)
center_fire <- as.matrix(center_fire)
threshold <- quantile(center_fire,0.05)
metric_CV <- function(emulation, X, threshold){
  tmp_where <- (X < threshold)
  return(sum((emulation[tmp_where]-X[tmp_where])^2))
}

###### ---------------------------------------------------------------------- ######
###### ---------------------------- Load NMF results ------------------------ ######
###### ---------------------------------------------------------------------- ######
ncomp <- 50
features <- read.csv(file = paste0("../../data/Vertical velocity/features",ncomp,".csv"), header = FALSE, row.names = NULL)
components <- read.csv(file = paste0("../../data/Vertical velocity/components",ncomp,".csv"), header = FALSE, row.names = NULL)


alpha = 0.5; tau <- 0.1; m <- 0.85
###### ---------------------------------------------------------------------- ######
###### ---------------------------- CROSS VALIDATION ------------------------ ######
###### ---------------------------------------------------------------------- ######
Metrix <- rep(NA, 10)
for(fold in 1:10){
  fold_ind <- seq(fold*10-9,fold*10,by=1)
  
  ###### ---------------------------------------------------------------------- ######
  ###### ------------------------ Marginal transformation --------------------- ######
  ###### ---------------------------------------------------------------------- ######
  X <- center_fire[, -fold_ind]
  
  X_holdout <- center_fire[, fold_ind]
  
  
  ###### ---------------------------------------------------------------------- ######
  ###### ----------------------- First initial guess -------------------------- ######
  ###### ---------------------------------------------------------------------- ######
  ## -------------------- Initial guess for the latent Z variables --------------------
  W_alpha <- data.matrix(features)
  Z_approx <- data.matrix(components[,-fold_ind]) + 1e-7
  Z_holdout <- data.matrix(components[,fold_ind]) + 1e-7
  
  Y_star <- (W_alpha)%*%(Z_approx)
  
  ## -------------------- Initial values for the weights --------------------
  # Z = w_1 %*% X or Z^T = X^T %*% w_1^T
  X_tensor <- torch_tensor(X,dtype=torch_float())
  W_alpha_tensor <- torch_tensor(W_alpha,dtype=torch_float())
  k <- ncol(W_alpha)
  tmp <- qr.solve(a=t(X), b=t(Z_approx))
  w_1 <- t(tmp)
  w_1 <- torch_tensor(w_1,dtype=torch_float(),requires_grad = TRUE)
  b_1 <- matrix(rep(1e-6, k), ncol=1)
  b_1 <- torch_tensor(b_1,dtype=torch_float(),requires_grad = TRUE)
  
  w_2 <- diag(k)
  w_2 <- torch_tensor(w_2,dtype=torch_float(),requires_grad = TRUE)
  b_2 <- matrix(rep(1e-6, k), ncol=1)
  b_2 <- torch_tensor(b_2,dtype=torch_float(),requires_grad = TRUE)
  
  w_4 <- diag(k)
  w_4 <- torch_tensor(w_4,dtype=torch_float(),requires_grad = TRUE)
  b_4 <- matrix(rep(1e-6, k), ncol=1)
  b_4 <- torch_tensor(b_4,dtype=torch_float(),requires_grad = TRUE)
  
  w_3 <- 0*diag(k)
  w_3 <- torch_tensor(w_3,dtype=torch_float(),requires_grad = TRUE)
  b_3 <- matrix(rep(-15,k), ncol=1)
  b_3 <- torch_tensor(b_3,dtype=torch_float(),requires_grad = TRUE)
  
  
  h <- w_1$mm(X_tensor)$add(b_1)$relu()
  h_1 <- w_2$mm(h)$add(b_2)$relu()
  sigma_sq_vec <- w_3$mm(h_1)$add(b_3)$exp()
  mu <- w_4$mm(h_1)$add(b_4)$relu()
  
  ## -------------------- Re-parameterization trick --------------------
  n.t <- ncol(X)
  n.s <- nrow(X)
  Epsilon <- matrix(abs(rnorm(k*n.t))+0.1, nrow=k)
  Epsilon <- torch_tensor(Epsilon,dtype=torch_float())
  v_t <- mu + sqrt(sigma_sq_vec)*Epsilon
  b_8 <- as_array((W_alpha_tensor))%*%as_array(v_t)-X 
  b_8 <- array(-relu(b_8), dim=dim(b_8))
  b_8 <- torch_tensor(b_8,dtype=torch_float(), requires_grad = TRUE)
  
  ## -- Ensure Y > 0
  y_approx <-  (W_alpha_tensor)$mm(v_t) + b_8
  which(as_array(y_approx)<0, arr.ind=TRUE)
  
  ## -- Ensure X/Y > 1
  indices_tmp <- which(as_array(X_tensor$divide(y_approx)< m), arr.ind = TRUE)
  if(nrow(indices_tmp)>0) {
    for(iter in 1:nrow(indices_tmp)){
      X_tensor[indices_tmp[iter,1], indices_tmp[iter,2]] <- y_approx[indices_tmp[iter,1], indices_tmp[iter,2]]*(m+0.1)
    }
  }
  
  
  
  ## -------------------- Add in the Frechet noise --------------------
  library(VGAM)
  Epsilon_frechet_indep <- matrix(VGAM::rfrechet(n.s*n.t, location = m, shape=3, scale = tau), nrow=n.s)
  X_approx <- matrix(NA, nrow=n.s, ncol=n.t)
  for (iter in 1:n.t){
    X_approx[,iter] <- as_array(Epsilon_frechet_indep[,iter]* y_approx[,iter])
  }
  
  
  ## -------------------- Visualize the X_approx --------------------
  ind <- 3
  tmp_range <- range(c(log(as_array(y_approx[,ind]))), log(X[,ind]))
  ggplot(stations) + geom_raster(aes(x=x, y=z, fill=log(X[,ind]))) +
    scale_fill_gradientn(colours = topo.colors(100), name = paste0('Original replicate #', ind), na.value = NA, limits=tmp_range) 
  
  ggplot(stations) + geom_raster(aes(x=x, y=z, fill=log(as_array(y_approx[,ind])))) +
    scale_fill_gradientn(colours = topo.colors(100), name = paste0('Smooth replicate #', ind), na.value = NA, limits=tmp_range) 
  
  
  
  
  ## -------------------- w_prime and b_prime for theta --------------------
  w_1_prime <- as_array(w_1) #matrix(rnorm(k*n.s,0,0.001), nrow=k)
  w_1_prime <- torch_tensor(w_1_prime,dtype=torch_float(),requires_grad = TRUE)
  w_2_prime <- matrix(diag(k), nrow=k)
  w_2_prime <- torch_tensor(w_2_prime,dtype=torch_float(),requires_grad = TRUE)
  w_3_prime <- matrix(rep(0,k*k), nrow=k)
  w_3_prime <- torch_tensor(w_3_prime,dtype=torch_float(),requires_grad = TRUE)
  w_4_prime <- matrix(diag(k), nrow=k)
  w_4_prime <- torch_tensor(w_4_prime,dtype=torch_float(),requires_grad = TRUE)
  b_1_prime <- matrix(rep(1e-6, k), ncol=1)
  b_1_prime <- torch_tensor(b_1_prime,dtype=torch_float(),requires_grad = TRUE)
  b_2_prime <- matrix(rep(1e-6, k), ncol=1)
  b_2_prime <- torch_tensor(b_2_prime,dtype=torch_float(),requires_grad = TRUE)
  b_3_prime <- matrix(rep(-10,k), ncol=1)
  b_3_prime <- torch_tensor(b_3_prime,dtype=torch_float(),requires_grad = TRUE)
  b_4_prime <- matrix(rep(1e-6, k), ncol=1)
  b_4_prime <- torch_tensor(b_4_prime,dtype=torch_float(),requires_grad = TRUE)
  
  w_5 <- diag(k) #matrix(rnorm(k*k,0,0.001), nrow=k)
  w_5 <- torch_tensor(w_5,dtype=torch_float(),requires_grad = TRUE)
  w_6 <- diag(k) #matrix(rnorm(k*k,0,0.001), nrow=k)
  w_6 <- torch_tensor(w_6,dtype=torch_float(),requires_grad = TRUE)
  w_7 <- diag(k) #matrix(rnorm(k*k,0,0.0001), nrow=k)
  w_7 <- torch_tensor(w_7,dtype=torch_float(),requires_grad = TRUE)
  b_5 <- matrix(rep(1e-6,k), ncol=1) #matrix(rnorm(k), ncol=1)
  b_5 <- torch_tensor(b_5,dtype=torch_float(),requires_grad = TRUE)
  b_6 <- matrix(rep(1e-6,k), ncol=1) #matrix(rnorm(k), ncol=1)
  b_6 <- torch_tensor(b_6,dtype=torch_float(),requires_grad = TRUE)
  b_7 <- matrix(rep(1e-6,k), ncol=1) #matrix(rnorm(k,0,0.0001), ncol=1)
  b_7 <- torch_tensor(b_7,dtype=torch_float(),requires_grad = TRUE)
  
  w_1_velocity <- torch_zeros(w_1$size())
  w_2_velocity <- torch_zeros(w_2$size())
  w_3_velocity <- torch_zeros(w_3$size())
  w_4_velocity <- torch_zeros(w_4$size())
  b_1_velocity <- torch_zeros(b_1$size())
  b_2_velocity <- torch_zeros(b_2$size())
  b_3_velocity <- torch_zeros(b_3$size())
  b_4_velocity <- torch_zeros(b_4$size())
  w_1_prime_velocity <- torch_zeros(w_1_prime$size())
  w_2_prime_velocity <- torch_zeros(w_2_prime$size())
  w_3_prime_velocity <- torch_zeros(w_3_prime$size())
  w_4_prime_velocity <- torch_zeros(w_4_prime$size())
  b_1_prime_velocity <- torch_zeros(b_1_prime$size())
  b_2_prime_velocity <- torch_zeros(b_2_prime$size())
  b_3_prime_velocity <- torch_zeros(b_3_prime$size())
  b_4_prime_velocity <- torch_zeros(b_4_prime$size())
  
  w_5_velocity <- torch_zeros(w_5$size())
  w_6_velocity <- torch_zeros(w_6$size())
  w_7_velocity <- torch_zeros(w_7$size())
  b_5_velocity <- torch_zeros(b_5$size())
  b_6_velocity <- torch_zeros(b_6$size())
  b_7_velocity <- torch_zeros(b_7$size())
  
  b_8_velocity <- torch_zeros(b_8$size())
  
  Epsilon_prime <- t(mvtnorm::rmvnorm(n.t, mean=rep(0, k), sigma = diag(rep(1, k))))
  Epsilon_prime <- torch_tensor(Epsilon_prime,dtype=torch_float())
  
  ###### ---------------------------------------------------------------------- ######
  ###### --------------------   Loss and Optimization ------------------------- ######
  ###### ---------------------------------------------------------------------- ######
  
  ###### ----------------------  network parameters --------------------------- ######
  ###### ---------------------- Training hyperparameters ---------------------- ######
  mu_momentum <- 0.9              # constant momentum
  eta_max <- 1e-20
  eta_min <- 1e-30
  T_w <- 5000                     # warmup epochs
  max_epochs <- 10000            # large cap; early stopping will usually stop sooner
  delta_stop <- 1e-7             # stopping threshold based on rolling ELBO means
  window_elbo <- 100             # compare latest 100 vs previous 100 epochs
  
  lrelu <- nn_leaky_relu(-0.01)
  
  # store ELBO history
  elbo_hist <- rep(NA_real_, max_epochs)
  
  # learning-rate schedule: linear warmup, then half-cosine decay
  get_lr <- function(epoch, T_w, max_epochs, eta_max, eta_min) {
    if (epoch <= T_w) {
      return(eta_max * epoch / T_w)
    } else {
      progress <- (epoch - T_w) / (max_epochs - T_w)
      progress <- min(max(progress, 0), 1)
      return(eta_min + 0.5 * (eta_max - eta_min) * (1 + cos(pi * progress)))
    }
  }
  
  ###### ---------------------------------------------------------------------- ######
  ###### --------------------   Loss and Optimization ------------------------- ######
  ###### ---------------------------------------------------------------------- ######
  
  n <- 1e3
  vec <- Zolo_A(pi*seq(1/2, n - 1/2, 1) / n, alpha)
  Zolo_vec <- torch_tensor(matrix(vec, nrow = 1, ncol = n), dtype = torch_float(), requires_grad = FALSE)
  Zolo_vec_double <- torch_tensor(Zolo_vec, dtype = torch_float64(), requires_grad = FALSE)
  
  const  <- 1 / (1 - alpha)
  const1 <- 1 / (1 - alpha) - 1
  const3 <- log(const1)
  
  for (t in 1:max_epochs) {
    
    ## current learning rate from warmup + cosine decay schedule
    learning_rate <- get_lr(
      epoch = t,
      T_w = T_w,
      max_epochs = max_epochs,
      eta_max = eta_max,
      eta_min = eta_min
    )
    
    ### -------- Forward pass ---------
    Epsilon <- matrix(abs(rnorm(k * n.t)) + 0.1, nrow = k)
    Epsilon_prime <- t(mvtnorm::rmvnorm(n.t, mean = rep(0, k), sigma = diag(rep(1, k))))
    
    Epsilon <- torch_tensor(Epsilon, dtype = torch_float())
    Epsilon_prime <- torch_tensor(Epsilon_prime, dtype = torch_float())
    
    ### -------- Encoder for v_t --------
    h <- w_1$mm(X_tensor)$add(b_1)$relu()
    h_1 <- w_2$mm(h)$add(b_2)$relu()
    sigma_sq_vec <- w_3$mm(h_1)$add(b_3)$exp()
    mu <- w_4$mm(h_1)$add(b_4)$relu()
    
    ### -------- Encoder for v_t_prime --------
    h_prime <- w_1_prime$mm(X_tensor)$add(b_1_prime)$relu()
    h_1_prime <- w_2_prime$mm(h_prime)$add(b_2_prime)$relu()
    
    ### -------- Activation via Laplace transformation --------
    h_1_prime_laplace <- h_1_prime$multiply(-0.2)$exp()$mean(dim = 2)
    h_1_prime_t <- h_1_prime_laplace$log()$multiply(-1)
    h_1_prime_to_theta <- (0.2 - h_1_prime_t$pow(2))$pow(2)$divide(4 * h_1_prime_t$pow(2))$view(c(k, 1))
    theta_propagate <- h_1_prime_to_theta$expand(c(k, n.t))
    
    sigma_sq_vec_prime <- w_3_prime$mm(theta_propagate)$add(b_3_prime)$exp()
    mu_prime <- w_4_prime$mm(theta_propagate)$add(b_4_prime)
    
    ### -------- Re-parameterization trick --------
    v_t <- mu + torch_sqrt(sigma_sq_vec) * Epsilon
    v_t_prime <- mu_prime + torch_sqrt(sigma_sq_vec_prime) * Epsilon_prime
    
    ### -------- Decoder --------
    l <- w_5$mm(v_t_prime)$add(b_5)$relu()
    l_1 <- w_6$mm(l)$add(b_6)$relu()
    theta_t <- w_7$mm(l_1)$add(b_7)$relu()
    
    y_star <- lrelu(W_alpha_tensor$mm(v_t)$add(b_8)) + 1.73266e-13
    
    ### -------- ELBO -------- 
    ## Part 1
    standardized <- X_tensor$divide(y_star)$sub(m)
    leak <- as_array(sum(standardized < 0))
    if (leak > 0 && leak <= 1e4) standardized$abs_()
    
    part1 <- -3 * standardized$log()$sum() -
      y_star$log()$sum() -
      tau * standardized$pow(-2)$sum()
    
    ## Part 2
    V_t <- v_t$view(c(k * n.t, 1))
    Theta_t <- theta_t$view(c(k * n.t, 1))
    
    part_log_v1 <- V_t$pow(-const)$mm(Zolo_vec)
    part_log_v2 <- (-V_t$pow(-const1)$mm(Zolo_vec))$exp()
    part_log_v3 <- Theta_t$pow(alpha) - Theta_t$mul(V_t)
    
    part2 <- (part_log_v1$mul(part_log_v2)$mean(dim = 2)$log() +
                part_log_v3$view(k * n.t)$add(const3))$sum()
    
    if (as_array(part2) == -Inf) {
      part2_tmp <- part_log_v1$log()[, 1] + (-V_t$pow(-const1)$mm(Zolo_vec))[, 1]
      part2 <- (part2_tmp + part_log_v3$view(k * n.t)$add(const3))$sum()
    }
    
    ## Part 3
    part3 <- Epsilon$pow(2)$sum() / 2 +
      Epsilon_prime$pow(2)$sum() / 2 +
      sigma_sq_vec$log()$sum() +
      sigma_sq_vec_prime$log()$sum()
    
    elbo <- (part1 + part2 + part3) / (n.s * n.t)   # maximize this
    
    if (!is.finite(elbo$item())) {
      cat("Stopping: non-finite ELBO at epoch", t, "\n")
      break
    }
    if (as.numeric(torch_isnan(elbo)) == 1) {
      cat("Stopping: NaN ELBO at epoch", t, "\n")
      break
    }
    
    elbo_hist[t] <- elbo$item()
    
    if (t %% 100 == 0) {
      cat(sprintf("Epoch: %d   ELBO: %.10f   lr: %.3e\n", t, elbo$item(), learning_rate))
    }
    
    ### -------- Backpropagation --------
    # We maximize ELBO, so backprop through -ELBO
    loss_to_minimize <- -elbo
    loss_to_minimize$backward()
    
    ### -------- Update weights with momentum --------
    with_no_grad({
      w_1_velocity <<- mu_momentum * w_1_velocity - learning_rate * w_1$grad
      w_1$add_(w_1_velocity)
      
      w_2_velocity <<- mu_momentum * w_2_velocity - learning_rate * w_2$grad
      w_2$add_(w_2_velocity)
      
      w_3_velocity <<- mu_momentum * w_3_velocity - learning_rate * w_3$grad
      w_3$add_(w_3_velocity)
      
      w_4_velocity <<- mu_momentum * w_4_velocity - learning_rate * w_4$grad
      w_4$add_(w_4_velocity)
      
      b_1_velocity <<- mu_momentum * b_1_velocity - learning_rate * b_1$grad
      b_1$add_(b_1_velocity)
      
      b_2_velocity <<- mu_momentum * b_2_velocity - learning_rate * b_2$grad
      b_2$add_(b_2_velocity)
      
      b_3_velocity <<- mu_momentum * b_3_velocity - learning_rate * b_3$grad
      b_3$add_(b_3_velocity)
      
      b_4_velocity <<- mu_momentum * b_4_velocity - learning_rate * b_4$grad
      b_4$add_(b_4_velocity)
      
      w_1_prime_velocity <<- mu_momentum * w_1_prime_velocity - learning_rate * w_1_prime$grad
      w_1_prime$add_(w_1_prime_velocity)
      
      w_2_prime_velocity <<- mu_momentum * w_2_prime_velocity - learning_rate * w_2_prime$grad
      w_2_prime$add_(w_2_prime_velocity)
      
      w_3_prime_velocity <<- mu_momentum * w_3_prime_velocity - learning_rate * w_3_prime$grad
      w_3_prime$add_(w_3_prime_velocity)
      
      w_4_prime_velocity <<- mu_momentum * w_4_prime_velocity - learning_rate * w_4_prime$grad
      w_4_prime$add_(w_4_prime_velocity)
      
      b_1_prime_velocity <<- mu_momentum * b_1_prime_velocity - learning_rate * b_1_prime$grad
      b_1_prime$add_(b_1_prime_velocity)
      
      b_2_prime_velocity <<- mu_momentum * b_2_prime_velocity - learning_rate * b_2_prime$grad
      b_2_prime$add_(b_2_prime_velocity)
      
      b_3_prime_velocity <<- mu_momentum * b_3_prime_velocity - learning_rate * b_3_prime$grad
      b_3_prime$add_(b_3_prime_velocity)
      
      b_4_prime_velocity <<- mu_momentum * b_4_prime_velocity - learning_rate * b_4_prime$grad
      b_4_prime$add_(b_4_prime_velocity)
      
      w_5_velocity <<- mu_momentum * w_5_velocity - learning_rate * w_5$grad
      w_5$add_(w_5_velocity)
      
      w_6_velocity <<- mu_momentum * w_6_velocity - learning_rate * w_6$grad
      w_6$add_(w_6_velocity)
      
      w_7_velocity <<- mu_momentum * w_7_velocity - learning_rate * w_7$grad
      w_7$add_(w_7_velocity)
      
      b_5_velocity <<- mu_momentum * b_5_velocity - learning_rate * b_5$grad
      b_5$add_(b_5_velocity)
      
      b_6_velocity <<- mu_momentum * b_6_velocity - learning_rate * b_6$grad
      b_6$add_(b_6_velocity)
      
      b_7_velocity <<- mu_momentum * b_7_velocity - learning_rate * b_7$grad
      b_7$add_(b_7_velocity)
      
      b_8_velocity <<- mu_momentum * b_8_velocity - learning_rate * b_8$grad
      b_8$add_(b_8_velocity)
      
      ## zero gradients
      w_1$grad$zero_();  w_2$grad$zero_();  w_3$grad$zero_();  w_4$grad$zero_()
      b_1$grad$zero_();  b_2$grad$zero_();  b_3$grad$zero_();  b_4$grad$zero_()
      
      w_1_prime$grad$zero_(); w_2_prime$grad$zero_(); w_3_prime$grad$zero_(); w_4_prime$grad$zero_()
      b_1_prime$grad$zero_(); b_2_prime$grad$zero_(); b_3_prime$grad$zero_(); b_4_prime$grad$zero_()
      
      w_5$grad$zero_();  w_6$grad$zero_();  w_7$grad$zero_()
      b_5$grad$zero_();  b_6$grad$zero_();  b_7$grad$zero_();  b_8$grad$zero_()
    })
    
    ### -------- Early stopping based on rolling ELBO means --------
    if (t >= 2 * window_elbo) {
      recent_mean <- mean(elbo_hist[(t - window_elbo + 1):t], na.rm = TRUE)
      prev_mean <- mean(elbo_hist[(t - 2 * window_elbo + 1):(t - window_elbo)], na.rm = TRUE)
      elbo_diff <- abs(recent_mean - prev_mean)
      
      if (t %% 10 == 0) {
        cat(sprintf("  rolling mean diff (latest %d vs previous %d): %.3e\n",
                    window_elbo, window_elbo, elbo_diff))
      }
      
      if (elbo_diff < delta_stop) {
        cat(sprintf("Early stopping at epoch %d: |mean_recent - mean_prev| = %.3e < %.1e\n",
                    t, elbo_diff, delta_stop))
        break
      }
    }
  }
  
  
  
  
  ###### ---------------------------------------------------------------------- ######
  ###### ----------------------------- Tail RMSE ------------------------------ ######
  ###### ---------------------------------------------------------------------- ######
  X_pred_tensor <- torch_tensor(X_holdout,dtype=torch_float())
  
  ### -------- Encoder for v_t --------
  h <- w_1_holdout$mm(X_pred_tensor)$add(b_1)$relu()
  h_1 <- w_2$mm(h)$add(b_2)$relu()
  sigma_sq_vec <- w_3$mm(h_1)$add(b_3)$exp()
  mu <- w_4$mm(h_1)$add(b_4)$relu()
  
  ### -------- Encoder for v_t_prime --------
  h_prime <- w_1_prime$mm(X_tensor)$add(b_1_prime)$relu()
  h_1_prime <- w_2_prime$mm(h_prime)$add(b_2_prime)$relu()
  
  ### -------- Activation via Laplace transformation --------
  h_1_prime_laplace <- h_1_prime$multiply(-0.2)$exp()$mean(dim=2)
  h_1_prime_t <- h_1_prime_laplace$log()$multiply(-1)
  h_1_prime_to_theta <- (0.2-h_1_prime_t$pow(2))$pow(2)$divide(4*h_1_prime_t$pow(2))$view(c(k,1))
  theta_propagate <- h_1_prime_to_theta$expand(c(k,length(fold_ind)))
  
  sigma_sq_vec_prime <- w_3_prime$mm(theta_propagate)$add(b_3_prime)$exp() #w_3_prime$mm(h_1_prime)$add(b_3_prime)$exp()
  mu_prime <- w_4_prime$mm(theta_propagate)$add(b_4_prime) #w_4_prime$mm(h_1_prime)$add(b_4_prime)
  
  
  Epsilon <- t(abs(mvtnorm::rmvnorm(length(fold_ind), mean=rep(0, k), sigma = diag(rep(1, k)))))
  Epsilon_prime <- t(mvtnorm::rmvnorm(length(fold_ind), mean=rep(0, k), sigma = diag(rep(1, k))))
  
  Epsilon <- torch_tensor(Epsilon, dtype=torch_float())
  Epsilon_prime <- torch_tensor(Epsilon_prime, dtype=torch_float())
  
  ### -------- Re-parameterization trick --------
  v_t <- mu + sqrt(sigma_sq_vec)*Epsilon
  v_t_prime <- mu_prime + sqrt(sigma_sq_vec_prime)*Epsilon_prime
  
  ### -------- Decoder --------
  l <- w_5$mm(v_t_prime)$add(b_5)$relu()
  l_1 <- w_6$mm(l)$add(b_6)$relu()
  theta_t <- w_7$mm(l_1)$add(b_7)$relu()
  
  y_star <- lrelu(W_alpha_tensor$mm(v_t)) + 1.73266e-13
  
  
  
  # ## -------------------- Convert back to original scale --------------------
  fire_approx <- y_star
  Metrix[fold] <- metric_CV(as_array(fire_approx), X_holdout, threshold)
  # fire_approx <- as_array(fire_approx)
  # save(fire_approx, file="./emulation_fire_FOLD_1.RData")
  # save(fire_approx, file="./emulation_fire_1_ncomp_7.RData")
}


rm(features)
rm(components)
fire_approx <- as_array(fire_approx)
save(fire_approx, Metrix, file=paste0("./CV_metrix_k_", k, ".RData"))
save(fire_approx, Metrix, file=paste0("./CV_metrix_k_", k, "fold1.RData"))
