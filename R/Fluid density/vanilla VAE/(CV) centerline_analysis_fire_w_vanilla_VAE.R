library(torch)


source("./utils.R")

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

load("../../data/Fluid density/center_fire_input.RData")



###### ---------------------------------------------------------------------- ######
###### ------------------------ Marginal transformation --------------------- ######
###### ---------------------------------------------------------------------- ######
hist(center_fire[1,])
hist(center_fire[49,])
range(center_fire)




# ----------------------------
# Vanilla Gaussian VAE module
# ----------------------------
vae_module <- nn_module(
  classname = "VAE",
  initialize = function(input_dim, latent_dim = 50, hidden_dim = 50,
                        enforce_positive_output = FALSE) {
    
    self$enforce_positive_output <- enforce_positive_output
    
    self$enc_fc1 <- nn_linear(input_dim, hidden_dim)
    self$enc_fc2 <- nn_linear(hidden_dim, hidden_dim)
    self$enc_mu  <- nn_linear(hidden_dim, latent_dim)
    self$enc_lv  <- nn_linear(hidden_dim, latent_dim)
    
    self$dec_fc1 <- nn_linear(latent_dim, hidden_dim)
    self$dec_fc2 <- nn_linear(hidden_dim, hidden_dim)
    self$dec_out <- nn_linear(hidden_dim, input_dim)
    
    # Common VAE trick: start with small posterior variance
    with_no_grad({
      if (!is.null(self$enc_lv$bias)) self$enc_lv$bias$fill_(-5)
    })
  },
  
  encode = function(x) {
    h  <- nnf_relu(self$enc_fc1(x))
    h1 <- nnf_relu(self$enc_fc2(h))
    mu <- self$enc_mu(h1)
    lv <- self$enc_lv(h1)  # log(sigma^2)
    list(mu = mu, logvar = lv)
  },
  
  reparameterize = function(mu, logvar) {
    # z = mu + sigma * eps, eps ~ N(0, I)
    std <- torch_exp(0.5 * logvar)
    eps <- torch_randn_like(std)
    mu + std * eps
  },
  
  decode = function(z) {
    g  <- nnf_relu(self$dec_fc1(z))
    g1 <- nnf_relu(self$dec_fc2(g))
    xhat <- self$dec_out(g1)
    
    # If your data must be strictly positive, you can enforce it:
    if (self$enforce_positive_output) {
      xhat <- nnf_softplus(xhat) + 1e-12
    }
    xhat
  },
  
  forward = function(x) {
    enc <- self$encode(x)
    z <- self$reparameterize(enc$mu, enc$logvar)
    xhat <- self$decode(z)
    list(xhat = xhat, mu = enc$mu, logvar = enc$logvar, z = z)
  }
)

# ----------------------------
# ELBO loss: recon + KL
# ----------------------------
vae_loss <- function(xhat, x, mu, logvar, beta = 1.0,
                     recon = c("mse", "gaussian_nll")) {
  recon <- match.arg(recon)
  
  if (recon == "mse") {
    # Equivalent (up to constants) to Gaussian likelihood with fixed variance
    recon_loss <- nnf_mse_loss(xhat, x, reduction = "sum")
  } else {
    # Gaussian NLL with fixed observation variance sigma_x^2 = 1
    # NLL = 0.5*||x-xhat||^2 + const; kept here for clarity
    recon_loss <- 0.5 * torch_sum((x - xhat)^2)
  }
  
  # KL(q(z|x) || p(z)), with q = N(mu, diag(exp(logvar))), p = N(0, I)
  kl <- -0.5 * torch_sum(1 + logvar - mu$pow(2) - torch_exp(logvar))
  
  # Return average per sample (more stable across batch sizes)
  batch_size <- x$size(1)
  list(
    total = (recon_loss + beta * kl) / batch_size,
    recon = recon_loss / batch_size,
    kl    = kl / batch_size
  )
}

lr_warmup_cosine <- function(step, total_steps, warmup_steps,
                             lr_max, lr_min = 0) {
  if (warmup_steps < 1) warmup_steps <- 1
  if (warmup_steps >= total_steps) warmup_steps <- total_steps - 1
  
  if (step <= warmup_steps) {
    # linear warmup from 0 -> lr_max
    return(lr_max * step / warmup_steps)
  }
  
  # cosine decay from lr_max -> lr_min
  progress <- (step - warmup_steps) / (total_steps - warmup_steps)  # in (0,1]
  lr <- lr_min + 0.5 * (lr_max - lr_min) * (1 + cos(pi * progress))
  return(lr)
}

set_optimizer_lr <- function(opt, lr_new) {
  # optim_adam stores param groups; update each group's lr
  for (g in opt$param_groups) g$lr <- lr_new
  invisible(NULL)
}


# ----------------------------
# Training wrapper
# ----------------------------
train_vae <- function(X, latent_dim = 50, hidden_dim = 50,
                      batch_size = 32, epochs = 200,
                      lr = 1e-3, lr_min = 1e-6, warmup_frac = 0.05,
                      beta = 1.0,
                      enforce_positive_output = FALSE,
                      device = c("cpu", "cuda")) {
  
  device <- match.arg(device)
  dev <- if (device == "cuda" && cuda_is_available()) torch_device("cuda") else torch_device("cpu")
  
  X_t <- torch_tensor(t(X), dtype = torch_float())$to(device = dev)
  input_dim <- X_t$size(2)
  
  ds <- tensor_dataset(X_t)
  dl <- dataloader(ds, batch_size = batch_size, shuffle = TRUE)
  
  model <- vae_module(input_dim = input_dim, latent_dim = latent_dim,
                      hidden_dim = hidden_dim,
                      enforce_positive_output = enforce_positive_output)$to(device = dev)
  
  opt <- optim_adam(model$parameters, lr = lr)
  
  # ---- scheduler bookkeeping (per minibatch step) ----
  steps_per_epoch <- ceiling(X_t$size(1) / batch_size)
  total_steps <- as.integer(epochs * steps_per_epoch)
  warmup_steps <- as.integer(max(1, floor(warmup_frac * total_steps)))
  global_step <- 0L
  
  for (ep in 1:epochs) {
    model$train()
    total_ep <- 0
    recon_ep <- 0
    kl_ep <- 0
    nb <- 0
    
    coro::loop(for (b in dl) {
      global_step <- global_step + 1L
      
      # update LR for this step
      lr_now <- lr_warmup_cosine(global_step, total_steps, warmup_steps,
                                 lr_max = lr, lr_min = lr_min)
      set_optimizer_lr(opt, lr_now)
      
      x <- b[[1]]
      opt$zero_grad()
      
      out <- model(x)
      L <- vae_loss(out$xhat, x, out$mu, out$logvar, beta = beta, recon = "mse")
      L$total$backward()
      opt$step()
      
      total_ep <- total_ep + L$total$item()
      recon_ep <- recon_ep + L$recon$item()
      kl_ep    <- kl_ep + L$kl$item()
      nb <- nb + 1
    })
    
    # report LR at end of epoch
    cat(sprintf("Epoch %4d | lr %.3e | loss %.6f | recon %.6f | KL %.6f\n",
                ep, lr_now, total_ep/nb, recon_ep/nb, kl_ep/nb))
  }
  
  invisible(model)
}

emulate_time <- function(model, X, t_idx, n_draws = 1) {
  model$eval()
  
  # one sample: [1, n_s]
  x_row <- t(X)[t_idx, , drop = FALSE]
  x_row <- matrix(as.numeric(x_row), nrow = 1)  # ensure numeric/double
  
  x_t <- torch_tensor(x_row, dtype = torch_float(), device = torch_device("cpu"))
  
  with_no_grad({
    enc <- model$encode(x_t)
    mu <- enc$mu
    logvar <- enc$logvar
    
    xhat_mean <- model$decode(mu)
    
    draws <- lapply(seq_len(n_draws), function(i) {
      z <- model$reparameterize(mu, logvar)
      model$decode(z)
    })
  })
  
  list(
    xhat_mean  = as_array(xhat_mean)[1, ],
    xhat_draws = lapply(draws, function(a) as_array(a)[1, ])
  )
}





###### ---------------------------------------------------------------------- ######
###### ---------------------------- CROSS VALIDATION ------------------------ ######
###### ---------------------------------------------------------------------- ######
metric_CV <- function(emulation, X, threshold){
  tmp_where <- (X < threshold)
  return(sum((emulation[tmp_where]-X[tmp_where])^2))
}

for(ncomp in seq(7,90, by=6)){
  threshold <- quantile(center_fire,0.05)
  Metrix <- rep(NA, 10)
  for(fold in 1:10){
    fold_ind <- seq(fold*10-9,fold*10,by=1)
    
    ###### ---------------------------------------------------------------------- ######
    ###### ------------------------ Marginal transformation --------------------- ######
    ###### ---------------------------------------------------------------------- ######
    X <- (-center_fire[, -fold_ind])*1e5 + 1
    
    X_holdout <- (-center_fire[, fold_ind])*1e5 + 1
    
    model <- train_vae(
      log(X),
      latent_dim = ncomp, hidden_dim = 200,
      batch_size = 90, epochs = 2000,
      lr = 1e-4, lr_min = 1e-7, warmup_frac = 0.5,
      beta = 1,
      enforce_positive_output = TRUE,
      device = "cpu"
    )
    
    # Example:
    X_emu <- array(NA, dim=dim(X_holdout))
    for(ind in 1:ncol(X_holdout)){
      out <- emulate_time(model, log(X_holdout), t_idx = ind, n_draws = 10)
      X_emu[,ind] <- out$xhat_mean
      # xhat_1   <- out$xhat_draws[[1]]
    }
    
    
    
    
    ###### ---------------------------------------------------------------------- ######
    ###### ----------------------------- Tail RMSE ------------------------------ ######
    ###### ---------------------------------------------------------------------- ######
    # ## -------------------- Convert back to original scale --------------------
    fire_vanilla_VAE <- -(exp(X_emu) - 1)/1e5
    Metrix[fold] <- metric_CV(fire_vanilla_VAE, center_fire[,fold_ind], threshold)
    
    # fire_approx <- as_array(fire_approx)
    # save(fire_approx, file="./emulation_fire_FOLD_1.RData")
    # save(fire_approx, file="./emulation_fire_1_ncomp_7.RData")
    if(ncomp %in% c(7,49,90) & fold==1) save(fire_vanilla_VAE, Metrix, file=paste0("./CV_metrix_k_vanilla_VAE", ncomp, "fold1.RData"))
  }
  save(fire_vanilla_VAE, Metrix, file=paste0("./CV_metrix_k_vanilla_VAE", ncomp, ".RData"))
}

