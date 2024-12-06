library(torch)

setwd("~/Desktop/Turbulence/")
source("../GEV-GP_VAE/extCVAE/utils.R")

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



###### ---------------------------------------------------------------------- ######
###### ------------------------ Marginal transformation --------------------- ######
###### ---------------------------------------------------------------------- ######
hist(center_fire[1,])
hist(center_fire[49,])
range(center_fire)

X <- (-center_fire)*1e5 + 1


###### ---------------------------------------------------------------------- ######
###### ---------------------------- Load NMF results ------------------------ ######
###### ---------------------------------------------------------------------- ######
ncomp <- 50
features <- read.csv(file = paste0("./features",ncomp,".csv"), header = FALSE, row.names = NULL)
components <- read.csv(file = paste0("./components",ncomp,".csv"), header = FALSE, row.names = NULL)



alpha = 0.5; tau <- 0.1; m <- 0.85
###### ---------------------------------------------------------------------- ######
###### ----------------------- First initial guess -------------------------- ######
###### ---------------------------------------------------------------------- ######

## -------------------- Initial guess for the latent Z variables --------------------
W_alpha <- data.matrix(features)
Z_approx <- data.matrix(components) + 1e-7
rm(features)
rm(components)

hist(log(Z_approx))
Y_star <- (W_alpha)%*%(Z_approx)


which.knot <- 1 #20
plot(Z_approx[which.knot,], type='l', main=paste('Coefficients for Feature #', which.knot, sep=''), ylab=expression(Z[t]), xlab="t")
grid()

which.knot <- 2
plot(Z_approx[which.knot,], type='l', main=paste('Coefficients for Feature #', which.knot, sep=''), ylab=expression(Z[t]), xlab="t")
grid()

which.knot <- 3
plot(Z_approx[which.knot,], type='l', main=paste('Coefficients for Feature #', which.knot, sep=''), ylab=expression(Z[t]), xlab="t")
grid()

library(ggplot2)
ind <-1
ggplot(stations) + geom_raster(aes(x=x, y=z, fill=X[,ind])) +
  scale_fill_gradientn(colours = topo.colors(100), name = paste0('Original replicate #', ind), na.value = "transparent")

ggplot(stations) + geom_raster(aes(x=x, y=z, fill=(Y_star[,ind]))) +
  # geom_point(data=knots, aes(x=x,y=z),shape='+',size=4) +
  scale_fill_gradientn(colours = topo.colors(100), name = paste0('Original replicate #', ind), na.value = "transparent")


# ## -------------------- Convert back to original scale --------------------
# fire_approx <- -Y_star/1e5
# 
# library(ggplot2)
# library(ggh4x)
# pal <- RColorBrewer::brewer.pal(11,"Spectral")
# for(tt in 1:ncol(center_fire)){
#   tmp_dat <- data.frame(x=stations$x, z=stations$z, fire = center_fire[,tt], model="Data observation")
#   tmp_dat <- rbind(tmp_dat, data.frame(x=stations$x, z=stations$z, fire = fire_approx[,tt], model="Emulation"))
#   
#   if(abs(min(center_fire[,tt]))< 1e-3) {
#     brks <- c(-2e-04, -1e-04); labels = c('-2e-04', '-1e-04')}else{
#       brks <- seq(-4e-03, 0, length.out = 5); labels = c('-4e-03', '-3e-03', '-2e-03', '-1e-03',"0")
#     }
#   ggplot(tmp_dat) + geom_raster(aes(x = x, y = z, fill = fire)) +
#     scale_fill_gradientn(colours = pal, name=paste("Time", tt), breaks = brks, labels = labels, limits = c(min(c(center_fire[,tt], fire_approx[,tt])),-min(center_fire[,tt])/40)) +
#     scale_x_continuous(breaks=seq(0, 200,length.out = 5),labels=seq(0, 8000,length.out = 5)) +
#     scale_y_continuous(breaks=seq(0, 500,length.out = 6),labels=seq(0, 5000,length.out = 6)) + labs(x='x (m)', y='z (m)') +
#     force_panelsizes(rows = unit(4.5, "in"),
#                      cols = unit(4, "in")) +
#     facet_grid(cols = vars(model))
#   
#   ggsave(paste0("./Figures/fire",tt+100,".png"), width = 9, height = 4)
# }
# 
# 
# library(magick)
# library(magrittr)
# list.files(path='./Figures/', pattern = '*.png', full.names = TRUE) %>%
#   image_read() %>% # reads each path file
#   image_join() %>% # joins image
#   image_animate(fps=2) %>% # animates, can opt for number of loops
#   image_write("./fire.gif") # write to current dir


## -------------------- Initial values for the weights --------------------
# Z = w_1 %*% X or Z^T = X^T %*% w_1^T
X_tensor <- torch_tensor(X,dtype=torch_float())
W_alpha_tensor <- torch_tensor(W_alpha,dtype=torch_float())
k <- ncol(W_alpha)
tmp <- qr.solve(a=t(X), b=t(Z_approx))
w_1 <- t(tmp)
w_1 <- torch_tensor(w_1,dtype=torch_float(),requires_grad = TRUE)
b_1 <- matrix(rep(0,k), ncol=1) 
b_1 <- torch_tensor(b_1,dtype=torch_float(),requires_grad = TRUE)

w_2 <- diag(k)
w_2 <- torch_tensor(w_2,dtype=torch_float(),requires_grad = TRUE)
b_2 <- matrix(rep(0.00001,k), ncol=1)
b_2 <- torch_tensor(b_2,dtype=torch_float(),requires_grad = TRUE)

w_4 <- diag(k)
w_4 <- torch_tensor(w_4,dtype=torch_float(),requires_grad = TRUE)
b_4 <- matrix(rep(0,k), ncol=1) 
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
b_1_prime <- matrix(rep(0,k), ncol=1) #matrix(rnorm(k), ncol=1)
b_1_prime <- torch_tensor(b_1_prime,dtype=torch_float(),requires_grad = TRUE)
b_2_prime <- matrix(rep(0.05,k), ncol=1)
b_2_prime <- torch_tensor(b_2_prime,dtype=torch_float(),requires_grad = TRUE)
b_3_prime <- matrix(rep(-10,k), ncol=1)
b_3_prime <- torch_tensor(b_3_prime,dtype=torch_float(),requires_grad = TRUE)
b_4_prime <- matrix(rep(0,k), ncol=1)
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
learning_rate <- -1e-19
alpha_v <- 0.95
lrelu <- nn_leaky_relu(-0.01)

niter = 5000
n <- 1e3
# only depends on alpha; we are not updating alpha for now.
vec <- Zolo_A(pi*seq(1/2,n-1/2,1)/n, alpha)
Zolo_vec <- torch_tensor(matrix(vec, nrow=1,ncol=n), dtype=torch_float(), requires_grad = FALSE) 
Zolo_vec_double <- torch_tensor(Zolo_vec, dtype = torch_float64(), requires_grad = FALSE)
const <- 1/(1-alpha); const1 <- 1/(1-alpha)-1; const3 <- log(const1)
old_loss <- -Inf

# start_time <- Sys.time()
for (t in 1:niter) {
  if(t==1000) { learning_rate <- -1e-17; alpha_v <- 0.9}
  if(t==5000) { learning_rate <- -1e-14; alpha_v <- 0.7}
  if(t==10000) { learning_rate <- -5e-14; alpha_v <- 0.5}
  if(t==20000) { learning_rate <- -1e-16; alpha_v <- 0.5}
  if(t==80000) { learning_rate <- -1e-15; alpha_v <- 0.4}
  if(t==100000) { learning_rate <- -5e-15; alpha_v <- 0.4}
  if(t==120000) { learning_rate <- -1e-14; alpha_v <- 0.4}
  if(t==140000) { learning_rate <- -1e-13; alpha_v <- 0.4}
  if(t==300000) { learning_rate <- -1e-12; alpha_v <- 0.4}
  
  ### -------- Forward pass ---------
  Epsilon <- matrix(abs(rnorm(k*n.t))+0.1, nrow=k)
  Epsilon_prime <- t(mvtnorm::rmvnorm(n.t, mean=rep(0, k), sigma = diag(rep(1, k))))
  
  Epsilon <- torch_tensor(Epsilon,dtype=torch_float())
  Epsilon_prime <- torch_tensor(Epsilon_prime,dtype=torch_float())
  
  ### -------- Encoder for v_t --------
  h <- w_1$mm(X_tensor)$add(b_1)$relu()
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
  theta_propagate <- h_1_prime_to_theta$expand(c(k,n.t))
  
  sigma_sq_vec_prime <- w_3_prime$mm(theta_propagate)$add(b_3_prime)$exp() #w_3_prime$mm(h_1_prime)$add(b_3_prime)$exp()
  mu_prime <- w_4_prime$mm(theta_propagate)$add(b_4_prime) #w_4_prime$mm(h_1_prime)$add(b_4_prime)
  
  ### -------- Re-parameterization trick --------
  v_t <- mu + sqrt(sigma_sq_vec)*Epsilon
  v_t_prime <- mu_prime + sqrt(sigma_sq_vec_prime)*Epsilon_prime
  
  
  ### -------- Decoder --------
  l <- w_5$mm(v_t_prime)$add(b_5)$relu()
  l_1 <- w_6$mm(l)$add(b_6)$relu()
  theta_t <- w_7$mm(l_1)$add(b_7)$relu()
  
  y_star <- lrelu(W_alpha_tensor$mm(v_t)$add(b_8)) + 1.73266e-13
  
  
  ### -------- ELBO -------- 
  ### Part 1
  standardized <- X_tensor$divide(y_star)$sub(m)
  leak <- as_array(sum(standardized<0))
  if(leak>0 & leak<=1e4) standardized$abs_()
  part1 <- -3 * standardized$log()$sum() - y_star$log()$sum() - tau*standardized$pow(-2)$sum() # + n.s*n.t*log(tau)
  
  ### Part 2
  V_t <- v_t$view(c(k*n.t,1))
  Theta_t <- theta_t$view(c(k*n.t,1))
  part_log_v1  <- V_t$pow(-const)$mm(Zolo_vec) 
  part_log_v2  <- (-V_t$pow(-const1)$mm(Zolo_vec))$exp()
  part_log_v3 <- Theta_t$pow(alpha)-Theta_t$mul(V_t)
  part2 <- (part_log_v1$mul(part_log_v2)$mean(dim=2)$log()+part_log_v3$view(k*n.t)$add(const3))$sum()
  if(as_array(part2) == -Inf) {
    part2_tmp <- part_log_v1$log()[,1] + (-V_t$pow(-const1)$mm(Zolo_vec))[,1]
    part2 <- (part2_tmp+part_log_v3$view(k*n.t)$add(const3))$sum()
  }
  
  ### Part 3
  part3 = Epsilon$pow(2)$sum()/2 + Epsilon_prime$pow(2)$sum()/2 + sigma_sq_vec$log()$sum() + sigma_sq_vec_prime$log()$sum()
  res <- part1 + part2 + part3
  
  ### -------- compute loss -------- 
  loss <- (res/(n.s*n.t))
  if(!is.finite(loss$item())) break
  if (t %% 1 == 0)
    cat("Epoch: ", t, "   ELBO: ", loss$item(), "\n") # we want to maximize
  if (as.numeric(torch_isnan(loss))==1) break
  ### -------- Backpropagation --------
  # compute gradient of loss w.r.t. all tensors with requires_grad = TRUE
  loss$backward()
  # w_1$grad$argmax(dim=2)
  
  ### -------- Update weights -------- 
  
  # Wrap in with_no_grad() because this is a part we DON'T 
  # want to record for automatic gradient computation
  with_no_grad({
    old_loss <- loss$item()
    w_1_velocity <- alpha_v*w_1_velocity - learning_rate*w_1$grad
    w_1$add_(w_1_velocity)
    w_2_velocity <- alpha_v*w_2_velocity - learning_rate*w_2$grad
    w_2$add_(w_2_velocity)
    w_3_velocity <- alpha_v*w_3_velocity - learning_rate*w_3$grad
    w_3$add_(w_3_velocity)
    w_4_velocity <- alpha_v*w_4_velocity - learning_rate*w_4$grad
    w_4$add_(w_4_velocity)
    b_1_velocity <- alpha_v*b_1_velocity - learning_rate*b_1$grad
    b_1$add_(b_1_velocity)
    b_2_velocity <- alpha_v*b_2_velocity - learning_rate*b_2$grad
    b_2$add_(b_2_velocity)
    b_3_velocity <- alpha_v*b_3_velocity - learning_rate*b_3$grad
    b_3$add_(b_3_velocity)
    b_4_velocity <- alpha_v*b_4_velocity - learning_rate*b_4$grad
    b_4$add_(b_4_velocity)
    w_1_prime_velocity <- alpha_v*w_1_prime_velocity - learning_rate*w_1_prime$grad
    w_1_prime$add_(w_1_prime_velocity)
    w_2_prime_velocity <- alpha_v*w_2_prime_velocity - learning_rate*w_2_prime$grad
    w_2_prime$add_(w_2_prime_velocity)
    w_3_prime_velocity <- alpha_v*w_3_prime_velocity - learning_rate*w_3_prime$grad
    w_3_prime$add_(w_3_prime_velocity)
    w_4_prime_velocity <- alpha_v*w_4_prime_velocity - learning_rate*w_4_prime$grad
    w_4_prime$add_(w_4_prime_velocity)
    b_1_prime_velocity <- alpha_v*b_1_prime_velocity - learning_rate*b_1_prime$grad
    b_1_prime$add_(b_1_prime_velocity)
    b_2_prime_velocity <- alpha_v*b_2_prime_velocity - learning_rate*b_2_prime$grad
    b_2_prime$add_(b_2_prime_velocity)
    b_3_prime_velocity <- alpha_v*b_3_prime_velocity - learning_rate*b_3_prime$grad
    b_3_prime$add_(b_3_prime_velocity)
    b_4_prime_velocity <- alpha_v*b_4_prime_velocity - learning_rate*b_4_prime$grad
    b_4_prime$add_(b_4_prime_velocity)
    
    w_5_velocity <- alpha_v*w_5_velocity - learning_rate*w_5$grad
    w_5$add_(w_5_velocity)
    w_6_velocity <- alpha_v*w_6_velocity - learning_rate*w_6$grad
    w_6$add_(w_6_velocity)
    w_7_velocity <- alpha_v*w_7_velocity - learning_rate*w_7$grad
    w_7$add_(w_7_velocity)
    b_5_velocity <- alpha_v*b_5_velocity - learning_rate*b_5$grad
    b_5$add_(b_5_velocity)
    b_6_velocity <- alpha_v*b_6_velocity - learning_rate*b_6$grad
    b_6$add_(b_6_velocity)
    b_7_velocity <- alpha_v*b_7_velocity - learning_rate*b_7$grad
    b_7$add_(b_7_velocity)
    b_8_velocity <- alpha_v*b_8_velocity - learning_rate*b_8$grad
    b_8$add_(b_8_velocity)
    
    # Zero gradients after every pass, as they'd accumulate otherwise
    w_1$grad$zero_()
    w_2$grad$zero_()
    w_3$grad$zero_()
    w_4$grad$zero_()
    b_1$grad$zero_()
    b_2$grad$zero_()
    b_3$grad$zero_()
    b_4$grad$zero_()
    w_1_prime$grad$zero_()
    w_2_prime$grad$zero_()
    w_3_prime$grad$zero_()
    w_4_prime$grad$zero_()
    b_1_prime$grad$zero_()
    b_2_prime$grad$zero_()
    b_3_prime$grad$zero_()
    b_4_prime$grad$zero_()
    w_5$grad$zero_()
    w_6$grad$zero_()
    w_7$grad$zero_()
    b_5$grad$zero_()
    b_6$grad$zero_()
    b_7$grad$zero_()
    b_8$grad$zero_()
  })
  
}
# end_time <- Sys.time()
# end_time - start_time

# loss - part3/(n.s*n.t)
# plot(theta_sim[,1], as_array(theta_t[,3]), xlab=expression(paste('true ',theta)), ylab=expression(paste('CVAE ',theta)))
# abline(a = 0, b = 1, col='orange', lty=2)

# Z_approx <- as_array(v_t)
# save(Z_approx, file="./Z_approx_NMF_50.RData")


###### ---------------------------------------------------------------------- ######
###### ------------------------ Animation + Emulation ----------------------- ######
###### ---------------------------------------------------------------------- ######

### -------- Encoder for v_t --------
h <- w_1$mm(X_tensor)$add(b_1)$relu()
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
theta_propagate <- h_1_prime_to_theta$expand(c(k,n.t))

sigma_sq_vec_prime <- w_3_prime$mm(theta_propagate)$add(b_3_prime)$exp() #w_3_prime$mm(h_1_prime)$add(b_3_prime)$exp()
mu_prime <- w_4_prime$mm(theta_propagate)$add(b_4_prime) #w_4_prime$mm(h_1_prime)$add(b_4_prime)


Epsilon <- t(abs(mvtnorm::rmvnorm(n.t, mean=rep(0, k), sigma = diag(rep(1, k)))))
Epsilon_prime <- t(mvtnorm::rmvnorm(n.t, mean=rep(0, k), sigma = diag(rep(1, k))))

Epsilon <- torch_tensor(Epsilon, dtype=torch_float())
Epsilon_prime <- torch_tensor(Epsilon_prime, dtype=torch_float())

### -------- Re-parameterization trick --------
v_t <- mu + sqrt(sigma_sq_vec)*Epsilon
v_t_prime <- mu_prime + sqrt(sigma_sq_vec_prime)*Epsilon_prime

### -------- Decoder --------
l <- w_5$mm(v_t_prime)$add(b_5)$relu()
l_1 <- w_6$mm(l)$add(b_6)$relu()
theta_t <- w_7$mm(l_1)$add(b_7)$relu()

y_star <- lrelu(W_alpha_tensor$mm(v_t)$add(b_8)) + 1.73266e-13

## Decoder
station_Simulations_All <- matrix(VGAM::rfrechet(n.s*n.t,shape=2, location = m, scale = tau), nrow=n.s) * as_array((y_star))

## For tail RMSE
fire_approx <- -(y_star-1)/1e5
threshold <- quantile(center_fire,0.05)
metric_CV <- function(emulation, X, threshold){
  tmp_where <- (X < threshold)
  return(sum((emulation[tmp_where]-X[tmp_where])^2))
}
fire_approx <- as_array(fire_approx)
save(fire_approx, file="./emulation_fire_1.RData")
sqrt(metric_CV(fire_approx, center_fire, threshold)/(n.s*n.t))



# ## -------------------- Convert back to original scale --------------------



library(ggplot2)
library(ggh4x)
pal <- RColorBrewer::brewer.pal(11,"Spectral")
for(tt in 1:ncol(center_fire)){
  tmp_dat <- data.frame(x=stations$x, z=stations$z, fire = center_fire[,tt], model="Data observation")
  tmp_dat <- rbind(tmp_dat, data.frame(x=stations$x, z=stations$z, fire = fire_approx[,tt], model="Emulation"))
  
  if(abs(min(center_fire[,tt]))< 1e-3) {
    brks <- c(-2e-04, -1e-04); labels = c('-2e-04', '-1e-04')}else{
      brks <- seq(-4e-03, 0, length.out = 5); labels = c('-4e-03', '-3e-03', '-2e-03', '-1e-03',"0")
    }
  ggplot(tmp_dat) + geom_raster(aes(x = x, y = z, fill = fire)) +
    scale_fill_gradientn(colours = pal, name=paste("Time", tt), breaks = brks, labels = labels, limits = c(min(c(center_fire[,tt], fire_approx[,tt])),-min(center_fire[,tt])/40)) +
    scale_x_continuous(breaks=seq(0, 200,length.out = 5),labels=seq(0, 8000,length.out = 5)) +
    scale_y_continuous(breaks=seq(0, 500,length.out = 6),labels=seq(0, 5000,length.out = 6)) + labs(x='x (m)', y='z (m)') +
    force_panelsizes(rows = unit(4.5, "in"),
                     cols = unit(4, "in")) +
    facet_grid(cols = vars(model))
  
  ggsave(paste0("./Figures/fire",tt+100,".png"), width = 9, height = 4)
}


library(magick)
library(magrittr)
list.files(path='./Figures/', pattern = '*.png', full.names = TRUE) %>%
  image_read() %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=2) %>% # animates, can opt for number of loops
  image_write("./fire.gif") # write to current dir








###### ---------------------------------------------------------------------- ######
###### ----------------------------   QQ plots ------------------------------ ######
###### ---------------------------------------------------------------------- ######

###### ---------  Run decoder to get the simulated weather processes -------- ######
n.sim<-1000
station1_Simulations <- matrix(NA, nrow=n.s, ncol=n.sim)
station50_Simulations <- matrix(NA, nrow=n.s, ncol=n.sim)
station55_Simulations <- matrix(NA, nrow=n.s, ncol=n.sim)
station90_Simulations <- matrix(NA, nrow=n.s, ncol=n.sim)
station100_Simulations <- matrix(NA, nrow=n.s, ncol=n.sim)

### -------- Encoder for v_t --------
h <- w_1$mm(X_tensor)$add(b_1)$relu()
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
theta_propagate <- h_1_prime_to_theta$expand(c(k,n.t))

sigma_sq_vec_prime <- w_3_prime$mm(theta_propagate)$add(b_3_prime)$exp() #w_3_prime$mm(h_1_prime)$add(b_3_prime)$exp()
mu_prime <- w_4_prime$mm(theta_propagate)$add(b_4_prime) #w_4_prime$mm(h_1_prime)$add(b_4_prime)



for(iter in 1:n.sim){
  if(iter %% 100==0) cat('iter=', iter, '\n')
  Epsilon <- t(abs(mvtnorm::rmvnorm(n.t, mean=rep(0, k), sigma = diag(rep(10000, k)))))
  Epsilon_prime <- t(mvtnorm::rmvnorm(n.t, mean=rep(0, k), sigma = diag(rep(10000, k))))
  
  Epsilon <- torch_tensor(Epsilon, dtype=torch_float())
  Epsilon_prime <- torch_tensor(Epsilon_prime, dtype=torch_float())
  
  ### -------- Re-parameterization trick --------
  v_t <- mu + sqrt(sigma_sq_vec)*Epsilon
  v_t_prime <- mu_prime + sqrt(sigma_sq_vec_prime)*Epsilon_prime
  
  ### -------- Decoder --------
  l <- w_5$mm(v_t_prime)$add(b_5)$relu()
  l_1 <- w_6$mm(l)$add(b_6)$relu()
  theta_t <- w_7$mm(l_1)$add(b_7)$relu()
  
  y_star <- lrelu(W_alpha_tensor$mm(v_t)$add(b_8)) + 1.73266e-13
  
  ##Decoder
  station1_Simulations[, iter] <- VGAM::rfrechet(n.s,shape=2, location = m, scale = tau) * as_array((y_star[,1]))
  station50_Simulations[, iter] <- VGAM::rfrechet(n.s,shape=2, location = m, scale = tau) * as_array((y_star[,50]))
  station55_Simulations[, iter] <- VGAM::rfrechet(n.s,shape=2, location = m, scale = tau) * as_array((y_star[,55]))
  station90_Simulations[, iter] <- VGAM::rfrechet(n.s,shape=2, location = m, scale = tau) * as_array((y_star[,90]))
  station100_Simulations[, iter] <- VGAM::rfrechet(n.s,shape=2, location = m, scale = tau) * as_array((y_star[,100]))
}

ind <- 1
ggplot(stations) + geom_raster(aes(x=x, y=z, fill=X[,ind])) +
  scale_fill_gradientn(colours = topo.colors(100), name = paste0('Original replicate #', ind), na.value = "transparent")

ggplot(stations) + geom_raster(aes(x=x, y=z, fill= as_array(y_star[,1]))) +
  scale_fill_gradientn(colours = topo.colors(100), name = paste0('Emulated replicate #', ind), na.value = "transparent")



## -------- time 1 --------
ind <- 1
png(filename = "./Figures/qqplot_1.png", width = 360, height = 380)
extRemes::qqplot(center_fire[,ind], -(station1_Simulations[,floor(n.sim/2)]-1)/1e5, regress = FALSE,  
                 xlab=expression(paste('Observed ', X[1])), 
                 ylab=expression(paste('Emulated ', X[1])),
                 main=paste0("Fold # = ", 1, ", XVAE optimum setting"))
dev.off()


## -------- time 55 --------
ind <- 55
png(filename = "./Figures/qqplot_55.png", width = 360, height = 380)
extRemes::qqplot(X[,ind], station55_Simulations[,floor(n.sim/2)], regress = FALSE,  
                 xlab=expression(paste('Observed ', X[55])), 
                 ylab=expression(paste('Emulated ', X[55])))
dev.off()

## -------- time 90 --------
ind <- 90
png(filename = "./Figures/qqplot_90.png", width = 360, height = 380)
extRemes::qqplot(center_fire[,ind], -(station90_Simulations[,floor(n.sim/2)]-1)/1e5, regress = FALSE,  
                 xlab=expression(paste('Observed ', X[90])), ylim=c(-4e-04,0),
                 ylab=expression(paste('Emulated ', X[90])),
                 main=paste0("Fold # = ", 10, ", XVAE optimum setting"))
dev.off()




library(ggh4x)
pal <- rev(RColorBrewer::brewer.pal(11, "PiYG"))
## ----------------------------------
## ------------- Time 1 -------------
## ----------------------------------
plot_dat <- data.frame(logX=c(log(station1_Simulations[500,]),log(station1_Simulations[3861,]),log(station1_Simulations[17323,]),log(station1_Simulations[80017,]),log(station1_Simulations[83014,])),
                        ind =factor(rep(c("X[500,1]", "X[3861,1]", "X[17323,1]", "X[80017,1]", "X[83014,1]"), each=1000), levels = c("X[500,1]", "X[3861,1]", "X[17323,1]", "X[80017,1]", "X[83014,1]")))
ggplot(plot_dat, aes(x=logX, fill=ind, color=ind))  +
  geom_histogram(aes(y=ifelse(after_stat(density) > 0, after_stat(density), NA)), position="identity", alpha=0.5, binwidth=0.05)  + theme_bw() +
  xlab("log(X)") + ylab("Density") + guides(fill=guide_legend(title="Time 1"), color=guide_legend(title="Time 1"))+
  scale_fill_manual(values = c('#030fab', '#9c7605', '#fa860a', '#03fcec', '#b3079f')) +
  scale_color_manual(values = c('#030fab', '#9c7605', '#fa860a', '#03fcec', '#b3079f')) +
  scale_x_continuous(limits=c(0,7.4), expand = c(0, 0)) +
  scale_y_continuous(limits=c(0,6.1), expand = c(0, 0)) +
  theme(legend.title = element_text(face='bold')) + 
  force_panelsizes(cols = unit(5, "in"),
                   rows = unit(2.4, "in"))
ggsave(filename="./Figures/histogram_time1.pdf", width=7, height=3)

pts <- data.frame(x = stations$x[c(500,3861,17323,80017,83014)], z = stations$z[c(500,3861,17323,80017,83014)], ind = factor(1:5))
tmp_dat <- data.frame(x=stations$x, z=stations$z, fire = log(X[,1]), model="Data observation")
ggplot(tmp_dat) + geom_raster(aes(x = x, y = z, fill = fire), alpha=0.9) +
  geom_point(data=pts, aes(x=x,y=z, col=ind), size=3, shape=17) +
  scale_fill_gradientn(colours = pal, name=paste("Time", 1)) +
  scale_color_manual(values = c('#030fab', '#9c7605', '#fa860a', '#03fcec', '#b3079f')) +
  scale_x_continuous(expand = c(0, 0), breaks=seq(0, 200,length.out = 5),labels=seq(0, 8000,length.out = 5)) +
  scale_y_continuous(expand = c(0, 0), breaks=seq(0, 500,length.out = 6),labels=seq(0, 5000,length.out = 6)) + labs(x='x (m)', y='z (m)') +
  theme(legend.title = element_text(face='bold'), panel.grid=element_blank()) + guides(color="none") +
  force_panelsizes(rows = unit(2.4, "in"),
                   cols = unit(2.44, "in"))
ggsave(filename="./Figures/X_time1.pdf", width=3.9, height=3)



## -----------------------------------
## ------------- Time 50 -------------
## -----------------------------------
plot_dat <- data.frame(logX=c(log(station50_Simulations[500,]),log(station50_Simulations[3861,]),log(station50_Simulations[17323,]),log(station50_Simulations[80017,]),log(station50_Simulations[83014,])),
                       ind =factor(rep(c("X[500,50]", "X[3861,50]", "X[17323,50]", "X[80017,50]", "X[83014,50]"), each=1000), levels = c("X[500,50]", "X[3861,50]", "X[17323,50]", "X[80017,50]", "X[83014,50]")))
ggplot(plot_dat, aes(x=logX, fill=ind, color=ind))  +
  geom_histogram(aes(y=ifelse(after_stat(density) > 0, after_stat(density), NA)), position="identity", alpha=0.5, binwidth=0.05) + theme_bw() +
  xlab("log(X)") + ylab("Density") + guides(fill=guide_legend(title="Time 50"), color=guide_legend(title="Time 50"))+
  scale_fill_manual(values = c('#030fab', '#9c7605', '#fa860a', '#03fcec', '#b3079f')) +
  scale_color_manual(values = c('#030fab', '#9c7605', '#fa860a', '#03fcec', '#b3079f')) +
  scale_x_continuous(limits=c(0,7.4), expand = c(0, 0)) +
  scale_y_continuous(limits=c(0,5.3), expand = c(0, 0)) +
  theme(legend.title = element_text(face='bold')) + 
  force_panelsizes(cols = unit(5, "in"),
                   rows = unit(2.4, "in"))
ggsave(filename="./Figures/histogram_time50.pdf", width=7, height=3)

pts <- data.frame(x = stations$x[c(500,3861,17323,80017,83014)], z = stations$z[c(500,3861,17323,80017,83014)], ind = factor(1:5))
tmp_dat <- data.frame(x=stations$x, z=stations$z, fire = log(X[,50]), model="Data observation")
ggplot(tmp_dat) + geom_raster(aes(x = x, y = z, fill = fire), alpha=0.9) +
  geom_point(data=pts, aes(x=x,y=z, col=ind), size=3, shape=17) +
  scale_fill_gradientn(colours = pal, name=paste("Time", 50)) +
  scale_color_manual(values = c('#030fab', '#9c7605', '#fa860a', '#03fcec', '#b3079f')) +
  scale_x_continuous(expand = c(0, 0), breaks=seq(0, 200,length.out = 5),labels=seq(0, 8000,length.out = 5)) +
  scale_y_continuous(expand = c(0, 0), breaks=seq(0, 500,length.out = 6),labels=seq(0, 5000,length.out = 6)) + labs(x='x (m)', y='z (m)') +
  theme(legend.title = element_text(face='bold'), panel.grid=element_blank()) + guides(color="none") +
  force_panelsizes(rows = unit(2.4, "in"),
                   cols = unit(2.44, "in"))
ggsave(filename="./Figures/X_time50.pdf", width=3.9, height=3)




## ------------------------------------
## ------------- Time 100 -------------
## ------------------------------------
plot_dat <- data.frame(logX=c(log(station100_Simulations[500,]),log(station100_Simulations[3861,]),log(station100_Simulations[17323,]),log(station100_Simulations[80017,]),log(station100_Simulations[83014,])),
                       ind =factor(rep(c("X[500,100]", "X[3861,100]", "X[17323,100]", "X[80017,100]", "X[83014,100]"), each=1000), levels = c("X[500,100]", "X[3861,100]", "X[17323,100]", "X[80017,100]", "X[83014,100]")))
ggplot(plot_dat, aes(x=logX, fill=ind, color=ind))  +
  geom_histogram(aes(y=ifelse(after_stat(density) > 0, after_stat(density), NA)), position="identity", alpha=0.5, binwidth=0.05) + theme_bw() +
  xlab("log(X)") + ylab("Density") + guides(fill=guide_legend(title="Time 100"), color=guide_legend(title="Time 100"))+
  scale_fill_manual(values = c('#030fab', '#9c7605', '#fa860a', '#03fcec', '#b3079f')) +
  scale_color_manual(values = c('#030fab', '#9c7605', '#fa860a', '#03fcec', '#b3079f')) +
  scale_x_continuous(limits=c(0,7.4), expand = c(0, 0)) +
  scale_y_continuous(limits=c(0,6.5), expand = c(0, 0)) +
  theme(legend.title = element_text(face='bold')) + 
  force_panelsizes(cols = unit(5, "in"),
                   rows = unit(2.4, "in"))
ggsave(filename="./Figures/histogram_time100.pdf", width=7, height=3)

pts <- data.frame(x = stations$x[c(500,3861,17323,80017,83014)], z = stations$z[c(500,3861,17323,80017,83014)], ind = factor(1:5))
tmp_dat <- data.frame(x=stations$x, z=stations$z, fire = log(X[,100]), model="Data observation")
ggplot(tmp_dat) + geom_raster(aes(x = x, y = z, fill = fire), alpha=0.9) +
  geom_point(data=pts, aes(x=x,y=z, col=ind), size=3, shape=17) +
  scale_fill_gradientn(colours = pal, name=paste("Time", 100)) +
  scale_color_manual(values = c('#030fab', '#9c7605', '#fa860a', '#03fcec', '#b3079f')) +
  scale_x_continuous(expand = c(0, 0), breaks=seq(0, 200,length.out = 5),labels=seq(0, 8000,length.out = 5)) +
  scale_y_continuous(expand = c(0, 0), breaks=seq(0, 500,length.out = 6),labels=seq(0, 5000,length.out = 6)) + labs(x='x (m)', y='z (m)') +
  theme(legend.title = element_text(face='bold'), panel.grid=element_blank()) + guides(color="none") +
  force_panelsizes(rows = unit(2.4, "in"),
                   cols = unit(2.44, "in"))
ggsave(filename="./Figures/X_time100.pdf", width=3.9, height=3)





## --------------------------------------------------------------------------------------------
## ----------------------------------- Joint Density ------------------------------------------
## --------------------------------------------------------------------------------------------

min_tmp1 <- range(station1_Simulations[500,])[1]
min_tmp2 <- range(station1_Simulations[3861,])[1]
f1 <- MASS::kde2d(x=log(station1_Simulations[500,]), y=log(station1_Simulations[3861,]), h=0.2, n = 200, 
                  lims = c(log(min_tmp1)-0.1, log(min_tmp1)+0.4, log(min_tmp2)-0.1, log(min_tmp2)+0.5))
min_tmp1 <- range(station50_Simulations[500,])[1]
min_tmp2 <- range(station50_Simulations[3861,])[1]
f2 <- MASS::kde2d(log(station50_Simulations[500,]), log(station50_Simulations[3861,]), h=0.3, n = 200,
                  lims = c(log(min_tmp1)-0.1, log(min_tmp1)+1.4, log(min_tmp2)-0.1, log(min_tmp2)+1.2))
min_tmp1 <- range(station100_Simulations[500,])[1]
min_tmp2 <- range(station100_Simulations[3861,])[1]
f3 <- MASS::kde2d(log(station100_Simulations[500,]), log(station100_Simulations[3861,]), h=0.5, n = 200, 
                  lims = c(log(min_tmp1)-0.1, log(min_tmp1)+2, log(min_tmp2)-0.1, log(min_tmp2)+2))
upper <- max(f1$z, f2$z, f3$z)
lower_x <- min(f1$x, f2$x, f3$x); upper_x <- max(f1$x, f2$x, f3$x)
lower_y <- min(f1$y, f2$y, f3$y); upper_y <- max(f1$y, f2$y, f3$y) 


pal <- hcl.colors(10, "Spectral")[10:1]

plot_dat1 <- expand.grid(x=f1$x, y=f1$y)
plot_dat1$den <- as.vector(f1$z)
plot1 <- ggplot(plot_dat1) + geom_raster(aes(x=x, y=y, fill=den)) +
  geom_contour(col="black", aes(x=x, y=y, z=den)) +
  labs(x= expression(log(X[500])), y = expression(log(X[3861]))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + ggtitle("Time 1")+
  scale_fill_gradientn(colours = pal, name = "Density", na.value = NA, limits = c(0, upper)) +
  theme(plot.title = element_text(hjust = 0.5, face='bold'))


plot_dat2 <- expand.grid(x=f2$x, y=f2$y)
plot_dat2$den <- as.vector(f2$z)
plot2 <- ggplot(plot_dat2) + geom_raster(aes(x=x, y=y, fill=den)) +
  geom_contour(col="black", aes(x=x, y=y, z=den)) +
  labs(x= expression(log(X[500])), y = expression(log(X[3861]))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + ggtitle("Time 50")+
  scale_fill_gradientn(colours = pal, name = "Density", na.value = NA, limits = c(0, upper)) +
  theme(plot.title = element_text(hjust = 0.5, face='bold'))

plot_dat3 <- expand.grid(x=f3$x, y=f3$y)
plot_dat3$den <- as.vector(f3$z)
plot3 <- ggplot(plot_dat3) + geom_raster(aes(x=x, y=y, fill=den)) +
  geom_contour(col="black", aes(x=x, y=y, z=den)) +
  labs(x= expression(log(X[500])), y = expression(log(X[3861]))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + ggtitle("Time 100")+
  scale_fill_gradientn(colours = pal, name = "Density", na.value = NA, limits = c(0, upper)) +
  theme(plot.title = element_text(hjust = 0.5, face='bold'))

library(cowplot)
prow <- plot_grid(plot1+ theme(legend.position="none"), plot2+ theme(legend.position="none"), plot3+ theme(legend.position="none"),
                  align = "h", nrow = 1)
legend <- get_legend(
  # create some space to the left of the legend
  plot3 + theme(legend.box.margin = margin(0, 0, 0, 0))
)
plot_grid(prow, legend, rel_widths = c(3, .4))
ggsave(file = paste0("./Figures/bivariate_pdf.pdf"), width = 8, height =2.5)






min_tmp1 <- range(station1_Simulations[500,])[1]
min_tmp2 <- range(station1_Simulations[17323,])[1]
f1 <- MASS::kde2d(log(station1_Simulations[500,]), log(station1_Simulations[17323,]), h=0.4, n = 200, 
                  lims = c(log(min_tmp1)-0.2, log(min_tmp1)+0.5, log(min_tmp2)-0.1, log(min_tmp2)+1.8))
min_tmp1 <- range(station50_Simulations[500,])[1]
min_tmp2 <- range(station50_Simulations[17323,])[1]
f2 <- MASS::kde2d(log(station50_Simulations[500,]), log(station50_Simulations[17323,]), h=0.5, n = 200, 
                  lims = c(log(min_tmp1)-0.1, log(min_tmp1)+1.4, log(min_tmp2)-0.2, log(min_tmp2)+1))
min_tmp1 <- range(station100_Simulations[500,])[1]
min_tmp2 <- range(station100_Simulations[17323,])[1]
f3 <- MASS::kde2d(log(station100_Simulations[500,]), log(station100_Simulations[17323,]), h=0.5, n = 200, 
                  lims = c(log(min_tmp1)-0.1, log(min_tmp1)+2, log(min_tmp2)-0.1, log(min_tmp2)+0.9))
upper <- max(f1$z, f2$z, f3$z)
lower_x <- min(f1$x, f2$x, f3$x); upper_x <- max(f1$x, f2$x, f3$x)
lower_y <- min(f1$y, f2$y, f3$y); upper_y <- max(f1$y, f2$y, f3$y) 


pal <- hcl.colors(10, "Spectral")[10:1]

plot_dat1 <- expand.grid(x=f1$x, y=f1$y)
plot_dat1$den <- as.vector(f1$z)
plot1 <- ggplot(plot_dat1) + geom_raster(aes(x=x, y=y, fill=den)) +
  geom_contour(col="black", aes(x=x, y=y, z=den)) +
  labs(x= expression(log(X[500])), y = expression(log(X[17323]))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + ggtitle("Time 1")+
  scale_fill_gradientn(colours = pal, name = "Density", na.value = NA, limits = c(0, upper)) +
  theme(plot.title = element_text(hjust = 0.5, face='bold'))


plot_dat2 <- expand.grid(x=f2$x, y=f2$y)
plot_dat2$den <- as.vector(f2$z)
plot2 <- ggplot(plot_dat2) + geom_raster(aes(x=x, y=y, fill=den)) +
  geom_contour(col="black", aes(x=x, y=y, z=den)) +
  labs(x= expression(log(X[500])), y = expression(log(X[17323]))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + ggtitle("Time 50")+
  scale_fill_gradientn(colours = pal, name = "Density", na.value = NA, limits = c(0, upper)) +
  theme(plot.title = element_text(hjust = 0.5, face='bold'))

plot_dat3 <- expand.grid(x=f3$x, y=f3$y)
plot_dat3$den <- as.vector(f3$z)
plot3 <- ggplot(plot_dat3) + geom_raster(aes(x=x, y=y, fill=den)) +
  geom_contour(col="black", aes(x=x, y=y, z=den)) +
  labs(x= expression(log(X[500])), y = expression(log(X[17323]))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + ggtitle("Time 100")+
  scale_fill_gradientn(colours = pal, name = "Density", na.value = NA, limits = c(0, upper)) +
  theme(plot.title = element_text(hjust = 0.5, face='bold'))

library(cowplot)
prow <- plot_grid(plot1+ theme(legend.position="none"), plot2+ theme(legend.position="none"), plot3+ theme(legend.position="none"),
                  align = "h", nrow = 1)
legend <- get_legend(
  # create some space to the left of the legend
  plot3 + theme(legend.box.margin = margin(0, 0, 0, 0))
)
plot_grid(prow, legend, rel_widths = c(3, .4))
ggsave(file = paste0("./Figures/bivariate_pdf2.pdf"), width = 8, height =2.5)




min_tmp1 <- range(station1_Simulations[80017,])[1]
min_tmp2 <- range(station1_Simulations[83014,])[1]
f1 <- MASS::kde2d(log(station1_Simulations[80017,]), log(station1_Simulations[83014,]), h=0.4, n = 200, 
                  lims = c(log(min_tmp1)-0.25, log(min_tmp1)+0.75, log(min_tmp2)-0.25, log(min_tmp2)+0.75))
min_tmp1 <- range(station50_Simulations[80017,])[1]
min_tmp2 <- range(station50_Simulations[83014,])[1]
f2 <- MASS::kde2d(log(station50_Simulations[80017,]), log(station50_Simulations[83014,]), h=0.5, n = 200, 
                  lims = c(log(min_tmp1)-0.25, log(min_tmp1)+0.85, log(min_tmp2)-0.25, log(min_tmp2)+0.85))
min_tmp1 <- range(station100_Simulations[80017,])[1]
min_tmp2 <- range(station100_Simulations[83014,])[1]
f3 <- MASS::kde2d(log(station100_Simulations[80017,]), log(station100_Simulations[83014,]), h=0.4, n = 200, 
                  lims = c(log(min_tmp1)-0.25, log(min_tmp1)+0.65, log(min_tmp2)-0.25, log(min_tmp2)+0.65))
upper <- max(f1$z, f2$z, f3$z)
lower_x <- min(f1$x, f2$x, f3$x); upper_x <- max(f1$x, f2$x, f3$x)
lower_y <- min(f1$y, f2$y, f3$y); upper_y <- max(f1$y, f2$y, f3$y) 


pal <- hcl.colors(10, "Spectral")[10:1]

plot_dat1 <- expand.grid(x=f1$x, y=f1$y)
plot_dat1$den <- as.vector(f1$z)
plot1 <- ggplot(plot_dat1) + geom_raster(aes(x=x, y=y, fill=den)) +
  geom_contour(col="black", aes(x=x, y=y, z=den)) +
  labs(x= expression(log(X[80017])), y = expression(log(X[83014]))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + ggtitle("Time 1")+
  scale_fill_gradientn(colours = pal, name = "Density", na.value = NA, limits = c(0, upper)) +
  theme(plot.title = element_text(hjust = 0.5, face='bold'))


plot_dat2 <- expand.grid(x=f2$x, y=f2$y)
plot_dat2$den <- as.vector(f2$z)
plot2 <- ggplot(plot_dat2) + geom_raster(aes(x=x, y=y, fill=den)) +
  geom_contour(col="black", aes(x=x, y=y, z=den)) +
  labs(x= expression(log(X[80017])), y = expression(log(X[83014]))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + ggtitle("Time 50")+
  scale_fill_gradientn(colours = pal, name = "Density", na.value = NA, limits = c(0, upper)) +
  theme(plot.title = element_text(hjust = 0.5, face='bold'))

plot_dat3 <- expand.grid(x=f3$x, y=f3$y)
plot_dat3$den <- as.vector(f3$z)
plot3 <- ggplot(plot_dat3) + geom_raster(aes(x=x, y=y, fill=den)) +
  geom_contour(col="black", aes(x=x, y=y, z=den)) +
  labs(x= expression(log(X[80017])), y = expression(log(X[83014]))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + ggtitle("Time 100")+
  scale_fill_gradientn(colours = pal, name = "Density", na.value = NA, limits = c(0, upper)) +
  theme(plot.title = element_text(hjust = 0.5, face='bold'))

library(cowplot)
prow <- plot_grid(plot1+ theme(legend.position="none"), plot2+ theme(legend.position="none"), plot3+ theme(legend.position="none"),
                  align = "h", nrow = 1)
legend <- get_legend(
  # create some space to the left of the legend
  plot3 + theme(legend.box.margin = margin(0, 0, 0, 0))
)
plot_grid(prow, legend, rel_widths = c(3, .4))
ggsave(file = paste0("./Figures/bivariate_pdf3.pdf"), width = 8, height =2.5)
