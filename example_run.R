### Simple example to run NIFTY
##### Two-dim non-Gaussian example 
setwd("~/Documents/NIFTY_package")
library(GGally)
library(ggplot2)
library(infinitefactor)
source("Rfuncs/sampling_funcs.R") 
source("Rfuncs/run_NIFTY+.R")
n = 400
p = 2
y = matrix(0, n, p) 
set.seed(1)

normalize <- function(x){
  (x - mean(x)) / sd(x)
}
z1 = normalize(rbeta(n, .4, .4))   
z2 = normalize(rgamma(n, 1, 1))   
y[, 1:2] = cbind(z1, z2) 
y <- apply(y, 2, normalize)
    
### Run diffusion map for pre-train
q = 3
y_dm <- dm(y, 1e1, q, 1.5) 
K0 = 2
rep_per_k <- 1
y_augmented <- y_dm$y[,(q-K0+1):q] 
plot(y[,1],y_augmented[,1])
plot(y[,2],y_augmented[,2])
W0 <- t(matrix(rep(diag(K0), rep_per_k), nrow=K0))  
H0 = rep_per_k*K0 
 

### Train nifty and obtain mcmc samples
param_tr <- train_nifty(step_a = .01, step_u = .01, y = y, y_augmented = y_augmented, W = W0, H = H0, K = K0, n_iter = 5000, burn = 2000, L = 20, sigma_init = 1e-2, lam = 0.1)
u = param_tr$u
u_list=unlist(param_tr$u_list)
plot(array(u_list,c(100,5000))[5,] ) 
plot(unlist(lapply(param_tr$u_list, `[[`, 20)))
plot(unlist(lapply(param_tr$eta_list, `[[`, 1))[2000:4000]) 
plot(unlist(lapply(param_tr$sig_list, `[[`, 1)))
plot(unlist(lapply(param_tr$sig_list, `[[`, 2)))

u = param_tr$u
eta_tr = param_tr$eta
Gamma_tr = param_tr$Gamma
intercept_tr = param_tr$intercept
y_fit = t(t(eta_tr%*%t(Gamma_tr))+intercept_tr) 
par(mfrow=c(1,2))
plot(y_fit) 
plot(y)

plot(y[,1]-y_fit[,1])
plot(y[,2]-y_fit[,2])

