### Simple example to run NIFTY on a curve-shaped data
##### Two-dim curve example
setwd("~/Documents/NIFTY_package")
library(GGally)
library(ggplot2)
library(infinitefactor)
source("Rfuncs/sampling_funcs.R") 
source("Rfuncs/run_NIFTY+.R")
###2-factors
L1_list = matrix(0, 100, 4)
L2_list = matrix(0, 100, 4)
W_list = matrix(0, 100, 4)
n <- 500 
m_list = c(100,200,500,800)
mid = 400
seed = 1
p = 10
y = matrix(0,n,p)
#Gamma0 = matrix(rnorm(p,4),4,p)
z1 = (rbeta(n, .5,.5)-.5) * 2 # rnorm(n)#
z2 = (rbeta(n, .5,.5)-.5) * 2 # rnorm(n)#
x1 = z1
x2 = z1^2
x3 = z2
x4 = z2^2
y[,1:4] = cbind(x1,x2,x3,x4)  
y = y + matrix(rnorm(n*p),n,p)*.01
set.seed(seed)
train = sample(1:n, mid, replace = FALSE)
y_tr = y[train,]
y_te = y[-train,]

#contour3d(y, color = U,  grid=F, box=F, pch = '.')
data_ <- y  
q = 4
eps_band = c(1, 2, 3)
C_band = c(1.5)
r = choose_band(y_tr, eps_band, q, C_band)  
eps = eps_band[r[1]]
C = C_band[r[2]]
print(c(eps,C))
y_dm <- dm(y_tr, .6, q, 1.2) 
y_embed <- y_dm$y   
plot(x1[train],y_embed[,4])

K0 = 2
y_augmented <-  y_embed[,(q-K0+1):q]
y_augmented_tr =cbind(y_tr[,1], y_tr[,3]) 
par(mfrow=c(2,2))
plot(y_tr)
plot(y_augmented[,1],y_tr[,1],main=seed)
plot(y_augmented[,1],y_tr[,3])
rep_per_k <- 2
W0 <- t(matrix(rep(diag(K0), rep_per_k), nrow=K0))  
H0 = rep_per_k*K0 
### Train nifty and obtain mcmc samples
param_tr <- train_nifty(step_a = .01, step_u = .01, y = y_tr, y_augmented = y_augmented_tr, W = W0, H = H0, K = K0, n_iter = 5000, burn = 2000, L = 15, sigma_init = 1e-2, lam = 0.1)
u = param_tr$u  
eta_tr = param_tr$eta
Gamma_tr = param_tr$Gamma
intercept_tr = param_tr$intercept
y_fit = t(t(eta_tr%*%t(Gamma_tr))+intercept_tr) 
par(mfrow=c(1,2))
plot(y_fit) 
plot(y) 
