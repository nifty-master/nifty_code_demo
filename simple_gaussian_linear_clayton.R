### Figure 1
### Generate data with clayton copula
par(mfrow=c(1,1))
n = 400
library(copula)
set.seed(100)
myCop <-  claytonCopula(param=3,dim=2)#normalCopula(param=c(0.4,0.2,-0.8), dim = 3, dispstr = "un")
myMvd <- mvdc(copula=myCop, margins=c("norm", "norm"), paramMargins=list(list(mean = 0, sd =1), list(mean = 0, sd =1)))
y <- rMvdc(n, myMvd)
plot(y, pch=16, col='steelblue')

def_break <- function(y){
  seq(range(y)[1],range(y)[2],(range(y)[2]-range(y)[1])/40)
}
par(mfrow=c(2,1))
hist(y[,1], xaxt='n', yaxt='n', xlim=range(y[,1]), breaks=def_break(y[,1]), col='steelblue')
hist(y[,2], xaxt='n', yaxt='n', xlim=range(y[,2]), breaks=def_break(y[,2]), col='steelblue')

### ppca sampler
par(mfrow=c(1,1))
sample_gamma <- function(sig_G, Gamma, sigma2, eta, y){
  H = dim(eta)[2]
  p = dim(y)[2]
  D_1 <- diag(H) / sig_G
  for (j in 1:p){
    v = solve(D_1 + t(eta)%*%eta/sigma2[j])
    vchol = chol(v)
    mean = v %*% t(eta) %*% y[,j]/ sigma2[j]
    Gamma[j,] <-  mean +vchol %*% rnorm(H)
  } 
  return(Gamma)
}

sample_sigma <- function(y, Gamma, sigma2, eta, a_sigma, b_sigma){  
  p = dim(y)[2]
  n = dim(y)[1]
  residual <- y - t(Gamma %*% t(eta))
  for (j in 1:p){
    sigma2[j] = 1 / rgamma(1, shape = a_sigma + n/2, rate = b_sigma + sum((residual[,j])^2/2))
  }
  return(sigma2)
}

sample_eta <- function(y, eta, Gamma, sigma2){
  for (i in 1:n){
    H = dim(Gamma)[2]
    v = solve(diag(H)+t(Gamma)%*%diag(1/sigma2)%*%Gamma)
    vchol=chol(v)
    mean = v %*% t(Gamma)%*%diag(1/sigma2)%*%y[i,]
    eta[i,] <- mean + vchol %*% rnorm(H)
  }
  return(eta)
}


### setting hyperparamters
p= 2
H = 2
Gamma = matrix(1,p,H)
sigma2 = as.vector(matrix(.01,p,1))
eta = y
a_sigma = 100
b_sigma = 1
sig_G = 1

threshold <- function(matrix, lam){
  return(sign(matrix) * (abs(matrix) - lam) * ((abs(matrix) - lam) > 0) )
}
pb = txtProgressBar(min = 0, max = 1500, initial = 0) 

### start sampling
for (steps in 1:5000){
  setTxtProgressBar(pb,steps)
  Gamma <- sample_gamma(sig_G, Gamma, sigma2, eta, y)
  theta <- threshold(Gamma, 1)
  theta_svd <- svd(theta)
  d <- theta_svd$d
  d1 <- threshold(d, 2)
  Gamma <- theta_svd$u %*% diag(sqrt(d1))
  sigma2 <- sample_sigma(y, Gamma, sigma2, eta, a_sigma, b_sigma)
  eta <- sample_eta(y, eta, Gamma, sigma2)
}
close(pb)
#plot(eta%*%t(Gamma))

par(mfrow=c(1,1))
plot(eta, pch=16, col='cadetblue4')

par(mfrow=c(2,1))
hist(eta[,1], xaxt='n', yaxt='n', xlim=c(-5,5), breaks=seq(-5,5,1/2), col='cadetblue4')
hist(eta[,2], xaxt='n', yaxt='n', xlim=c(-5,5), breaks=seq(-4,4,1/2), col='cadetblue4')
 

par(mfrow=c(1,1))
new_eta = matrix(rnorm(n*p),n,p)
new_y = new_eta%*%t(Gamma) + cbind(rnorm(n,sd=sigma2[1]),rnorm(n,sd=sigma2[2]))
plot(new_y, pch=16, col='chocolate')

par(mfrow=c(2,1))
hist(new_y[,1], xaxt='n', yaxt='n', xlim=c(-5,5), breaks=seq(-5,5,1/2), col='chocolate')
hist(new_y[,2], xaxt='n', yaxt='n', xlim=c(-5,5), breaks=seq(-4,4,1/2), col='chocolate')
