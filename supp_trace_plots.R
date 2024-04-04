#repeat_trace_plots
#curve
# 
param_tr <- train_nifty(step_a = .001, step_u = .001, y = y_tr, y_augmented = y_augmented_tr, W = W0, H = H0, K = K0, n_iter = 30000, L = 5, sigma_init = 1e-3, lam = 0.1)

par(mfrow=c(1,1))
y_list = unlist(param_tr$y_list)
ind = seq(1, 1000*29999+1, 1000)

ylist = y_list[ind]
ind1 = seq(20000,30000,1)
plot(detrend(ylist[ind1])*2+1.2,type='l',lwd=.5)

acf(ylist[6000:11000])

library(pracma)
llh_list = param_tr$llh_list
ind <- seq(20000,30000,1)
llh <- llh_list[ind] 
#llh[llh< -1e6 ] = median(llh) + rnorm(1)
plot(llh/300,type='l',lwd=.5)

llh<-detrend(llh, tt = 'linear', bp = c())
plot(llh/1000-1000,type='l')
acf(llh)
 

#gaussian
param_tr <- train_nifty(step_a = .001, step_u = .01, y = y_tr, y_augmented = y_augmented_tr, W = W0, H = H0, K = K0, n_iter = 20000, L = 10, lam = 0.05)

par(mfrow=c(1,1))
y_list = unlist(param_tr$y_list)
ind = seq(0000*2000, 2000*19999+1, 2000)
ylist = y_list[ind]
ind1 = seq(10000,19000,2)
plot(detrend(ylist[ind1])*10-3,type='l',lwd=.5)

library(pracma)
llh_list = param_tr$llh_list
ind <- seq(10000,20000,2)
llh <- llh_list[ind] 
#llh[llh< -1e6 ] = median(llh) + rnorm(1)
plot(detrend(llh/100)-400,type='l',lwd=.5)
acf(llh)
#mariginal
param_tr <- train_nifty(step_a = .01, step_u = .01, y = y_tr, y_augmented = y_augmented_tr, W = W0, H = H0, K = K0, n_iter = 20000, L = 8, lam = 0.2, sigma_init = 1e-2)


par(mfrow=c(1,1))
y_list = unlist(param_tr$y_list)
ind = seq(0000*2000, 2000*19999+1, 2000)
ylist = y_list[ind]
ind1 = seq(10000,19000,2)
plot((ylist[ind1]),type='l',lwd=1)

library(pracma)
llh_list = param_tr$llh_list
ind <- seq(10000,20000,5)
llh <- llh_list[ind] 
#llh[llh< -1e6 ] = median(llh) + rnorm(1)
plot(detrend(llh/100)-500,type='l',lwd=1)



##bird
par(mfrow=c(1,1))
y_list = unlist(param_tr$y_list)
ind = seq(2, 73000*9999+2, 73000)

ylist = y_list[ind]
ind1 = seq(5000,10000,2)
plot((ylist[ind1])*10,type='l',lwd=.5)

acf(ylist[6000:10000])

library(pracma)
llh_list = param_tr$llh_list
ind <- seq(2000,10000,1)
llh <- llh_list[ind] 
#llh[llh< -1e6 ] = median(llh) + rnorm(1)
plot(llh/100,type='l',lwd=.5)

llh<-detrend(llh, tt = 'linear', bp = c())
plot(llh/1000-1000,type='l')
acf(llh)
