train_nifty <- function(y, y_augmented, n_iter, burn, lam, W, H, K, L, step_a = 1e-1, step_u = 1e-1, sigma_init = 1e-3){
  n = length(y)
  y_star <- cbind(y, y_augmented)      
  ###Initialization   
  n <- dim(y_star)[1]
  p <- dim(y)[2] 
  u <- apply(y_augmented, 2, function(x)((x - min(x))/(max(x)-min(x)))*(1-1e-4)) 
  alpha <- array(runif(L*H)*5, c(L, H))
  alpha_aug <- array(runif(L*K)*5, c(L, K))
  
  intercept <- runif(p) 
  intercept_aug <- runif(K)
  Gamma <- matrix(1, p, H)  
  
  #####param for horseshoe prior 
  nu = 1
  xi = 1
  global = 1e-3
  local = Gamma
  fix_sigma = FALSE
  sig_G <- 1
  sig_I <- 1
  sig_a <- 1
  sigma2 <- as.vector(rep(sigma_init, p)) 
  sigma2_aug <- as.vector(rep(1e-3, K)) 
  if (K==1){sigma2_aug = as.matrix(sigma_init)}
  u_piece <- compute_upiece(u, W, L)
  eta <- compute_eta(alpha, u_piece, intercept, n, p, H, K, L)
  u_piece_aug <- compute_upiece(u, diag(K), L)
  eta_aug <- compute_eta(alpha_aug, u_piece_aug, intercept_aug, n, K, K, K, L)
  
  ###Initialize return list
  pb = txtProgressBar(1, n_iter, style= 3)
  count_a = 0
  count_u = 0
  u_list <- vector("list",n_iter)
  a_list <- numeric(n_iter) 
  sig_list <- vector("list",n_iter)
  p_star <- p + K
  log_prob_list <- numeric(n_iter)
  
  Omega <- diag(L)*2
  Omega[1,1] <- 1
  Omega[L,L] <- 1
  for (i in 1:(L-1)){
    Omega[i,i+1] <- -1
    Omega[i+1,i] <- -1
  }
  ###hyper-parameters
  a_sigma = n * 100
  b_sigma = 1 
  eps_u = 1e-4 
  eps_a = 1e-4
  Omega <- (Omega + diag(L)*.1)*10
  Gamma_aug <- (runif(K))
  if (K==1){Gamma_aug = as.matrix(Gamma_aug)}
  Gamma_list = vector("list", n_iter)
  eta_list = vector("list", n_iter)
  alpha_list = vector("list", n_iter)
  int_list = vector("list", n_iter) 
  y_list <- vector("list", n_iter)
  u_list <- vector("list", n_iter)
  
  for (step in 1:n_iter){
    Gamma <- sample_Gamma(t(t(y)-intercept), Gamma, eta, sigma2, sig_G, n, p, H, L, K, global, local)
    ln <- sample_local(y, Gamma, eta, sigma2, n, p, H, L, K, global, local, nu, xi)
    local <- ln$local 
    nu <- ln$nu  
    Gamma_aug <- pre_sample_Gamma(t(t(y_augmented)-intercept_aug), Gamma_aug, eta_aug, sigma2_aug, W, sig_G, n, K, K, L, K)
    intercept <- sample_intercept_p(y, intercept, Gamma, eta, sigma2, sig_I, n, p, L, K)
    intercept_aug<- sample_intercept_p(y_augmented, intercept_aug, diag(Gamma_aug), eta_aug, sigma2_aug, sig_I, n, K, L, K)
    eta <- compute_eta(alpha, u_piece, intercept, n, p, H, K, L)
    a_e <- sample_a(Omega,t(t(y)-intercept), sig_a, count_a, alpha, eta, eps_a, u_piece, intercept, W, Gamma, sigma2, n, p, H, K, L) 
    a_a <- sample_a(Omega,t(t(y_augmented)-intercept_aug), sig_a, count_a, alpha_aug, eta_aug, eps_a, u_piece_aug, intercept_aug, diag(K), diag(Gamma_aug), sigma2_aug, n, K, K, K, L)
    alpha_aug <- a_a[[1]]
    eta_aug <- a_a[[2]]
    alpha <- a_e[[1]]
    eta <- a_e[[2]]
    count_a = a_e[[3]] 
    u_e <-sample_u_star(y, y_augmented, count_u, u, eta, eta_aug, eps_u, u_piece, u_piece_aug, Gamma, diag(Gamma_aug), alpha, alpha_aug, sigma2, sigma2_aug, lam, intercept, intercept_aug, W, n, p, H, L, K)
    u <- u_e[[1]]
    u_piece <- u_e[[2]]
    eps_u = step_u/max(abs(u_e[[5]]))
    eps_a = step_a/max(abs(a_e[[4]]))
    u_piece_aug <- compute_upiece(u, diag(K), L)
    eta <- u_e[[3]]
    eta_aug <- compute_eta(alpha_aug, u_piece_aug, intercept_aug, n, K, K, K, L) 
    if (step > burn/2){sigma2 <- sample_sigma(y, Gamma, eta, sigma2, fix_sigma, a_sigma, b_sigma, n, p, H, L, K)}
    count_u <- u_e[[4]] 
    u_list[[step]] <- u[,1]
    a_list[step] <- alpha[1]
    ypre = (t(eta%*%t(Gamma))+intercept)  
    int_list[[step]]<- intercept
    y_list[[step]] <- ypre
    sig_list[[step]] <- sigma2 
    Gamma_list[[step]] <- Gamma 
    eta_list[[step]] <- eta
    alpha_list[[step]] <- alpha
    int_list[[step]] <- intercept
    log_prob_list[step] <- compute_logprob(y, K, H, Gamma, sigma2, intercept, W, alpha, u, eta)
    setTxtProgressBar(pb, step)
  }
  close(pb)  
  
  aligned = jointRot(Gamma_list, eta_list)
  Gamma_list = aligned$lambda
  eta_list = aligned$eta 
  Gamma_list <- lapply(Gamma_list, msf, pivot=Gamma_list[[n_iter]])
  return(list("Gamma_list"=Gamma_list, "y_list"=y_list," a_list"=a_list, "u_list"= u_list, "sig_list"=sig_list,"int_list"=int_list, "eta_list"=eta_list, "alpha_list"=alpha_list,
              "global"=global,"local"=local, 'llh_list' = log_prob_list,
              'u' = u, 'eta'=eta, 'alpha'=alpha, 'alpha_aug'=alpha_aug,  
              'intercept'=intercept, 'intercept_aug'=intercept_aug, 'Gamma'=Gamma, 'Gamma_aug'=Gamma_aug, 'sigma2'=sigma2, 'sigma2_aug'=sigma2_aug,
              'lam'=lam, 'W'=W, 'H'=H, 'K'=K, 'L'=L))
}
 
