##############################################################################
### Gibbs sampler for Gaussian MTD
##############################################################################
gibbsGMTD <- function(obs, mtdorder, weight_param, prior, tuning, starting, updateWeight, niter, nburn, nthin){
  
  # prior
  mu_0 <- prior$mu_0
  sigma2_0 <- prior$sigma2_0
  u_0 <- prior$u_0
  v_0 <- prior$v_0
  
  # tuning parameters
  step_size <- tuning$step_size
  
  #################
  # initialization
  #################
  ## initial values of the parameters
  rho <- starting$rho
  sigma2 <- starting$sigma2
  lag_mat <- as.matrix(lagmat(obs, mtdorder))
  dat <- obs[-(1:mtdorder)]
  trun_size <- length(obs) - mtdorder
  weight <- extraDistr::rdirichlet(1, rep(1 / mtdorder, mtdorder))
  data_label <- sample(1:mtdorder, trun_size, replace = TRUE,  weight)
  unique_label <- as.numeric(names(table(data_label)))
  n_j <- rep(0, mtdorder)
  n_j[unique_label] <- as.numeric(table(data_label))
  dep_lag <- as.numeric(array(NA, dim = trun_size))
  for (t in 1:trun_size){
    dep_lag[t] <- lag_mat[t, data_label[t]]
  }
  dat_rho <- rho[data_label]
  logprho <- function(x, dd, lags, mu, sigma2) {
    logp <- dnorm(dd, (1 - x) * mu + x * lags, sqrt(sigma2 * (1 - x^2)), log = TRUE)
    sum(logp)
  }
  
  ## empty stuffs to save samples
  weight_save <- array(NA, dim = c(mtdorder, niter - nburn))
  rho_save <- array(NA, dim = c(mtdorder, niter - nburn))
  mu_save <- array(NA, dim = niter - nburn)
  sigma2_save <- array(NA, dim = niter - nburn)
  
  #######
  # mcmc
  #######
  for (iter in 1:niter){
    
    ## update mu
    mu <- updateGMTDmu(dat, dep_lag, dat_rho, sigma2, mu_0, sigma2_0)

    ## update sigma2
    sigma2 <- updateGMTDsigma2(dat, dep_lag, dat_rho, mu, trun_size, u_0, v_0)
    
    ## update rho
    rho <- updateGMTDrho(rho, dat, dep_lag, data_label, mtdorder, step_size, mu, sigma2, logprho)
    
    ## update labels
    labels <- updateGMTDlabel(dat, lag_mat, mtdorder, weight, rep(mu, mtdorder), rep(sigma2, mtdorder), rho)
    data_label <- as.numeric(labels$data_label)
    dep_lag <- as.numeric(labels$dep_lag)
    unique_label <- as.numeric(names(table(data_label)))
    n_j <- rep(0, mtdorder)
    n_j[unique_label] <- as.numeric(table(data_label))
    dat_rho <- rho[data_label]
    
    ## sample weights
    weight <- updateWeight(n_j, weight_param)
    
    ## save samples
    if (iter > nburn){
      weight_save[, iter - nburn] <- weight
      rho_save[, iter - nburn] <- rho
      mu_save[iter - nburn] <- mu
      sigma2_save[iter - nburn] <- sigma2
    }
  }
  
  ######################
  # thinning and output
  ######################
  selc_index <- seq(1, niter - nburn, by = nthin)
  list(weight = weight_save[, selc_index], 
       rho = rho_save[, selc_index],
       mu = mu_save[selc_index],
       sigma2 = sigma2_save[selc_index])
  
}

##############################################################################
### Gibbs sampler for Poisson MTD
##############################################################################
gibbsPMTD <- function(obs, mtdorder, weight_param, prior, tuning, starting, updateWeight, niter, nburn, nthin){
  
  # Prior
  u_la <- prior$u_la
  v_la<- prior$v_la
  u_th <- prior$u_th
  v_th <- prior$v_th
  
  #################
  # initialization
  #################
  ## initial values of the parameters
  lag_mat <- as.matrix(lagmat(obs, mtdorder))
  dat <- obs[-(1:mtdorder)]
  trun_size <- length(obs) - mtdorder
  weight <- extraDistr::rdirichlet(1, rep(1 / mtdorder, mtdorder))
  data_label <- sample(1:mtdorder, trun_size, replace = TRUE,  weight)
  unique_label <- as.numeric(names(table(data_label)))
  n_j <- as.numeric(table(data_label))
  dep_lag <- array(NA, dim = trun_size)
  for (t in 1:trun_size) dep_lag[t] <- lag_mat[t, data_label[t]]
  temp_ut <- rpois(trun_size, 1)
  temp_ut2 <- apply(cbind(temp_ut, dat - dep_lag), 1, max)
  ut <- apply(cbind(temp_ut2, dat), 1, min)
  
  ## empty stuffs to save samples
  weight_save <- array(NA, dim = c(mtdorder, niter - nburn))
  ut_save <- array(NA, dim = c(trun_size, niter - nburn))
  la_save <- array(NA, dim = niter - nburn)
  th_save <- array(NA, dim = niter - nburn)
  
  #######
  # mcmc
  #######
  for (iter in 1:niter){
    
    ## update th
    th <- updatePMTDth(dat, dep_lag, ut, u_th, v_th)
 
    ## update lambda
    la <- updatePMTDla(ut, trun_size, u_la, v_la)

    ## update ut
    ut <- as.numeric(updatePMTDut(ut, dat, dep_lag, la, la/(1-th)-la))
    
    ## update labels
    labels <- updatePMTDlabel(dat, lag_mat, mtdorder, weight, th, ut)
    data_label <- as.numeric(labels$data_label)
    dep_lag <- as.numeric(labels$dep_lag)
    unique_label <- as.numeric(names(table(data_label)))
    n_j <- rep(0, mtdorder)
    n_j[unique_label] <- as.numeric(table(data_label))
    
    ## update weights
    weight <- updateWeight(n_j, weight_param)
    
    ## save samples
    if (iter > nburn){
      weight_save[, iter - nburn] <- weight
      ut_save[, iter - nburn] <- ut
      la_save[iter - nburn] <- la
      th_save[iter - nburn] <- th
    }
  }
  
  ####################
  # thinning and output
  ######################  
  selc_index <- seq(1, niter - nburn, by = nthin)
  list(weight = weight_save[, selc_index],
       ut = ut_save[, selc_index],
       la = la_save[selc_index],
       th = th_save[selc_index])
}

##############################################################################
### Gibbs sampler for negative binomial MTD
#############################################################################
gibbsNBMTD <- function(obs, mtdorder, weight_param, prior, tuning, starting, updateWeight, niter, nburn, nthin) {
  
  # prior
  u_th <- prior$u_th
  v_th <- prior$v_th
  u_psi <- prior$u_psi
  v_psi <- prior$v_psi
  u_kap <- prior$u_kap
  v_kap <- prior$v_kap
  
  # tuning parameters
  se_kap <- tuning$se_kap
  
  #################
  # initialization
  #################
  ## initial values of the parameters
  kap <- starting$kap
  lag_mat <- as.matrix(lagmat(obs, mtdorder))
  dat <- obs[-(1:mtdorder)]
  trun_size <- length(obs) - mtdorder
  weight <- extraDistr::rdirichlet(1, rep(1 / mtdorder, mtdorder))
  data_label <- sample(1:mtdorder, trun_size, replace = TRUE,  weight)
  unique_label <- as.numeric(names(table(data_label)))
  n_j <- as.numeric(table(data_label))
  dep_lag <- array(NA, dim = trun_size)
  for (t in 1:trun_size) dep_lag[t] <- lag_mat[t, data_label[t]]
  temp_ut <- rpois(trun_size, 1)
  temp_ut2 <- apply(cbind(temp_ut, dat - dep_lag), 1, max)
  ut <- apply(cbind(temp_ut2, dat), 1, min)
  
  ## empty stuffs to save samples
  weight_save <- array(NA, dim = c(mtdorder, niter - nburn))
  ut_save <- array(NA, dim = c(trun_size, niter - nburn))
  psi_save <- array(NA, dim = niter - nburn)
  th_save <- array(NA, dim = niter - nburn)
  kap_save <- array(NA, dim = niter - nburn)
  kap_accept_count <- array(0, dim = niter)
  
  #######
  # mcmc
  #######
  for (iter in 1:niter){
    
    ## update th
    th <- updateNBMTDth(dat, dep_lag, ut, u_th, v_th)

    ## update psi
    psi <- updateNBMTDpsi(dep_lag, ut, trun_size, kap, u_psi, v_psi)

    ## update ut
    ut <- as.numeric(updateNBMTDut(ut, dat, dep_lag, th, psi, kap))
    
    ## update kap
    kap_res <- updateNBMTDkap(kap, se_kap, ut, dep_lag, psi, u_kap, v_kap)
    kap <- kap_res$kap
    if (kap_res$accept) kap_accept_count[iter] <- 1
    
    ## update labels
    labels <- updateNBMTDlabel(dat, lag_mat, mtdorder, weight, th, ut, psi, kap)
    data_label <- as.numeric(labels$data_label)
    dep_lag <- as.numeric(labels$dep_lag)
    unique_label <- as.numeric(names(table(data_label)))
    n_j <- rep(0, mtdorder)
    n_j[unique_label] <- as.numeric(table(data_label))
    
    ## update weights
    weight <- updateWeight(n_j, weight_param)
    
    ## save samples
    if (iter > nburn){
      weight_save[, iter - nburn] <- weight
      ut_save[, iter - nburn] <- ut
      psi_save[iter - nburn] <- psi
      th_save[iter - nburn] <- th
      kap_save[iter - nburn] <- kap
    }
  }
  
  ####################
  # thinning and output
  ######################  
  selc_index <- seq(1, niter - nburn, by = nthin)
  list(weight = weight_save[, selc_index],
       ut = ut_save[, selc_index],
       psi = psi_save[selc_index],
       th = th_save[selc_index],
       kap = kap_save[selc_index],
       accept_cout = kap_accept_count)
}

##############################################################################
### Gibbs sampler for Lomax MTD regression model 
#############################################################################
gibbsRegLMTD <- function(obs, XX, mtdorder, weight_param, prior, tuning, starting, updateWeight, niter, nburn, nthin){
  
  # prior
  u_alpha <- prior$u_alpha
  v_alpha <- prior$v_alpha
  u_phi <- prior$u_phi
  v_phi <- prior$v_phi
  
  # tuning parameter
  se_phi <- tuning$se_phi
  se_bb <- tuning$se_bb
  
  #################
  # initialization
  ################
  ## initial values of the parameters
  phi <- starting$phi
  bb <- starting$bb
  mu_obs <- exp(XX %*% bb)
  eps_lag_mat <- as.matrix(lagmat(obs / mu_obs, mtdorder))
  eps_dat <- (obs / mu_obs)[-(1:mtdorder)]
  trun_size <- length(obs) - mtdorder
  weight <- extraDistr::rdirichlet(1, rep(1 / mtdorder, mtdorder))
  data_label <- sample(1:mtdorder, trun_size, replace = TRUE,  weight)
  unique_label <- as.numeric(names(table(data_label)))
  n_j <- as.numeric(table(data_label))
  eps_dep_lag <- array(NA, dim = trun_size)
  for (t in 1:trun_size){
    eps_dep_lag[t] <- eps_lag_mat[t, data_label[t]]
  }
  
  ## empty stuffs to save samples
  weight_save <- array(NA, dim = c(mtdorder, niter - nburn))
  alpha_save <- array(NA, dim = niter - nburn)
  phi_save <- array(NA, dim = niter - nburn)
  bb_save <- array(NA, dim = c(ncol(XX), niter - nburn))
  phi_accept_count <- array(0, dim = niter)
  bb_accept_count <- array(0, dim = c(length(bb), niter))
  
  # mcmc
  
  for (iter in 1:niter){
    
    ## update alpha
    alpha <- updateRegLMTDalpha(eps_dat, eps_dep_lag, trun_size, phi, u_alpha, v_alpha)
    
    ## update phi
    phi_res <- updateRegLMTDphi(phi, se_phi, eps_dat, eps_dep_lag, trun_size, u_alpha, v_alpha, u_phi, v_phi)
    phi <- phi_res$phi
    if (phi_res$accept) phi_accept_count[iter] <- 1
    
    ## update beta
    bb_res <- updateRegLMTDbeta(bb, se_bb, bb_accept_count, eps_dat, eps_dep_lag, eps_lag_mat,
                                mtdorder, obs, XX, data_label, trun_size, phi, u_alpha, v_alpha, iter)
    
    bb <- bb_res$bb
    bb_accept_count <- bb_res$bb_accept_count
    eps_dat <- bb_res$eps_dat
    eps_dep_lag <- bb_res$eps_dep_lag
    eps_lag_mat <- bb_res$eps_lag_mat
    
    ## update labels
    labels <- updateLMTDlabel(eps_dat, eps_lag_mat, mtdorder, weight, rep(alpha, mtdorder), rep(phi, mtdorder))
    data_label <- as.numeric(labels$data_label)
    eps_dep_lag <- as.numeric(labels$dep_lag)
    unique_label <- as.numeric(names(table(data_label)))
    n_j <- rep(0, mtdorder)
    n_j[unique_label] <- as.numeric(table(data_label))
    
    ## update weights
    weight <- updateWeight(n_j, weight_param)
    
    ## save samples
    if (iter > nburn){
      weight_save[, iter - nburn] <- weight
      alpha_save[iter - nburn] <- alpha
      phi_save[iter - nburn] <- phi
      bb_save[, iter - nburn] <- bb
    }
  }
  
  ######################
  # thinning and output
  ######################
  selc_index <- seq(1, niter - nburn, by = nthin)
  list(weight = weight_save[, selc_index],
       alpha = alpha_save[selc_index],
       phi = phi_save[selc_index],
       phi_accept_count = phi_accept_count,
       bb = bb_save[, selc_index])
}