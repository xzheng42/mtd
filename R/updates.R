##############################################################################
### updates for the weights
##############################################################################
updateDIRweight <- function(m, weight_param) {
  dir_shape <- weight_param$dir_shape
  extraDistr::rdirichlet(1, m + dir_shape)
}

updateSBweight <- function(m, weight_param) {
  alpha <- weight_param$alpha
  mtdorder <- weight_param$mtdorder
  m_cumsum <- cumsum(rev(m[-1]))
  stick_break <- c(rbeta(mtdorder - 1, 1 + m[-mtdorder], alpha + rev(m_cumsum)), 1)
  log_weight <- log(stick_break) + c(0, cumsum(log(1 - stick_break[-mtdorder])))
  exp(log_weight)
}

updateCDPweight <- function(m, weight_param) {
  cutoff <- weight_param$cutoff
  alpha_0 <- weight_param$alpha_0
  a_G0 <- weight_param$a_G0
  b_G0 <- weight_param$b_G0
  pp <- pbeta(cutoff, a_G0, b_G0)
  extraDistr::rdirichlet(1, m + alpha_0 * diff(pp))
}

##############################################################################
### updates for the Gaussian MTD
##############################################################################
updateGMTDmu <- function(dat, dep_lag, dat_rho, sigma2, mu_0, sigma2_0) {
  b_mu <- sum((1 - dat_rho)^2 / (1 - dat_rho^2))
  c_mu <- sum((1 - dat_rho) * (dat - dat_rho * dep_lag) / (1 - dat_rho^2))
  sigma2_1 <- 1 / (1 / sigma2_0 + b_mu / sigma2)
  mu_1 <- sigma2_1 * (mu_0 / sigma2_0 + c_mu / sigma2)
  rnorm(1, mu_1, sqrt(sigma2_1))
}

updateGMTDsigma2 <- function(dat, dep_lag, dat_rho, mu, trun_size, u_0, v_0) {
  u_1 <- u_0 + trun_size / 2
  ee <- (dat - mu - dat_rho * (dep_lag - mu))^2 / (1 - dat_rho^2)
  v_1 <- v_0 + sum(ee) / 2
  1 / rgamma(1, u_1, v_1)
}

updateGMTDrho <- function(rho, dat, dep_lag, data_label, mtdorder, step_size, mu, sigma2, logf) {
  for (l in 1:mtdorder) {
    dd_l <- dat[which(data_label==l)]
    lags_l <- dep_lag[which(data_label==l)]
    rho[l] <- uniSlice(x0 = rho[l], g = logf, w = step_size, m = Inf,
                       lower = -1, upper = 1, dd_l, lags_l, mu, sigma2)
  }
  rho
}

##############################################################################
### updates for the Poisson MTD
##############################################################################
updatePMTDth <- function(dat, dep_lag, ut, u_th, v_th) {
  rbeta(1, u_th + sum(dat-ut), v_th + sum(dep_lag - dat + ut))
}

updatePMTDla <- function(ut, trun_size, u_la, v_la) {
  rgamma(1, u_la + sum(ut), v_la + trun_size)
}

##############################################################################
### updates for the Negative binomial MTD
##############################################################################
updateNBMTDth <- function(dat, dep_lag, ut, u_th, v_th) {
  rbeta(1, u_th + sum(dat - ut), v_th + sum(dep_lag - dat + ut))
}

updateNBMTDpsi <- function(dep_lag, ut, trun_size, kap, u_psi, v_psi) {
  rbeta(1, u_psi + trun_size * kap + sum(dep_lag), v_psi + sum(ut))
}

updateNBMTDkap <- function(kap, se_kap, ut, dep_lag, alp, u_kap, v_kap) {
  log_kap <- log(kap)
  prop_log_kap <- rnorm(1, log_kap, se_kap)
  prop_kap <- exp(prop_log_kap)
  prop_loglik <- 
    sum(dnbinom(ut, prop_kap + dep_lag, alp, log = TRUE)) + 
    dgamma(prop_kap, u_kap, v_kap, log = TRUE) + prop_log_kap
  cur_loglik <- 
    sum(dnbinom(ut, kap + dep_lag, alp, log = TRUE)) + 
    dgamma(kap, u_kap, v_kap, log = TRUE) + log_kap
  diff_loglik <- prop_loglik - cur_loglik
  if (diff_loglik > log(runif(1))) {
    return(list(kap = prop_kap, accept = TRUE))
  } else {
    return(list(kap = kap, accept = FALSE))
  }
}

#############################################################################
### updates for the Lomax MTD regression model 
#############################################################################
updateRegLMTDalpha <- function(eps_dat, eps_dep_lag, trun_size, phi, u_alpha, v_alpha) {
  u_alpha_1 <- u_alpha + trun_size
  v_alpha_1 <- v_alpha + sum(log(phi + eps_dat + eps_dep_lag) -  log(phi + eps_dep_lag))
  rgamma(1, u_alpha_1, v_alpha_1)
}

updateRegLMTDphi <- function(phi, se_phi, eps_dat, eps_dep_lag, trun_size, u_alpha, v_alpha, u_phi, v_phi) {
  log_phi <- log(phi)
  prop_log_phi <- rnorm(1, log_phi, se_phi)
  prop_phi <- exp(prop_log_phi)
  prop_loglik <-
    - (trun_size + u_alpha) * 
    log(sum(log(1 + eps_dat / (prop_phi + eps_dep_lag))) + v_alpha) -
    sum(log(prop_phi + eps_dep_lag)) - 
    sum(log(1 + eps_dat / (prop_phi + eps_dep_lag))) +
    dgamma(1 / prop_phi, u_phi, v_phi, log = TRUE) + 2 * log(1 / prop_phi) + prop_log_phi 
  cur_loglik <-
    - (trun_size + u_alpha) * 
    log(sum(log(1 + eps_dat / (phi + eps_dep_lag))) + v_alpha) -
    sum(log(phi + eps_dep_lag)) - 
    sum(log(1 + eps_dat / (phi + eps_dep_lag))) +
    dgamma(1 / phi, u_phi, v_phi, log = TRUE) + 2 * log(1 / phi) + log_phi
  diff_loglik <- prop_loglik  - cur_loglik
  if (diff_loglik > log(runif(1))) {
    return(list(phi = prop_phi, accept = TRUE))
  } else {
    return(list(phi = phi, accept = FALSE))
  }
}

updateRegLMTDbeta <- function(bb, se_bb, bb_accept_count, eps_dat, eps_dep_lag, eps_lag_mat,
                              mtdorder, obs, XX, data_label, trun_size, phi, u_alpha, v_alpha, iter) {
  for (j in 1:length(bb)) {
    prop_bb <- bb
    prop_bb[j] <- rnorm(1, bb[j], se_bb[j])
    prop_mu_obs <- exp(XX %*% prop_bb)
    prop_eps_lag_mat <- as.matrix(lagmat(obs / prop_mu_obs, mtdorder))
    prop_eps_dat <- (obs / prop_mu_obs)[-(1:mtdorder)]
    prop_eps_dep_lag <- array(NA, dim = trun_size)
    for (t in 1:trun_size) {
      prop_eps_dep_lag[t] <- prop_eps_lag_mat[t, data_label[t]]
    }
    prop_loglik <- 
      - (trun_size + u_alpha) *
      log(sum(log(1 + prop_eps_dat / (phi + prop_eps_dep_lag))) + v_alpha) -
      sum(log(phi + prop_eps_dep_lag)) -
      sum(log(1 + prop_eps_dat / (phi + prop_eps_dep_lag))) -
      sum((XX %*% prop_bb)[(-(1:mtdorder))])
    cur_loglik <-
      - (trun_size + u_alpha) *
      log(sum(log(1 + eps_dat / (phi + eps_dep_lag))) + v_alpha) -
      sum(log(phi + eps_dep_lag)) -
      sum(log(1 + eps_dat / (phi + eps_dep_lag))) -
      sum((XX %*% bb)[-(1:mtdorder)])
    diff_loglik <- prop_loglik - cur_loglik
    if (diff_loglik > log(runif(1))) {
      bb <- prop_bb
      eps_dat <- prop_eps_dat
      eps_dep_lag <- prop_eps_dep_lag
      eps_lag_mat <- prop_eps_lag_mat
      bb_accept_count[j, iter] <- 1
    }
  }
  list(bb = bb, bb_accept_count = bb_accept_count, eps_dat = eps_dat, eps_dep_lag = eps_dep_lag, eps_lag_mat = eps_lag_mat)
}

#############################################################################
### slice sampler
### the code is obtained from the https://www.cs.toronto.edu/~radford/ftp/slice-R-prog.
### the code was developed by Radford Neal, and is available for free use;
### see https://www.cs.toronto.edu/~radford/software-online.html.
#############################################################################
uniSlice <- function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, ...)
{
  # Check the validity of the arguments.
  
  if (!is.numeric(x0) || length(x0)!=1
      || !is.function(g) 
      || !is.numeric(w) || length(w)!=1 || w<=0 
      || !is.numeric(m) || !is.infinite(m) && (m<=0 || m>1e9 || floor(m)!=m)
      || !is.numeric(lower) || length(lower)!=1 || x0<lower
      || !is.numeric(upper) || length(upper)!=1 || x0>upper
      || upper<=lower 
      # || !is.null(gx0) && (!is.numeric(gx0) || length(gx0)!=1)
  )
  { 
    stop ("Invalid slice sampling argument")
  }
  
  # Keep track of the number of calls made to this function.
  
  #uni.slice.calls <<- uni.slice.calls + 1
  
  # Find the log density at the initial point, if not already known.
  
  # if (is.null(gx0))
  # { uni.slice.evals <<- uni.slice.evals + 1
  # gx0 <- g(x0)
  # }
  gx0 <- g(x0, ...)
  
  # Determine the slice level, in log terms.
  
  logy <- gx0 - rexp(1)
  
  # Find the initial interval to sample from.
  
  u <- runif(1, 0, w)
  L <- x0 - u
  R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff
  
  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.
  
  if (is.infinite(m))  # no limit on number of steps
  { 
    repeat
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L, ...) <= logy) break
      L <- L - w
    }
    
    repeat
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R, ...) <= logy) break
      R <- R + w
    }
  }
  
  else if (m>1)  # limit on steps, bigger than one
  { 
    J <- floor(runif(1,0,m))
    K <- (m-1) - J
    
    while (J>0)
    { if (L<=lower) break
      uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
      J <- J - 1
    }
    
    while (K>0)
    { if (R>=upper) break
      uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
      K <- K - 1
    }
  }
  
  # Shrink interval to lower and upper bounds.
  
  if (L<lower) 
  { L <- lower
  }
  if (R>upper)
  { R <- upper
  }
  
  # Sample from the interval, shrinking it on each rejection.
  
  repeat
  { 
    x1 <- runif(1, L, R)
    
    #uni.slice.evals <<- uni.slice.evals + 1
    gx1 <- g(x1, ...)
    
    if (gx1 >= logy) break
    
    if (x1 > x0) 
    { R <- x1
    }
    else 
    { L <- x1
    }
  }
  
  # Return the point sampled, with its log density attached as an attribute.
  
  #attr(x1,"log.density") <- gx1
  return (x1)
}