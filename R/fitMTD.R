#' Fit a time series MTD model
#'
#' This function fits an MTD model with a pre-specified stationary marginal distribution.
#
#' @param obs a vector of complete data (not truncated).
#' @param mtdorder the model order. Note mtdorder > 1.
#' @param family a quoted keyword that specifies the stationary marginal distribution.
#'        Supported keywords are: "gaussian", "poisson", "negative binomial".
#' @param weight a quoted keyword that specifies the prior for the weights.
#'        Supported keywords are: "dir", "sb", "cdp". 
#' @param prior a list of priors for model parameters.
#' @param tuning a list of tuning parameters  for model parameters.
#' @param starting a list of starting values for model parameters.
#' @param mcmc_param a list of MCMC parameters that contains, 
#'        niter: the number of iterations;
#'        nburn: the number of burn-in samples;
#'        nthin: the number of steps to obtain thinned samples.
#'
#' @references 
#' Zheng, X., Kottas, A., and Sansó, B. (2021),
#' "On construction and estimation of stationary mixture transition distribution models,"
#' \emph{arXiv:2010.12696}.
#' 
#' @export
#'
tsMTD <- function(obs, mtdorder, family, weight, prior, tuning, starting, mcmc_param) {
  
  ####################################################
  ### observation
  ####################################################
  if (length(obs) <= mtdorder) {
    stop("error: the data size must be greater than the model order.")
  }
  
  ####################################################
  ### model order
  ####################################################
  if (mtdorder%%1 !=0) {
    stop("error: the model order must be an integer.")
  }

  if (mtdorder < 2) {
    stop("error: the model order must be greater than two.")
  }
  
  ####################################################
  ### family
  ####################################################
  if (missing(family)) {
    stop("error: family must be specified")
  }
  
  family <- tolower(family)
  family_names <- c("gaussian", "poisson", "negative binomial")
  
  if (!family %in% family_names) {
    stop("error: specified family '", family, "' is not a valid option;
         available families are ", paste(family_names, collapse = ", ", sep = "") ,".")
  }
  
  ####################################################
  ### priors for the weights
  ####################################################  
  if (missing(weight)) {
    stop("error: weight must be specified")
  }
  
  weight <- tolower(weight)
  weight_names <- c("dir", "sb", "cdp")
  
  if (!weight %in% weight_names) {
    stop("error: specified prior for the weight '", weight, "' is not a valid option;
         available priors are ", 
         paste(weight_names, collapse = ", ", sep = "") ,".")
  }
  
  if (weight == "dir") {
    
    g <- updateDIRweight
    dir_shape <- prior$dir_shape
    
    if (is.null(dir_shape)) {
      stop("error: shape parameter for the Dirichlet prior must be specified.")
    }

    if (any(dir_shape <= 0)) {
      stop("error: dir_shap must be a positive vector of length mtdorder")
    }    
    
    weight_param <- list(dir_shape = dir_shape)
    
  } 
  else if (weight == "sb") {
    
    g <- updateSBweight
    alpha <- prior$alpha
    
    if (is.null(alpha)) {
      stop("error: concentration parameter for the truncated stick-breaking prior must be specified.")
    }
    
    if (alpha <= 0) stop("error: alpha must be positive.")

    weight_param <- list(alpha = alpha, mtdorder = mtdorder)
    
  } 
  else if (weight == "cdp") {
    
    g <- updateCDPweight
    alpha_0 <- prior$alpha_0
    a_G0 <- prior$a_G0
    b_G0 <- prior$b_G0
    
    if (is.null(alpha_0)) {
      stop("error: concentration parameter for the cdf-based prior must be specified.")
    }
    
    if (is.null(a_G0) || is.null(b_G0)) {
      stop("error: base distribution parameter for the cdf-based prior must be specified.")
    }

    if (alpha_0 <= 0) stop("error: alpha_0 must be positive.")
    if (a_G0 <= 0) stop("error: a_G0 must be positive")
    if (b_G0 <= 0) stop("error: b_G0 must be positive")
    
    weight_param <- list(alpha_0 = alpha_0, a_G0 = a_G0, b_G0 = b_G0, 
                         cutoff = seq(0, 1, length = mtdorder + 1))
    
  } 
  
  ####################################################
  ### mcmc
  #################################################### 
  niter <- mcmc_param$niter
  nburn <- mcmc_param$nburn
  nthin <- mcmc_param$nthin
  
  if (family == "gaussian") {
    
    ### check prior
    if (is.null(prior$mu_0)) stop("error: hyperparameter mu_0 must be specified.")
    if (is.null(prior$sigma2_0)) stop("error: hyperparameter sigma2_0 must be specified.")
    if (is.null(prior$u_0)) stop("error: hyperparameter u_0 must be specified.")
    if (is.null(prior$v_0)) stop("error: hyperparameter v_0 must be specified.")
    
    if (prior$sigma2_0 <= 0) stop("error: hyperparameter sigma2_0 must be positive.")
    if (prior$u_0 <= 0) stop("error: hyperparameter u_0 must be positive,")
    if (prior$v_0 <= 0) stop("error: hyperparameter v_0 must be positive.")
    
    ### check tuning parameter
    if (is.null(tuning$step_size)) stop("error: tuning parameter step_size must be specified.")
    if (tuning$step_size <= 0) stop("error: tuning parameter step_size must be positive.")
    
    ### check starting values
    if (is.null(starting$rho)) stop("error: starting value of rho must be specified.")
    if (any(starting$rho < 0, starting$rho >= 1)) {
      stop("error: starting value of rho must be a positive vector and each element must be between 0 and 1.")
    }
    if (is.null(starting$sigma2)) stop("error: starting value of sigma2 must be specified.")
    
    gibbsGMTD(obs, mtdorder, weight_param, prior, tuning, starting, g, niter, nburn, nthin)
    
  } 
  else if (family == "poisson") {
    
    ### check prior
    if (is.null(prior$u_la)) stop("error: hyperparameter u_la must be specified.")
    if (is.null(prior$v_la)) stop("error: hyperparameter v_la must be specified.")
    if (is.null(prior$u_th)) stop("error: hyperparameter u_th must be specified.")
    if (is.null(prior$v_th)) stop("error: hyperparameter v_th must be specified.")
    
    gibbsPMTD(obs, mtdorder, weight_param, prior, tuning, starting, g, niter, nburn, nthin)
    
  }
  else if (family == "negative binomial") {
    
    ### check prior
    if (is.null(prior$u_th)) stop("error: hyperparameter u_th must be specified.")
    if (is.null(prior$v_th)) stop("error: hyperparameter v_th must be specified.")
    if (is.null(prior$u_psi)) stop("error: hyperparameter u_psi must be specified.")
    if (is.null(prior$v_psi)) stop("error: hyperparameter v_psi must be specified.")
    if (is.null(prior$u_kap)) stop("error: hyperparameter u_kap must be specified.")
    if (is.null(prior$v_kap)) stop("error: hyperparameter v_kap must be specified.")
    
    ### check tuning parameter
    if (is.null(tuning$se_kap)) stop("error: tuning parameter se_kap must be specified.")
    if (tuning$se_kap <= 0) stop("error: tuning parameter step_size must be positive.")
    
    ### check starting values
    if (is.null(starting$kap)) stop("error: starting value of kap must be specified")
    if (starting$kap <= 0) stop("error: starting value kap must be positive.")
    
    gibbsNBMTD(obs, mtdorder, weight_param, prior, tuning, starting, g, niter, nburn, nthin)
    
  }
  
}

#' Fit a time series regression MTD model
#' 
#' This function fits a regression model of which the error term is an 
#' MTD model with a pre-specified stationary marginal distribution.
#
#' @param formula a symbolic description of the regression model. 
#'                Depending on the marginal distribution of the MTD, the regression
#'                model can be additive or multiplicative.
#' @param mtdorder the model order. Note mtdorder > 1.
#' @param family a quoted keyword that specifies the stationary marginal distribution.
#'        Supported keyword is "lomax".
#' @param weight a quoted keyword that specifies the prior for the weights.
#'        Supported keywords are: "dir", "sb", "cdp". 
#' @param prior a list of priors for model parameters.
#' @param tuning a list of tuning parameters  for model parameters.
#' @param starting a list of starting values for model parameters.
#' @param mcmc_param a list of MCMC parameters that contains, 
#'        niter: the number of iterations;
#'        nburn: the number of burn-in samples;
#'        nthin: the number of steps to obtain thinned samples.
#'
#' @references 
#' Zheng, X., Kottas, A., and Sansó, B. (2021),
#' "On construction and estimation of stationary mixture transition distribution models,"
#' \emph{arXiv:2010.12696}.
#'
#' @export
#'
tsRegMTD <- function(formula, mtdorder, family, weight, prior, tuning, starting, mcmc_param) {
  
  ####################################################
  ### formula
  ####################################################
  if (missing(formula)) {
    stop("error: the formula must be specified.")
  }
  
  if (class(formula) == "formula") {
    ff <- model.frame(formula)
    yy <- as.vector(ff[, 1])
    XX <- as.matrix(ff[,-1])
  } else {
    stop("error: the formula is misspecified.")
  }
  
  if (length(yy) <= mtdorder) {
    stop("error: the data size must be greater than the model order.")
  }

  
  ####################################################
  ### model order
  ####################################################
  if (mtdorder%%1 !=0) {
    stop("error: the model order must be an integer.")
  }

  if (mtdorder < 2) {
    stop("error: the model order must be greater than two.")
  }

  ####################################################
  ### family
  ####################################################
  if (missing(family)) {
    stop("error: family must be specified")
  }
  
  family <- tolower(family)
  family_names <- c("lomax")
  
  if (!family %in% family_names) {
    stop("error: specified family '", family, "' is not a valid option;
         available families are ", paste(family_names, collapse = ", ", sep = "") ,".")
  }  
  
  ####################################################
  ### priors for the weights
  ####################################################  
  if (missing(weight)) {
    stop("error: weight must be specified")
  }
  
  weight <- tolower(weight)
  weight_names <- c("dir", "sb", "cdp")
  
  if (!weight %in% weight_names) {
    stop("error: specified prior for the weight '", weight, "' is not a valid option;
         available priors are ", 
         paste(weight_names, collapse = ", ", sep = "") ,".")
  }
  
  if (weight == "dir") {
    
    g <- updateDIRweight
    dir_shape <- prior$dir_shape
    
    if (is.null(dir_shape)) {
      stop("error: dir_shape for the Dirichlet prior must be specified.")
    }
    
    if (any(dir_shape <= 0)) {
      stop("error: dir_shap must be a positive vector of length mtdorder")
    }
    
    weight_param <- list(dir_shape = dir_shape)
    
  } 
  else if (weight == "sb") {
    
    g <- updateSBweight
    alpha <- prior$alpha
    
    if (is.null(alpha)) {
      stop("error: alpha for the truncated stick-breaking prior must be specified.")
    }
    
    if (alpha <= 0) stop("error: alpha must be positive.")
    
    weight_param <- list(alpha = alpha, mtdorder = mtdorder)
    
  } 
  else if (weight == "cdp") {
    
    g <- updateCDPweight
    alpha_0 <- prior$alpha_0
    a_G0 <- prior$a_G0
    b_G0 <- prior$b_G0
    
    if (is.null(alpha_0)) {
      stop("error: alpha_0for the cdf-based prior must be specified.")
    }
    
    if (is.null(a_G0) || is.null(b_G0)) {
      stop("error: base distribution parameter a_G0 and b_G0 for the cdf-based prior must be specified.")
    }
    
    if (alpha_0 <= 0) stop("error: alpha_0 must be positive.")
    if (a_G0 <= 0) stop("error: a_G0 must be positive")
    if (b_G0 <= 0) stop("error: b_G0 must be positive")
    
    weight_param <- list(alpha_0 = alpha_0, a_G0 = a_G0, b_G0 = b_G0, 
                         cutoff = seq(0, 1, length = mtdorder + 1))
    
  } 
  
  ####################################################
  ### mcmc
  ####################################################   
  niter <- mcmc_param$niter
  nburn <- mcmc_param$nburn
  nthin <- mcmc_param$nthin
  
  if (family == "lomax") {
    
    ### check prior
    if (is.null(prior$u_alpha)) stop("error: hyperparameter u_alpha must be specified.")
    if (is.null(prior$v_alpha)) stop("error: hyperparameter v_alpha must be specified.")
    if (is.null(prior$u_phi)) stop("error: hyperparameter u_phi must be specified.")
    if (is.null(prior$v_phi)) stop("error: hyperparameter v_phi must be specified.")
    
    if (prior$u_alpha <= 0) stop("error: hyperparameter u_alpha must be positive.")
    if (prior$v_alpha <= 0) stop("error: hyperparameter v_alpha must be positive.")
    if (prior$u_phi <= 0) stop("error: hyperparameter u_phi must be positive,")
    if (prior$v_phi <= 0) stop("error: hyperparameter v_phi must be positive.")
    
    ### check tuning parameter
    if (is.null(tuning$se_phi)) stop("error: tuning paramter se_phi must be specified.")
    if (is.null(tuning$se_bb)) stop("error: tuning parameter se_bb must be specified.")
    
    if (tuning$se_phi <= 0) stop("error: tuning parameter se_phi must be positive.")
    if (any(tuning$se_bb <= 0)) stop("error: tuning parameter se_bb must be a positive vector.")
    
    ### check starting values
    if (is.null(starting$phi)) stop("error: starting value of phi must be specified.")
    if (is.null(starting$bb)) stop("error: starting value of regression paramter bb must be specified.")
    
    if (starting$phi <= 0) stop("error: starting value of phi must be positive")
    if (any(starting$bb < 0)) stop("error: starting value of regression parameter bb must be a nonnegative vector.")

    gibbsRegLMTD(yy, XX, mtdorder, weight_param, prior, tuning, starting, g, niter, nburn, nthin)

  } 
  
}
