#' Simulation from an MTD model
#' 
#' This function generates data from an MTD model with a given family of stationary marginal distribution.
#' 
#' @param mtdorder the model order. Note mtdorder > 1.
#' @param weight a numeric vector of weights of the MTD model.
#' @param family a quoted keyword that specifies the stationary marginal distribution.
#'        Supported keywords are: "gaussian", "poisson, negative binomial and lomax".
#' @param param a list of pre-specified parameter values of the model.
#' @param trun a logical value. If TRUE, the first mtdorder elements need to be specified.
#' @param init_data a numeric vector of the initial values of the data to be generated. 
#'                  The length of the vector is mtdorder.
#' @param size the number of data to be generated conditional on the initial values.
#' 
#' @return 
#' The return object is a list comprising:
#' \describe{
#'   \item{y}{a numeric vector of simulated data.}
#'   \item{label}{a numeric vector of configuration variables that identify which component of the mtd was used to generate the data.}
#'   \item{lags}{a numeric vector that contains the specific lags were used to generate the data.}
#'   \item{u}{if there is, a numeric vector of the auxiliary variables used to generate the data.}
#' }
#' 
#' @references 
#' Zheng, X., Kottas, A., and Sansó, B. (2021),
#' "On construction and estimation of stationary mixture transition distribution models,"
#' \emph{arXiv:2010.12696}. 
#' 
#' @export
rMTD <- function(mtdorder, weight, family, param, size, trun = FALSE, init_data = NULL) {
  
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
  ### weights
  ####################################################
  if (length(weight) < 2) {
    stop("error: the length of the weight must be at least two.")
  }
  
  if (mtdorder != length(weight)) {
    stop("error: the length of the weight must be the same as the mtdorder.")
  }

  if (!all.equal(sum(weight), 1)) stop("error: sum of the weight must be 1.")

  if (any(weight <= 0, weight >= 1)) {
    stop("error: weight must be a positive vector and each element must be between 0 and 1.")
  }

  ####################################################
  ### family
  ####################################################
  if (missing(family)) {
    stop("error: family must be specified")
  }
  
  family <- tolower(family)
  family_names <- c("gaussian", "poisson", "negative binomial", "lomax")
  
  if (!family %in% family_names) {
    stop("error: specified family '", family, "' is not a valid option;
         available families are ", paste(family_names, collapse = ", ", sep = "") ,".")
  }  
  
  ####################################################
  ### size
  ####################################################
  if (trun && size <= mtdorder) {
    stop("error: the data size must be greater than the model order ")
  }
  
  ####################################################
  ### init_data
  ####################################################
  if (trun && is.null(init_data)) {
    stop("error: initial values must be specified.")
  }

  if (trun && length(init_data) != mtdorder) {
    stop("error: the length of init_data must be the same as mtdorder.")
  }
  
  ####################################################
  ### family and parameters
  ####################################################
  if (family == "gaussian") {
    
    rho <- param$rho
    mu <- param$mu
    sigma2 <- param$sigma2
    
    if (is.null(rho)) stop("error: rho must be specified.")
    if (is.null(mu)) stop("error: mu must be specified.")
    if (is.null(sigma2)) stop("error: sigma2 must be specified.")

    if (any(rho <= 0, rho >= 1)) {
      stop("error: rho must be a positive vector and each element must be between 0 and 1.")
    }
    if (mu < 0) stop("error: mu must be non-negative.")
    if (sigma2 <= 0) stop("error: sigma2 must be positive.")
    
    rGMTD(mtdorder, weight, rho, mu, sigma2, size, trun, init_data)
    
  } 
  else if (family == "poisson") {
    
    la <- param$la
    ga <- param$ga

    if (is.null(la)) stop("error: la must be specified.")
    if (is.null(ga)) stop("error: ga must be specified.")

    if (la <= 0) stop("error: la must be positive.")
    if (ga <= 0) stop("error: ga must be positive.")
    
    rPMTD(mtdorder, weight, la, ga, size, trun, init_data)
    
  } 
  else if (family == "negative binomial") {
    
    la <- param$la
    ga <- param$ga
    kap <- param$kap
    eta <- param$eta
    
    if (is.null(la)) stop("error: la must be specified.")
    if (is.null(ga)) stop("error: ga must be specified.")
    if (is.null(kap)) stop("error: kap must be specified.")
    if (is.null(eta)) stop("error: eta must be specified.")

    if (la <= 0) stop("error: la must be positive.")
    if (ga <= 0) stop("error: ga must be positive.")
    if (kap <= 0) stop("error: kap must be positive.")
    if (eta <= 0) stop("error: eta must be positive.")
    
    rNBMTD(mtdorder, weight, la, ga, kap, eta, size, trun, init_data)
    
  } 
  else if (family == "lomax") {
    
    alpha <- param$alpha
    phi <- param$phi
    
    if (is.null(alpha)) stop("error: alpha must be specified.")
    if (is.null(phi)) stop("error: phi must be specified.")

    if (alpha <= 1) stop("error: alpha must be greater than 1.")
    if (phi <= 0) stop("error: phi must be positive.")
    
    rLMTD(mtdorder, weight, alpha, phi, size, init_data)
    
  } 
  
}

rGMTD <- function(mtdorder, weight, rho, mu, sigma2, size, trun, init_data){
  
  y <- array(NA, dim = size)
  
  if (trun) {
    
    y[1:mtdorder] <- init_data
    label <- sample(1:mtdorder, size = size - mtdorder, replace = TRUE, prob = weight)
    lags <- array(NA, length(label))
    
    for (i in (mtdorder + 1):size) {
      
      ilabel <- label[i - mtdorder]
      ilag <- y[i - ilabel]
      lags[i - mtdorder] <- ilag
      y[i] <- rnorm(1, mu + rho[ilabel] * (ilag - mu), sqrt(sigma2 * (1 - rho[ilabel]^2)))
      
    }

  } else {
    
    label <- array(NA, dim = size)
    lags <- array(NA, dim = size)
    y[1] <- rnorm(1, mu, sqrt(sigma2))  
    lags[2] <- y[1]
    y[2] <- rnorm(1, mu + rho[1] * (lags[2] - mu), sqrt(sigma2 * (1 - rho[1]^2)))
    
    if (mtdorder > 2) {
      
      for (i in 3:mtdorder) {
        
        iweight <- weight[1:(i-1)]
        ilabel <- sample(1:(i-1), size = 1, prob = iweight)
        label[i] <- ilabel
        ilag <- y[i - ilabel]
        lags[i] <- ilag
        y[i] <- rnorm(1, mu + rho[ilabel] * (ilag - mu), sqrt(sigma2 * (1 - rho[ilabel]^2)))
        
      }
      
    }
    
    label[(mtdorder + 1):size] <- sample(1:mtdorder, size = size - mtdorder, replace = TRUE, prob = weight)

    for (i in (mtdorder + 1):size){
      
      ilabel <- label[i]
      ilag <- y[i - ilabel]
      lags[i] <- ilag
      y[i] <- rnorm(1, mu + rho[ilabel] * (ilag - mu), sqrt(sigma2 * (1 - rho[ilabel]^2)))
      
    }

  }
  
  list(y = y, label = label, lags = lags)

}

rPMTD <- function(mtdorder, weight, la, ga, size, trun, init_data){

  y <- array(NA, dim = size)
  p <- ga / (la + ga)
  
  if (trun) {

    y[1:mtdorder] <- init_data
    label <- sample(1:mtdorder, size = size - mtdorder, replace = TRUE, prob = weight)
    lags <- array(NA, length(label))
    u <- array(NA, length(label))

    for (i in (mtdorder + 1):size){
      
      ilabel <- label[i - mtdorder]
      ilag <- y[i - ilabel]
      lags[i - mtdorder] <- ilag
      ut <- rpois(1, la)
      u[i - mtdorder] <- ut
      y[i] <- ut + rbinom(1, ilag, p)
      
    }

  } else {

    label <- array(NA, dim = size)
    lags <- array(NA, dim = size)
    u <- array(NA, dim = size)
    y[1] <- rpois(1, la + ga)
    lags[2] <- y[1]
    ut <- rpois(1, la)
    u[2] <- ut
    y[2] <- ut + rbinom(1, lags[2], p)
    
    if (mtdorder > 2) {
      
      for (i in 3:mtdorder) {
        
        iweight <- weight[1:(i-1)]
        ilabel <- sample(1:(i-1), size = 1, prob = iweight)
        label[i] <- ilabel
        ilag <- y[i - ilabel]
        lags[i] <- ilag
        ut <- rpois(1, la)
        u[i] <- ut
        y[i] <- ut + rbinom(1, ilag, p)
        
      }
      
    }
    
    label[(mtdorder + 1):size] <- sample(1:mtdorder, size = size - mtdorder, replace = TRUE, prob = weight)

    for (i in (mtdorder + 1):size){
      
      ilabel <- label[i]
      ilag <- y[i - ilabel]
      lags[i] <- ilag
      ut <- rpois(1, la)
      u[i] <- ut
      y[i] <- ut + rbinom(1, ilag, p)
      
    }

  }
  
  list(y = y, label = label, lags = lags, u = u)

}

rNBMTD <- function(mtdorder, weight, la, ga, kap, eta, size, trun, init_data){

  y <- array(NA, dim = size)
  th <- ga / (la + ga)
  psi <- 1 - la / (2 * la + ga + eta)

  if (trun) {

    y[1:mtdorder] <- init_data
    label <- sample(1:mtdorder, size = size - mtdorder, replace = TRUE, prob = weight)
    lags <- array(NA, length(label))
    u <- array(NA, length(label))
    
    for (i in (mtdorder + 1):(size)){
      
      ilabel <- label[i - mtdorder]
      ilag <- y[i - ilabel]
      lags[i - mtdorder] <- ilag
      ut <- rnbinom(1, kap + ilag, psi)
      u[i - mtdorder] <- ut
      y[i] <- ut + rbinom(1, ilag, th)
      
    }

  } else {

    label <- array(NA, dim = size)
    lags <- array(NA, dim = size)
    u <- array(NA, dim = size)
    y[1] <- rnbinom(1, kap, eta / (la + ga + eta))
    lags[2] <- y[1]
    ut <- rnbinom(1, kap + lags[2], psi)
    u[2] <- ut
    y[2] <- ut + rbinom(1, lags[2], th)
    
    if (mtdorder > 2) {
      
      for (i in 3:mtdorder) {
        
        iweight <- weight[1:(i-1)]
        ilabel <- sample(1:(i-1), size = 1, prob = iweight)
        label[i] <- ilabel
        ilag <- y[i - ilabel]
        lags[i] <- ilag
        ut <- rnbinom(1, kap + ilag, psi)
        u[i] <- ut
        y[i] <- ut + rbinom(1, ilag, th)
        
      }
      
    }
    
    label[(mtdorder + 1):size] <- sample(1:mtdorder, size = size - mtdorder, replace = TRUE, prob = weight)

    for (i in (mtdorder + 1):size){
      
      ilabel <- label[i]
      ilag <- y[i - ilabel]
      lags[i] <- ilag
      ut <- rnbinom(1, kap + ilag, psi)
      u[i] <- ut
      y[i] <- ut + rbinom(1, ilag, th)
      
    }
    
  }

  list(y = y, label = label, lags = lags, u = u)
  
}

 
rLMTD <- function(mtdorder, weight, alpha, phi, size, init_data){

  y <- array(NA, dim = size)

  if (trun) {

    y[1:mtdorder] <- init_data
    label <- sample(1:mtdorder, size = size - mtdorder, replace = TRUE, prob = weight)
    lags <- array(NA, length(label))
    
    for (i in (mtdorder + 1):size) {
      
      ilabel <- label[i - mtdorder]
      ilag <- y[i - ilabel]
      lags[i - mtdorder] <- ilag
      y[i] <- extraDistr::rlomax(1, 1 / (phi + ilag), alpha)
      
    }

  } else {

    label <- array(NA, dim = size)
    lags <- array(NA, dim = size)
    y[1] <- extraDistr::rlomax(1, 1 / phi, alpha - 1)
    lags[2] <- y[1]
    y[2] <- extraDistr::rlomax(1, 1 / (phi + lags[2]), alpha)
    
    if (mtdorder > 2) {
      
      for (i in 3:mtdorder) {
        
        iweight <- weight[1:(i-1)]
        ilabel <- sample(1:(i-1), size = 1, prob = iweight)
        label[i] <- ilabel
        ilag <- y[i - ilabel]
        lags[i] <- ilag
        y[i] <- extraDistr::rlomax(1, 1 / (phi + ilag), alpha)
        
      }
      
    }
    
    label[(mtdorder + 1):size] <- sample(1:mtdorder, size = size - mtdorder, replace = TRUE, prob = weight)

    for (i in (mtdorder + 1):size){
      
      ilabel <- label[i]
      ilag <- y[i - ilabel]
      lags[i] <- ilag
      y[i] <- extraDistr::rlomax(1, 1 / (phi + ilag), alpha)
      
    }

  }

  list(y = y, label = label, lags = lags)
  
}