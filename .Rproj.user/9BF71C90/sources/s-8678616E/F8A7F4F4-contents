#' One-step ahead prediction
#'
#' This function generates one-step-ahead posterior predictive distribution 
#' given a fitted MTD model.
#
#' @param res a list of MCMC results from the MTD model; an output from functions 
#'        either \code{\link{tsMTD}} or \code{\link{tsRegMTD}}.
#' @param family a quoted keyword that specifies the stationary marginal distribution.
#'        Supported keywords are: "gaussian", "poisson", "negative binomial".
#' @param lags a numeric matrix of the lags of the observation. Each row 
#'        corresponds to an observation. The dimension of the matrix is 
#'        truncated data size x mtdorder.
#' @param probs a numeric vector of probabilities with values in \eqn{[0, 1]}.
#'
#' @references 
#' Zheng, X., Kottas, A., and Sansó, B. (2021),
#' "On construction and estimation of stationary mixture transition distribution models,"
#' \emph{arXiv:2010.12696}.
#' 
#' @export
#'
predMTD <- function(res, family, lags, probs){
  
  if (missing(family)) {
    
    stop("error: family must be specified")
    
  }
  
  family <- tolower(family)
  family_names <- c("poisson", "negative binomial")
  
  if (!family %in% family_names) {
    
    stop("error: specified family '", family, "' is not a valid option;
         available families are ", 
         paste(family_names, collapse = ", ", sep = "") ,".")
    
  }
  
  if (any(probs <= 0, probs >= 1)) {
    
    stop("error: probs must be a positive vector and each element must be between 0 and 1.")
    
  }
  
  if (family == "poisson") {
    
    weight <- res$weight
    la <- res$la
    th <- res$th
    
    if (nrow(weight) != ncol(lags)) {
      
      stop("error: the model order used to fit the MTD is different from the number of lags of each observaton.")
      
    }
    
    if (is.null(la) || is.null(th) || is.null(weight)) {
      
      stop("error: lack of posterior samples.")
      
    }
    
    mtdorder <- nrow(weight)
    
    out <- predPMTD(mtdorder, weight, lags, la, th, probs) 
      
  }
  else if (family == "negative binomial") {
    
    weight <- res$weight
    th <- res$th
    psi <- res$psi
    kap <- res$kap
    
    if (nrow(weight) != ncol(lags)) {
      
      stop("error: the model order used to fit the MTD is different from the number of lags of each observaton.")
      
    }
    
    
    if (is.null(th) || is.null(psi) || is.null(kap) || is.null(weight)) {
      
      stop("error: lack of posterior samples.")
      
    }
    
    mtdorder <- nrow(weight)
    
    out <- predNBMTD(mtdorder, weight, lags, th, psi, kap, probs) 
    
  }
  
  out
  
}