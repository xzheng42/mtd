#' Randomize quantile residual analysis for MTD models
#' 
#' This function performs model checking for MTD models based on the computed randomized quantile residuals.
#
#' @param res a list of MCMC results from the MTD model; an output from functions either \code{\link{tsMTD}} or \code{\link{tsRegMTD}}.
#' @param obs a numeric vector of complete data (not truncated).
#' @param family a quoted keyword that specifies the stationary marginal distribution.
#'        Supported keywords are: "gaussian", "poisson" and "negative binomial".
#' @param covar covariates. The default is NULL. The argument is specified when a regression MTD model is desired.
#'
#' @return 
#' The return object is a matrix. Each column corresponds to posterior samples of an iteration.
#' The dimension of the matrix is truncated data size x posterior sample size.
#' 
#' @references 
#' Dunn, P. K. and Smyth, G. K. (1996), 
#' “Randomized quantile residuals,”
#' \emph{Journal of Computational and Graphical Statistics}.
#' 
#' Zheng, X., Kottas, A., and Sansó, B. (2021),
#' "On construction and estimation of stationary mixture transition distribution models,"
#' \emph{arXiv:2010.12696}.
#'
#' @export
#'
rqrMTD <- function(res, obs, family, covar = NULL) {
  
  if (missing(family)) {
    
    stop("error: family must be specified")
    
  }
  
  family <- tolower(family)
  family_names <- c("poisson", "lomax")
  
  if (!family %in% family_names) {
    
    stop("error: specified family '", family, "' is not a valid option;
         available families are ", 
         paste(family_names, collapse = ", ", sep = "") ,".")
    
  }
  
  if (family == "poisson") {
    
    weight <- res$weight
    la <- res$la
    th <- res$th
    
    if (is.null(la) || is.null(th) || is.null(weight)) {
      
      stop("error: lack of posterior samples.")
      
    }
    
    rqr <- rqrPMTD(obs, weight, la, la / (1 - th) - la)
    
    if (!is.null(covar)) {
      
      stop("error: we don't support Poisson MTD model with covariates now.")
      
    }
    
  } 
  else if (family == "lomax") {
    
    if (is.null(covar)) {
      
      stop("error: need to include covariates.")
      
    } else {
      
      weight <- res$weight
      bb <- res$bb
      alpha <- res$alpha
      phi <- res$phi
      
      if (is.null(bb) || is.null(alpha) || is.null(phi) || is.null(weight)) {
        
        stop("error: lack of posterior samples.")
        
      }
      
      rqr <- rqrRegLMTD(obs, covar, weight, bb, alpha, phi)
      
    }
  } 

  rqr

}