#' Lag matrix of a time series
#' 
#' This function computes a given number of lags for each observation of a time series.
# 
#' @param obs a numeric vector of observations
#' @param k an integer indicating the number of lags. 
#' 
#' @return 
#' The return object is a matrix. Each row corresponds to k lags of the observation.
#' The dimension of the matrix is truncated data size x k.
#' 
#' @export
#'
lagmat <- function(obs, k){
  if (k == 1){
    lag_mat <- rev(rev(obs)[-1])
  } else {
    grid <- seq(1, k - 1, by = 1)
    lag_dat <- sapply(grid, function(x) rev(rev(obs)[-(1:x)])[-(1:(k - x))])
    lag_mat <- cbind(lag_dat, rev(rev(obs)[-(1:k)]))
  }
  lag_mat
}

#' Obtain posterior estimates
#' 
#' @param x an object for which "summarize" is desired.
#' @param point_est a quoted keyword that specifies the point estimate.
#'                  Supported keywords are "mean", "median", "both".
#' @param lower_prob Lower probability.
#' @param upper_prob Upper probability.
#' @param k an integer indicating the number of decimal places to be used.
#' 
#' @export
#' 
summarize <- function(x, point_est, lower_prob = 0.025, upper_prob = 0.975, k = 2) {
  UseMethod("summarize", x)
}

#' Obtain posterior estimates
#' 
#' @param x a numeric vector.
#' @param point_est a quoted keyword that specifies the point estimate.
#'                  Supported keywords are "mean", "median", "both".
#' @param lower_prob Lower probability.
#' @param upper_prob Upper probability.
#' @param k an integer indicating the number of decimal places to be used.
#' 
#' @export
#'  
summarize.numeric <- function(x, point_est = "mean", lower_prob = 0.025, upper_prob = 0.975, k = 2){

  point_est <- tolower(point_est)
  point_est_names <- c("mean", "median", "both")
  
  if (!point_est %in% point_est_names) {
    stop("error: specified type of point estimate '", point_est, "' is not a valid option;
         available types are ", paste(point_est_names, collapse = ", ", sep = "") ,".")
  }

  if (lower_prob <=0 || lower_prob >= 1) stop("error: lower_prob must be between 0 and 1.")
  if (upper_prob <=0 || upper_prob >= 1) stop("error: upper_prob must be between 0 and 1.")
  if (lower_prob >= upper_prob) stop("error: lower_prob must be smaller than upper_prob.")

  if (point_est == "mean"){
    pe <- specRound(mean(x), k)
  } 
  else if (point_est == "median"){
    pe <- specRound(median(x), k)
  } 
  else if (point_est == "both")
    pe <- c(specRound(mean(x), k), specRound(median(x), k))

  ci <- specCI(x, lower_prob, upper_prob, k)

  paste(pe, ci)

}

#' Obtain posterior estimates
#' 
#' @param x a matrix. Each row corresponds to a parameter.
#' @param point_est a quoted keyword that specifies the point estimate.
#'                  Supported keywords are "mean", "median", "both".
#' @param lower_prob Lower probability.
#' @param upper_prob Upper probability.
#' @param k an integer indicating the number of decimal places to be used.
#' 
#' @export
#' 
summarize.matrix <- function(x, point_est = "mean", lower_prob = 0.025, upper_prob = 0.975, k = 2){

  point_est <- tolower(point_est)
  point_est_names <- c("mean", "median", "both")
  
  if (!point_est %in% point_est_names) {
    stop("error: specified type of point estimate '", point_est, "' is not a valid option;
         available types are ", paste(point_est_names, collapse = ", ", sep = "") ,".")
  }

  if (lower_prob <=0 || lower_prob >= 1) stop("error: lower_prob must be between 0 and 1.")
  if (upper_prob <=0 || upper_prob >= 1) stop("error: upper_prob must be between 0 and 1.")
  if (lower_prob >= upper_prob) stop("error: lower_prob must be smaller than upper_prob.")

  if (point_est == "mean"){
    pe <- specRound(rowMeans(x), k)
  } 
  else if (point_est == "median"){
    pe <- specRound(apply(x, 1, median), k)
  } 
  else if (point_est == "both"){
    pe <- cbind(specRound(rowMeans(x), k), specRound(apply(x, 1, median), k))
  } 

  ci <- apply(x, 1, specCI, lower_prob, upper_prob, k)
  
  paste(pe, ci)
}

# round numbers
specRound <- function(x, k = 3) trimws(format(round(x, k), nsmall = k))

# obtain credible intervals
specCI <- function(x, lower_prob = 0.025, upper_prob = 0.975, k = 2){
  q <- quantile(x, probs = c(lower_prob, upper_prob))
  ci <- paste0('(', specRound(q[1],k), ', ', specRound(q[2],k), ')')
  ci
}