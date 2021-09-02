#include <iostream>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// update auxiliary variables of the Poisson MTD model
// [[Rcpp::export]]
arma::colvec updatePMTDut(const arma::colvec& ut,
                          const arma::colvec& dat,
                          const arma::colvec& dep_lag,
                          const double& la,
                          const double& ga) {
  
  int trun_size = dat.n_rows;
  double th = ga / (la + ga);
  double ut_lb, val_len, prop_ut, cur_ut, idat, ilag_dep, prop_loglik, cur_loglik, diff_loglik;
  arma::Col<int> prop_ut_vals;
  arma::colvec idistr;
  arma::colvec new_ut(trun_size);

  for (int i = 0; i < trun_size; ++i) {
    
    idat = dat(i);
    ilag_dep = dep_lag(i);
    cur_ut = ut(i);
    ut_lb = max(0.0, idat - ilag_dep);
    val_len = idat - ut_lb + 1.0;
    idistr.set_size(val_len);
    idistr.fill(1.0 / val_len);
    
    if (ut_lb == idat) {
      prop_ut = idat;
    } else {
      prop_ut_vals = arma::linspace<arma::Col<int>>(ut_lb, idat, val_len);
      prop_ut = as_scalar(Rcpp::RcppArmadillo::sample(prop_ut_vals, 1, 1, idistr));
    }
    
    prop_loglik = R::dbinom(idat - prop_ut, ilag_dep, th, true) + R::dpois(prop_ut, la, true);
    cur_loglik = R::dbinom(idat - cur_ut, ilag_dep, th, true) + R::dpois(cur_ut, la, true);
    diff_loglik = prop_loglik - cur_loglik;
    if (diff_loglik > log(arma::randu())) {
      new_ut(i) = prop_ut;
    } else {
      new_ut(i) = cur_ut;
    }
  }

  return new_ut;
  
}

// update auxiliary variables of the negative binomial MTD model
// [[Rcpp::export]]
arma::colvec updateNBMTDut(const arma::colvec& ut,
                           const arma::colvec& dat,
                           const arma::colvec& dep_lag,
                           const double& th,
                           const double& alp,
                           const double& kap) {
  
  int trun_size = dat.n_rows;
  double ut_lb, val_len, prop_ut, cur_ut, idat, ilag_dep, prop_loglik, cur_loglik, diff_loglik;
  arma::Col<int> prop_ut_vals;
  arma::colvec idistr;
  arma::colvec new_ut(trun_size);
  
  for (int i = 0; i < trun_size; ++i) {
    idat = dat(i);
    ilag_dep = dep_lag(i);
    cur_ut = ut(i);
    ut_lb = max(0.0, idat - ilag_dep);
    if (ut_lb == idat) {
      prop_ut = idat;
    } else {
      val_len = idat - ut_lb + 1.0;
      idistr.set_size(val_len);
      idistr.fill(1.0 / val_len);
      prop_ut_vals = arma::linspace<arma::Col<int>>(ut_lb, idat, val_len);
      prop_ut = as_scalar(Rcpp::RcppArmadillo::sample(prop_ut_vals, 1, 1, idistr));
    }
    prop_loglik = R::dbinom(idat - prop_ut, ilag_dep, th, true) + R::dnbinom(prop_ut, kap + ilag_dep, alp, true);
    cur_loglik = R::dbinom(idat - cur_ut, ilag_dep, th, true) + R::dnbinom(cur_ut, kap + ilag_dep, alp, true);
    diff_loglik = prop_loglik - cur_loglik;
    if (diff_loglik > log(arma::randu())) {
      new_ut(i) = prop_ut;
    } else {
      new_ut(i) = cur_ut;
    }
  }
  
  return new_ut;
  
}

// update labels of the Gaussian MTD model
// [[Rcpp::export]]
List updateGMTDlabel(const arma::colvec& data, 
                     const arma::mat& lags, 
                     const int& mtdorder, 
                     const arma::colvec& weight, 
                     const arma::colvec& mu, 
                     const arma::colvec& sigma2,
                     const arma::colvec& rho) {

  int data_size = data.n_rows;
  int label;
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(1, mtdorder, mtdorder); 
  arma::Col<int> labels(data_size);
  arma::colvec dep_lag(data_size), idata(mtdorder);
  arma::colvec ilag, imu, isigma2, iloglik, idistr;
  
  for (int i = 0; i < data_size; ++i){
    idata.fill(data(i));
    ilag = lags.row(i).t();
    imu = mu + rho % (ilag - mu);
    isigma2 = sigma2 % (1.0 - arma::pow(rho, 2));
    iloglik = log(arma::normpdf(idata, imu, arma::sqrt(isigma2)));
    idistr = log(weight) + iloglik;
    idistr = exp(idistr - idistr.max());
    idistr = idistr / sum(idistr);
    label = as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, idistr));
    labels(i) = label;
    dep_lag(i) = ilag(label - 1);
  }

  return List::create(Named("data_label") = labels, Named("dep_lag") = dep_lag);

}

// update labels of the Poisson MTD model
// [[Rcpp::export]]
List updatePMTDlabel(const arma::colvec& data, 
                     const arma::mat& lags, 
                     const int& mtdorder, 
                     const arma::colvec& weight, 
                     const double& th, 
                     const arma::colvec& ut) {

  int data_size = data.n_rows;
  int label, idata;
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(1, mtdorder, mtdorder); 
  arma::Col<int> labels(data_size);
  arma::colvec dep_lag(data_size), iloglik(mtdorder);
  arma::colvec ilag, idistr;
  
  for (int i = 0; i < data_size; ++i){
    idata = data(i) - ut(i);
    ilag = lags.row(i).t();
    iloglik.fill(NA_REAL);
    for (int j = 0; j < mtdorder; ++j){
      iloglik(j) = R::dbinom(idata, ilag(j), th, true);
    }
    idistr = log(weight) + iloglik;
    idistr = exp(idistr - idistr.max());
    idistr = idistr / sum(idistr);
    label = as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, idistr));
    labels(i) = label;
    dep_lag(i) = ilag(label - 1);
  }

  return List::create(Named("data_label") = labels, Named("dep_lag") = dep_lag);

}

// update labels of the negative binomial MTD model
// [[Rcpp::export]]
List updateNBMTDlabel(const arma::colvec& data, 
                      const arma::mat& lags, 
                      const int& mtdorder, 
                      const arma::colvec& weight, 
                      const double& th, 
                      const arma::colvec& ut,
                      const double& alp,
                      const double& kap) {

  int data_size = data.n_rows;
  int label, idata, iut;
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(1, mtdorder, mtdorder); 
  arma::Col<int> labels(data_size);
  arma::colvec dep_lag(data_size), iloglik(mtdorder);
  arma::colvec ilag, idistr;
  
  for (int i = 0; i < data_size; ++i){
    iut = ut(i);
    idata = data(i) - iut;
    ilag = lags.row(i).t();
    iloglik.fill(NA_REAL);
    for (int j = 0; j < mtdorder; ++j){
      iloglik(j) = R::dbinom(idata, ilag(j), th, true) + R::dnbinom(iut, kap + ilag(j), alp, true);
    }
    idistr = log(weight) + iloglik;
    idistr = exp(idistr - idistr.max());
    idistr = idistr / sum(idistr);
    label = as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, idistr));
    labels(i) = label;
    dep_lag(i) = ilag(label - 1);
  }

  return List::create(Named("data_label") = labels, Named("dep_lag") = dep_lag);

}

// update labels of the Lomax MTD regression model
// [[Rcpp::export]]
List updateLMTDlabel(const arma::colvec& data, 
                     const arma::mat& lags, 
                     const int& mtdorder, 
                     const arma::colvec& weight, 
                     const arma::colvec& alpha, 
                     const arma::colvec& phi) {

  int data_size = data.n_rows;
  int label;
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(1, mtdorder, mtdorder); 
  arma::Col<int> labels(data_size);
  arma::colvec dep_lag(data_size), idata(mtdorder);
  arma::colvec ilag, iloglik, idistr;

  for (int i = 0; i < data_size; ++i){
    idata.fill(data(i));
    ilag = lags.row(i).t();
    iloglik = log(alpha) + alpha % log(phi + ilag) - (alpha + 1.0) % log(phi + ilag + idata);
    idistr = log(weight) + iloglik;
    idistr = exp(idistr - idistr.max());
    idistr = idistr / sum(idistr);
    label = as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, idistr));
    labels(i) = label;
    dep_lag(i) = ilag(label - 1);
  }

  return List::create(Named("data_label") = labels, Named("dep_lag") = dep_lag);

}

// one-step ahead posterior predictive distribution of the Poisson MTD model
// [[Rcpp::export]]
arma::mat predPMTD(const int& mtdorder, 
                   const arma::mat& weight, 
                   const arma::mat& lags,
                   const arma::colvec& la, 
                   const arma::colvec& th,
                   const arma::colvec& probs) {

  int sample_size = weight.n_cols;
  int trun_size = lags.n_rows;
  int probs_size = probs.n_rows;
  int label;
  double iter_la, iter_th, iter_ut, dep_lag;
  arma::mat preds(probs_size, trun_size);
  arma::colvec ipred(sample_size);
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(0, mtdorder - 1, mtdorder);
  arma::colvec ilags, iter_weight;

  for (int i = 0; i < trun_size; ++ i) {
    
    ilags = lags.row(i).t();
    
    for (int iter = 0; iter < sample_size; ++iter) {
      iter_la = la(iter);
      iter_th = th(iter);
      iter_weight = weight.col(iter);
      label = arma::as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, iter_weight));
      dep_lag = ilags(label);
      iter_ut = R::rpois(iter_la);
      ipred(iter) = iter_ut + R::rbinom(dep_lag, iter_th);
    }
    
    preds.col(i) = arma::quantile(ipred, probs);
    
  }

  return preds;

}

// one-step ahead posterior predictive distribution of the negative binomial MTD model
// [[Rcpp::export]]
arma::mat predNBMTD(const int& mtdorder, 
                    const arma::mat& weight, 
                    const arma::mat& lags,
                    const arma::colvec& th, 
                    const arma::colvec& psi,
                    const arma::colvec& kap,
                    const arma::colvec& probs) {
  
  int sample_size = weight.n_cols;
  int trun_size = lags.n_rows;
  int probs_size = probs.n_rows;
  int label;
  double iter_th, iter_psi, iter_kap, iter_ut, dep_lag;
  arma::mat preds(probs_size, trun_size);
  arma::colvec ipred(sample_size);
  arma::Col<int> label_vals = arma::linspace<arma::Col<int>>(0, mtdorder - 1, mtdorder);
  arma::colvec ilags, iter_weight;
  
  for (int i = 0; i < trun_size; ++ i) {
    
    ilags = lags.row(i).t();
    
    for (int iter = 0; iter < sample_size; ++iter) {
      
      iter_th = th(iter);
      iter_psi = psi(iter);
      iter_kap = kap(iter);
      iter_weight = weight.col(iter);
      label = arma::as_scalar(Rcpp::RcppArmadillo::sample(label_vals, 1, 1, iter_weight));
      dep_lag = ilags(label);
      iter_ut = R::rnbinom(iter_kap + dep_lag, iter_psi);
      ipred(iter) = iter_ut + R::rbinom(dep_lag, iter_th);
      
    }
    
    preds.col(i) = arma::quantile(ipred, probs);
    
  }
  
  return preds;
  
}


