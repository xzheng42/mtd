#include <iostream>
#include <RcppArmadillo.h>
#include <math.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double dCondBiPois(const int& uu,
                   const int& vv,
                   const double& th1,
                   const double& th2,
                   const double& th12) {
  
  int zz = min(uu, vv);                    
  double logp1 = log(th12) - log(th2 + th12);
  double logp2 = log(th2) - log(th2 + th12);
  double cc1, cc2, dens;
  arma::Col<int> grid = arma::linspace<arma::Col<int>>(0, zz, zz + 1); 
  arma::colvec yy(zz + 1);
  
  for (int j = 0; j < zz + 1; ++j) {
    cc1 = lgamma(vv + 1) - lgamma(j + 1) - lgamma(vv - j + 1);
    cc2 = j * logp1 + (vv - j) * logp2 + (uu - j) * log(th1) - lgamma(uu - j + 1);
    yy(j) = exp(cc1 + cc2);
  }
  
  dens = exp(-th1) * sum(yy);
  
  return dens;
  
}

// [[Rcpp::export]]
double pCondBiPois(const int& uu,
                   const int& vv,
                   const double& th1,
                   const double& th2,
                   const double& th12) {
  
  arma::colvec yy(uu + 1);  
  
  for (int i = 0; i < uu + 1; ++i) {
    yy(i) = dCondBiPois(i, vv, th1, th2, th12);
  }
  
  return sum(yy);
  
}

// [[Rcpp::export]]
double pPMTD(const int& obs,
             const arma::Col<int>& lags,
             const arma::colvec& weight,
             const double& la,
             const double& ga) {
  
  int mtdorder = lags.n_rows;
  arma::colvec yy(mtdorder);
  
  for (int j = 0; j < mtdorder; ++j) {
    yy(j) = weight(j) * pCondBiPois(obs, lags(j), la, la, ga);
  }
  
  return sum(yy);
  
}

// [[Rcpp::export]]
arma::mat rqrPMTD(const arma::Col<int>& obs,
                  const arma::mat& weight,
                  const arma::colvec& la,
                  const arma::colvec& ga) {
  
  int mtdorder = weight.n_rows;
  int trun_size = obs.n_rows - mtdorder;
  int sample_size = la.n_rows; 
  int iobs;
  double iter_la, iter_ga, at, bt, pt;
  arma::colvec iter_weight;
  arma::Col<int> ilags_rev, ilags;
  arma::mat rt(trun_size, sample_size);
  
  for (int iter = 0; iter < sample_size; ++iter) {
    
    iter_la = la(iter);
    iter_ga = ga(iter);
    iter_weight = weight.col(iter);
    
    for (int i = 0; i < trun_size; ++i) {
      
      iobs = obs(i + mtdorder);
      ilags_rev = obs.rows(i, i + mtdorder - 1);
      ilags = arma::reverse(ilags_rev);
      bt = pPMTD(iobs, ilags, iter_weight, iter_la, iter_ga);
      at = pPMTD(iobs - 1, ilags, iter_weight, iter_la, iter_ga);
      pt = R::runif(at, bt);
      rt(i, iter) = R::qnorm(pt, 0.0, 1.0, true, false);
      
    }
    
  }
  
  return rt;
  
}

// [[Rcpp::export]]
double pLMTD(const double& obs,
             const arma::colvec& lags,
             const arma::colvec& weight,
             const double& alpha,
             const double& phi) {
  
  int mtdorder = lags.n_rows;
  arma::colvec yy(mtdorder);
  
  for (int j = 0; j < mtdorder; ++j) {
    yy(j) = weight(j) * (1.0 - pow((1.0 + obs / (lags(j) + phi)), -alpha));
  }
  
  return sum(yy);
  
}

// [[Rcpp::export]]
arma::mat rqrRegLMTD(const arma::colvec& obs,
                     const arma::mat& xx,
                     const arma::mat& weight,
                     const arma::mat& bb, 
                     const arma::colvec& alpha,
                     const arma::colvec& phi) {
  
  int mtdorder = weight.n_rows;
  int trun_size = obs.n_rows - mtdorder;
  int sample_size = alpha.n_rows;
  arma::rowvec ixx;
  arma::colvec eps, ilags_rev, ilags, ipp, iter_bb, iter_weight;
  double iEps, iter_alpha, iter_phi, pt;
  arma::mat rt(trun_size, sample_size);
  
  for (int iter = 0; iter < sample_size; ++iter) {
    
    iter_bb = bb.col(iter);
    iter_alpha = alpha(iter);
    iter_phi = phi(iter);
    iter_weight = weight.col(iter);
    eps = obs % arma::exp(-xx * iter_bb);
    
    for (int i = 0; i < trun_size; ++i) {
      iEps = eps(i + mtdorder);
      ilags_rev = eps.rows(i, i + mtdorder - 1);
      ilags = arma::reverse(ilags_rev);
      pt = pLMTD(iEps, ilags, iter_weight, iter_alpha, iter_phi);
      rt(i, iter) = R::qnorm(pt, 0.0, 1.0, true, false);
    }
    
  }
  
  return rt;
  
}
