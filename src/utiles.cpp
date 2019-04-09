#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//' @importFrom Rcpp sourceCpp
//' @useDynLib fda2d
//' @export
// [[Rcpp::export]]
arma::mat c_hat(arma::mat E0,arma::mat F0,arma::colvec cy){
  int m = E0.n_rows;
  int n = F0.n_rows;
  int k;
  int p = E0.n_cols;
  int q = F0.n_cols;
  arma::mat Phi(m*n,p*q);
  for (int j =0; j < n; j++ )
  {
    for (int i =0; i < m; i++ )
    {
      k = j*m+i;
      Phi.row(k) =arma::trans(arma::vectorise(arma::trans(E0.row(i))*F0.row(j)));
        //arma::trans(cy(arma::span(0,p*q-1)));
    }
  }
  arma::colvec coef = arma::solve(Phi, cy);
  return coef;
}
