#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadilloExtensions/sample.h>
using namespace arma;
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector compute_upiece_cpp(NumericMatrix u, NumericMatrix W, int L) {
  int n = u.nrow();
  int H = W.nrow();
  int K = u.ncol();

  NumericVector u_piece(L * n * K * H);

  for (int i = 0; i < L; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < K; k++) {
        double u_value = u(j, k);
        u_piece[i + L * (j + n * (k + K * H))] = u_value >= i * 1.0 / L ? 1.0 / L : 0.0;
        u_piece[i + L * (j + n * (k + K * H))] += (u_value >= (i - 1) * 1.0 / L) & (u_value < i * 1.0 / L) ? (u_value - (i - 1) * 1.0 / L) : 0.0;
      }
    }
  }

  NumericVector weighted_u_piece(L * n * K * H);
  for (int i = 0; i < L; i++) {
    NumericMatrix u_piece_mat = NumericMatrix::create(u_piece[Range(i, i + L * n * K * H - 1)]);
    NumericMatrix W_mat = NumericMatrix::create(W[Range(i, i + L * n - 1)]);
    NumericMatrix product = u_piece_mat * W_mat;
    NumericVector sum_result = rowSums(product);
    NumericMatrix result_mat = NumericMatrix::create(sum_result);
    NumericMatrix transposed_result = transpose(result_mat);
    NumericMatrix reshaped_result = transposed_result.reshape(H, n * K);
    NumericMatrix permuted_result = reshaped_result.permute(Dimension(0, 2, 1, 3));
    NumericVector final_result = as<NumericVector>(permuted_result);
    weighted_u_piece[Range(i, i + L * n * K * H - 1)] = final_result;
  }

  return weighted_u_piece;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix msf(arma::mat lambda, arma::mat pivot) {
  arma::mat refr = join_rows(lambda, -lambda);
  int k = lambda.n_cols;
  uvec ind = regspace<uvec> (0, k-1);
  uvec perm(k);
  arma::vec signs(k);
  rowvec norms(2*k);
  unsigned int w, c, wc;
  arma::mat diff, diffsq;
  
  for(int i=0; i<k; i++){
    diff = refr.each_col() - pivot.col(i);
    diffsq = square(diff);
    norms = sum(diffsq);
    w = index_min(norms);
    c = refr.n_cols / 2;
    if(w>=c){
      wc = w-c;
      signs(i) = -1;
      perm(i) = ind(wc);
      refr.shed_col(w);
      refr.shed_col(wc); 
      ind.shed_row(wc);} 
    else {
      wc = w+c;
      signs(i) = 1;
      perm(i) = ind(w);
      refr.shed_col(wc); 
      refr.shed_col(w);
      ind.shed_row(w);}
    }
  
  arma::mat permmat = zeros<arma::mat>(k,k);
  for(int i=0; i<k; i++){
    permmat(perm(i), i) = signs(i);
  }
  
  lambda *= permmat;
  return Rcpp::wrap(lambda);
}

// [[Rcpp::export]]
Rcpp::NumericVector msfOUT(arma::mat lambda, arma::mat pivot) {
  arma::mat refr = join_rows(lambda, -lambda);
  int k = lambda.n_cols;
  uvec ind = regspace<uvec> (0, k-1);
  uvec perm(k);
  arma::vec signs(k);
  rowvec norms(2*k);
  unsigned int w, c, wc;
  arma::mat diff, diffsq;
  
  for(int i=0; i<k; i++){
    diff = refr.each_col() - pivot.col(i);
    diffsq = square(diff);
    norms = sum(diffsq);
    w = index_min(norms);
    c = refr.n_cols / 2;
    if(w>=c){
      wc = w-c;
      signs(i) = -1;
      perm(i) = ind(wc);
      refr.shed_col(w);
      refr.shed_col(wc); 
      ind.shed_row(wc);} 
    else {
      wc = w+c;
      signs(i) = 1;
      perm(i) = ind(w);
      refr.shed_col(wc); 
      refr.shed_col(w);
      ind.shed_row(w);}
  }
  
  arma::vec out = (perm + ones<arma::vec>(k)) % signs;
  
  return Rcpp::wrap(out);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix aplr(arma::mat matr, arma::vec perm){
  int k = matr.n_cols;
  arma::mat permmat = zeros<arma::mat>(k,k);
  arma::vec perms = abs(perm) - ones<arma::vec>(k);
  arma::vec signs = sign(perm);
  
  for(int i=0; i<k; i++){
    permmat(perms(i), i) = signs(i);
  }
  
  matr *= permmat;
  return(Rcpp::wrap(matr));
}
