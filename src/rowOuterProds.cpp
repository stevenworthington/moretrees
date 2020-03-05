#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rowOuterProds(NumericMatrix X) {
  int n = X.rows();
  int K = X.cols();
  int J = K * (K + 1) / 2;
  NumericMatrix out(n , J);
  
  int col = 0;
  
  for(int i = 0; i < n; i++) {
    col = 0;
    for(int j = 0; j < K; j++) {
      for(int k = j; k < K; k++) {
        out(i, col) = X(i, j) * X(i, k);
        col = col + 1;
      }
    }
  }
  
  return out;
}