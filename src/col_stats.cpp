#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix col_stats(NumericMatrix obj) {
  NumericMatrix result(2, obj.ncol());

  for (int j = 0; j < obj.ncol(); ++j) {
    double sum = 0.0;
    double n_valid = 0.0;
    for (int i = 0; i < obj.nrow(); ++i) {
      double x = obj(i, j);
      if (R_finite(x)) {
        sum += x;
        n_valid += 1.0;
      }
    }
    result(0, j) = sum;
    result(1, j) = n_valid;
  }
  return result;
}
