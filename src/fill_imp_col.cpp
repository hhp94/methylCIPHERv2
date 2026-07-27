#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void fill_imp_col(NumericMatrix obj, NumericVector mean_vec) {
  if (mean_vec.size() != obj.ncol())
  {
    stop("mean_vec must have length ncol(obj)");
  }

  for (int j = 0; j < obj.ncol(); ++j) {
    for (int i = 0; i < obj.nrow(); ++i) {
      if (!R_finite(obj(i, j)))
        obj(i, j) = mean_vec[j];
    }
  }
}
