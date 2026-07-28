#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// In-place cohort-mean fill, one column pointer at a time.
//
// Written as a select rather than `if (missing) assign`: the test per element
// is unavoidable, but the BRANCH is not, and a branch is what stops the loop
// vectorizing. Every lane does the same work and the comparison becomes a
// blend. Measured 1.6 -> 1.2 ms on a 100 x 39,025 slice at -O2; against the
// original sugar-indexed version, 5.5 -> 1.2 ms.
//
// `std::isnan` and not `!std::isfinite`: this fills MISSING values, and an
// Inf is not one. It cannot arrive here anyway -- `check_col_values()` aborts
// on the first Inf in the panel before any cache is built -- but if that ever
// changed, leaving the Inf in place carries it visibly into the score, where
// quietly swapping it for a cohort mean would not.

// [[Rcpp::export]]
void fill_imp_col(NumericMatrix obj, NumericVector mean_vec)
{
  if (mean_vec.size() != obj.ncol())
  {
    stop("mean_vec must have length ncol(obj)");
  }

  const R_xlen_t nc = obj.ncol(), nr = obj.nrow();
  double *x = obj.begin();
  const double *mu = mean_vec.begin();

  for (R_xlen_t j = 0; j < nc; ++j)
  {
    double *col = x + j * nr;
    const double m = mu[j];
    for (R_xlen_t i = 0; i < nr; ++i)
    {
      const double v = col[i];
      col[i] = std::isnan(v) ? m : v;
    }
  }
}
