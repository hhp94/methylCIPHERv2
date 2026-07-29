#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// fill every non-finite entry (NA, NaN, +/-Inf) with the column mean.
// same predicate as col_stats(), which computed that mean over the same set
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
      col[i] = std::isfinite(v) ? v : m;
    }
  }
}
