#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// cols is 1-based into obj (NULL = all); stats/overflow_col positions are relative to cols.
// Non-finite values (NA, NaN, +/-Inf) are missing; fill_imp_col() fills the same set.
// [[Rcpp::export]]
List col_stats(NumericMatrix obj, Nullable<IntegerVector> cols = R_NilValue)
{
  const R_xlen_t nr = obj.nrow();
  IntegerVector idx = cols.isNull()
                          ? seq_len(obj.ncol())
                          : IntegerVector(cols.get());
  const R_xlen_t nc = idx.size();
  const int nc_obj = obj.ncol();

  const double *x = obj.begin();
  const int *ix = idx.begin();
  for (R_xlen_t j = 0; j < nc; ++j)
  {
    if (ix[j] == NA_INTEGER || ix[j] < 1 || ix[j] > nc_obj)
    {
      stop("cols must be column indices in 1..ncol(obj)");
    }
  }

  NumericMatrix stats(2, nc);
  IntegerVector row_obs(nr);
  double *st = stats.begin();
  int *robs = row_obs.begin();

  bool any_lt0 = false, any_gt1 = false, any_inf = false;
  for (R_xlen_t j = 0; j < nc; ++j)
  {
    const double *col = x + (static_cast<R_xlen_t>(ix[j]) - 1) * nr;
    double sum = 0.0, n_obs = 0.0;
    for (R_xlen_t i = 0; i < nr; ++i)
    {
      const double v = col[i];
      if (std::isfinite(v))
      {
        sum += v;
        n_obs += 1.0;
        robs[i] += 1;
        if (v < 0.0)
          any_lt0 = true;
        else if (v > 1.0)
          any_gt1 = true;
      }
      else if (!std::isnan(v))
      {
        // missing like the rest, but the caller is told it was an Inf
        any_inf = true;
      }
    }

    // Inf never enters sum now, so only a finite-value overflow can poison it
    if (!std::isfinite(sum))
    {
      return List::create(
          _["stats"] = R_NilValue,
          _["row_obs"] = R_NilValue,
          _["any_lt0"] = any_lt0,
          _["any_gt1"] = any_gt1,
          _["any_inf"] = any_inf,
          _["overflow_col"] = IntegerVector::create(static_cast<int>(j + 1)));
    }

    st[2 * j] = sum;
    st[2 * j + 1] = n_obs;
  }

  stats.attr("dimnames") = List::create(
      CharacterVector::create("sum", "n_obs"),
      R_NilValue);
  return List::create(
      _["stats"] = stats,
      _["row_obs"] = row_obs,
      _["any_lt0"] = any_lt0,
      _["any_gt1"] = any_gt1,
      _["any_inf"] = any_inf,
      _["overflow_col"] = R_NilValue);
}
