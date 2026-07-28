#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
List col_stats(NumericMatrix obj)
{
  const R_xlen_t nc = obj.ncol(), nr = obj.nrow();
  NumericMatrix stats(2, nc);
  IntegerVector row_obs(nr);

  const double *x = obj.begin();
  double *st = stats.begin();
  int *robs = row_obs.begin();

  bool any_lt0 = false, any_gt1 = false;
  for (R_xlen_t j = 0; j < nc; ++j)
  {
    const double *col = x + j * nr;
    double sum = 0.0, n_obs = 0.0;
    for (R_xlen_t i = 0; i < nr; ++i)
    {
      const double v = col[i];
      if (!std::isnan(v))
      {
        sum += v;
        n_obs += 1.0;
        robs[i] += 1;
        if (v < 0.0)
          any_lt0 = true;
        else if (v > 1.0)
          any_gt1 = true;
      }
    }

    if (!std::isfinite(sum))
    {
      R_xlen_t at = -1;
      for (R_xlen_t i = 0; i < nr; ++i)
      {
        if (std::isinf(col[i]))
        {
          at = i;
          break;
        }
      }
      return List::create(
          _["stats"] = R_NilValue,
          _["row_obs"] = R_NilValue,
          _["any_lt0"] = any_lt0,
          _["any_gt1"] = any_gt1,
          _["inf_at"] = IntegerVector::create(
              at < 0 ? NA_INTEGER : static_cast<int>(at + 1),
              static_cast<int>(j + 1)));
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
      _["inf_at"] = R_NilValue);
}
