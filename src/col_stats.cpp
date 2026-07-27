#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List col_stats(NumericMatrix obj)
{
  const int nc = obj.ncol(), nr = obj.nrow();
  NumericMatrix stats(2, nc);
  bool any_lt0 = false, any_gt1 = false;
  for (int j = 0; j < nc; ++j)
  {
    double sum = 0.0;
    double n_obs = 0.0;
    for (int i = 0; i < nr; ++i)
    {
      double x = obj(i, j);
      if (R_finite(x))
      {
        sum += x;
        n_obs += 1.0;
        if (x < 0.0)
          any_lt0 = true;
        else if (x > 1.0)
          any_gt1 = true;
      }
      else if (!ISNAN(x))
      {
        // not finite and not NA/NaN, so +/-Inf: bail out and report the
        // position; R aborts through cli, this kernel never throws
        return List::create(
            _["stats"] = R_NilValue,
            _["any_lt0"] = any_lt0,
            _["any_gt1"] = any_gt1,
            _["inf_at"] = IntegerVector::create(i + 1, j + 1));
      }
    }
    stats(0, j) = sum;
    stats(1, j) = n_obs;
  }
  stats.attr("dimnames") = List::create(
      CharacterVector::create("sum", "n_obs"),
      R_NilValue);
  return List::create(
      _["stats"] = stats,
      _["any_lt0"] = any_lt0,
      _["any_gt1"] = any_gt1,
      _["inf_at"] = R_NilValue);
}