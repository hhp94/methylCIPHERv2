#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// One column-major sweep over the panel slice, which is the widest thing this
// package touches. Three choices make the inner loop cheap, in order of what
// they were worth (measured, 100 x 39,025, -O2):
//
//   * a raw column pointer instead of `obj(i, j)`, which recomputes i + j * nr
//   * `std::isnan` instead of `ISNAN` / `R_finite` -- the libR predicates are
//     out-of-line calls, so nothing across them vectorizes
//   * ONE predicate in the loop, and no early exit out of it
//
// The last needs explaining. There are two questions to ask of a value -- "is
// it missing?" (NA/NaN) and "is it usable?" (also excludes +/-Inf) -- and the
// loop only asks the first. An Inf is let through into `sum`, which poisons it
// to +/-Inf or NaN, and one `isfinite(sum)` per COLUMN catches that: nr * nc
// tests become nc. Locating the offending row is then a cold rescan of the one
// bad column, which happens at most once per call and on the way to an abort.
//
// So the classification the caller depends on still holds exactly --
//   NA_real_, NaN -> skipped, not counted    (missing; the fill handles it)
//   +/-Inf        -> poisons the column sum  (bad data; report the position)
//   anything else -> summed and counted
// -- it is just established per column rather than per element.
//
// `!isfinite(sum)` has one other cause: a finite column whose sum overflows.
// That needs values around 1e306, which no beta matrix has, but the value gate
// only WARNS on out-of-range input so nothing stops one arriving. The rescan
// tells the two apart -- no Inf found means overflow -- and reports row NA,
// which `check_col_values()` words differently. Both are hard stops either way.

// [[Rcpp::export]]
List col_stats(NumericMatrix obj)
{
  const R_xlen_t nc = obj.ncol(), nr = obj.nrow();
  NumericMatrix stats(2, nc);
  // per-row observed count over this submatrix. The row gate reads it, so the
  // all-NA-sample check costs nothing beyond the sweep already being made.
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
      // cold path: this column holds an Inf, or it overflowed. Find out which,
      // and report a 1-based position; R aborts through cli, this kernel never
      // throws. `stats` is not filled in, so the caller must not read it.
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

    // stats is 2 x nc and column-major: (sum, n_obs) of column j sit adjacent
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
