#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// row-moment-only pass over columns not in the subset
static void sweep_rows_complement(const double *x, R_xlen_t nr, R_xlen_t nc_obj,
                                  const std::vector<char> &in_subset,
                                  const int *robs, int *rcomp,
                                  double *rmean, double *rm2,
                                  const double *inv)
{
  for (R_xlen_t j = 0; j < nc_obj; ++j)
  {
    if (in_subset[j])
      continue;
    const double *col = x + j * nr;
    for (R_xlen_t i = 0; i < nr; ++i)
    {
      const double v = col[i];
      if (std::isfinite(v))
      {
        const int n = robs[i] + (++rcomp[i]);
        const double delta = v - rmean[i];
        rmean[i] += delta * inv[n];
        rm2[i] += delta * (v - rmean[i]);
      }
    }
  }
}

// column stats over the subset. RM=true fuses row moments into the same pass.
// with the complement pass, each matrix column is read exactly once
template <bool RM>
static R_xlen_t sweep_cols(const double *x, R_xlen_t nr, R_xlen_t nc,
                           const int *ix, // nullptr => identity (all cols)
                           double *st, int *robs,
                           double *rmean, double *rm2,
                           const double *inv,
                           bool &any_lt0, bool &any_gt1, bool &any_inf)
{
  for (R_xlen_t j = 0; j < nc; ++j)
  {
    const R_xlen_t col_off = ix ? (static_cast<R_xlen_t>(ix[j]) - 1) : j;
    const double *col = x + col_off * nr;
    double sum = 0.0, n_obs = 0.0;
    for (R_xlen_t i = 0; i < nr; ++i)
    {
      const double v = col[i];
      if (std::isfinite(v))
      {
        sum += v;
        n_obs += 1.0;
        const int n = ++robs[i];
        if constexpr (RM)
        {
          // Welford with the division replaced by a cached reciprocal:
          // n is bounded by ncol, so inv[] covers every possible count
          const double delta = v - rmean[i];
          rmean[i] += delta * inv[n];
          rm2[i] += delta * (v - rmean[i]);
        }
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
    // Inf never enters sum, so only a finite-value overflow can poison it
    if (!std::isfinite(sum))
    {
      return j + 1;
    }
    st[2 * j] = sum;
    st[2 * j + 1] = n_obs;
  }
  return 0;
}

// [[Rcpp::export]]
List col_stats(NumericMatrix obj,
               Nullable<IntegerVector> cols = R_NilValue,
               bool row_moments = false)
{
  const R_xlen_t nr = obj.nrow();
  const int nc_obj = obj.ncol();
  const double *x = obj.begin();

  // NULL-cols fast path: no seq_len(ncol) materialization, no validation
  const bool all_cols = cols.isNull();
  IntegerVector idx = all_cols ? IntegerVector(0) : IntegerVector(cols.get());
  const R_xlen_t nc = all_cols ? static_cast<R_xlen_t>(nc_obj) : idx.size();
  const int *ix = all_cols ? nullptr : idx.begin();

  // marker for the complement pass. only needed when the domains differ
  const bool split = row_moments && !all_cols;
  std::vector<char> in_subset;
  if (split)
    in_subset.assign(static_cast<size_t>(nc_obj), 0);

  if (!all_cols)
  {
    for (R_xlen_t j = 0; j < nc; ++j)
    {
      if (ix[j] == NA_INTEGER || ix[j] < 1 || ix[j] > nc_obj)
      {
        stop("cols must be column indices in 1..ncol(obj)");
      }
      if (split)
      {
        // a duplicate would feed the same column into the row moments
        // twice while the complement marker excludes it only once
        if (in_subset[static_cast<size_t>(ix[j]) - 1])
        {
          stop("cols must not contain duplicate indices when row_moments = TRUE");
        }
        in_subset[static_cast<size_t>(ix[j]) - 1] = 1;
      }
    }
  }

  NumericMatrix stats(2, nc);
  IntegerVector row_obs(nr);
  double *st = stats.begin();
  int *robs = row_obs.begin();

  // initialize
  IntegerVector row_obs_comp_v = row_moments ? IntegerVector(nr) : IntegerVector(0);
  int *rcomp = row_moments ? row_obs_comp_v.begin() : nullptr;

  // accumulators allocated only when requested
  NumericVector row_mean_v = row_moments ? NumericVector(nr) : NumericVector(0);
  NumericVector row_m2_v = row_moments ? NumericVector(nr) : NumericVector(0);
  double *rmean = row_moments ? row_mean_v.begin() : nullptr;
  double *rm2 = row_moments ? row_m2_v.begin() : nullptr;

  // reciprocal table: one division per possible count instead of one per
  // finite cell
  std::vector<double> inv_v;
  const double *inv = nullptr;
  if (row_moments)
  {
    inv_v.resize(static_cast<size_t>(nc_obj) + 1);
    inv_v[0] = 0.0; // never indexed; n >= 1 on every update
    for (int n = 1; n <= nc_obj; ++n)
      inv_v[static_cast<size_t>(n)] = 1.0 / static_cast<double>(n);
    inv = inv_v.data();
  }

  bool any_lt0 = false, any_gt1 = false, any_inf = false;

  // pass 1: subset columns, fused with row moments when requested
  const R_xlen_t overflow =
      row_moments
          ? sweep_cols<true>(x, nr, nc, ix, st, robs, rmean, rm2, inv,
                             any_lt0, any_gt1, any_inf)
          : sweep_cols<false>(x, nr, nc, ix, st, robs, rmean, rm2, inv,
                              any_lt0, any_gt1, any_inf);

  if (overflow != 0)
  {
    return List::create(
        _["stats"] = R_NilValue,
        _["row_obs"] = R_NilValue,
        _["row_obs_complement"] = R_NilValue,
        _["row_mean"] = R_NilValue,
        _["row_m2"] = R_NilValue,
        _["any_lt0"] = any_lt0,
        _["any_gt1"] = any_gt1,
        _["any_inf"] = any_inf,
        _["overflow_col"] = IntegerVector::create(static_cast<int>(overflow)));
  }

  // pass 2: remaining columns, row moments only. each column is read
  // exactly once across the two passes
  if (split)
    sweep_rows_complement(x, nr, nc_obj, in_subset, robs, rcomp, rmean, rm2, inv);

  stats.attr("dimnames") = List::create(
      CharacterVector::create("sum", "n_obs"),
      R_NilValue);
  return List::create(
      _["stats"] = stats,
      _["row_obs"] = row_obs,
      _["row_obs_complement"] = row_moments ? (SEXP)row_obs_comp_v : R_NilValue,
      _["row_mean"] = row_moments ? (SEXP)row_mean_v : R_NilValue,
      _["row_m2"] = row_moments ? (SEXP)row_m2_v : R_NilValue,
      _["any_lt0"] = any_lt0,
      _["any_gt1"] = any_gt1,
      _["any_inf"] = any_inf,
      _["overflow_col"] = R_NilValue);
}
