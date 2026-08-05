#include <Rcpp.h>
#include <climits>
#include <cmath>
#include <cstdint>
#include <vector>
using namespace Rcpp;

// internal kernel. caller validates moment_sets (see check_moment_sets() in R).
// moments on iff moment_sets non-NULL. K is 1..8. zero-count rows: mean 0 / m2 0 / count 0.

// observed value range over scanned columns, seeded at beta bounds. `col` is 1-based within `cols`.
struct value_range
{
  double min_val = 0.0;
  double max_val = 1.0;
  int min_col = NA_INTEGER;
  int max_col = NA_INTEGER;
  bool any_inf = false;
};

// one welford step with a cached reciprocal (n bounded by ncol).
static inline void welford_update(int &n, double &mean, double &m2,
                                  double v, const double *inv)
{
  const int nn = ++n;
  const double delta = v - mean;
  mean += delta * inv[nn];
  m2 += delta * (v - mean);
}

// single-column scan. RM=true fuses row moments. overlapping sets share one signature block.
template <bool RM>
static bool scan_col(const double *col, R_xlen_t nr,
                     double *st_pair, int *robs,
                     int *gn, double *gmean, double *gm2,
                     const double *inv,
                     value_range &vr, int col_pos)
{
  double sum = 0.0, n_obs = 0.0;
  for (R_xlen_t i = 0; i < nr; ++i)
  {
    const double v = col[i];
    if (std::isfinite(v))
    {
      sum += v;
      n_obs += 1.0;
      ++robs[i];
      if constexpr (RM)
        welford_update(gn[i], gmean[i], gm2[i], v, inv);
      // min_val <= 0 <= 1 <= max_val always, so the two are exclusive
      if (v < vr.min_val)
      {
        vr.min_val = v;
        vr.min_col = col_pos;
      }
      else if (v > vr.max_val)
      {
        vr.max_val = v;
        vr.max_col = col_pos;
      }
    }
    else if (!std::isnan(v))
    {
      // missing like the rest, but the caller is told it was an Inf
      vr.any_inf = true;
    }
  }
  // inf never enters sum, so only a finite-value overflow can poison it
  if (!std::isfinite(sum))
    return false;
  st_pair[0] = sum;
  st_pair[1] = n_obs;
  return true;
}

// column stats over the subset. fused scan consumes marked columns.
static R_xlen_t sweep_cols(const double *x, R_xlen_t nr, R_xlen_t nc,
                           const int *ix,              // nullptr => identity (all cols)
                           std::vector<char> *pending, // nullptr => no moments
                           const int *slot,
                           int *an, double *am, double *a2,
                           double *st, int *robs,
                           const double *inv,
                           value_range &vr)
{
  for (R_xlen_t j = 0; j < nc; ++j)
  {
    const R_xlen_t col_off = ix ? (static_cast<R_xlen_t>(ix[j]) - 1) : j;
    const double *col = x + col_off * nr;
    const bool fuse = pending && (*pending)[static_cast<size_t>(col_off)];
    bool ok;
    if (fuse)
    {
      (*pending)[static_cast<size_t>(col_off)] = 0;
      const R_xlen_t off = static_cast<R_xlen_t>(slot[col_off]) * nr;
      ok = scan_col<true>(col, nr, st + 2 * j, robs,
                          an + off, am + off, a2 + off, inv,
                          vr, static_cast<int>(j) + 1);
    }
    else
    {
      ok = scan_col<false>(col, nr, st + 2 * j, robs,
                           nullptr, nullptr, nullptr, inv,
                           vr, static_cast<int>(j) + 1);
    }
    if (!ok)
      return j + 1;
  }
  return 0;
}

// moment-only pass over columns whose pending marker survived pass 1
static void sweep_moments_remaining(const double *x, R_xlen_t nr, R_xlen_t nc_obj,
                                    const std::vector<char> &pending,
                                    const int *slot,
                                    int *an, double *am, double *a2,
                                    const double *inv)
{
  for (R_xlen_t j = 0; j < nc_obj; ++j)
  {
    if (!pending[static_cast<size_t>(j)])
      continue;
    const R_xlen_t off = static_cast<R_xlen_t>(slot[j]) * nr;
    const double *col = x + j * nr;
    for (R_xlen_t i = 0; i < nr; ++i)
    {
      const double v = col[i];
      if (std::isfinite(v))
        welford_update(an[off + i], am[off + i], a2[off + i], v, inv);
    }
  }
}

// [[Rcpp::export]]
List col_stats(NumericMatrix obj,
               Nullable<IntegerVector> cols = R_NilValue,
               Nullable<List> moment_sets = R_NilValue)
{
  const R_xlen_t nr = obj.nrow();
  const int nc_obj = obj.ncol();
  const double *x = obj.begin();

  // null-cols fast path: no seq_len(ncol) materialization, no validation
  const bool all_cols = cols.isNull();
  IntegerVector idx = all_cols ? IntegerVector(0) : IntegerVector(cols.get());
  const R_xlen_t nc = all_cols ? static_cast<R_xlen_t>(nc_obj) : idx.size();
  const int *ix = all_cols ? nullptr : idx.begin();

  if (!all_cols)
  {
    // stats is an R matrix, whose column dimension is int
    if (nc > static_cast<R_xlen_t>(INT_MAX))
      stop("length(cols) exceeds the maximum matrix dimension");
    for (R_xlen_t j = 0; j < nc; ++j)
    {
      if (ix[j] == NA_INTEGER || ix[j] < 1 || ix[j] > nc_obj)
        stop("cols must be column indices in 1..ncol(obj)");
    }
  }

  // moments are on iff sets were given. NULL means none at all
  const bool moments = moment_sets.isNotNull();

  // per-column membership mask over the moment sets (1-based indices). bits encode cross-set overlap.
  int K = 0;
  CharacterVector set_names(0);
  std::vector<uint8_t> mask;
  if (moments)
  {
    List msets(moment_sets.get());
    K = msets.size();
    if (K < 1)
      stop("moment_sets must contain at least one set (or be NULL)");
    if (K > 8)
      stop("at most 8 moment sets are supported (uint8_t mask)");
    if (msets.hasAttribute("names"))
      set_names = msets.names();
    mask.assign(static_cast<size_t>(nc_obj), 0);
    for (int k = 0; k < K; ++k)
    {
      IntegerVector s = as<IntegerVector>(msets[k]);
      const R_xlen_t ns = s.size();
      for (R_xlen_t j = 0; j < ns; ++j)
      {
        if (s[j] == NA_INTEGER || s[j] < 1 || s[j] > nc_obj)
          stop("moment_sets must contain column indices in 1..ncol(obj)");
        mask[static_cast<size_t>(s[j]) - 1] |= static_cast<uint8_t>(1u << k);
      }
    }
  }

  // distinct signatures -> accumulator slots (flat 256-entry table).
  int NS = 0;
  std::vector<uint8_t> sigs;
  std::vector<int> slot;
  std::vector<char> pending;
  if (moments)
  {
    int slot_of[256];
    for (int i = 0; i < 256; ++i)
      slot_of[i] = -1;
    slot.assign(static_cast<size_t>(nc_obj), -1);
    pending.assign(static_cast<size_t>(nc_obj), 0);
    for (int j = 0; j < nc_obj; ++j)
    {
      const uint8_t mj = mask[static_cast<size_t>(j)];
      if (!mj)
        continue; // column in no set: moments never touch it
      int s = slot_of[mj];
      if (s < 0)
      {
        s = NS++;
        slot_of[mj] = s;
        sigs.push_back(mj);
      }
      slot[static_cast<size_t>(j)] = s;
      pending[static_cast<size_t>(j)] = 1;
    }
  }

  NumericMatrix stats(2, static_cast<int>(nc));
  IntegerVector row_obs(nr);
  double *st = stats.begin();
  int *robs = row_obs.begin();

  // one nr-length welford block per distinct signature.
  std::vector<double> a_mean, a_m2;
  std::vector<int> a_n;
  if (moments)
  {
    a_mean.assign(static_cast<size_t>(NS) * nr, 0.0);
    a_m2.assign(static_cast<size_t>(NS) * nr, 0.0);
    a_n.assign(static_cast<size_t>(NS) * nr, 0);
  }

  // reciprocal table: one division per possible count.
  std::vector<double> inv_v;
  const double *inv = nullptr;
  if (moments)
  {
    inv_v.resize(static_cast<size_t>(nc_obj) + 1);
    inv_v[0] = 0.0; // defensive only. every real update/merge has n >= 1
    for (int n = 1; n <= nc_obj; ++n)
      inv_v[static_cast<size_t>(n)] = 1.0 / static_cast<double>(n);
    inv = inv_v.data();
  }

  value_range vr;

  // nr == 0 skips the pointer-based passes (empty vector::data may be null).
  if (nr > 0)
  {
    // pass 1: stats subset. marked columns inside it are fused and consumed
    const R_xlen_t overflow =
        sweep_cols(x, nr, nc, ix,
                   moments ? &pending : nullptr,
                   moments ? slot.data() : nullptr,
                   a_n.data(), a_mean.data(), a_m2.data(),
                   st, robs, inv, vr);

    if (overflow != 0)
    {
      return List::create(
          _["stats"] = R_NilValue,
          _["row_obs"] = R_NilValue,
          _["row_moment_obs"] = R_NilValue,
          _["row_mean"] = R_NilValue,
          _["row_m2"] = R_NilValue,
          _["min_val"] = vr.min_val,
          _["max_val"] = vr.max_val,
          _["min_col"] = vr.min_col,
          _["max_col"] = vr.max_col,
          _["any_inf"] = vr.any_inf,
          _["overflow_col"] = IntegerVector::create(static_cast<int>(overflow)));
    }

    // pass 2: marked columns not covered by pass 1. each column read at most once.
    if (moments && !all_cols)
      sweep_moments_remaining(x, nr, nc_obj, pending, slot.data(),
                              a_n.data(), a_mean.data(), a_m2.data(), inv);
  }

  // outputs stay at function scope until List::create protects the return.
  IntegerMatrix Nout(moments ? static_cast<int>(nr) : 0, K);
  NumericMatrix Mout(moments ? static_cast<int>(nr) : 0, K);
  NumericMatrix M2out(moments ? static_cast<int>(nr) : 0, K);

  // expand signatures into per-set outputs via Chan's pairwise combine.
  if (moments && nr > 0)
  {
    for (int k = 0; k < K; ++k)
    {
      int *on = Nout.begin() + static_cast<R_xlen_t>(k) * nr;
      double *om = Mout.begin() + static_cast<R_xlen_t>(k) * nr;
      double *o2 = M2out.begin() + static_cast<R_xlen_t>(k) * nr;
      const uint8_t bit = static_cast<uint8_t>(1u << k);
      for (int s = 0; s < NS; ++s)
      {
        if (!(sigs[static_cast<size_t>(s)] & bit))
          continue;
        const R_xlen_t off = static_cast<R_xlen_t>(s) * nr;
        const int *gn = a_n.data() + off;
        const double *gm = a_mean.data() + off;
        const double *g2 = a_m2.data() + off;
        for (R_xlen_t i = 0; i < nr; ++i)
        {
          const int nb = gn[i];
          if (!nb)
            continue;
          const int na = on[i];
          if (na == 0)
          {
            on[i] = nb;
            om[i] = gm[i];
            o2[i] = g2[i];
            continue;
          }
          const int nt = na + nb;
          const double wb = static_cast<double>(nb) * inv[nt];
          const double delta = gm[i] - om[i];
          om[i] += delta * wb;
          o2[i] += g2[i] + delta * (delta * (static_cast<double>(na) * wb));
          on[i] = nt;
        }
      }
    }
  }
  if (moments && set_names.size() == K)
  {
    List dn = List::create(R_NilValue, set_names);
    Nout.attr("dimnames") = dn;
    Mout.attr("dimnames") = dn;
    M2out.attr("dimnames") = dn;
  }

  stats.attr("dimnames") = List::create(
      CharacterVector::create("sum", "n_obs"),
      R_NilValue);
  return List::create(
      _["stats"] = stats,
      _["row_obs"] = row_obs,
      _["row_moment_obs"] = moments ? (SEXP)Nout : R_NilValue,
      _["row_mean"] = moments ? (SEXP)Mout : R_NilValue,
      _["row_m2"] = moments ? (SEXP)M2out : R_NilValue,
      _["min_val"] = vr.min_val,
      _["max_val"] = vr.max_val,
      _["min_col"] = vr.min_col,
      _["max_col"] = vr.max_col,
      _["any_inf"] = vr.any_inf,
      _["overflow_col"] = R_NilValue);
}
