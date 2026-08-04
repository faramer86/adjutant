#include "adjutant_pgs_eval.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cerrno>

#include "adjutant_cmdline.h"
#include "adjutant_decompress.h"
#include "include/plink2_bits.h"
#include "include/plink2_htable.h"
#include "include/plink2_string.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Continued-fraction expansion for the incomplete beta function.
static double BetaCf(double aa, double bb, double xx) {
  const int32_t kMaxIter = 300;
  const double kEps = 3.0e-15;
  const double kFpMin = 1.0e-300;
  const double qab = aa + bb;
  const double qap = aa + 1.0;
  const double qam = aa - 1.0;
  double cc = 1.0;
  double dd = 1.0 - qab * xx / qap;
  if (fabs(dd) < kFpMin) {
    dd = kFpMin;
  }
  dd = 1.0 / dd;
  double hh = dd;
  for (int32_t mm = 1; mm <= kMaxIter; ++mm) {
    const int32_t m2 = 2 * mm;
    double numer = mm * (bb - mm) * xx / ((qam + m2) * (aa + m2));
    dd = 1.0 + numer * dd;
    if (fabs(dd) < kFpMin) { dd = kFpMin; }
    cc = 1.0 + numer / cc;
    if (fabs(cc) < kFpMin) { cc = kFpMin; }
    dd = 1.0 / dd;
    hh *= dd * cc;
    numer = -(aa + mm) * (qab + mm) * xx / ((aa + m2) * (qap + m2));
    dd = 1.0 + numer * dd;
    if (fabs(dd) < kFpMin) { dd = kFpMin; }
    cc = 1.0 + numer / cc;
    if (fabs(cc) < kFpMin) { cc = kFpMin; }
    dd = 1.0 / dd;
    const double del = dd * cc;
    hh *= del;
    if (fabs(del - 1.0) < kEps) {
      break;
    }
  }
  return hh;
}

// Regularized incomplete beta I_x(a, b).
static double BetaI(double aa, double bb, double xx) {
  if (xx <= 0.0) {
    return 0.0;
  }
  if (xx >= 1.0) {
    return 1.0;
  }
  const double front = exp(lgamma(aa + bb) - lgamma(aa) - lgamma(bb) +
                           aa * log(xx) + bb * log1p(-xx));
  if (xx < (aa + 1.0) / (aa + bb + 2.0)) {
    return front * BetaCf(aa, bb, xx) / aa;
  }
  return 1.0 - front * BetaCf(bb, aa, 1.0 - xx) / bb;
}

double PgsEvalTPval(double tt, double df) {
  if (!(df > 0.0)) {
    return 1.0;
  }
  // Two-sided: P(|T| > |t|) = I_{df/(df+t^2)}(df/2, 1/2).
  return BetaI(0.5 * df, 0.5, df / (df + tt * tt));
}

// dest = X' X, packed predictor_ct x predictor_ct, column-major.
static void ComputeXtx(const double* xx, uint32_t sample_ct,
                       uint32_t predictor_ct, double* dest) {
  for (uint32_t cc = 0; cc != predictor_ct; ++cc) {
    for (uint32_t rr = 0; rr != predictor_ct; ++rr) {
      const double* col_r = &(xx[(uintptr_t)rr * sample_ct]);
      const double* col_c = &(xx[(uintptr_t)cc * sample_ct]);
      double acc = 0.0;
      for (uint32_t ii = 0; ii != sample_ct; ++ii) {
        acc += col_r[ii] * col_c[ii];
      }
      dest[(uintptr_t)cc * predictor_ct + rr] = acc;
    }
  }
}

// resid = vv - X (X'X)^-1 X' vv
static void Residualize(const PgsEvalOlsCtx* ctx, const double* vv,
                        double* resid) {
  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t predictor_ct = ctx->predictor_ct;
  double* xtv = ctx->work;             // predictor_ct
  double* coef = &(ctx->work[predictor_ct]);  // predictor_ct
  for (uint32_t pp = 0; pp != predictor_ct; ++pp) {
    const double* col = &(ctx->xx[(uintptr_t)pp * sample_ct]);
    double acc = 0.0;
    for (uint32_t ii = 0; ii != sample_ct; ++ii) {
      acc += col[ii] * vv[ii];
    }
    xtv[pp] = acc;
  }
  for (uint32_t rr = 0; rr != predictor_ct; ++rr) {
    double acc = 0.0;
    for (uint32_t cc = 0; cc != predictor_ct; ++cc) {
      acc += ctx->xtx_inv[(uintptr_t)cc * predictor_ct + rr] * xtv[cc];
    }
    coef[rr] = acc;
  }
  // PgsEvalOlsOne residualises in place, passing the same pointer for vv and
  // resid.  memcpy requires non-overlapping regions, so identical pointers are
  // undefined behaviour -- it would appear to work and then break under a
  // different compiler or a sanitiser.
  if (resid != vv) {
    memcpy(resid, vv, sample_ct * sizeof(double));
  }
  for (uint32_t pp = 0; pp != predictor_ct; ++pp) {
    const double* col = &(ctx->xx[(uintptr_t)pp * sample_ct]);
    const double beta_p = coef[pp];
    for (uint32_t ii = 0; ii != sample_ct; ++ii) {
      resid[ii] -= beta_p * col[ii];
    }
  }
}

// Null every pointer PgsEvalOlsCleanup frees, and set fold_ct to the "no
// cross-validation" value PgsEvalOlsSetFolds overwrites.
//
// Both Inits need this, and PgsEvalLogisticInit needs it on its embedded
// context before its own failure returns, which happen before PgsEvalOlsInit
// runs.  Keeping it in one function is what makes adding an owned pointer to
// PgsEvalOlsCtx safe: miss one of two open-coded copies and cleanup frees
// uninitialised memory on a path nothing routinely exercises.
static void PgsEvalOlsPreinit(PgsEvalOlsCtx* ctx) {
  ctx->xtx_inv = nullptr;
  ctx->e_y = nullptr;
  ctx->work = nullptr;
  ctx->fold_assign = nullptr;
  ctx->fold_xtx_inv = nullptr;
  ctx->fold_beta = nullptr;
  ctx->cv_work = nullptr;
  ctx->fold_ct = 1;
  ctx->r2_null_oof = 0.0;
}

BoolErr PgsEvalOlsInit(const double* yy, const double* xx, uint32_t sample_ct,
                       uint32_t predictor_ct, PgsEvalOlsCtx* ctx) {
  ctx->xx = xx;
  ctx->yy = yy;
  ctx->sample_ct = sample_ct;
  ctx->predictor_ct = predictor_ct;
  PgsEvalOlsPreinit(ctx);
  ctx->xtx_inv = S_CAST(double*, malloc((uintptr_t)predictor_ct * predictor_ct * sizeof(double)));
  ctx->e_y = S_CAST(double*, malloc(sample_ct * sizeof(double)));
  ctx->work = S_CAST(double*, malloc((sample_ct + 2 * predictor_ct) * sizeof(double)));
  if ((!ctx->xtx_inv) || (!ctx->e_y) || (!ctx->work)) {
    PgsEvalOlsCleanup(ctx);
    return 1;
  }
  ComputeXtx(xx, sample_ct, predictor_ct, ctx->xtx_inv);

  MatrixInvertBuf1* inv_buf = S_CAST(MatrixInvertBuf1*,
      malloc(kMatrixInvertBuf1CheckedAlloc * predictor_ct));
  double* dbl_2d_buf = S_CAST(double*, malloc((uintptr_t)predictor_ct * MAXV(predictor_ct, 3) * sizeof(double)));
  if ((!inv_buf) || (!dbl_2d_buf) ||
      InvertSymmdefMatrixChecked(predictor_ct, ctx->xtx_inv, inv_buf, dbl_2d_buf)) {
    free(inv_buf);
    free(dbl_2d_buf);
    PgsEvalOlsCleanup(ctx);
    return 1;  // rank-deficient: e.g. a constant covariate vs the intercept
  }
  ReflectMatrix(predictor_ct, ctx->xtx_inv);
  free(inv_buf);
  free(dbl_2d_buf);

  Residualize(ctx, yy, ctx->e_y);

  double ysum = 0.0;
  for (uint32_t ii = 0; ii != sample_ct; ++ii) {
    ysum += yy[ii];
  }
  const double ymean = ysum / S_CAST(double, sample_ct);
  double tss = 0.0;
  for (uint32_t ii = 0; ii != sample_ct; ++ii) {
    const double dev = yy[ii] - ymean;
    tss += dev * dev;
  }
  double rss = 0.0;
  for (uint32_t ii = 0; ii != sample_ct; ++ii) {
    rss += ctx->e_y[ii] * ctx->e_y[ii];
  }
  ctx->tss = tss;
  ctx->rss_null = rss;
  ctx->r2_null = (tss > 0.0) ? (1.0 - rss / tss) : 0.0;
  return 0;
}

void PgsEvalOlsCleanup(PgsEvalOlsCtx* ctx) {
  free(ctx->xtx_inv);
  free(ctx->e_y);
  free(ctx->work);
  free(ctx->fold_assign);
  free(ctx->fold_xtx_inv);
  free(ctx->fold_beta);
  free(ctx->cv_work);
  ctx->xtx_inv = nullptr;
  ctx->e_y = nullptr;
  ctx->work = nullptr;
  ctx->fold_assign = nullptr;
  ctx->fold_xtx_inv = nullptr;
  ctx->fold_beta = nullptr;
  ctx->cv_work = nullptr;
  ctx->fold_ct = 1;
}

// Analytic 95% CI for the incremental R^2, via the Fisher z transform of the
// partial correlation between y and the score, both residualised on the
// covariates.  O(1) on top of dot products PgsEvalOlsOne has already computed.
//
//   r  = (e_g . e_y) / sqrt((e_g . e_g)(e_y . e_y))
//   z  = atanh(r),  SE(z) = 1 / sqrt(n - predictor_ct - 2)
//   CI(r)    = tanh(z +/- 1.96 * SE)
//   CI(dR2)  = CI(r)^2 * (1 - R2_null)
//
// Checked against a 2000-replicate percentile bootstrap on both fixtures (see
// Tests/pgs_eval_fixtures/generate_ci.R).  The two agree closely, and neither
// is uniformly the wider: measured widths are 0.15170 analytic against 0.13832
// bootstrap on ols_reference.txt, and 0.11261 against 0.11387 on
// logistic_reference.txt.  So it runs ~10% wide on one fixture and ~1% narrow
// on the other -- agreement, not conservatism.  Two fixtures cannot establish a
// direction anyway, and the tests assert overlap, never a bound.
//
// eyey_dot is ctx->rss_null, which is TSS * (1 - R2_null).  Substituting that
// into the last line gives r^2 * (1 - R2_null) = gy^2 / (gg * TSS), which is
// exactly r2_incr -- so the interval brackets its own point estimate
// everywhere except in the clamped corner below, which ClampCiToPoint then
// takes care of.  The tests assert containment on every fixture row.
//
// Split in two so the cross-validated path can reach the Fisher z machinery
// with its own partial correlation (see PgsEvalOofR2) without duplicating it.
static void IncrR2CiFromR(double rr, double r2_null, uint32_t sample_ct,
                          uint32_t predictor_ct, double* lo, double* hi) {
  if (sample_ct <= predictor_ct + 2) {
    *lo = 0.0 / 0.0; *hi = 0.0 / 0.0; return;
  }
  if (rr > 0.9999999) { rr = 0.9999999; }
  if (rr < -0.9999999) { rr = -0.9999999; }
  const double zz = atanh(rr);
  const double se = 1.0 / sqrt(S_CAST(double, sample_ct - predictor_ct - 2));
  const double lo_r = tanh(zz - 1.959963984540054 * se);
  const double hi_r = tanh(zz + 1.959963984540054 * se);
  const double scale = 1.0 - r2_null;
  // Squaring an interval that straddles zero: the minimum of r^2 is 0, NOT
  // lo_r^2.  Without this the interval can exclude its own point estimate
  // whenever the effect is near null -- precisely the case the CI exists for.
  const double lo_sq = (lo_r <= 0.0 && hi_r >= 0.0)
                     ? 0.0 : MINV(lo_r * lo_r, hi_r * hi_r);
  const double hi_sq = MAXV(lo_r * lo_r, hi_r * hi_r);
  *lo = lo_sq * scale;
  *hi = hi_sq * scale;
}

// Force the emitted bounds to contain the point estimate they are printed next
// to.  A row whose interval excludes its own estimate contradicts itself, and
// that invariant is the reason the interval is worth having at all, so it is
// enforced here rather than argued about.
//
// The Fisher z path gives containment for free except in one corner: the
// |r| < 0.9999999 clamp in IncrR2CiFromR caps the upper bound at
// tanh(atanh(0.9999999) + 1.96/sqrt(n-p-2))^2 * (1 - R2_null), which a
// near-perfect fit exceeds (see PgsEvalOofR2 for the exact condition).  Both
// paths reach that corner; both call this.
//
// NaN-safe by construction: both comparisons are false when either side is
// NaN, so a NaN interval -- too few degrees of freedom, or a negative
// out-of-fold r2_incr -- stays NaN rather than collapsing onto the estimate.
static void ClampCiToPoint(double point, double* lo, double* hi) {
  if (*lo > point) { *lo = point; }
  if (*hi < point) { *hi = point; }
}

// The in-sample entry point: r straight from the dot products PgsEvalOlsOne
// already has.
static void IncrR2Ci(double gy_dot, double gg_dot, double eyey_dot,
                     double r2_null, uint32_t sample_ct,
                     uint32_t predictor_ct, double* lo, double* hi) {
  const double denom = sqrt(gg_dot * eyey_dot);
  if (!(denom > 0.0)) {
    *lo = 0.0 / 0.0; *hi = 0.0 / 0.0; return;
  }
  IncrR2CiFromR(gy_dot / denom, r2_null, sample_ct, predictor_ct, lo, hi);
}

// ---------------------------------------------------------------------------
// k-fold cross-validation
// ---------------------------------------------------------------------------

// SplitMix64 (Steele, Lea & Flood 2014).  Self-contained and reproducible from
// a seed, which is all the fold assignment needs; deliberately not plink2's
// shared sfmt stream, so drawing folds cannot shift anything else that draws
// from it.
static uint64_t SplitMix64(uint64_t* state) {
  uint64_t zz = (*state += 0x9e3779b97f4a7c15ULL);
  zz = (zz ^ (zz >> 30)) * 0xbf58476d1ce4e5b9ULL;
  zz = (zz ^ (zz >> 27)) * 0x94d049bb133111ebULL;
  return zz ^ (zz >> 31);
}

// Shuffle idxs[0, ct) (Fisher-Yates), then deal folds round-robin, so this
// stratum contributes floor(ct / fold_ct) or ceil(ct / fold_ct) samples to
// every fold.  Called once per stratum off a single PRNG stream.
static void AssignFoldsInStratum(uint32_t* idxs, uint32_t ct, uint32_t fold_ct,
                                 uint64_t* rng, uint32_t* fold_assign) {
  for (uint32_t ii = ct; ii > 1; --ii) {
    // Modulo on a full 64-bit draw.  Any one outcome's probability differs
    // from 1/ii by at most ii / 2^64, so at the largest ii a uint32_t sample
    // count can reach the bias is under 2^-32 -- far below the granularity of
    // the fold sizes it would have to move.
    const uint32_t jj = S_CAST(uint32_t, SplitMix64(rng) % ii);
    const uint32_t tmp = idxs[ii - 1];
    idxs[ii - 1] = idxs[jj];
    idxs[jj] = tmp;
  }
  for (uint32_t ii = 0; ii != ct; ++ii) {
    fold_assign[idxs[ii]] = ii % fold_ct;
  }
}

BoolErr PgsEvalOlsSetFolds(PgsEvalOlsCtx* ctx, uint32_t fold_ct,
                           uint32_t stratify_binary, uint64_t seed) {
  if (fold_ct <= 1) {
    return 0;  // no cross-validation; ctx stays exactly as PgsEvalOlsInit left it
  }
  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t predictor_ct = ctx->predictor_ct;
  const double* yy = ctx->yy;
  const double* xx = ctx->xx;

  ctx->fold_assign = S_CAST(uint32_t*, malloc(sample_ct * sizeof(uint32_t)));
  ctx->fold_xtx_inv = S_CAST(double*, malloc(
      (uintptr_t)fold_ct * predictor_ct * predictor_ct * sizeof(double)));
  ctx->fold_beta = S_CAST(double*, malloc(
      (uintptr_t)fold_ct * predictor_ct * sizeof(double)));
  ctx->cv_work = S_CAST(double*, malloc(
      (3 * (uintptr_t)sample_ct + 2 * predictor_ct) * sizeof(double)));
  uint32_t* strat_idx = S_CAST(uint32_t*, malloc(sample_ct * sizeof(uint32_t)));
  // Sized for predictor_ct; the MAXV(., 3) is load-bearing because dpocon
  // inside InvertSymmdefMatrixChecked needs a 3*N work vector.
  MatrixInvertBuf1* inv_buf = S_CAST(MatrixInvertBuf1*,
      malloc(kMatrixInvertBuf1CheckedAlloc * predictor_ct));
  double* dbl_2d_buf = S_CAST(double*, malloc(
      (uintptr_t)predictor_ct * MAXV(predictor_ct, 3) * sizeof(double)));
  if ((!ctx->fold_assign) || (!ctx->fold_xtx_inv) || (!ctx->fold_beta) ||
      (!ctx->cv_work) || (!strat_idx) || (!inv_buf) || (!dbl_2d_buf)) {
    free(strat_idx); free(inv_buf); free(dbl_2d_buf);
    return 1;  // PgsEvalOlsCleanup frees whatever did get allocated
  }

  {
    uint64_t rng = seed;
    if (stratify_binary) {
      // Controls, then cases, off one stream.  Without the split, a fold of an
      // imbalanced phenotype can come out with almost no cases, and the R^2 of
      // a fit on that fold means nothing.
      for (uint32_t pass = 0; pass != 2; ++pass) {
        uint32_t ct = 0;
        for (uint32_t ii = 0; ii != sample_ct; ++ii) {
          if ((yy[ii] > 0.5) == (pass == 1)) {
            strat_idx[ct++] = ii;
          }
        }
        AssignFoldsInStratum(strat_idx, ct, fold_ct, &rng, ctx->fold_assign);
      }
    } else {
      for (uint32_t ii = 0; ii != sample_ct; ++ii) {
        strat_idx[ii] = ii;
      }
      AssignFoldsInStratum(strat_idx, sample_ct, fold_ct, &rng, ctx->fold_assign);
    }
  }

  // Per training set: factor X'X, fit the null model, and predict into the
  // held-out fold.  None of this depends on the PGS column, so it is paid once
  // for the whole model list rather than once per model.
  double* resid_null = ctx->cv_work;
  double* xty = &(ctx->cv_work[3 * (uintptr_t)sample_ct]);
  for (uint32_t ff = 0; ff != fold_ct; ++ff) {
    double* xtx = &(ctx->fold_xtx_inv[(uintptr_t)ff * predictor_ct * predictor_ct]);
    double* beta = &(ctx->fold_beta[(uintptr_t)ff * predictor_ct]);
    uint32_t train_ct = 0;
    for (uint32_t ii = 0; ii != sample_ct; ++ii) {
      if (ctx->fold_assign[ii] != ff) {
        ++train_ct;
      }
    }
    if (train_ct <= predictor_ct + 1) {
      free(strat_idx); free(inv_buf); free(dbl_2d_buf);
      return 1;
    }
    for (uint32_t cc = 0; cc != predictor_ct; ++cc) {
      const double* col_c = &(xx[(uintptr_t)cc * sample_ct]);
      for (uint32_t rr = 0; rr <= cc; ++rr) {
        const double* col_r = &(xx[(uintptr_t)rr * sample_ct]);
        double acc = 0.0;
        for (uint32_t ii = 0; ii != sample_ct; ++ii) {
          if (ctx->fold_assign[ii] != ff) {
            acc += col_r[ii] * col_c[ii];
          }
        }
        xtx[(uintptr_t)cc * predictor_ct + rr] = acc;
      }
      double acc = 0.0;
      for (uint32_t ii = 0; ii != sample_ct; ++ii) {
        if (ctx->fold_assign[ii] != ff) {
          acc += col_c[ii] * yy[ii];
        }
      }
      xty[cc] = acc;
    }
    if (InvertSymmdefMatrixChecked(predictor_ct, xtx, inv_buf, dbl_2d_buf)) {
      free(strat_idx); free(inv_buf); free(dbl_2d_buf);
      return 1;  // e.g. a covariate that happens to be constant on this fold
    }
    ReflectMatrix(predictor_ct, xtx);
    for (uint32_t rr = 0; rr != predictor_ct; ++rr) {
      double acc = 0.0;
      for (uint32_t cc = 0; cc != predictor_ct; ++cc) {
        acc += xtx[(uintptr_t)cc * predictor_ct + rr] * xty[cc];
      }
      beta[rr] = acc;
    }
    for (uint32_t ii = 0; ii != sample_ct; ++ii) {
      if (ctx->fold_assign[ii] != ff) {
        continue;
      }
      double eta = 0.0;
      for (uint32_t kk = 0; kk != predictor_ct; ++kk) {
        eta += xx[(uintptr_t)kk * sample_ct + ii] * beta[kk];
      }
      resid_null[ii] = yy[ii] - eta;
    }
  }
  free(strat_idx);
  free(inv_buf);
  free(dbl_2d_buf);

  double ss = 0.0;
  for (uint32_t ii = 0; ii != sample_ct; ++ii) {
    ss += resid_null[ii] * resid_null[ii];
  }
  ctx->r2_null_oof = (ctx->tss > 0.0) ? (1.0 - ss / ctx->tss) : 0.0;
  ctx->fold_ct = fold_ct;  // last: a failure above leaves the context in-sample
  return 0;
}

// Out-of-fold R^2 family for one PGS column.  Overwrites r2_null / r2_full /
// r2_incr and the interval on r2_incr; beta / se / pval are left exactly as the
// whole-sample fit set them.
//
// All three R^2 values come from the same pooled set -- every sample predicted
// once, by the fit that did not see it -- against the same TSS, so
// r2_incr == r2_full - r2_null holds by construction and not by luck.
//
// gmean / gsd are the ones PgsEvalOlsOne already computed, so the standardised
// column here is the same one the in-sample fit used.
//
// Returns 1 if a training set leaves the PGS collinear with the covariates,
// which the whole-sample fit can pass while a fold does not.
static BoolErr PgsEvalOofR2(const PgsEvalOlsCtx* ctx, const double* gg,
                            double gmean, double gsd, PgsEvalResult* result) {
  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t predictor_ct = ctx->predictor_ct;
  const double* xx = ctx->xx;
  const double* yy = ctx->yy;
  const uint32_t* fold_assign = ctx->fold_assign;
  const double* resid_null = ctx->cv_work;
  double* incr  = &(ctx->cv_work[sample_ct]);
  double* g_std = &(ctx->cv_work[2 * (uintptr_t)sample_ct]);
  double* xtg   = &(ctx->cv_work[3 * (uintptr_t)sample_ct]);
  double* gamma = &(ctx->cv_work[3 * (uintptr_t)sample_ct + predictor_ct]);

  for (uint32_t ii = 0; ii != sample_ct; ++ii) {
    g_std[ii] = (gg[ii] - gmean) / gsd;
  }

  for (uint32_t ff = 0; ff != ctx->fold_ct; ++ff) {
    const double* xtx_inv = &(ctx->fold_xtx_inv[(uintptr_t)ff * predictor_ct * predictor_ct]);
    const double* beta_null = &(ctx->fold_beta[(uintptr_t)ff * predictor_ct]);
    // Frisch-Waugh-Lovell inside the training set: regress the PGS on the
    // covariates there, then the held-out prediction of the full model is the
    // held-out prediction of the null model plus beta_g times the same
    // residualisation evaluated on the held-out rows.  Fitting the full model
    // outright would give the identical predictions.
    for (uint32_t kk = 0; kk != predictor_ct; ++kk) {
      const double* col = &(xx[(uintptr_t)kk * sample_ct]);
      double acc = 0.0;
      for (uint32_t ii = 0; ii != sample_ct; ++ii) {
        if (fold_assign[ii] != ff) {
          acc += col[ii] * g_std[ii];
        }
      }
      xtg[kk] = acc;
    }
    for (uint32_t rr = 0; rr != predictor_ct; ++rr) {
      double acc = 0.0;
      for (uint32_t cc = 0; cc != predictor_ct; ++cc) {
        acc += xtx_inv[(uintptr_t)cc * predictor_ct + rr] * xtg[cc];
      }
      gamma[rr] = acc;
    }
    double gg_dot = 0.0;
    double gy_dot = 0.0;
    uint32_t train_ct = 0;
    for (uint32_t ii = 0; ii != sample_ct; ++ii) {
      double e_g = g_std[ii];
      double eta = 0.0;
      for (uint32_t kk = 0; kk != predictor_ct; ++kk) {
        const double xval = xx[(uintptr_t)kk * sample_ct + ii];
        e_g -= xval * gamma[kk];
        eta += xval * beta_null[kk];
      }
      if (fold_assign[ii] == ff) {
        incr[ii] = e_g;  // scaled by beta_g below, once it is known
        continue;
      }
      gg_dot += e_g * e_g;
      gy_dot += e_g * (yy[ii] - eta);
      ++train_ct;
    }
    if (!(gg_dot > 1e-10 * S_CAST(double, train_ct - 1))) {
      return 1;
    }
    const double beta_g = gy_dot / gg_dot;
    for (uint32_t ii = 0; ii != sample_ct; ++ii) {
      if (fold_assign[ii] == ff) {
        incr[ii] *= beta_g;
      }
    }
  }

  double ss_res_full = 0.0;
  for (uint32_t ii = 0; ii != sample_ct; ++ii) {
    const double dev = resid_null[ii] - incr[ii];
    ss_res_full += dev * dev;
  }
  result->r2_null = ctx->r2_null_oof;
  result->r2_full = (ctx->tss > 0.0) ? (1.0 - ss_res_full / ctx->tss) : 0.0;
  result->r2_incr = result->r2_full - result->r2_null;

  // The interval, from the same Fisher z path as the in-sample one.
  //
  // The partial correlation it needs is the one implied by the reported
  // r2_incr: r = sqrt(r2_incr / (1 - r2_null)).  In-sample that is an identity
  // -- rearranging it gives r^2 * (1 - R2_null) = r2_incr, the line IncrR2Ci's
  // containment guarantee rests on -- so this is the same definition, not a
  // second one, and at fold_ct == 1 it returns the same interval.  (Only |r|
  // survives, but the interval depends on r only through r^2, so the sign is
  // not information the output ever carried.)
  //
  // Taking r instead as corr(y - yhat_null, yhat_full - yhat_null) directly,
  // which is the same thing in-sample, does NOT work here.  Write a and b for
  // those two pooled vectors, s = a.b and q = b.b.  In-sample s == q, because
  // beta was chosen to make it so; out of fold it is not, and
  //
  //     s^2 / (q * TSS)  -  r2_incr  =  (s - q)^2 / (q * TSS)  >=  0
  //
  // so the interval would be centred on a quantity strictly above the r2_incr
  // the row reports, by an amount that grows with how badly the model failed to
  // generalise.  Containment would then hold only by luck, and the invariant
  // the whole interval is worth having for is that it never does.
  //
  // Here it holds by construction, for the same reason as in-sample: CI(r)
  // contains r, squaring maps it to an interval bracketing r^2 (with the
  // straddle-zero rule), and scaling by 1 - r2_null > 0 preserves the order.
  //
  // The one crack is the |r| < 0.9999999 clamp inside IncrR2CiFromR, which the
  // in-sample path shares.  Once r is clamped the upper bound stops tracking
  // it, so the bound falls below the estimate as soon as
  //
  //     r > tanh(atanh(0.9999999) + 1.96 / sqrt(n - p - 2))
  //
  // which is 0.99999995 at n = 30 and falls toward 0.9999999 as n grows -- a
  // residual sum of squares below roughly 1e-7 to 2e-7 of the null's, not the
  // zero an earlier version of this comment claimed.  A score column equal to
  // the phenotype plus 1e-5 noise reaches it, which is what handing the tool
  // the phenotype itself as a model by mistake looks like.  ClampCiToPoint
  // below closes it on both paths; the unit tests assert containment there at
  // full precision, which pgs_eval_check.sh cannot (dtoa_g prints ~6 digits).
  //
  // A negative r2_incr, which out of fold is a real and meaningful answer (the
  // PGS predicts held-out samples worse than nothing), has no such r.  Rather
  // than report an interval on a correlation that does not exist -- necessarily
  // non-negative, hence excluding the point estimate -- both bounds are NaN and
  // the writers emit NA.
  //
  // APPROXIMATE: out-of-fold predictions from different folds share training
  // samples, so they are not independent and the nominal 95% is nominal.  The
  // dependence falls off with fold_ct; it is weak by fold_ct = 5.
  const double scale = 1.0 - result->r2_null;
  if ((result->r2_incr < 0.0) || (!(scale > 0.0))) {
    result->r2_incr_l95 = 0.0 / 0.0;
    result->r2_incr_u95 = 0.0 / 0.0;
  } else {
    IncrR2CiFromR(sqrt(result->r2_incr / scale), result->r2_null,
                  sample_ct, predictor_ct,
                  &result->r2_incr_l95, &result->r2_incr_u95);
    ClampCiToPoint(result->r2_incr, &result->r2_incr_l95,
                   &result->r2_incr_u95);
  }
  return 0;
}

void PgsEvalOlsOne(const PgsEvalOlsCtx* ctx, const double* gg,
                   PgsEvalResult* result) {
  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t predictor_ct = ctx->predictor_ct;
  const double nan_val = 0.0 / 0.0;
  // Every field, including the binary-trait ones this path never computes:
  // *result is typically a row of a malloc'd array, so a field left unwritten
  // is indeterminate rather than reliably wrong.  See PgsEvalResult.
  result->r2_null = ctx->r2_null;
  result->beta = nan_val;
  result->se = nan_val;
  result->pval = nan_val;
  result->r2_full = nan_val;
  result->r2_incr = nan_val;
  result->r2_incr_l95 = nan_val;
  result->r2_incr_u95 = nan_val;
  result->r2_liability = nan_val;
  result->r2_liability_l95 = nan_val;
  result->r2_liability_u95 = nan_val;
  result->auc_null = nan_val;
  result->auc_full = nan_val;
  result->converged = 0;

  // Standardise so beta is per standard deviation.
  double gsum = 0.0;
  for (uint32_t ii = 0; ii != sample_ct; ++ii) {
    gsum += gg[ii];
  }
  const double gmean = gsum / S_CAST(double, sample_ct);
  double gvar = 0.0;
  for (uint32_t ii = 0; ii != sample_ct; ++ii) {
    const double dev = gg[ii] - gmean;
    gvar += dev * dev;
  }
  if (!(gvar > 0.0)) {
    return;  // degenerate PGS column
  }
  const double gsd = sqrt(gvar / S_CAST(double, sample_ct - 1));

  double* e_g = &(ctx->work[2 * predictor_ct]);  // sample_ct scratch
  for (uint32_t ii = 0; ii != sample_ct; ++ii) {
    e_g[ii] = (gg[ii] - gmean) / gsd;
  }
  // Residualise in place.  Safe because Residualize reads vv fully into xtv
  // before touching resid, and skips the self-copy when the pointers match.
  Residualize(ctx, e_g, e_g);

  double gg_dot = 0.0;
  double gy_dot = 0.0;
  for (uint32_t ii = 0; ii != sample_ct; ++ii) {
    gg_dot += e_g[ii] * e_g[ii];
    gy_dot += e_g[ii] * ctx->e_y[ii];
  }
  if (!(gg_dot > 1e-10 * S_CAST(double, sample_ct - 1))) {
    return;  // PGS is collinear with the covariates
  }
  const double beta = gy_dot / gg_dot;
  const double rss_full = ctx->rss_null - gy_dot * gy_dot / gg_dot;
  if (!(rss_full > 0.0)) {
    return;  // degenerate: the PGS explains all remaining variance
  }
  const double df = S_CAST(double, sample_ct) - S_CAST(double, predictor_ct) - 1.0;
  if (!(df > 0.0)) {
    return;
  }
  const double se = sqrt((rss_full / df) / gg_dot);
  result->beta = beta;
  result->se = se;
  result->pval = PgsEvalTPval(beta / se, df);
  result->r2_full = (ctx->tss > 0.0) ? (1.0 - rss_full / ctx->tss) : 0.0;
  result->r2_incr = result->r2_full - ctx->r2_null;
  // Cross-validation replaces the R^2 family and its interval outright, so the
  // in-sample interval is not computed at all under --folds -- it would be
  // overwritten unread.  beta, se and pval above stay in-sample either way,
  // which is the point (see PgsEvalOlsSetFolds).
  if (ctx->fold_ct > 1) {
    if (PgsEvalOofR2(ctx, gg, gmean, gsd, result)) {
      return;  // collinear inside a training set; converged stays 0
    }
  } else {
    IncrR2Ci(gy_dot, gg_dot, ctx->rss_null, ctx->r2_null, sample_ct,
             predictor_ct, &result->r2_incr_l95, &result->r2_incr_u95);
    ClampCiToPoint(result->r2_incr, &result->r2_incr_l95,
                   &result->r2_incr_u95);
  }
  result->converged = 1;
}

// ---------------------------------------------------------------------------
// Binary-trait metrics (AUC, liability-scale R², qnorm)
// ---------------------------------------------------------------------------

// AUC via Mann-Whitney U statistic.
//
// Rank predicted probabilities ascending, break ties with average ranks, then
//   AUC = (sum_of_case_ranks - n_case*(n_case+1)/2) / (n_case * n_ctrl)
//
// Average-rank tie handling is required: a constant predictor assigns the same
// probability to every sample, so all ranks must average to (n+1)/2, giving
// AUC = 0.5 exactly.  Skipping the average and using raw ranks instead
// produces an incorrect value for the constant-predictor case.
//
// pp:     fitted probabilities (length sample_ct)
// yy:     0/1 labels (length sample_ct)
// buf:    scratch of length sample_ct (filled with indices, sorted in place)
// Returns AUC in [0, 1], or 0.5 when case_ct == 0 or ctrl_ct == 0.
static double AucMannWhitney(const double* pp, const double* yy,
                             uint32_t sample_ct, uint32_t case_ct,
                             uint32_t* idx_buf) {
  const uint32_t ctrl_ct = sample_ct - case_ct;
  if (case_ct == 0 || ctrl_ct == 0) {
    return 0.5;
  }

  // Fill index buffer and sort by probability ascending.
  for (uint32_t ii = 0; ii != sample_ct; ++ii) {
    idx_buf[ii] = ii;
  }
  // Insertion sort is fine for small arrays; logistic eval runs a few hundred
  // models each with at most a few thousand samples.
  for (uint32_t ii = 1; ii != sample_ct; ++ii) {
    uint32_t tmp = idx_buf[ii];
    uint32_t jj = ii;
    while (jj > 0 && pp[idx_buf[jj - 1]] > pp[tmp]) {
      idx_buf[jj] = idx_buf[jj - 1];
      --jj;
    }
    idx_buf[jj] = tmp;
  }

  // Assign average ranks within ties.  We accumulate case rank sums.
  // Rank is 1-based.
  double case_rank_sum = 0.0;
  uint32_t ii = 0;
  while (ii < sample_ct) {
    // Find the run of equal probabilities.
    const double pval = pp[idx_buf[ii]];
    uint32_t jj = ii + 1;
    while (jj < sample_ct && pp[idx_buf[jj]] == pval) {
      ++jj;
    }
    // Positions ii..jj-1 (0-based), ranks ii+1..jj (1-based).
    // Average rank = (ii+1 + jj) / 2.0 (sum / count = ((ii+1)+jj)/2).
    const double avg_rank = 0.5 * S_CAST(double, ii + 1 + jj);
    for (uint32_t kk = ii; kk != jj; ++kk) {
      if (yy[idx_buf[kk]] > 0.5) {
        case_rank_sum += avg_rank;
      }
    }
    ii = jj;
  }

  const double n_case = S_CAST(double, case_ct);
  const double n_ctrl = S_CAST(double, ctrl_ct);
  const double u_stat = case_rank_sum - n_case * (n_case + 1.0) * 0.5;
  return u_stat / (n_case * n_ctrl);
}

// ---------------------------------------------------------------------------
// qnorm: inverse normal CDF, Wichura's AS241 algorithm.
//
// Reference: Wichura MJ (1988) Algorithm AS241: The Percentage Points of the
// Normal Distribution.  Appl. Statist. 37(3): 477-484.
//
// Maximum absolute error is approximately 1.15e-9 over (0, 1).
// Verified values: qnorm(0.5) = 0, qnorm(0.975) = 1.959963984540054,
//                 qnorm(0.99) = 2.3263478740408408.
// ---------------------------------------------------------------------------
static double Qnorm(double pp) {
  static const double kA[8] = {
     3.3871328727963666080e0,  1.3314166789178437745e2,
     1.9715909503065514427e3,  1.3731693765509461125e4,
     4.5921953931549871457e4,  6.7265770927008700853e4,
     3.3430575583588128105e4,  2.5090809287301226727e3
  };
  static const double kB[8] = {
     1.0,                      4.2313330701600911252e1,
     6.8718700749205790830e2,  5.3941960214247511077e3,
     2.1213794301586595867e4,  3.9307895800092710610e4,
     2.8729085735721942674e4,  5.2264952788528545610e3
  };
  static const double kC[8] = {
     1.42343711074721209144e0,  4.63033784615654529590e0,
     5.76949722146864628717e0,  3.64784832476320460504e0,
     1.27045825245236838258e0,  2.41780725177450611770e-1,
     2.27001187347206136437e-2, 7.74545433631369706865e-4
  };
  static const double kD[8] = {
     1.0,                       2.05319162663775882187e0,
     1.67638483018615382528e0,  6.89767334985100004550e-1,
     1.48103976427480074590e-1, 1.51986665636164571966e-2,
     5.47593808499534494600e-4, 1.05075007164441684324e-9
  };
  static const double kE[8] = {
     6.65790464350110377720e0,  5.46378491116411436990e0,
     1.78482653991729133580e0,  2.96560571828504891230e-1,
     2.65321895265761230930e-2, 1.24266094738807843860e-3,
     2.71155556874348757815e-5, 2.01033439929228813265e-7
  };
  static const double kF[8] = {
     1.0,                       5.99832206555887937690e-1,
     1.36929880922735805310e-1, 1.48753612908506508940e-2,
     7.86869131145613259100e-4, 1.84631831751005468180e-5,
     1.42151175831644588870e-7, 2.04426310338993978564e-15
  };

  if (pp <= 0.0) return -1.0 / 0.0;  // -inf
  if (pp >= 1.0) return  1.0 / 0.0;  // +inf
  if (pp == 0.5) return 0.0;

  const double qq = pp - 0.5;
  double rr, num, den;
  // AS241 coefficients are listed with index 0 = constant term (lowest degree).
  // Horner evaluation must start from the highest-degree end (index 7) and
  // work down to index 0 to correctly give poly(r) = A[0] + A[1]*r + ... + A[7]*r^7.
  if (fabs(qq) <= 0.425) {
    rr = 0.180625 - qq * qq;
    num = kA[7]; for (int i = 6; i >= 0; --i) { num = num * rr + kA[i]; }
    den = kB[7]; for (int i = 6; i >= 0; --i) { den = den * rr + kB[i]; }
    return qq * num / den;
  }

  rr = sqrt(-log(qq < 0.0 ? pp : 1.0 - pp));
  double result;
  if (rr <= 5.0) {
    rr -= 1.6;
    num = kC[7]; for (int i = 6; i >= 0; --i) { num = num * rr + kC[i]; }
    den = kD[7]; for (int i = 6; i >= 0; --i) { den = den * rr + kD[i]; }
    result = num / den;
  } else {
    rr -= 5.0;
    num = kE[7]; for (int i = 6; i >= 0; --i) { num = num * rr + kE[i]; }
    den = kF[7]; for (int i = 6; i >= 0; --i) { den = den * rr + kF[i]; }
    result = num / den;
  }
  result = (qq < 0.0) ? -result : result;

  // One Newton-Raphson refinement step (Halley's method):
  //   result -= (Phi(result) - pp) / phi(result)
  // Phi(x) = 0.5 * erfc(-x / sqrt(2)), phi(x) = exp(-x^2/2) / sqrt(2*pi).
  // This brings accuracy from ~1e-8 (AS241 alone) to <1e-12.
  const double phi_r = exp(-0.5 * result * result) / sqrt(2.0 * M_PI);
  if (phi_r > 0.0) {
    const double Phi_r = 0.5 * erfc(-result / sqrt(2.0));
    result -= (Phi_r - pp) / phi_r;
  }
  return result;
}

// Liability-scale R² transform (Lee et al. 2012).
//
//   R2_liability = R2_observed * (K*(1-K))^2 / (z^2 * P*(1-P))
//
// where K = population prevalence, P = sample case proportion,
// t = qnorm(1-K), z = N(t) = exp(-t^2/2)/sqrt(2*pi).
//
// The factor (K*(1-K))^2 / (z^2 * P*(1-P)) adjusts for the difference in
// ascertainment between the study sample and the population:
//   - at P == K (population sample):  factor > 1 when K is small (inflates)
//   - at K=0.01, P=0.50 (case-control 50/50): factor < 1 (shrinks)
//
// Returns NaN if prevalence < 0 (caller signals "not given") or if K is
// outside (0,1) or if P is 0 or 1.
static double LiabilityScaleR2(double r2_observed, double prevalence,
                                double case_proportion) {
  const double kk = prevalence;
  const double pp = case_proportion;
  if (!(kk > 0.0 && kk < 1.0)) return 0.0 / 0.0;
  if (!(pp > 0.0 && pp < 1.0)) return 0.0 / 0.0;
  const double tt = Qnorm(1.0 - kk);
  // z = dnorm(t) = exp(-t^2/2) / sqrt(2*pi)
  const double z = exp(-0.5 * tt * tt) / sqrt(2.0 * M_PI);
  if (!(z > 0.0)) return 0.0 / 0.0;
  const double k1mk = kk * (1.0 - kk);
  const double p1mp = pp * (1.0 - pp);
  const double factor = (k1mk * k1mk) / (z * z * p1mp);
  return r2_observed * factor;
}

// Test-visible non-static wrappers for the metric helpers.
double AucMannWhitneyTest(const double* pp, const double* yy,
                          uint32_t sample_ct, uint32_t case_ct,
                          uint32_t* idx_buf) {
  return AucMannWhitney(pp, yy, sample_ct, case_ct, idx_buf);
}
double LiabilityScaleR2Test(double r2_observed, double prevalence,
                             double case_proportion) {
  return LiabilityScaleR2(r2_observed, prevalence, case_proportion);
}
double QnormTest(double pp) {
  return Qnorm(pp);
}
void IncrR2CiFromRTest(double rr, double r2_null, uint32_t sample_ct,
                       uint32_t predictor_ct, double* lo, double* hi) {
  IncrR2CiFromR(rr, r2_null, sample_ct, predictor_ct, lo, hi);
}

// ---------------------------------------------------------------------------
// Logistic regression: plain Newton-Raphson IRLS (binomial GLM)
// ---------------------------------------------------------------------------
//
// Uses 1/(1+exp(-x)) rather than plink2's logistic_v_unsafe, which would pull
// in all 435 lines of plink2_fmath.h.  That vectorised approximate-exp was
// built for genome-wide scans; for a few hundred small regressions the true
// exp is more accurate and costs nothing measurable.
//
// dbl_2d_buf must be sized predictor_ct * MAXV(predictor_ct, 3) doubles.
// That sizing is load-bearing: LAPACK dpocon needs a 3*N work vector inside
// InvertSymmdefMatrixChecked, and undersizing caused a heap overflow in the
// OLS path.
static BoolErr IrlsFitLogistic(const double* xx, const double* yy,
                                uint32_t sample_ct, uint32_t predictor_ct,
                                double* beta,
                                MatrixInvertBuf1* inv_buf,
                                double* dbl_2d_buf,
                                double* xtwx,    // predictor_ct^2 scratch
                                double* score,   // predictor_ct scratch
                                double* pp,      // sample_ct scratch
                                double* ww) {    // sample_ct scratch
  // Initialise coefficients to 0.
  for (uint32_t pp_idx = 0; pp_idx != predictor_ct; ++pp_idx) {
    beta[pp_idx] = 0.0;
  }

  const uint32_t kMaxIter = 25;
  const double kConvTol = 1e-9;

  for (uint32_t iter = 0; iter != kMaxIter; ++iter) {
    // Compute p = sigmoid(X beta) and W = diag(p*(1-p)).
    for (uint32_t ii = 0; ii != sample_ct; ++ii) {
      double eta = 0.0;
      for (uint32_t kk = 0; kk != predictor_ct; ++kk) {
        eta += xx[(uintptr_t)kk * sample_ct + ii] * beta[kk];
      }
      const double pii = 1.0 / (1.0 + exp(-eta));
      pp[ii] = pii;
      ww[ii] = pii * (1.0 - pii);
    }

    // Compute X'WX (symmetric, only lower triangle needed for LAPACK).
    for (uint32_t cc = 0; cc != predictor_ct; ++cc) {
      for (uint32_t rr = 0; rr <= cc; ++rr) {
        const double* col_r = &(xx[(uintptr_t)rr * sample_ct]);
        const double* col_c = &(xx[(uintptr_t)cc * sample_ct]);
        double acc = 0.0;
        for (uint32_t ii = 0; ii != sample_ct; ++ii) {
          acc += col_r[ii] * ww[ii] * col_c[ii];
        }
        xtwx[(uintptr_t)cc * predictor_ct + rr] = acc;
      }
    }

    // score = X'(y - p)
    for (uint32_t kk = 0; kk != predictor_ct; ++kk) {
      const double* col = &(xx[(uintptr_t)kk * sample_ct]);
      double acc = 0.0;
      for (uint32_t ii = 0; ii != sample_ct; ++ii) {
        acc += col[ii] * (yy[ii] - pp[ii]);
      }
      score[kk] = acc;
    }

    // Invert X'WX.  InvertSymmdefMatrixChecked fills upper triangle.
    if (InvertSymmdefMatrixChecked(predictor_ct, xtwx, inv_buf, dbl_2d_buf)) {
      return 1;  // singular or near-singular
    }
    ReflectMatrix(predictor_ct, xtwx);

    // delta = (X'WX)^-1 * score
    double max_delta = 0.0;
    for (uint32_t rr = 0; rr != predictor_ct; ++rr) {
      double acc = 0.0;
      for (uint32_t cc = 0; cc != predictor_ct; ++cc) {
        acc += xtwx[(uintptr_t)cc * predictor_ct + rr] * score[cc];
      }
      beta[rr] += acc;
      const double abs_delta = fabs(acc);
      if (abs_delta > max_delta) {
        max_delta = abs_delta;
      }
    }

    if (max_delta < kConvTol) {
      break;
    }
    if (iter == kMaxIter - 1) {
      return 1;  // did not converge
    }
  }

  return 0;
}

BoolErr PgsEvalLogisticInit(const double* yy, const double* xx,
                            uint32_t sample_ct, uint32_t predictor_ct,
                            PgsEvalLogisticCtx* ctx) {
  ctx->xx = xx;
  ctx->yy = yy;
  ctx->sample_ct = sample_ct;
  ctx->predictor_ct = predictor_ct;
  ctx->beta = nullptr;
  ctx->work = nullptr;
  // So PgsEvalLogisticCleanup, which cleans up ctx->ols unconditionally, is
  // safe on every failure return below -- including the ones that happen
  // before PgsEvalOlsInit runs.
  PgsEvalOlsPreinit(&ctx->ols);

  uint32_t case_ct = 0;
  for (uint32_t ii = 0; ii != sample_ct; ++ii) {
    if (yy[ii] > 0.5) {
      ++case_ct;
    }
  }
  ctx->case_ct = case_ct;

  const uint32_t full_predictor_ct = predictor_ct + 1;

  // The null model's coefficients, kept for the lifetime of the context.  The
  // full model is fitted per PGS column into a buffer PgsEvalLogisticOne owns,
  // so predictor_ct is the whole requirement here.
  ctx->beta = S_CAST(double*, malloc(predictor_ct * sizeof(double)));
  if (!ctx->beta) {
    PgsEvalLogisticCleanup(ctx);
    return 1;
  }

  // work layout (sized for the full model throughout):
  //   [0,                         full_predictor_ct^2)  xtwx
  //   [full_predictor_ct^2,       +full_predictor_ct)   score
  //   [+full_predictor_ct,        +sample_ct)           pp
  //   [+sample_ct,                +sample_ct)           ww
  const uintptr_t work_size =
      (uintptr_t)full_predictor_ct * full_predictor_ct  // xtwx
      + full_predictor_ct                               // score
      + 2 * sample_ct;                                  // pp + ww
  ctx->work = S_CAST(double*, malloc(work_size * sizeof(double)));
  if (!ctx->work) {
    PgsEvalLogisticCleanup(ctx);
    return 1;
  }

  // Temporary inversion buffers sized for the full model (full_predictor_ct).
  MatrixInvertBuf1* inv_buf = S_CAST(MatrixInvertBuf1*,
      malloc(kMatrixInvertBuf1CheckedAlloc * full_predictor_ct));
  // dbl_2d_buf: full_predictor_ct * MAXV(full_predictor_ct, 3) --
  // load-bearing sizing: LAPACK dpocon needs a 3*N work vector inside
  // InvertSymmdefMatrixChecked, and undersizing it is a live heap overflow
  // for predictor_ct < 3.
  double* dbl_2d_buf = S_CAST(double*, malloc(
      (uintptr_t)full_predictor_ct * MAXV(full_predictor_ct, 3) * sizeof(double)));
  if (!inv_buf || !dbl_2d_buf) {
    free(inv_buf);
    free(dbl_2d_buf);
    PgsEvalLogisticCleanup(ctx);
    return 1;
  }

  // Fit the null model (covariates only, predictor_ct columns).
  double* xtwx  = ctx->work;
  double* score = xtwx  + (uintptr_t)full_predictor_ct * full_predictor_ct;
  double* pp    = score + full_predictor_ct;
  double* ww    = pp    + sample_ct;

  if (IrlsFitLogistic(xx, yy, sample_ct, predictor_ct,
                      ctx->beta, inv_buf, dbl_2d_buf,
                      xtwx, score, pp, ww)) {
    free(inv_buf);
    free(dbl_2d_buf);
    PgsEvalLogisticCleanup(ctx);
    return 1;
  }

  // AUC of the null (covariates-only) model fitted probabilities.
  // After IrlsFitLogistic, pp[] holds the final fitted probabilities.
  // We need a temporary index buffer for AucMannWhitney.
  {
    uint32_t* idx_buf = S_CAST(uint32_t*, malloc(sample_ct * sizeof(uint32_t)));
    if (!idx_buf) {
      free(inv_buf);
      free(dbl_2d_buf);
      PgsEvalLogisticCleanup(ctx);
      return 1;
    }
    ctx->auc_null = AucMannWhitney(pp, yy, sample_ct, case_ct, idx_buf);
    free(idx_buf);
  }

  free(inv_buf);
  free(dbl_2d_buf);

  // The OLS companion, on the same 0/1 outcome and the same design matrix,
  // supplying the observed-scale R² family for every PGS column.
  //
  // Its failure is checked rather than assumed away, but not because it is
  // likely to disagree about singularity: IRLS starts at beta = 0, where
  // W = 0.25 I makes X'WX a scalar multiple of X'X, and the reciprocal-
  // condition test inside the InvertSymmdefMatrixChecked both of them call is
  // scale-invariant, so the fit above has already rejected a singular X.  What
  // is left is allocation failure.
  if (PgsEvalOlsInit(yy, xx, sample_ct, predictor_ct, &ctx->ols)) {
    PgsEvalLogisticCleanup(ctx);
    return 1;
  }
  return 0;
}

void PgsEvalLogisticCleanup(PgsEvalLogisticCtx* ctx) {
  free(ctx->beta);
  free(ctx->work);
  ctx->beta = nullptr;
  ctx->work = nullptr;
  PgsEvalOlsCleanup(&ctx->ols);
}

void PgsEvalLogisticOne(PgsEvalLogisticCtx* ctx, const double* gg,
                        double prevalence, PgsEvalResult* result) {
  const uint32_t sample_ct    = ctx->sample_ct;
  const uint32_t predictor_ct = ctx->predictor_ct;  // null-model predictor count
  const uint32_t full_predictor_ct = predictor_ct + 1;  // adds the PGS column
  const double nan_val = 0.0 / 0.0;

  result->beta         = nan_val;
  result->se           = nan_val;
  result->pval         = nan_val;
  result->r2_null      = nan_val;
  result->r2_full      = nan_val;
  result->r2_incr      = nan_val;
  result->r2_incr_l95  = nan_val;
  result->r2_incr_u95  = nan_val;
  result->r2_liability = nan_val;
  result->r2_liability_l95 = nan_val;
  result->r2_liability_u95 = nan_val;
  result->auc_null     = nan_val;
  result->auc_full     = nan_val;
  result->converged    = 0;

  if (sample_ct <= full_predictor_ct + 1) {
    return;  // not enough degrees of freedom
  }

  // ---------------------------------------------------------------------------
  // Observed-scale R², from OLS on the 0/1 outcome.  Frisch-Waugh-Lovell on an
  // X'X that PgsEvalOlsInit already factored, so it is cheap; run it first
  // because r2_incr is the sort key and a row that cannot produce one is
  // useless however well the logistic fit goes.
  //
  // *result goes straight in rather than into a local whose fields are then
  // copied across: PgsEvalOlsOne writes the struct exhaustively, so this stays
  // right when a field is added to PgsEvalResult.
  // ---------------------------------------------------------------------------
  PgsEvalOlsOne(&ctx->ols, gg, result);
  if (!result->converged) {
    return;
  }
  // OLS converged, the logistic fit has not been attempted yet.  Every early
  // return below must still report failure, so drop the flag back to 0 -- the
  // OLS beta/se/pval it leaves behind are never read while converged is 0.
  result->converged = 0;

  // Standardise the PGS so beta is per standard deviation.  PgsEvalOlsOne
  // rejected a zero-variance column with this identical computation, so gsd > 0.
  double gsum = 0.0;
  for (uint32_t ii = 0; ii != sample_ct; ++ii) {
    gsum += gg[ii];
  }
  const double gmean = gsum / S_CAST(double, sample_ct);
  double gvar = 0.0;
  for (uint32_t ii = 0; ii != sample_ct; ++ii) {
    const double dev = gg[ii] - gmean;
    gvar += dev * dev;
  }
  const double gsd = sqrt(gvar / S_CAST(double, sample_ct - 1));

  // Build the full design matrix: [null-model columns | standardised PGS].
  // Column-major, sample_ct * full_predictor_ct.
  double* xx_full = S_CAST(double*, malloc(
      (uintptr_t)sample_ct * full_predictor_ct * sizeof(double)));
  if (!xx_full) {
    return;
  }
  // Copy null-model columns.
  memcpy(xx_full, ctx->xx, (uintptr_t)sample_ct * predictor_ct * sizeof(double));
  // Append standardised PGS as the last column.
  double* pgs_col = xx_full + (uintptr_t)predictor_ct * sample_ct;
  for (uint32_t ii = 0; ii != sample_ct; ++ii) {
    pgs_col[ii] = (gg[ii] - gmean) / gsd;
  }

  // Per-call inversion buffers sized for the full model.
  // dbl_2d_buf must be full_predictor_ct * MAXV(full_predictor_ct, 3): the
  // MAXV(.,3) is load-bearing because dpocon inside InvertSymmdefMatrixChecked
  // needs a 3*N work vector.  Getting this wrong is a heap overflow.
  MatrixInvertBuf1* inv_buf = S_CAST(MatrixInvertBuf1*,
      malloc(kMatrixInvertBuf1CheckedAlloc * full_predictor_ct));
  double* dbl_2d_buf = S_CAST(double*, malloc(
      (uintptr_t)full_predictor_ct * MAXV(full_predictor_ct, 3) * sizeof(double)));
  // beta_full: full_predictor_ct coefficients for the full model.
  double* beta_full = S_CAST(double*, malloc(full_predictor_ct * sizeof(double)));

  if (!inv_buf || !dbl_2d_buf || !beta_full) {
    free(xx_full); free(inv_buf); free(dbl_2d_buf); free(beta_full);
    return;
  }

  // Reuse ctx->work for xtwx/score/pp/ww scratch (sized for full model in Init).
  double* xtwx  = ctx->work;
  double* score = xtwx  + (uintptr_t)full_predictor_ct * full_predictor_ct;
  double* pp    = score + full_predictor_ct;
  double* ww    = pp    + sample_ct;

  // Under complete separation, fitted probabilities saturate at 0 and 1,
  // X'WX becomes singular (W -> 0), and InvertSymmdefMatrixChecked returns 1,
  // or IRLS fails to meet the convergence criterion within 25 iterations.
  // Either outcome makes IrlsFitLogistic return 1.  We report converged = 0
  // rather than returning a garbage coefficient.
  if (IrlsFitLogistic(xx_full, ctx->yy, sample_ct, full_predictor_ct,
                      beta_full, inv_buf, dbl_2d_buf,
                      xtwx, score, pp, ww)) {
    free(xx_full); free(inv_buf); free(dbl_2d_buf); free(beta_full);
    return;  // non-convergence (e.g. separation); converged stays 0
  }

  // SE of the PGS coefficient = sqrt( (X'WX)^-1 [pgs, pgs] ).
  // After the final inversion inside IrlsFitLogistic + ReflectMatrix,
  // xtwx holds (X'WX)^-1.
  const uint32_t pgs_idx = full_predictor_ct - 1;
  const double var_beta = xtwx[(uintptr_t)pgs_idx * full_predictor_ct + pgs_idx];
  if (!(var_beta > 0.0)) {
    free(xx_full); free(inv_buf); free(dbl_2d_buf); free(beta_full);
    return;
  }
  const double se      = sqrt(var_beta);
  const double beta_pgs = beta_full[pgs_idx];

  // Two-sided Wald p-value: erfc(|z|/sqrt(2)) == 2*Phi(-|z|) for a standard
  // normal.  This is correct for a plain binomial GLM -- it is exactly what
  // R's summary.glm reports for family=binomial.  Do not replace with a
  // likelihood-ratio test: LRT is only needed when Firth is active (Wald
  // degrades under separation via the Hauck-Donner effect), and Firth has been
  // removed from this path (recoverable from git eeae1e2 if ever re-added).
  const double zz   = beta_pgs / se;
  const double pval = erfc(fabs(zz) / sqrt(2.0));

  // ---------------------------------------------------------------------------
  // AUC of the full model predicted probabilities.
  // pp[] now holds the fitted probabilities from IrlsFitLogistic (full model).
  // ---------------------------------------------------------------------------
  double auc_full;
  {
    uint32_t* idx_buf = S_CAST(uint32_t*, malloc(sample_ct * sizeof(uint32_t)));
    if (!idx_buf) {
      free(xx_full); free(inv_buf); free(dbl_2d_buf); free(beta_full);
      return;
    }
    auc_full = AucMannWhitney(pp, ctx->yy, sample_ct, ctx->case_ct, idx_buf);
    free(idx_buf);
  }

  // ---------------------------------------------------------------------------
  // Liability-scale R² (Lee et al. 2012), applied to the observed-scale
  // incremental R² -- the scale the transform is derived for.
  // Returns NaN when prevalence < 0 (not provided) or if inputs are degenerate.
  // ---------------------------------------------------------------------------
  const double case_proportion = S_CAST(double, ctx->case_ct) /
                                  S_CAST(double, sample_ct);
  const double r2_liability = LiabilityScaleR2(result->r2_incr, prevalence,
                                                case_proportion);
  // Given K and P the transform is multiplication by a positive constant, so it
  // carries the observed-scale interval across endpoint by endpoint and keeps
  // them in order.  No separate interval theory is needed on this scale.
  const double r2_liability_l95 = LiabilityScaleR2(result->r2_incr_l95,
                                                   prevalence, case_proportion);
  const double r2_liability_u95 = LiabilityScaleR2(result->r2_incr_u95,
                                                   prevalence, case_proportion);

  // Only the logistic-owned fields; the R² family is already in place from
  // PgsEvalOlsOne above.
  result->beta         = beta_pgs;
  result->se           = se;
  result->pval         = pval;
  result->r2_liability = r2_liability;
  result->r2_liability_l95 = r2_liability_l95;
  result->r2_liability_u95 = r2_liability_u95;
  result->auc_null     = ctx->auc_null;
  result->auc_full     = auc_full;
  result->converged    = 1;

  free(xx_full); free(inv_buf); free(dbl_2d_buf); free(beta_full);
}

// ---------------------------------------------------------------------------
// BestPgsReport: read back .sscore, rank by incremental R², write .best_pgs
// ---------------------------------------------------------------------------

// The .best_pgs column layout, in one place.  21 columns:
//   1 PGS_ID  2 PHENO  3 TYPE  4 N  5 N_CASE  6 FOLDS
//   7 BETA  8 OR  9 SE  10 P  11 R2_NULL  12 R2_FULL
//   13 R2_INCR  14 R2_INCR_L95  15 R2_INCR_U95
//   16 R2L_INCR  17 R2L_INCR_L95  18 R2L_INCR_U95
//   19 AUC_NULL  20 AUC_FULL  21 AUC_INCR
// pgs_eval_check.sh Case 0 asserts this string exactly, and every other case
// there indexes by column number, so the two move together or not at all.
static const char kBestPgsHeader[] =
    "#PGS_ID\tPHENO\tTYPE\tN\tN_CASE\tFOLDS\tBETA\tOR\tSE\tP\tR2_NULL\tR2_FULL\t"
    "R2_INCR\tR2_INCR_L95\tR2_INCR_U95\tR2L_INCR\tR2L_INCR_L95\tR2L_INCR_U95\t"
    "AUC_NULL\tAUC_FULL\tAUC_INCR\n";

// Columns 7..21, i.e. everything an unconverged row has nothing to say about.
// The literal run of "NA\t" this count replaces was one short for the whole
// life of the 16-column layout: an unconverged row emitted 15 fields against 16
// names and nothing caught it, because no case in the suite produced one.  One
// now does, and asserts the field count.
static const uint32_t kBestPgsNaFieldCt = 15;

static char* WriteNaRun(uint32_t field_ct, char* wp) {
  for (uint32_t ii = 0; ii != field_ct; ++ii) {
    if (ii) { *wp++ = '\t'; }
    wp = strcpya_k(wp, "NA");
  }
  return wp;
}

// dtoa_g renders a NaN as "nan"; every other absent value in this file is "NA".
static char* WriteDoubleOrNa(double dxx, char* wp) {
  if (isnan(dxx)) {
    return strcpya_k(wp, "NA");
  }
  return dtoa_g(dxx, wp);
}

// qsort comparator: sort by r2_incr descending; unconverged rows go last.
typedef struct BestPgsRowStruct {
  PgsEvalResult result;
  // model name length is bounded; store a pointer into a separately-managed
  // name buffer (valid for the lifetime of the call).
  char* name;
} BestPgsRow;

static int BestPgsRowCmp(const void* aa, const void* bb) {
  const BestPgsRow* ra = S_CAST(const BestPgsRow*, aa);
  const BestPgsRow* rb = S_CAST(const BestPgsRow*, bb);
  // unconverged rows always go last
  if (ra->result.converged != rb->result.converged) {
    return rb->result.converged - ra->result.converged;
  }
  // sort by r2_incr descending; NaN comparisons are false so they end up last
  if (ra->result.r2_incr > rb->result.r2_incr) return -1;
  if (ra->result.r2_incr < rb->result.r2_incr) return  1;
  return 0;
}

// PgsEvalOlsSetFolds fails identically on the QT and CC paths -- only the
// context to release first differs -- so the message lives here rather than
// once per path.  It was deduplicated once before and drifted back in.
static PglErr FoldNullFitFailed(const char* pheno_name, uint32_t fold_ct) {
  logerrprintfww("Error: --give-best-pgs: phenotype '%s': could not fit the null model on all %u training sets. Reduce --folds, or drop a covariate that is constant within a fold.\n",
                 pheno_name, fold_ct);
  return kPglRetInconsistentInput;
}

// ---------------------------------------------------------------------------
// FitWritePhenoBlock: shared inner loop — constant-covariate drop, null model
// fit, per-score fits, sort, and write one phenotype group to outfile.
//
// score_block: column-major [score_col_ct x complete_ct], already extracted
//   to the complete-case sample set.
// X_complete: column-major [complete_ct_max x (1+covar_ct)]; columns
//   0..covar_ct are filled with intercept + covar values for complete cases.
//   The function may compact this array in place (constant-covariate drop).
// covar_names_arr: array of covar_ct C strings (may be nullptr if covar_ct=0).
// fold_ct / fold_seed: --folds and --seed; fold_ct == 1 means no
//   cross-validation, and is the only value that reproduces the pre-CV output.
//
// keep_idxs / keep_ct: --pgs-pheno-pairs.  keep_idxs == nullptr means fit every
//   score column, which is the unfiltered behaviour; otherwise keep_idxs[0..
//   keep_ct) are the score-column indices this phenotype was paired with, in
//   score-column order.  keep_ct == 0 means this phenotype appears in no pair,
//   and the whole block -- null model included -- is skipped: fitting a null
//   model for a phenotype nobody asked about can fail on its own (too few
//   complete cases, collinear covariates) and take out a run that was never
//   going to report it.
//
// Returns kPglRetSuccess on success; opens *outfile_ptr and writes the header
// on the first block that has anything to write.  That is keyed on the file
// still being closed rather than on a phenotype index, because under a filter
// the first phenotype is not necessarily the first one to produce rows.
// pheno_rows[0..fit_ct) are overwritten with the sorted results.
static PglErr FitWritePhenoBlock(uint32_t complete_ct, uint32_t covar_ct,
                                  double* y_complete, double* X_complete,
                                  const double* score_block,
                                  char** score_col_names, uint32_t score_col_ct,
                                  const char* const* covar_names_arr,
                                  const char* pheno_name, PhenoDtype type_code,
                                  double prevalence, uint32_t fold_ct,
                                  uint64_t fold_seed,
                                  const uint32_t* keep_idxs, uint32_t keep_ct,
                                  FILE** outfile_ptr,
                                  char* outname, char* outname_end,
                                  BestPgsRow* pheno_rows) {
  PglErr reterr = kPglRetSuccess;
  const uint32_t fit_ct = keep_idxs? keep_ct : score_col_ct;
  if (!fit_ct) {
    return kPglRetSuccess;
  }

  if (prevalence >= 0.0 && type_code == kPhenoDtypeQt) {
    logprintfww("Warning: --prevalence has no effect on quantitative phenotype '%s' and will be ignored.\n",
                pheno_name);
  }

  uint32_t predictor_ct = 1 + covar_ct;  // intercept + covars; may shrink

  // --- Constant-covariate drop (complete-case dependent). ---
  {
    uint32_t dst_col = 1;
    for (uint32_t ci = 0; ci != covar_ct; ++ci) {
      const double* src = &X_complete[S_CAST(uint64_t, 1 + ci) * complete_ct];
      double sum = 0.0;
      for (uint32_t ii = 0; ii != complete_ct; ++ii) {
        sum += src[ii];
      }
      const double mean = sum / S_CAST(double, complete_ct);
      double ss = 0.0;
      for (uint32_t ii = 0; ii != complete_ct; ++ii) {
        const double dev = src[ii] - mean;
        ss += dev * dev;
      }
      if (!(ss > 0.0)) {
        const char* cname = covar_names_arr ? covar_names_arr[ci] : "(unknown)";
        logerrprintfww("Warning: --give-best-pgs: covariate '%s' is constant across the %u complete-case samples for phenotype '%s'; dropping it.\n",
                       cname, complete_ct, pheno_name);
        continue;
      }
      if (dst_col != 1 + ci) {
        double* dst = &X_complete[S_CAST(uint64_t, dst_col) * complete_ct];
        memcpy(dst, src, complete_ct * sizeof(double));
      }
      ++dst_col;
    }
    predictor_ct = dst_col;
  }

  // Re-check minimum sample-size after dropping constant covariates.
  if (complete_ct <= predictor_ct + 1) {
    logerrprintfww("Error: --give-best-pgs: phenotype '%s': too few complete-case observations (%u) for %u predictors (after dropping constant covariates).\n",
                   pheno_name, complete_ct, predictor_ct);
    return kPglRetInconsistentInput;
  }

  // --- Fit null model and all score columns. ---
  if (type_code == kPhenoDtypeQt) {
    // OLS path (Frisch-Waugh-Lovell).
    PgsEvalOlsCtx ctx;
    if (PgsEvalOlsInit(y_complete, X_complete, complete_ct, predictor_ct, &ctx)) {
      logerrprintfww("Error: --give-best-pgs: phenotype '%s': covariate matrix is singular even after dropping constant covariates. Two or more covariates are probably collinear.\n",
                     pheno_name);
      return kPglRetInconsistentInput;
    }
    // No-op at fold_ct == 1, so the k = 1 path is untouched by this call.
    if (PgsEvalOlsSetFolds(&ctx, fold_ct, 0, fold_seed)) {
      PgsEvalOlsCleanup(&ctx);
      return FoldNullFitFailed(pheno_name, fold_ct);
    }

    for (uint32_t ki = 0; ki != fit_ct; ++ki) {
      const uint32_t sc = keep_idxs? keep_idxs[ki] : ki;
      pheno_rows[ki].name = score_col_names[sc];
      // score_block is already the complete-case selection, and PgsEvalOlsOne
      // does not write through its gg pointer.
      PgsEvalOlsOne(&ctx, &score_block[S_CAST(uint64_t, sc) * complete_ct],
                    &pheno_rows[ki].result);
    }

    PgsEvalOlsCleanup(&ctx);

  } else {
    // Logistic path (case/control).  Two models are fitted on the same
    // complete-case data: logistic for BETA/OR/SE/P/AUC, and OLS on the 0/1
    // outcome for the R² family.  y_complete is already 0/1 doubles here.
    // PgsEvalLogisticInit builds both from this one design matrix, so they
    // cannot be handed different data.
    PgsEvalLogisticCtx lctx;
    if (PgsEvalLogisticInit(y_complete, X_complete, complete_ct, predictor_ct, &lctx)) {
      logerrprintfww("Error: --give-best-pgs: phenotype '%s': the covariates-only null model could not be fitted. Two or more covariates are probably collinear, which dropping constant ones does not catch; allocation failure is the other possibility.\n",
                     pheno_name);
      return kPglRetInconsistentInput;
    }
    // Cross-validation lives entirely in the OLS companion, which is where the
    // binary path's R^2 family already comes from; stratify on the same 0/1
    // outcome that companion was built from.  AUC is not cross-validated -- it
    // comes from the logistic fit, which is not refit per fold.
    if (PgsEvalOlsSetFolds(&lctx.ols, fold_ct, 1, fold_seed)) {
      PgsEvalLogisticCleanup(&lctx);
      return FoldNullFitFailed(pheno_name, fold_ct);
    }

    for (uint32_t ki = 0; ki != fit_ct; ++ki) {
      const uint32_t sc = keep_idxs? keep_idxs[ki] : ki;
      pheno_rows[ki].name = score_col_names[sc];
      PgsEvalLogisticOne(&lctx, &score_block[S_CAST(uint64_t, sc) * complete_ct],
                         prevalence, &pheno_rows[ki].result);
    }

    PgsEvalLogisticCleanup(&lctx);
  }

  // --- Sort within this phenotype group by R2_INCR desc, unconverged last. ---
  qsort(pheno_rows, fit_ct, sizeof(BestPgsRow), BestPgsRowCmp);

  // The output's shape -- one group per phenotype, best first -- says "this
  // one", and nothing in it can say "none of them".  Row 0 now holds the
  // largest R2_INCR, so one comparison covers the group; !(x > 0) rather than
  // x <= 0 so an all-unconverged group, whose R2_INCR is NaN, warns too.
  if (!(pheno_rows[0].result.r2_incr > 0.0)) {
    logerrprintfww("Warning: --give-best-pgs: no model has a positive R2_INCR for phenotype '%s'; the first row is the least bad of them, not evidence that it predicts.\n",
                   pheno_name);
  }

  // --- Open output file and write unified header on the first block. ---
  if (!*outfile_ptr) {
    OutnameZstSet(".best_pgs", 0, outname_end);
    if (unlikely(fopen_checked(outname, FOPEN_WB, outfile_ptr))) {
      return kPglRetOpenFail;
    }

    const uint32_t hdr_len = S_CAST(uint32_t, sizeof(kBestPgsHeader) - 1);
    if (unlikely(fwrite(kBestPgsHeader, hdr_len, 1, *outfile_ptr) != 1)) {
      return kPglRetWriteFail;
    }
  }

  // --- Write rows for this phenotype. ---
  // Count n_case from y_complete (needed for CC rows).
  uint32_t n_case_out = 0;
  if (type_code == kPhenoDtypeCc) {
    for (uint32_t ii = 0; ii != complete_ct; ++ii) {
      if (y_complete[ii] > 0.5) ++n_case_out;
    }
  }

  {
    char wbuf[4096];
    for (uint32_t ki = 0; ki != fit_ct; ++ki) {
      const BestPgsRow* rr = &pheno_rows[ki];
      char* wp = wbuf;

      // PGS_ID
      wp = strcpya(wp, rr->name);
      *wp++ = '\t';
      // PHENO
      wp = strcpya(wp, pheno_name);
      *wp++ = '\t';

      if (type_code == kPhenoDtypeCc) {
        wp = strcpya_k(wp, "CC\t");
        // N
        wp = u32toa(complete_ct, wp); *wp++ = '\t';
        // N_CASE
        wp = u32toa(n_case_out, wp); *wp++ = '\t';
        // FOLDS: 1 means the model was fitted and scored on the same samples.
        wp = u32toa(fold_ct, wp); *wp++ = '\t';

        if (!rr->result.converged) {
          wp = WriteNaRun(kBestPgsNaFieldCt, wp);
        } else {
          // BETA (log-odds)
          wp = dtoa_g(rr->result.beta, wp);         *wp++ = '\t';
          // OR = exp(BETA)
          wp = dtoa_g(exp(rr->result.beta), wp);    *wp++ = '\t';
          // SE (log-odds scale)
          wp = dtoa_g(rr->result.se, wp);           *wp++ = '\t';
          // P
          wp = dtoa_g(rr->result.pval, wp);         *wp++ = '\t';
          // R2_NULL (observed scale, OLS on the 0/1 outcome)
          wp = dtoa_g(rr->result.r2_null, wp);      *wp++ = '\t';
          // R2_FULL (observed scale, OLS on the 0/1 outcome)
          wp = dtoa_g(rr->result.r2_full, wp);      *wp++ = '\t';
          // R2_INCR, then its 95% CI.  The CI is NA in the
          // n == predictor_ct + 2 corner, where the fit has a degree of freedom
          // left but the Fisher z standard error does not, and under
          // cross-validation when R2_INCR came out negative.
          wp = dtoa_g(rr->result.r2_incr, wp);      *wp++ = '\t';
          wp = WriteDoubleOrNa(rr->result.r2_incr_l95, wp); *wp++ = '\t';
          wp = WriteDoubleOrNa(rr->result.r2_incr_u95, wp); *wp++ = '\t';
          // R2L_INCR and its CI; all three are NA together, without --prevalence.
          wp = WriteDoubleOrNa(rr->result.r2_liability, wp);     *wp++ = '\t';
          wp = WriteDoubleOrNa(rr->result.r2_liability_l95, wp); *wp++ = '\t';
          wp = WriteDoubleOrNa(rr->result.r2_liability_u95, wp); *wp++ = '\t';
          // AUC_NULL
          wp = dtoa_g(rr->result.auc_null, wp);     *wp++ = '\t';
          // AUC_FULL
          wp = dtoa_g(rr->result.auc_full, wp);     *wp++ = '\t';
          // AUC_INCR
          wp = dtoa_g(rr->result.auc_full - rr->result.auc_null, wp);
        }
      } else {
        // QT rows: NA in N_CASE, OR, the R2L family, and AUC_*.
        wp = strcpya_k(wp, "QT\t");
        // N
        wp = u32toa(complete_ct, wp); *wp++ = '\t';
        // N_CASE = NA
        wp = strcpya_k(wp, "NA\t");
        // FOLDS
        wp = u32toa(fold_ct, wp); *wp++ = '\t';

        if (!rr->result.converged) {
          wp = WriteNaRun(kBestPgsNaFieldCt, wp);
        } else {
          // BETA
          wp = dtoa_g(rr->result.beta, wp);    *wp++ = '\t';
          // OR = NA
          wp = strcpya_k(wp, "NA\t");
          // SE
          wp = dtoa_g(rr->result.se, wp);      *wp++ = '\t';
          // P
          wp = dtoa_g(rr->result.pval, wp);    *wp++ = '\t';
          // R2_NULL
          wp = dtoa_g(rr->result.r2_null, wp); *wp++ = '\t';
          // R2_FULL
          wp = dtoa_g(rr->result.r2_full, wp); *wp++ = '\t';
          // R2_INCR and its 95% CI
          wp = dtoa_g(rr->result.r2_incr, wp); *wp++ = '\t';
          wp = WriteDoubleOrNa(rr->result.r2_incr_l95, wp); *wp++ = '\t';
          wp = WriteDoubleOrNa(rr->result.r2_incr_u95, wp); *wp++ = '\t';
          // R2L_INCR, R2L_INCR_L95, R2L_INCR_U95, AUC_NULL, AUC_FULL, AUC_INCR
          wp = WriteNaRun(6, wp);
        }
      }

      *wp++ = '\n';

      const uint32_t line_len = S_CAST(uint32_t, wp - wbuf);
      if (unlikely(fwrite(wbuf, line_len, 1, *outfile_ptr) != 1)) {
        return kPglRetWriteFail;
      }
    }
  }

  return reterr;
}

// --folds against one phenotype's complete-case count: k must be at least 1,
// leave at least two samples per fold, and leave the smallest training set --
// complete_ct minus ceil(complete_ct / k) -- with more than predictor_ct + 1
// samples, which is what PgsEvalOlsSetFolds enforces once fitting has started.
//
// A hard error, not a skip: k is something the user asked for outright, so a k
// the data cannot honour is a misconfiguration rather than a property of that
// phenotype.  Both input paths call it, so there is one wording for both.
//
// Call it only on a phenotype already known to be usable.  Before it moved
// below the case/control and constant-phenotype checks, a phenotype that would
// have been skipped with a warning aborted the whole run instead.
//
// predictor_ct is the pre-drop count, intercept plus every covariate.  Constant
// covariates are dropped per phenotype inside FitWritePhenoBlock, which can only
// lower it, so this can be up to one covariate stricter than the fit it guards
// -- the same pre-drop count the degrees-of-freedom skip in BestPgsValidate
// already reasons with.  Strict is the right direction: the reason this
// function exists is to fail before the scoring pass rather than after it.
static PglErr BestPgsCheckFolds(uint32_t fold_ct, uint32_t complete_ct,
                                uint32_t predictor_ct, const char* pheno_name) {
  // The flag parser rejects 0 outright, so this is the unit's own precondition
  // rather than the path a user reaches.  Worth keeping because the failure it
  // prevents is silent: PgsEvalOlsSetFolds treats fold_ct <= 1 as "no
  // cross-validation", so a 0 would produce in-sample numbers under a FOLDS
  // column reading 0.
  if (fold_ct < 1) {
    logerrputs("Error: --give-best-pgs: --folds must be at least 1.\n");
    return kPglRetInconsistentInput;
  }
  // Two samples per fold is the floor at which a held-out fold has any spread
  // to be scored on at all.
  if (fold_ct > complete_ct / 2) {
    logerrprintfww("Error: --give-best-pgs: phenotype '%s': --folds %u needs at least %u complete-case observations, and there are %u.\n",
                   pheno_name, fold_ct, 2 * fold_ct, complete_ct);
    return kPglRetInconsistentInput;
  }
  // The same constraint PgsEvalOlsSetFolds enforces per training set, but
  // computed up front, where it is just as computable: the largest held-out
  // fold takes ceil(complete_ct / fold_ct) samples, so the smallest training
  // set gets the rest, and fitting needs more than predictor_ct + 1 of them.
  // Checking only the fold-size floor above let 10 samples with 3 covariates
  // and --folds 2 through, to fail later inside the fitter -- by which point
  // ScoreReport has written and discarded a whole .sscore, and on a
  // multi-phenotype run earlier phenotypes' rows are already in the .best_pgs.
  if (fold_ct > 1) {
    const uint32_t min_train_ct =
        complete_ct - (complete_ct + fold_ct - 1) / fold_ct;
    if (min_train_ct <= predictor_ct + 1) {
      // Naming predictor_ct here would name a count that is not what gets
      // fitted -- constant covariates are dropped later, per phenotype -- and
      // "drop a covariate" would be advice about work the tool already does
      // for itself.  State the check instead, caveat included.
      logerrprintfww("Error: --give-best-pgs: phenotype '%s': --folds %u leaves only %u samples in the smallest training set, and fitting the intercept, %u covariate(s) and the PGS there needs at least %u. That covariate count is taken before constant covariates are dropped, so this can refuse a k the fit would have managed; the check is deliberately strict because it runs before the scoring pass. Reduce --folds.\n",
                     pheno_name, fold_ct, min_train_ct, predictor_ct - 1,
                     predictor_ct + 2);
      return kPglRetInconsistentInput;
    }
  }
  return kPglRetSuccess;
}

// ---------------------------------------------------------------------------
// One phenotype's complete-case data, in the form every guard wants it: dense,
// route-independent, already selected down to the samples whose phenotype value
// and every covariate value are present.
//
// The two input routes arrive here from opposite directions -- --pfile walks
// PhenoCol bitvectors and is told the type, --sscore walks joined arrays of
// doubles and infers it -- and used to carry a private copy of every check
// downstream of that.  The copies drifted: --sscore shipped with no
// constant-phenotype guard, so a zero-variance phenotype got a confidence
// interval on an R2 that did not exist, and with missing values coerced to
// controls, which inflated N and later stratified folds on fabricated labels.
// Both were fixed one route at a time, after shipping.  One block and one
// validator is what stops the next one.
//
// There is deliberately no X here.  None of the guards reads the design matrix,
// and FitWritePhenoBlock -- which does -- is already shared by both routes, so
// carrying X would add a member nothing reads.
typedef struct {
  const double* y;    // [n]; 0/1 when type is kPhenoDtypeCc
  uint32_t n;         // complete-case count
  uint32_t covar_ct;
  PhenoDtype type;    // kPhenoDtypeQt or kPhenoDtypeCc
  const char* name;
  const char* route;  // "--give-best-pgs" or "--sscore": the flag the user
                      // typed.  The two routes' messages were identical apart
                      // from this prefix, so carrying it lets one validator
                      // keep both wordings verbatim.
} PgsPhenoBlock;

// Every check a phenotype must pass before it is worth fitting: enough complete
// cases for the predictors, a usable case/control split, a quantitative
// phenotype that actually varies, and a --folds the data can honour.  Both
// routes call this, so neither can acquire a guard the other lacks.
//
// *usable_ptr is 0 when the phenotype was skipped -- the warning is already
// printed and the caller should move to the next one.  The return value is
// non-zero only for the --folds failure, which is a misconfiguration the user
// asked for outright and takes the whole run down.
//
// The --folds check comes last on purpose.  Above the skips, a phenotype that
// would merely have been dropped with a warning aborted the run instead.
static PglErr ValidatePhenoBlock(const PgsPhenoBlock* pb, uint32_t fold_ct,
                                 uint32_t* usable_ptr) {
  *usable_ptr = 0;
  const uint32_t predictor_ct = 1 + pb->covar_ct;

  if (!pb->n) {
    logerrprintfww("Warning: %s: phenotype '%s': no samples have both a non-missing phenotype and non-missing values for every covariate; skipping.\n",
                   pb->route, pb->name);
    return kPglRetSuccess;
  }
  if (pb->n <= predictor_ct + 1) {
    logerrprintfww("Warning: %s: phenotype '%s': too few complete-case observations (%u) for %u predictors; skipping.\n",
                   pb->route, pb->name, pb->n, predictor_ct);
    return kPglRetSuccess;
  }

  if (pb->type == kPhenoDtypeCc) {
    // y is 0/1 here, so this counts exactly what the case bitvector did.
    uint32_t case_ct = 0;
    for (uint32_t ii = 0; ii != pb->n; ++ii) {
      if (pb->y[ii] > 0.5) {
        ++case_ct;
      }
    }
    if (!case_ct) {
      logerrprintfww("Warning: %s: phenotype '%s': no cases in the complete-case sample set; skipping.\n",
                     pb->route, pb->name);
      return kPglRetSuccess;
    }
    if (case_ct == pb->n) {
      logerrprintfww("Warning: %s: phenotype '%s': no controls in the complete-case sample set; skipping.\n",
                     pb->route, pb->name);
      return kPglRetSuccess;
    }
  } else {
    // A constant phenotype has no variance for a PGS to explain, and the
    // fitters do not refuse it: PgsEvalOlsInit's `tss > 0` guards make R2_NULL,
    // R2_FULL and R2_INCR all exactly 0, and the CI is then a Fisher z
    // half-width at r = 0 -- a plausible-looking bound on an R2 that does not
    // exist.  The --sscore route printed exactly that until bdbeaca.
    double ysum = 0.0;
    for (uint32_t ii = 0; ii != pb->n; ++ii) {
      ysum += pb->y[ii];
    }
    const double ymean = ysum / S_CAST(double, pb->n);
    double yss = 0.0;
    for (uint32_t ii = 0; ii != pb->n; ++ii) {
      const double dev = pb->y[ii] - ymean;
      yss += dev * dev;
    }
    if (!(yss > 0.0)) {
      logerrprintfww("Warning: %s: phenotype '%s' is constant; skipping.\n",
                     pb->route, pb->name);
      return kPglRetSuccess;
    }
  }

  const PglErr fold_err = BestPgsCheckFolds(fold_ct, pb->n, predictor_ct,
                                            pb->name);
  if (unlikely(fold_err)) {
    return fold_err;
  }

  *usable_ptr = 1;
  return kPglRetSuccess;
}

// The --pfile route's producer: select one phenotype's complete cases out of
// the PhenoCol bitvectors and write them densely.  Returns the count.
//
// y receives the outcome (0/1 for case/control) and must have room for
// sample_ct doubles.  keep_ri, when non-null, receives the .sscore row index of
// each complete case, which is how the caller selects the score matrix down to
// the same samples.  X, when non-null, must have room for
// sample_ct * (1 + covar_ct) doubles and is filled column-major with the
// intercept and the covariates.
//
// The validation pass passes null for both: its guards read only y, and it runs
// before the scores exist.  This one walk replaces the four identical
// bitvector walks BestPgsValidate used to make per phenotype.
static uint32_t PfileFillPhenoBlock(const uintptr_t* sample_include,
                                    uint32_t sample_ct,
                                    const PhenoCol* the_pheno,
                                    const PhenoCol* covar_cols,
                                    uint32_t covar_ct,
                                    double* y, double* X, uint32_t* keep_ri) {
  const uint32_t is_cc = (the_pheno->type_code == kPhenoDtypeCc);
  uint32_t complete_ct = 0;
  {
    uintptr_t si_base = 0;
    uintptr_t si_bits = sample_include[0];
    for (uint32_t ri = 0; ri != sample_ct; ++ri) {
      const uintptr_t suidx = BitIter1(sample_include, &si_base, &si_bits);
      if (!IsSet(the_pheno->nonmiss, suidx)) continue;
      uint32_t all_covar_ok = 1;
      for (uint32_t ci = 0; ci != covar_ct; ++ci) {
        if (!IsSet(covar_cols[ci].nonmiss, suidx)) {
          all_covar_ok = 0;
          break;
        }
      }
      if (!all_covar_ok) continue;
      y[complete_ct] = is_cc? (IsSet(the_pheno->data.cc, suidx)? 1.0 : 0.0)
                            : the_pheno->data.qt[suidx];
      if (keep_ri) {
        keep_ri[complete_ct] = ri;
      }
      ++complete_ct;
    }
  }

  // X's column stride is the complete-case count, so it cannot be filled until
  // that count is known -- hence the second walk rather than one fused pass.
  if (X) {
    for (uint32_t ii = 0; ii != complete_ct; ++ii) {
      X[ii] = 1.0;  // intercept
    }
    uint32_t out_idx = 0;
    uintptr_t si_base = 0;
    uintptr_t si_bits = sample_include[0];
    for (uint32_t ri = 0; ri != sample_ct; ++ri) {
      const uintptr_t suidx = BitIter1(sample_include, &si_base, &si_bits);
      if (!IsSet(the_pheno->nonmiss, suidx)) continue;
      uint32_t all_covar_ok = 1;
      for (uint32_t ci = 0; ci != covar_ct; ++ci) {
        if (!IsSet(covar_cols[ci].nonmiss, suidx)) {
          all_covar_ok = 0;
          break;
        }
      }
      if (!all_covar_ok) continue;
      for (uint32_t ci = 0; ci != covar_ct; ++ci) {
        const PhenoCol* cc = &covar_cols[ci];
        X[S_CAST(uint64_t, 1 + ci) * complete_ct + out_idx] =
            (cc->type_code == kPhenoDtypeQt)? cc->data.qt[suidx]
                                            : S_CAST(double, cc->data.cat[suidx]);
      }
      ++out_idx;
    }
  }
  return complete_ct;
}

PglErr BestPgsValidate(const uintptr_t* sample_include, uint32_t sample_ct,
                       const PhenoCol* pheno_cols, const char* pheno_names,
                       uintptr_t max_pheno_name_blen, uint32_t pheno_ct,
                       const PhenoCol* covar_cols, uint32_t covar_ct,
                       uint32_t fold_ct,
                       uint32_t* usable_pheno_idxs, uint32_t* usable_pheno_ct_ptr) {
  // Validate each phenotype individually; skip unusable ones with a warning.
  // Error only if no phenotype is usable.
  // usable_pheno_idxs: caller-allocated array of at least pheno_ct uint32_t.
  // *usable_pheno_ct_ptr: set to the count of entries written.

  if (pheno_ct == 0) {
    logerrputs("Error: --give-best-pgs requires a phenotype.\n");
    return kPglRetInconsistentInput;
  }

  // The shared guards read y, so the complete-case selection has to be
  // materialised.  One buffer, refilled per phenotype; X is not needed here.
  double* y_complete =
      S_CAST(double*, malloc(sample_ct? (sample_ct * sizeof(double))
                                      : sizeof(double)));
  if (unlikely(!y_complete)) {
    return kPglRetNomem;
  }

  uint32_t usable_ct = 0;
  for (uint32_t pi = 0; pi != pheno_ct; ++pi) {
    const PhenoCol* the_pheno = &pheno_cols[pi];
    const char* pheno_name = &pheno_names[pi * max_pheno_name_blen];

    // Categorical phenotypes have no fit to attempt.  This is the one guard the
    // --sscore route cannot share: there a type is inferred from numbers and
    // never comes out categorical.
    if (the_pheno->type_code == kPhenoDtypeCat) {
      logerrprintfww("Warning: --give-best-pgs: phenotype '%s' is categorical; skipping.\n",
                     pheno_name);
      continue;
    }

    PgsPhenoBlock pb;
    pb.y = y_complete;
    pb.n = PfileFillPhenoBlock(sample_include, sample_ct, the_pheno,
                               covar_cols, covar_ct, y_complete, nullptr,
                               nullptr);
    pb.covar_ct = covar_ct;
    pb.type = the_pheno->type_code;
    pb.name = pheno_name;
    pb.route = "--give-best-pgs";

    uint32_t usable;
    const PglErr reterr = ValidatePhenoBlock(&pb, fold_ct, &usable);
    if (unlikely(reterr)) {
      free(y_complete);
      return reterr;
    }
    if (usable) {
      usable_pheno_idxs[usable_ct++] = pi;
    }
  }
  free(y_complete);

  if (usable_ct == 0) {
    logerrputs("Error: --give-best-pgs: no usable phenotype found.\n");
    return kPglRetInconsistentInput;
  }
  *usable_pheno_ct_ptr = usable_ct;
  return kPglRetSuccess;
}

// ---------------------------------------------------------------------------
// --pgs-pheno-pairs
// ---------------------------------------------------------------------------
// A filter, not a selector.  Every phenotype in the --pheno file is loaded as
// usual; this restricts which (phenotype, PGS) combinations get fitted.
//
// That distinction is what keeps the feature to one implementation.  The two
// routes select phenotypes in different places -- --pfile through plink2's
// LoadPhenos far upstream, --sscore through its own header parse -- so a flag
// that drove phenotype loading would need a copy in each, which is how the
// --pfile and --sscore phenotype checks drifted apart before, until
// ValidatePhenoBlock gave them one set of guards.  Filtering at
// FitWritePhenoBlock,
// which both routes already share, needs one.

void PgsPhenoPairsPreinit(PgsPhenoPairs* pp) {
  pp->pheno_names = nullptr;
  pp->pgs_names = nullptr;
  pp->pair_ct = 0;
}

void PgsPhenoPairsFree(PgsPhenoPairs* pp) {
  if (pp->pheno_names) {
    for (uint32_t ii = 0; ii != pp->pair_ct; ++ii) {
      free(pp->pheno_names[ii]);
    }
    free(pp->pheno_names);
  }
  if (pp->pgs_names) {
    for (uint32_t ii = 0; ii != pp->pair_ct; ++ii) {
      free(pp->pgs_names[ii]);
    }
    free(pp->pgs_names);
  }
  PgsPhenoPairsPreinit(pp);
}

// Copy [start, end) into a fresh C string.  Returns nullptr on nomem.
static char* PairsDupToken(const char* start, const char* end) {
  const uintptr_t slen = S_CAST(uintptr_t, end - start);
  char* out = S_CAST(char*, malloc(slen + 1));
  if (!out) {
    return nullptr;
  }
  memcpy(out, start, slen);
  out[slen] = '\0';
  return out;
}

PglErr LoadPgsPhenoPairs(const char* fname, PgsPhenoPairs* pp) {
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  PgsPhenoPairsPreinit(pp);

  uint32_t cap = 16;
  pp->pheno_names = S_CAST(char**, malloc(cap * sizeof(char*)));
  pp->pgs_names = S_CAST(char**, malloc(cap * sizeof(char*)));
  if (unlikely((!pp->pheno_names) || (!pp->pgs_names))) {
    reterr = kPglRetNomem;
    goto LoadPgsPhenoPairs_ret_1;
  }

  reterr = SizeAndInitTextStream(fname, bigstack_left() / 8, 1, &txs);
  if (unlikely(reterr)) {
    logerrprintfww("Error: --pgs-pheno-pairs: failed to open %s .\n", fname);
    goto LoadPgsPhenoPairs_ret_1;
  }

  {
    uint32_t line_idx = 0;
    for (;;) {
      char* line = TextGet(&txs);
      if (!line) {
        break;
      }
      ++line_idx;
      // Skip leading whitespace, then blank lines and # comments.
      char* tok0 = line;
      while ((*tok0 == ' ') || (*tok0 == '\t')) {
        ++tok0;
      }
      if ((*tok0 == '\0') || (*tok0 == '\n') || (*tok0 == '\r') || (*tok0 == '#')) {
        continue;
      }
      char* tok0_end = tok0;
      while ((*tok0_end != ' ') && (*tok0_end != '\t') && (*tok0_end != '\n') &&
             (*tok0_end != '\r') && (*tok0_end != '\0')) {
        ++tok0_end;
      }
      char* tok1 = tok0_end;
      while ((*tok1 == ' ') || (*tok1 == '\t')) {
        ++tok1;
      }
      if ((*tok1 == '\0') || (*tok1 == '\n') || (*tok1 == '\r')) {
        logerrprintfww("Error: --pgs-pheno-pairs: line %u of %s has one column; a pair needs two (phenotype, then PGS).\n",
                       line_idx, fname);
        reterr = kPglRetMalformedInput;
        goto LoadPgsPhenoPairs_ret_1;
      }
      char* tok1_end = tok1;
      while ((*tok1_end != ' ') && (*tok1_end != '\t') && (*tok1_end != '\n') &&
             (*tok1_end != '\r') && (*tok1_end != '\0')) {
        ++tok1_end;
      }

      if (pp->pair_ct == cap) {
        cap *= 2;
        char** new_ph = S_CAST(char**, realloc(pp->pheno_names, cap * sizeof(char*)));
        if (new_ph) {
          pp->pheno_names = new_ph;
        }
        char** new_pg = S_CAST(char**, realloc(pp->pgs_names, cap * sizeof(char*)));
        if (new_pg) {
          pp->pgs_names = new_pg;
        }
        if (unlikely((!new_ph) || (!new_pg))) {
          reterr = kPglRetNomem;
          goto LoadPgsPhenoPairs_ret_1;
        }
      }
      // Both allocations happen before pair_ct is bumped, so a failure on the
      // second cannot leave a half-filled slot for PgsPhenoPairsFree to walk.
      char* ph = PairsDupToken(tok0, tok0_end);
      char* pg = PairsDupToken(tok1, tok1_end);
      if (unlikely((!ph) || (!pg))) {
        free(ph);
        free(pg);
        reterr = kPglRetNomem;
        goto LoadPgsPhenoPairs_ret_1;
      }
      pp->pheno_names[pp->pair_ct] = ph;
      pp->pgs_names[pp->pair_ct] = pg;
      ++pp->pair_ct;
    }

    if (TextStreamErrcode2(&txs, &reterr)) {
      // TextStream stops AT the first empty line and only afterwards reports
      // that there was content past it, so a blank line used to group the file
      // would otherwise truncate the pair list and fail with plink2's internal
      // wording.  Blank lines at the END are fine -- a file with a trailing
      // newline is the normal case -- and only interior ones are refused,
      // which is the rule every other input this tool reads already follows.
      if (!strcmp(TextStreamError(&txs), kShortErrInteriorEmptyLine)) {
        logerrprintfww("Error: --pgs-pheno-pairs: blank line inside %s. Blank lines are allowed only at the end of the file; comment a line out with '#' instead.\n",
                       fname);
        reterr = kPglRetMalformedInput;
      } else {
        TextStreamErrPrint(fname, &txs);
      }
      goto LoadPgsPhenoPairs_ret_1;
    }
    if (unlikely(!pp->pair_ct)) {
      logerrprintfww("Error: --pgs-pheno-pairs: %s contains no pairs.\n", fname);
      reterr = kPglRetMalformedInput;
      goto LoadPgsPhenoPairs_ret_1;
    }
  }

 LoadPgsPhenoPairs_ret_1:
  CleanupTextStream2(fname, &txs, &reterr);
  if (reterr) {
    PgsPhenoPairsFree(pp);
  }
  return reterr;
}

static uint32_t PairsNameInList(const char* name, const char* const* list,
                                uint32_t list_ct) {
  for (uint32_t ii = 0; ii != list_ct; ++ii) {
    if (!strcmp(name, list[ii])) {
      return 1;
    }
  }
  return 0;
}

PglErr PgsPairsValidate(const PgsPhenoPairs* pp,
                        const char* const* pheno_names, uint32_t pheno_ct,
                        const char* const* score_col_names, uint32_t score_col_ct) {
  // Fatal, not skipped.  A typo that silently dropped a pair would leave a
  // report that looks complete and is not -- and the whole point of the flag is
  // that the reader trusts the row list to be the one they asked for.
  for (uint32_t ii = 0; ii != pp->pair_ct; ++ii) {
    if (unlikely(!PairsNameInList(pp->pheno_names[ii], pheno_names, pheno_ct))) {
      logerrprintfww("Error: --pgs-pheno-pairs: '%s' is not a phenotype in this run.\n",
                     pp->pheno_names[ii]);
      return kPglRetInconsistentInput;
    }
    if (unlikely(!PairsNameInList(pp->pgs_names[ii], score_col_names, score_col_ct))) {
      logerrprintfww("Error: --pgs-pheno-pairs: '%s' is not a score column in this run.\n",
                     pp->pgs_names[ii]);
      return kPglRetInconsistentInput;
    }
  }
  return kPglRetSuccess;
}

uint32_t PgsPairsKeepIdxs(const PgsPhenoPairs* pp, const char* pheno_name,
                          const char* const* score_col_names,
                          uint32_t score_col_ct, uint32_t* keep_idxs) {
  // Walk the score columns in their own order, not the pairs file's, so the
  // rows a phenotype produces do not depend on the order pairs were listed in.
  // Duplicate pairs collapse for the same reason.
  uint32_t keep_ct = 0;
  for (uint32_t sc = 0; sc != score_col_ct; ++sc) {
    for (uint32_t ii = 0; ii != pp->pair_ct; ++ii) {
      if ((!strcmp(pheno_name, pp->pheno_names[ii])) &&
          (!strcmp(score_col_names[sc], pp->pgs_names[ii]))) {
        keep_idxs[keep_ct++] = sc;
        break;
      }
    }
  }
  return keep_ct;
}

// ---------------------------------------------------------------------------
// --give-corr-matrix
// ---------------------------------------------------------------------------
// Pearson correlation across the PGS columns, over every sample in the score
// matrix.  Deliberately independent of phenotypes and covariates: this is a
// property of the scores, so both routes call it with their own matrix and
// neither needs a phenotype to have been loaded.
//
// Consequence, intended and documented in --help: the correlations are NOT
// restricted to the complete-case set any particular R2 was fitted on, so with
// missing phenotype or covariate values they will not exactly reproduce a
// correlation computed from the rows behind a given .best_pgs line.
//
// A column is "invalid" if it has zero variance or holds a non-finite value.
// Its whole row and column are NA rather than 0: the correlation is undefined,
// and 0 would read as "uncorrelated", which is a different claim.
PglErr PgsCorrReport(const double* score_data, uint32_t model_ct,
                     uint32_t sample_ct, const char* const* model_names,
                     char* outname, char* outname_end) {
  PglErr reterr = kPglRetSuccess;
  FILE* outfile = nullptr;
  double* mean = nullptr;
  double* rss = nullptr;   // sum of squared deviations from the mean
  double* corr = nullptr;  // model_ct * model_ct, row-major

  if (unlikely(sample_ct < 2)) {
    logerrputs("Error: --give-corr-matrix: at least 2 samples are required; a correlation over\nfewer is undefined.\n");
    return kPglRetInconsistentInput;
  }

  {
    mean = S_CAST(double*, malloc(model_ct * sizeof(double)));
    rss = S_CAST(double*, malloc(model_ct * sizeof(double)));
    corr = S_CAST(double*, malloc(S_CAST(uint64_t, model_ct) * model_ct * sizeof(double)));
    if (unlikely((!mean) || (!rss) || (!corr))) {
      reterr = kPglRetNomem;
      goto PgsCorrReport_ret_1;
    }

    // Pass 1: per-column mean and sum of squared deviations.  rss is set to a
    // negative sentinel for an invalid column, which is the single condition
    // every later test reads -- so "non-finite" and "constant" cannot drift
    // into being handled differently.
    uint32_t invalid_ct = 0;
    for (uint32_t mi = 0; mi != model_ct; ++mi) {
      const double* col = &(score_data[S_CAST(uint64_t, mi) * sample_ct]);
      double sum = 0.0;
      uint32_t finite = 1;
      for (uint32_t si = 0; si != sample_ct; ++si) {
        const double val = col[si];
        if (!isfinite(val)) {
          finite = 0;
          break;
        }
        sum += val;
      }
      if (!finite) {
        mean[mi] = 0.0;
        rss[mi] = -1.0;
        ++invalid_ct;
        continue;
      }
      const double mu = sum / u31tod(sample_ct);
      double ss = 0.0;
      for (uint32_t si = 0; si != sample_ct; ++si) {
        const double dev = col[si] - mu;
        ss += dev * dev;
      }
      mean[mi] = mu;
      // A column that is constant to the last bit gives ss == 0 exactly.  Do
      // not widen this to a tolerance: a genuinely tiny-variance column still
      // has a defined correlation, and calling it NA would hide a real number.
      if (!(ss > 0.0)) {
        rss[mi] = -1.0;
        ++invalid_ct;
      } else {
        rss[mi] = ss;
      }
    }

    if (invalid_ct) {
      logerrprintf("Warning: --give-corr-matrix: %u score column%s constant or non-finite; reported as NA.\n",
                   invalid_ct, (invalid_ct == 1)? " is" : "s are");
    }

    // Pass 2: upper triangle, mirrored.  O(model_ct^2 * sample_ct), with no
    // copy of the score matrix -- the mean subtraction is recomputed rather
    // than materialised, which is one subtract against a second matrix in RAM.
    for (uint32_t ii = 0; ii != model_ct; ++ii) {
      double* corr_row_i = &(corr[S_CAST(uint64_t, ii) * model_ct]);
      if (rss[ii] < 0.0) {
        for (uint32_t jj = 0; jj != model_ct; ++jj) {
          corr_row_i[jj] = 0.0 / 0.0;
          corr[S_CAST(uint64_t, jj) * model_ct + ii] = 0.0 / 0.0;
        }
        continue;
      }
      // The pair loop below starts at ii + 1, so it never reaches the
      // diagonal: this assignment is the only thing that writes it, and 1 is
      // the value by definition rather than a rounding of one.  Running the
      // general formula on the diagonal instead was tried and also prints
      // exactly 1 on the suite's data, so this is a structural choice, not a
      // numerical one.
      corr_row_i[ii] = 1.0;
      const double* col_i = &(score_data[S_CAST(uint64_t, ii) * sample_ct]);
      const double mean_i = mean[ii];
      for (uint32_t jj = ii + 1; jj != model_ct; ++jj) {
        if (rss[jj] < 0.0) {
          continue;  // already filled with NaN by the invalid-column branch
        }
        const double* col_j = &(score_data[S_CAST(uint64_t, jj) * sample_ct]);
        const double mean_j = mean[jj];
        double cross = 0.0;
        for (uint32_t si = 0; si != sample_ct; ++si) {
          cross += (col_i[si] - mean_i) * (col_j[si] - mean_j);
        }
        double rr = cross / sqrt(rss[ii] * rss[jj]);
        // Rounding can push a perfectly collinear pair a few ulps past 1.
        if (rr > 1.0) {
          rr = 1.0;
        } else if (rr < -1.0) {
          rr = -1.0;
        }
        corr_row_i[jj] = rr;
        corr[S_CAST(uint64_t, jj) * model_ct + ii] = rr;
      }
    }

    // Write.  Full square matrix: the header row names the columns and every
    // data row starts with its own name, so a triangle never has to be
    // reconstructed and the reader never has to know which half it got.
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pgs_corr");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      reterr = kPglRetOpenFail;
      goto PgsCorrReport_ret_1;
    }
    fputs("#ID", outfile);
    for (uint32_t mi = 0; mi != model_ct; ++mi) {
      putc_unlocked('\t', outfile);
      fputs(model_names[mi], outfile);
    }
    putc_unlocked('\n', outfile);
    for (uint32_t ii = 0; ii != model_ct; ++ii) {
      fputs(model_names[ii], outfile);
      const double* corr_row_i = &(corr[S_CAST(uint64_t, ii) * model_ct]);
      for (uint32_t jj = 0; jj != model_ct; ++jj) {
        putc_unlocked('\t', outfile);
        const double rr = corr_row_i[jj];
        if (isnan(rr)) {
          fputs("NA", outfile);
        } else {
          // 12 significant digits, not the tool's usual dtoa_g 6.  Everything
          // in .best_pgs is an estimate carrying a standard error, so 6 digits
          // is already past the point of meaning; a correlation matrix is
          // routinely fed back into arithmetic -- pruning at a threshold,
          // decomposition -- where the extra digits are real.  It is also what
          // makes agreement with R's cor() checkable to 1e-10, which is the
          // only external evidence these numbers get.
          char buf[32];
          const int slen = snprintf(buf, sizeof(buf), "%.12g", rr);
          fwrite(buf, 1, S_CAST(uintptr_t, slen), outfile);
        }
      }
      putc_unlocked('\n', outfile);
    }
    if (unlikely(fclose_null(&outfile))) {
      reterr = kPglRetWriteFail;
      goto PgsCorrReport_ret_1;
    }
    *outname_end = '\0';
    logprintfww("--give-corr-matrix: %u x %u correlation matrix written to %s.pgs_corr .\n",
                model_ct, model_ct, outname);
  }

 PgsCorrReport_ret_1:
  if (outfile) {
    fclose(outfile);
  }
  free(corr);
  free(rss);
  free(mean);
  return reterr;
}

// Length of `nm` with plink2's _AVG / _SUM score suffix dropped, or slen when
// there is none to drop.  Case-sensitive, matching what plink2 writes: a
// lowercase `_avg` is left alone rather than half-recognised.
static uint32_t PgsIdSlen(const char* nm, uint32_t slen) {
  if (slen > 4) {
    const char* suf = &(nm[slen - 4]);
    if (memequal(suf, "_AVG", 4) || memequal(suf, "_SUM", 4)) {
      return slen - 4;
    }
  }
  return slen;
}

// Shorten each score column name in place to the PGS ID the reports print.
//
// Dropping the _AVG / _SUM suffix keeps IDs readable while a .sscore names one
// column per model.  cols=+scoresums names TWO -- "<model>_AVG" and
// "<model>_SUM" -- and dropping the suffix from both landed them on the same
// ID: .pgs_corr came out headed `#ID W1 W2 W1 W2` with those same names down
// its rows, which no reader indexing by name can tell apart (R's
// read.table(row.names=1) refuses it outright), and .best_pgs carried two rows
// per (phenotype, model) that looked identical and held different numbers.
//
// So a column keeps its full name when the shortened form would equal any
// other column's shortened form, or any other column's full name.  The second
// clause is not implied by the first: with W1_AVG, W1_SUM and W1_AVG_AVG
// selected, the first two keep their suffixes under the first clause, and
// W1_AVG_AVG would otherwise shorten straight onto the W1_AVG that just kept
// its own.  Together they leave every ID distinct unless the .sscore header
// repeats a column name outright, which no naming rule can repair.
//
// Both input routes call this.  They used to carry separate copies of the
// strip, and both had the collapse.
//
// O(n^2) in the score column count, each compare gated on a length test.  n is
// the number of models asked for; one regression per model costs far more.
//
// Truncation is in place -- every caller owns a malloc'd copy of the full
// column name, and the ID is always a prefix of it.  Every column is decided
// against the full set of unshortened names before any of them is cut, so the
// answer does not depend on the order they are visited in.
static PglErr PgsIdsFromScoreColNames(char** names, uint32_t col_ct) {
  // full_slens[] then the chosen ID length per column, in one allocation.
  uint32_t* slens = S_CAST(uint32_t*, malloc(S_CAST(uint64_t, col_ct) * 2 * sizeof(uint32_t)));
  if (!slens) {
    return kPglRetNomem;
  }
  uint32_t* full_slens = slens;
  uint32_t* id_slens = &(slens[col_ct]);
  for (uint32_t i = 0; i != col_ct; ++i) {
    full_slens[i] = strlen(names[i]);
  }

  for (uint32_t i = 0; i != col_ct; ++i) {
    const uint32_t id_slen = PgsIdSlen(names[i], full_slens[i]);
    id_slens[i] = full_slens[i];
    if (id_slen == full_slens[i]) {
      continue;  // no suffix to drop
    }
    uint32_t collides = 0;
    for (uint32_t j = 0; j != col_ct; ++j) {
      if (j == i) {
        continue;
      }
      // Against j's full name and against the ID j would take on its own.
      // Whichever of the two j ends up with, it is one of these, so j's own
      // decision cannot turn a miss here into a clash later.  Both branches
      // imply full_slens[j] >= id_slen, so the compare stays in bounds.
      if (((full_slens[j] == id_slen) ||
           (PgsIdSlen(names[j], full_slens[j]) == id_slen)) &&
          memequal(names[j], names[i], id_slen)) {
        collides = 1;
        break;
      }
    }
    if (!collides) {
      id_slens[i] = id_slen;
    }
  }

  for (uint32_t i = 0; i != col_ct; ++i) {
    names[i][id_slens[i]] = '\0';
  }
  free(slens);
  return kPglRetSuccess;
}

PglErr BestPgsReport(const char* sscore_fname,
                     const uintptr_t* sample_include, uint32_t sample_ct,
                     const PhenoCol* pheno_cols, const char* pheno_names,
                     uintptr_t max_pheno_name_blen, uint32_t pheno_ct,
                     const PhenoCol* covar_cols, const char* covar_names,
                     uintptr_t max_covar_name_blen, uint32_t covar_ct,
                     double prevalence, uint32_t fold_ct, uint64_t fold_seed,
                     uint32_t give_best_pgs, uint32_t corr_matrix,
                     const char* pairs_fname,
                     char* outname, char* outname_end) {
  PglErr reterr = kPglRetSuccess;
  PgsPhenoPairs pairs;
  PgsPhenoPairsPreinit(&pairs);
  uint32_t* pair_keep_idxs = nullptr;
  const char** pheno_names_arr = nullptr;
  TextStream txs;
  PreinitTextStream(&txs);
  FILE* outfile = nullptr;

  // --- allocations (plain malloc; freed in ret_1) ---
  // score column names (array of char*); populated during header parse
  uint32_t score_col_ct = 0;
  char**   score_col_names = nullptr;  // malloc'd array of malloc'd strings
  // flat arrays over sample_ct rows, score_col_ct columns (column-major)
  double*  score_data = nullptr;
  double*  y_complete = nullptr;
  double*  X_complete = nullptr;
  // .sscore row index of each complete case, per phenotype; how the score
  // matrix gets selected down to the sample set the phenotype was validated on.
  uint32_t* keep_ri = nullptr;
  // All rows across all phenotypes (pheno_usable_ct * score_col_ct entries).
  BestPgsRow* all_rows = nullptr;
  // Usable phenotype index list from BestPgsValidate.
  uint32_t*   usable_pheno_idxs = nullptr;

  {
    // -----------------------------------------------------------------------
    // 1. Validate phenotype/covariates (score-independent checks).
    //    Skipped entirely under --give-corr-matrix alone: a correlation matrix
    //    is a property of the scores, so demanding a usable phenotype for one
    //    would refuse a run that has everything it needs.
    // -----------------------------------------------------------------------
    uint32_t usable_pheno_ct = 0;
    if (give_best_pgs) {
      usable_pheno_idxs = S_CAST(uint32_t*, malloc(pheno_ct * sizeof(uint32_t)));
      if (!usable_pheno_idxs) {
        reterr = kPglRetNomem;
        goto BestPgsReport_ret_1;
      }
      reterr = BestPgsValidate(sample_include, sample_ct, pheno_cols, pheno_names,
                               max_pheno_name_blen, pheno_ct, covar_cols, covar_ct,
                               fold_ct, usable_pheno_idxs, &usable_pheno_ct);
      if (unlikely(reterr)) {
        goto BestPgsReport_ret_1;
      }
    }

    // Warn if --prevalence was given but ALL usable phenotypes are quantitative.
    // (Per-phenotype prevalence check is done in the loop below.)
    // The flag-parsing check in adjutant.cc ensures prevalence is only set when
    // --give-best-pgs is also set, so we only need the type-dependent check here.

    // -----------------------------------------------------------------------
    // 2. Parse the .sscore header: discover score columns
    // -----------------------------------------------------------------------
    reterr = SizeAndInitTextStream(sscore_fname, bigstack_left() / 8, 1, &txs);
    if (unlikely(reterr)) {
      logerrprintfww("Error: Failed to open %s .\n", sscore_fname);
      goto BestPgsReport_ret_1;
    }

    char* header_line = TextGet(&txs);
    if (!header_line) {
      if (!TextStreamErrcode2(&txs, &reterr)) {
        logerrputs("Error: --give-best-pgs: .sscore file is empty.\n");
        reterr = kPglRetMalformedInput;
      } else {
        TextStreamErrPrint(sscore_fname, &txs);
      }
      goto BestPgsReport_ret_1;
    }

    // Strip the leading '#' if present.
    char* col_iter = header_line;
    if (*col_iter == '#') {
      ++col_iter;
    }

    // Walk through all tab-delimited column names.
    // Score columns end in _AVG or _SUM but not NAMED_ALLELE_DOSAGE_SUM.
    // Allocate a generous initial array; we'll realloc as needed.
    uint32_t score_col_cap = 32;
    score_col_names = S_CAST(char**, malloc(score_col_cap * sizeof(char*)));
    if (!score_col_names) {
      reterr = kPglRetNomem;
      goto BestPgsReport_ret_1;
    }

    while (1) {
      // Find end of this column name token (tab or EOL).
      char* col_end = col_iter;
      while (*col_end != '\t' && *col_end != '\n' && *col_end != '\r' && *col_end != '\0') {
        ++col_end;
      }
      const uint32_t col_slen = S_CAST(uint32_t, col_end - col_iter);
      const int at_end = (*col_end == '\0' || *col_end == '\n' || *col_end == '\r');

      // Check for _AVG or _SUM suffix, excluding NAMED_ALLELE_DOSAGE_SUM.
      const char kDosageSumName[] = "NAMED_ALLELE_DOSAGE_SUM";
      const uint32_t is_dosage_sum = (col_slen == sizeof(kDosageSumName) - 1) &&
                                      memequal(col_iter, kDosageSumName, col_slen);
      if (!is_dosage_sum && col_slen > 4) {
        const char* suffix_ptr = col_end - 4;
        const int has_avg = (col_slen >= 4) && memequal(suffix_ptr, "_AVG", 4);
        const int has_sum = (col_slen >= 4) && memequal(suffix_ptr, "_SUM", 4);
        if (has_avg || has_sum) {
          // Store the full column name.  PgsIdsFromScoreColNames shortens it
          // to the PGS ID below, once the whole set is known -- the suffix can
          // only be dropped safely with the other columns' names in hand.
          char* nm = S_CAST(char*, malloc(col_slen + 1));
          if (!nm) {
            reterr = kPglRetNomem;
            goto BestPgsReport_ret_1;
          }
          memcpy(nm, col_iter, col_slen);
          nm[col_slen] = '\0';

          if (score_col_ct == score_col_cap) {
            score_col_cap *= 2;
            char** tmp = S_CAST(char**, realloc(score_col_names, score_col_cap * sizeof(char*)));
            if (!tmp) {
              free(nm);
              reterr = kPglRetNomem;
              goto BestPgsReport_ret_1;
            }
            score_col_names = tmp;
          }
          score_col_names[score_col_ct++] = nm;
        }
      }

      if (at_end) {
        break;
      }
      col_iter = col_end + 1;  // advance past the tab
    }

    if (score_col_ct == 0) {
      logerrputs("Error: --give-best-pgs: no score columns (_AVG or _SUM) found in .sscore header.\n");
      reterr = kPglRetMalformedInput;
      goto BestPgsReport_ret_1;
    }

    reterr = PgsIdsFromScoreColNames(score_col_names, score_col_ct);
    if (unlikely(reterr)) {
      goto BestPgsReport_ret_1;
    }

    // -----------------------------------------------------------------------
    // 3. Read data rows into score_data[score_col_ct * sample_ct] (col-major)
    // -----------------------------------------------------------------------
    // score_data: column-major [score_col_ct x sample_ct]
    score_data = S_CAST(double*, malloc(S_CAST(uint64_t, score_col_ct) * sample_ct * sizeof(double)));
    if (!score_data) {
      reterr = kPglRetNomem;
      goto BestPgsReport_ret_1;
    }

    // Walk .sscore rows.
    // Each data row: <IID> [pheno-cols] [ALLELE_CT] [DOSAGE_SUM] <score cols>
    // We don't care about the IID or non-score columns — we just need the
    // score columns at the end.  The count of non-score prefix columns is
    // (total_header_cols - score_col_ct).  We walk the header again to count.
    // Simpler: just count non-score prefix columns from the header we already
    // parsed by counting the total columns and subtracting score_col_ct.
    // But we've already consumed the header line.  Re-derive from score_col_ct:
    // We don't know the total columns without re-scanning.  Instead, we'll
    // skip past (total - score_col_ct) tokens at the start of each data row.
    // We need to count total columns in the header.
    // Simpler approach: count how many non-score-suffix columns there are by
    // difference from total.  But we only recorded score columns.
    // Best approach: just count the number of prefix columns to skip.
    // We'll re-derive it: we know the header has 1 IID + possibly FID/SID +
    // pheno cols + ALLELE_CT + DOSAGE_SUM + score_col_ct columns.
    // But rather than hardcoding, let's count the header tokens directly.
    // We already walked the header but didn't count total tokens.  We need to
    // do it differently.  Let's re-read the header line.
    // Actually, we can't rewind the TextStream.  We need to count during the
    // header scan.

    // Restart with a new text stream: close this one, reopen.
    CleanupTextStream2(sscore_fname, &txs, &reterr);
    PreinitTextStream(&txs);
    reterr = SizeAndInitTextStream(sscore_fname, bigstack_left() / 8, 1, &txs);
    if (unlikely(reterr)) {
      logerrprintfww("Error: Failed to reopen %s .\n", sscore_fname);
      goto BestPgsReport_ret_1;
    }

    header_line = TextGet(&txs);
    if (!header_line) {
      // can't happen: we read it earlier
      reterr = kPglRetMalformedInput;
      goto BestPgsReport_ret_1;
    }

    // Count total columns in the header AND the number of prefix (non-score)
    // columns to skip per data row.
    uint32_t total_header_cols = 0;
    uint32_t score_col_seen = 0;
    {
      char* ci = header_line;
      if (*ci == '#') ++ci;
      uint32_t col_i = 0;
      while (1) {
        char* ce = ci;
        while (*ce != '\t' && *ce != '\n' && *ce != '\r' && *ce != '\0') ++ce;
        const uint32_t slen = S_CAST(uint32_t, ce - ci);
        const int at_end2 = (*ce == '\0' || *ce == '\n' || *ce == '\r');

        const char kDosageSumName2[] = "NAMED_ALLELE_DOSAGE_SUM";
        const uint32_t is_ds2 = (slen == sizeof(kDosageSumName2) - 1) &&
                                  memequal(ci, kDosageSumName2, slen);
        if (!is_ds2 && slen > 4) {
          const char* sp = ce - 4;
          if (memequal(sp, "_AVG", 4) || memequal(sp, "_SUM", 4)) {
            ++score_col_seen;
          }
        }
        ++col_i;
        if (at_end2) { break; }
        ci = ce + 1;
      }
      total_header_cols = col_i;
    }

    const uint32_t prefix_col_ct = total_header_cols - score_col_ct;
    // sanity: score_col_seen should equal score_col_ct
    if (score_col_seen != score_col_ct || prefix_col_ct > total_header_cols) {
      logerrputs("Error: --give-best-pgs: internal header-parse inconsistency.\n");
      reterr = kPglRetMalformedInput;
      goto BestPgsReport_ret_1;
    }

    // Now read data rows.
    uintptr_t row_idx = 0;
    uintptr_t sample_uidx_base = 0;
    uintptr_t sample_include_bits = sample_include[0];
    for (; ; ++row_idx) {
      char* line = TextGet(&txs);
      if (!line) {
        break;
      }
      if (row_idx >= sample_ct) {
        logerrprintfww("Error: --give-best-pgs: .sscore file has more rows than expected (%u).\n", sample_ct);
        reterr = kPglRetMalformedInput;
        goto BestPgsReport_ret_1;
      }

      // Advance to the first score column (skip prefix_col_ct tokens).
      const char* tok = S_CAST(const char*, line);
      for (uint32_t skip = 0; skip != prefix_col_ct; ++skip) {
        tok = NextToken(tok);
        if (unlikely(!tok)) {
          logerrprintfww("Error: --give-best-pgs: too few columns on row %" PRIuPTR " of .sscore.\n", row_idx + 1);
          reterr = kPglRetMalformedInput;
          goto BestPgsReport_ret_1;
        }
      }

      // Parse score_col_ct doubles.
      for (uint32_t sc = 0; sc != score_col_ct; ++sc) {
        double val;
        const char* tok_end = ScanadvDouble(tok, &val);
        if (unlikely(!tok_end)) {
          logerrprintfww("Error: --give-best-pgs: invalid score value on row %" PRIuPTR ", score column %u.\n", row_idx + 1, sc);
          reterr = kPglRetMalformedInput;
          goto BestPgsReport_ret_1;
        }
        // column-major: score_data[sc * sample_ct + row_idx]
        score_data[S_CAST(uint64_t, sc) * sample_ct + row_idx] = val;
        tok = NextToken(tok_end);
      }

      // Advance the sample_uidx iterator (we don't use sample_uidx here
      // since phenotype values are looked up later during complete-case
      // selection, directly from the_pheno->data.qt[suidx]).
      BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
    }

    // Check end of file.
    if (TextStreamErrcode2(&txs, &reterr)) {
      TextStreamErrPrint(sscore_fname, &txs);
      goto BestPgsReport_ret_1;
    }

    if (row_idx != sample_ct) {
      logerrprintfww("Error: %s: .sscore file has %u rows but %u were expected.\n",
                     give_best_pgs? "--give-best-pgs" : "--give-corr-matrix",
                     S_CAST(uint32_t, row_idx), sample_ct);
      reterr = kPglRetMalformedInput;
      goto BestPgsReport_ret_1;
    }

    // -----------------------------------------------------------------------
    // 3b. Correlation matrix, if asked for.  Before the fitting loop and over
    //     the full sample_ct, not a per-phenotype complete-case subset.
    // -----------------------------------------------------------------------
    if (corr_matrix) {
      reterr = PgsCorrReport(score_data, score_col_ct, sample_ct,
                             S_CAST(const char* const*, score_col_names),
                             outname, outname_end);
      if (unlikely(reterr)) {
        goto BestPgsReport_ret_1;
      }
    }
    if (!give_best_pgs) {
      goto BestPgsReport_ret_1;  // reterr is kPglRetSuccess
    }

    // -----------------------------------------------------------------------
    // 3c. --pgs-pheno-pairs.  Both name sets are known here and nowhere
    //     earlier: the score columns come from the .sscore this function just
    //     read.  One validation site rather than two, at the cost of a typo
    //     surfacing after the scoring pass rather than before it.
    // -----------------------------------------------------------------------
    if (pairs_fname) {
      reterr = LoadPgsPhenoPairs(pairs_fname, &pairs);
      if (unlikely(reterr)) {
        goto BestPgsReport_ret_1;
      }
      // Validate against every LOADED phenotype, not just the usable ones: a
      // phenotype that exists but was skipped is not a typo, and reporting it
      // as one would be wrong.  It simply contributes no rows.
      pheno_names_arr = S_CAST(const char**, malloc(pheno_ct? pheno_ct * sizeof(char*) : sizeof(char*)));
      pair_keep_idxs = S_CAST(uint32_t*, malloc(score_col_ct * sizeof(uint32_t)));
      if (unlikely((!pheno_names_arr) || (!pair_keep_idxs))) {
        reterr = kPglRetNomem;
        goto BestPgsReport_ret_1;
      }
      for (uint32_t pi = 0; pi != pheno_ct; ++pi) {
        pheno_names_arr[pi] = &pheno_names[pi * max_pheno_name_blen];
      }
      reterr = PgsPairsValidate(&pairs, pheno_names_arr, pheno_ct,
                                S_CAST(const char* const*, score_col_names),
                                score_col_ct);
      if (unlikely(reterr)) {
        goto BestPgsReport_ret_1;
      }
    }

    // -----------------------------------------------------------------------
    // 4. Allocate row storage for ALL phenotypes x ALL score columns.
    //    Rows are grouped by phenotype (selection order) and sorted within
    //    each group; the loop below fills them in that order.
    // -----------------------------------------------------------------------
    all_rows = S_CAST(BestPgsRow*, malloc(S_CAST(uint64_t, usable_pheno_ct) * score_col_ct * sizeof(BestPgsRow)));
    if (!all_rows) {
      reterr = kPglRetNomem;
      goto BestPgsReport_ret_1;
    }

    // Largest possible complete-case set is sample_ct.  Allocate for that.
    y_complete = S_CAST(double*, malloc(sample_ct * sizeof(double)));
    // X_complete: column-major [sample_ct x (1 + covar_ct)] — worst case
    X_complete = S_CAST(double*, malloc(S_CAST(uint64_t, sample_ct) * (1 + covar_ct) * sizeof(double)));
    keep_ri = S_CAST(uint32_t*, malloc(sample_ct * sizeof(uint32_t)));
    if (!y_complete || !X_complete || !keep_ri) {
      reterr = kPglRetNomem;
      goto BestPgsReport_ret_1;
    }

    // -----------------------------------------------------------------------
    // 5. Loop over usable phenotypes: per-phenotype complete-case selection
    //    and call the shared fit+write helper.
    // -----------------------------------------------------------------------

    // Build a flat array of covar name pointers (for FitWritePhenoBlock).
    // These borrow from the stride-allocated covar_names buffer.
    const char** covar_names_arr = nullptr;
    if (covar_ct) {
      covar_names_arr = S_CAST(const char**, malloc(covar_ct * sizeof(char*)));
      if (!covar_names_arr) {
        reterr = kPglRetNomem;
        goto BestPgsReport_ret_1;
      }
      for (uint32_t ci = 0; ci != covar_ct; ++ci) {
        covar_names_arr[ci] = &covar_names[S_CAST(uint64_t, ci) * max_covar_name_blen];
      }
    }

    // score_block: complete-case-selected scores for one phenotype.
    // [score_col_ct x complete_ct_max] column-major, reused per phenotype.
    double* score_block = S_CAST(double*, malloc(
        S_CAST(uint64_t, score_col_ct) * sample_ct * sizeof(double)));
    if (!score_block) {
      free(covar_names_arr);
      reterr = kPglRetNomem;
      goto BestPgsReport_ret_1;
    }

    for (uint32_t pi_pos = 0; pi_pos != usable_pheno_ct; ++pi_pos) {
      const uint32_t pi = usable_pheno_idxs[pi_pos];
      const PhenoCol* the_pheno = &pheno_cols[pi];
      const char*     pheno_name = &pheno_names[pi * max_pheno_name_blen];
      const PhenoDtype type_code = the_pheno->type_code;

      // --- 5a. Complete-case selection: y, X, and the rows they came from. ---
      // The same producer BestPgsValidate ran its guards on, so the sample set
      // fitted here is the one that was validated, by construction.
      const uint32_t complete_ct =
          PfileFillPhenoBlock(sample_include, sample_ct, the_pheno,
                              covar_cols, covar_ct,
                              y_complete, X_complete, keep_ri);

      // --- 5b. Select the score matrix down to those same rows. ---
      for (uint32_t sc = 0; sc != score_col_ct; ++sc) {
        for (uint32_t ii = 0; ii != complete_ct; ++ii) {
          score_block[S_CAST(uint64_t, sc) * complete_ct + ii] =
              score_data[S_CAST(uint64_t, sc) * sample_ct + keep_ri[ii]];
        }
      }

      // --- 5c-5f. Constant-covariate drop, fit, sort, write. ---
      uint32_t keep_ct = 0;
      if (pairs_fname) {
        keep_ct = PgsPairsKeepIdxs(&pairs, pheno_name,
                                   S_CAST(const char* const*, score_col_names),
                                   score_col_ct, pair_keep_idxs);
      }
      BestPgsRow* pheno_rows = &all_rows[S_CAST(uint64_t, pi_pos) * score_col_ct];
      reterr = FitWritePhenoBlock(complete_ct, covar_ct, y_complete, X_complete,
                                  score_block,
                                  score_col_names, score_col_ct,
                                  covar_names_arr, pheno_name, type_code,
                                  prevalence, fold_ct, fold_seed,
                                  pair_keep_idxs, keep_ct, &outfile,
                                  outname, outname_end, pheno_rows);
      if (unlikely(reterr)) {
        free(score_block);
        free(covar_names_arr);
        goto BestPgsReport_ret_1;
      }
    }  // end phenotype loop

    free(score_block);
    free(covar_names_arr);

    // Under a filter every named phenotype can turn out to be one that was
    // skipped, leaving nothing written at all.  A header-only .best_pgs would
    // read as a complete report of nothing; say so instead.
    if (unlikely(!outfile)) {
      logerrputs("Error: --pgs-pheno-pairs: no pair could be evaluated; every phenotype named in\nthe pairs file was skipped.\n");
      reterr = kPglRetInconsistentInput;
      goto BestPgsReport_ret_1;
    }
    if (unlikely(fclose_null(&outfile))) {
      reterr = kPglRetWriteFail;
      goto BestPgsReport_ret_1;
    }

    // Null-terminate outname at the extension start for the log message.
    *outname_end = '\0';
    logprintfww("--give-best-pgs: Results written to %s.best_pgs .\n", outname);
  }

 BestPgsReport_ret_1:
  PgsPhenoPairsFree(&pairs);
  free(pair_keep_idxs);
  free(pheno_names_arr);
  if (outfile) fclose(outfile);
  CleanupTextStream2(sscore_fname, &txs, &reterr);
  free(all_rows);
  free(keep_ri);
  free(X_complete);
  free(y_complete);
  free(score_data);
  free(usable_pheno_idxs);
  if (score_col_names) {
    for (uint32_t sc = 0; sc != score_col_ct; ++sc) {
      free(score_col_names[sc]);
    }
    free(score_col_names);
  }
  return reterr;
}

// The complete-case sample set for one --sscore phenotype: the joined rows
// whose phenotype value and every covariate value are present.  Writes those
// rows' joined-array indices to keep_idx (capacity joined_ct) and returns how
// many there were.
//
// This is the --sscore half of the pair PfileFillPhenoBlock completes: one
// producer per input route, both handing the same PgsPhenoBlock to the same
// ValidatePhenoBlock.  This route used to have no equivalent and fitted on
// every joined row, which silently turned a missing case/control value into a
// control -- N and N_CASE then described a sample set that did not exist, and
// with --folds both the pre-flight k check and the stratified fold assignment
// ran on those fabricated labels.
static uint32_t SscoreCompleteCases(const double* pheno_data,
                                    uint32_t pheno_sample_ct,
                                    const uint32_t* join_pheno_idx, uint32_t ph,
                                    const double* covar_data,
                                    uint32_t covar_sample_ct,
                                    const uint32_t* join_covar_idx,
                                    uint32_t covar_col_ct, uint32_t joined_ct,
                                    uint32_t* keep_idx) {
  uint32_t keep_ct = 0;
  for (uint32_t si = 0; si != joined_ct; ++si) {
    if (isnan(pheno_data[S_CAST(uint64_t, ph) * pheno_sample_ct + join_pheno_idx[si]])) {
      continue;
    }
    uint32_t covars_ok = 1;
    for (uint32_t cv = 0; cv != covar_col_ct; ++cv) {
      if (isnan(covar_data[S_CAST(uint64_t, cv) * covar_sample_ct + join_covar_idx[si]])) {
        covars_ok = 0;
        break;
      }
    }
    if (covars_ok) {
      keep_idx[keep_ct++] = si;
    }
  }
  return keep_ct;
}

// -------------------------------------------------------------------
// How many columns precede the IID in a PLINK phenotype/covariate
// header?  PLINK writes both of these and both must be accepted:
//
//   #IID  <data...>        -> 0
//   #FID  IID  <data...>   -> 1
//
// Reading the second form as the first is not a cosmetic parse error.
// It makes the IID column itself a data column; every IID is a
// non-numeric string, so it scans as missing, every sample then has a
// "missing covariate", and the run dies reporting "no usable phenotype
// found" -- a message pointing nowhere near the header.  It also joins
// on FID, which matches nothing the moment FID differs from IID.
//
// Both readers call this, so the two cannot drift apart -- which is how
// this class of bug reached users before.
// -------------------------------------------------------------------
static uint32_t SscoreIdSkip(const char* hdr) {
  const char* t0e = hdr;
  while (*t0e != '\t' && *t0e != '\n' && *t0e != '\r' && *t0e != '\0') ++t0e;
  if ((t0e - hdr) != 3 || memcmp(hdr, "FID", 3) != 0) return 0;
  if (*t0e != '\t') return 0;
  const char* t1 = t0e + 1;
  const char* t1e = t1;
  while (*t1e != '\t' && *t1e != '\n' && *t1e != '\r' && *t1e != '\0') ++t1e;
  // A leading FID is an ID column only when IID follows it.  Otherwise it
  // is a data column that happens to be called FID, and skipping it would
  // silently discard the user's first phenotype.
  if ((t1e - t1) != 3 || memcmp(t1, "IID", 3) != 0) return 0;
  return 1;
}

// Column names that are never a score, whatever their position.  Two
// writers are in play -- this one emits
//   #IID SuperPop Population ALLELE_CT NAMED_ALLELE_DOSAGE_SUM W1_AVG ...
// and plink2 emits
//   #IID NMISS_ALLELE_CT DENOM NAMED_ALLELE_DOSAGE_SUM SCORE1_AVG ...
static const char* const kNonScoreCols[] = {
  "iid", "fid", "sid", "nmiss_allele_ct", "allele_ct", "denom",
  "named_allele_dosage_sum", "superpop", "population", nullptr
};

// Case-insensitive membership of kNonScoreCols.
static uint32_t SscoreIsNonScoreName(const char* nm, uint32_t slen) {
  for (uint32_t k = 0; kNonScoreCols[k]; ++k) {
    const char* ns = kNonScoreCols[k];
    if (strlen(ns) != slen) continue;
    uint32_t match = 1;
    for (uint32_t cc = 0; cc != slen; ++cc) {
      // Lowercase LETTERS only.  The obvious `c | 0x20` also maps '_'
      // (0x5F) to 0x7F, so no name here containing an underscore -- which
      // is most of them -- could ever match.  That stayed harmless while a
      // positional rule decided the outcome first, and became visible the
      // moment NAMED_ALLELE_DOSAGE_SUM had to be excluded on its name
      // rather than on its position.
      char c = nm[cc];
      if (c >= 'A' && c <= 'Z') c = S_CAST(char, c | 0x20);
      if (c != ns[cc]) { match = 0; break; }
    }
    if (match) return 1;
  }
  return 0;
}

// Does this column name carry plink2's score-column suffix?  The scoreavgs
// and scoresums column sets write "<name>_AVG" and "<name>_SUM", so the
// suffix names the score columns exactly wherever they sit in the header --
// which a positional rule cannot do once 'phenos' puts arbitrarily named
// phenotype columns in the middle of the fixed block.
//
// Callers must reject kNonScoreCols first: NAMED_ALLELE_DOSAGE_SUM ends in
// _SUM and is bookkeeping.
static uint32_t SscoreHasScoreSuffix(const char* nm, uint32_t slen) {
  if (slen < 5) return 0;
  const char* e = &(nm[slen - 4]);
  if (e[0] != '_') return 0;
  const char c1 = e[1] | 0x20, c2 = e[2] | 0x20, c3 = e[3] | 0x20;
  return ((c1 == 'a' && c2 == 'v' && c3 == 'g') ||
          (c1 == 's' && c2 == 'u' && c3 == 'm'));
}

// Advance `ci` past the FID column when there is one, so all the header
// arithmetic downstream counts from IID and needs no further adjustment.
// The data ROWS still contain the FID, so row-level token indices add
// id_skip separately.
static void SscoreSkipFidCol(char** ci_ptr, uint32_t id_skip) {
  if (!id_skip) return;
  char* ci = *ci_ptr;
  while (*ci != '\t' && *ci != '\n' && *ci != '\r' && *ci != '\0') ++ci;
  if (*ci == '\t') ++ci;
  *ci_ptr = ci;
}

// ---------------------------------------------------------------------------
// BestPgsFromSscore: --sscore standalone input mode.
//
// Reads a user-supplied .sscore plus phenotype and covariate files, joins
// them on IID, and feeds the plain arrays to FitWritePhenoBlock.  No .pgen,
// .psam or bigstack involved -- this path exists so a pre-computed score
// matrix can be evaluated without re-scoring.
//
// score_col_range: from --score-col-nums (may be an empty/zero RangeList to
// signal "use auto-detection").  Column numbers are 1-based and literal.
// When not supplied, every column after the leading fixed block is treated as
// a score column (heuristic -- see comment below).
//
// pheno_range_list: from --pheno-name or --pheno-col-nums (may have name_ct=0
// meaning "use all columns").
// ---------------------------------------------------------------------------
PglErr BestPgsFromSscore(const char* sscore_fname,
                         const RangeList* score_col_range,
                         const char* pheno_fname,
                         const RangeList* pheno_range_list,
                         const char* covar_fname,
                         double prevalence, uint32_t fold_ct, uint64_t fold_seed,
                         uint32_t corr_matrix, const char* pairs_fname,
                         char* outname, char* outname_end) {
  PglErr reterr = kPglRetSuccess;
  PgsPhenoPairs pairs;
  PgsPhenoPairsPreinit(&pairs);
  uint32_t* pair_keep_idxs = nullptr;
  TextStream txs;
  PreinitTextStream(&txs);
  FILE* outfile = nullptr;

  // All heap allocations; freed in ret_1.
  char**   score_col_names = nullptr;
  uint32_t score_col_ct = 0;
  // iid_buf: flat IID list from the .sscore (one C string per row).
  char**   sscore_iids = nullptr;
  uint32_t sscore_sample_ct = 0;
  // 1 when the .sscore leads with FID before IID, so the IID is token 1 rather
  // than token 0.  --bfile output always does: a .fam carries FID, so
  // FidColIsRequired is true under 'maybefid' and ScoreReport writes the
  // column.  Without this, adjutant could not read back its own --bfile
  // .sscore -- every row's ID came out as the FID, and 300 samples sharing
  // FID '0' were rejected as "duplicate IID '0'".
  uint32_t sscore_id_skip = 0;
  double*  score_data = nullptr;           // [score_col_ct x sscore_sample_ct] col-major
  double*  score_data_rowmajor = nullptr;  // [sscore_cap x score_col_ct] row-major (temp)

  char**   pheno_iids = nullptr;
  uint32_t pheno_sample_ct = 0;
  char**   pheno_col_names = nullptr;
  uint32_t pheno_col_ct = 0;
  double*  pheno_data = nullptr;   // [pheno_col_ct x pheno_sample_ct] col-major
  // PhenoDtype per selected phenotype (kPhenoDtypeQt or kPhenoDtypeCc).
  PhenoDtype* pheno_types = nullptr;

  char**   covar_iids = nullptr;
  uint32_t covar_sample_ct = 0;
  char**   covar_col_names = nullptr;
  uint32_t covar_col_ct = 0;
  double*  covar_data = nullptr;   // [covar_col_ct x covar_sample_ct] col-major

  // After the inner join: arrays over the intersection.
  uint32_t joined_ct = 0;
  uint32_t* keep_idx = nullptr;    // joined-array indices of one phenotype's
                                   // complete cases; refilled per phenotype
  double*  y_complete = nullptr;
  double*  X_complete = nullptr;
  double*  score_block = nullptr;
  BestPgsRow* all_rows = nullptr;

  // Helper: read a plain tab-delimited file (IID in col 0) into parallel
  // arrays.  Returns kPglRetMalformedInput on bad input.

  // -------------------------------------------------------------------------
  // 1. Parse the .sscore header to discover score columns.
  // -------------------------------------------------------------------------
  reterr = SizeAndInitTextStream(sscore_fname, bigstack_left() / 8, 1, &txs);
  if (unlikely(reterr)) {
    logerrprintfww("Error: --sscore: Failed to open %s .\n", sscore_fname);
    goto BestPgsFromSscore_ret_1;
  }

  {
    char* header_line = TextGet(&txs);
    if (!header_line) {
      if (!TextStreamErrcode2(&txs, &reterr)) {
        logerrputs("Error: --sscore: file is empty.\n");
        reterr = kPglRetMalformedInput;
      } else {
        TextStreamErrPrint(sscore_fname, &txs);
      }
      goto BestPgsFromSscore_ret_1;
    }

    // Strip the leading '#'.
    char* col_iter = header_line;
    if (*col_iter == '#') ++col_iter;
    // Same FID/IID probe the phenotype and covariate readers use.
    sscore_id_skip = SscoreIdSkip(col_iter);

    // Known non-score leading column names (case-insensitive match).
    // Without --score-col-nums the score columns are detected from their
    // names; see SscoreHasScoreSuffix.  --score-col-nums (literal 1-based)
    // remains the reliable form for a file that follows no convention.

    // Count total header columns.
    uint32_t total_header_cols = 0;
    {
      char* ci = col_iter;
      while (1) {
        char* ce = ci;
        while (*ce != '\t' && *ce != '\n' && *ce != '\r' && *ce != '\0') ++ce;
        ++total_header_cols;
        if (*ce == '\0' || *ce == '\n' || *ce == '\r') break;
        ci = ce + 1;
      }
    }
    if (total_header_cols == 0) {
      logerrputs("Error: --sscore: header line is empty.\n");
      reterr = kPglRetMalformedInput;
      goto BestPgsFromSscore_ret_1;
    }

    // Determine which columns are score columns.
    // score_col_bits[i] = 1 iff column (i+1) is a score column.
    uintptr_t* score_col_bits = nullptr;
    const uint32_t col_ctl = BitCtToWordCt(total_header_cols);

    if (score_col_range && score_col_range->name_ct) {
      // --score-col-nums given: use literal 1-based column numbers.
      score_col_bits = S_CAST(uintptr_t*, calloc(col_ctl, sizeof(uintptr_t)));
      if (!score_col_bits) {
        reterr = kPglRetNomem;
        goto BestPgsFromSscore_ret_1;
      }
      const uint32_t range_max = NumericRangeListMax(score_col_range);
      if (range_max > total_header_cols) {
        logerrprintfww("Error: --score-col-nums: column index %u exceeds the number of columns in %s (%u).\n",
                       range_max, sscore_fname, total_header_cols);
        free(score_col_bits);
        reterr = kPglRetInvalidCmdline;
        goto BestPgsFromSscore_ret_1;
      }
      // offset=1 means column 1 maps to bit 0.
      if (NumericRangeListToBitarr(score_col_range, total_header_cols, 1, 0, score_col_bits)) {
        logerrputs("Error: --score-col-nums: invalid range.\n");
        free(score_col_bits);
        reterr = kPglRetInvalidCmdline;
        goto BestPgsFromSscore_ret_1;
      }
    } else {
      // Auto-detect which columns hold scores.
      //
      // plink2 names them "<name>_AVG" and "<name>_SUM" (the scoreavgs and
      // scoresums column sets), so where that convention holds the suffix
      // identifies them exactly.  Preferring it is not a nicety: the DEFAULT
      // --score column set also writes 'phenos', which puts arbitrarily named
      // phenotype columns BETWEEN the IID and ALLELE_CT.  The old positional
      // rule -- "everything after the first unrecognised name is a score" --
      // assumed the non-score columns were a contiguous prefix, so a phenotype
      // swept itself, ALLELE_CT and NAMED_ALLELE_DOSAGE_SUM in as models.  When
      // the .sscore carries the phenotype under analysis, that column then
      // predicts itself and ranks first.
      //
      // NAMED_ALLELE_DOSAGE_SUM ends in _SUM but is bookkeeping, so the
      // non-score list is consulted before the suffix in both passes.
      score_col_bits = S_CAST(uintptr_t*, calloc(col_ctl, sizeof(uintptr_t)));
      if (!score_col_bits) {
        reterr = kPglRetNomem;
        goto BestPgsFromSscore_ret_1;
      }
      uint32_t first_score_col = total_header_cols;  // default: none
      {
        // Pass 1: columns named with plink2's score suffixes.
        char* ci = col_iter;
        uint32_t col_idx = 0;
        uint32_t suffix_ct = 0;
        while (1) {
          char* ce = ci;
          while (*ce != '\t' && *ce != '\n' && *ce != '\r' && *ce != '\0') ++ce;
          const uint32_t slen = S_CAST(uint32_t, ce - ci);
          if (!SscoreIsNonScoreName(ci, slen) && SscoreHasScoreSuffix(ci, slen)) {
            SetBit(col_idx, score_col_bits);
            if (first_score_col == total_header_cols) first_score_col = col_idx;
            ++suffix_ct;
          }
          if (*ce == '\0' || *ce == '\n' || *ce == '\r') break;
          ci = ce + 1;
          ++col_idx;
        }

        if (!suffix_ct) {
          // Pass 2: no column uses the suffixes -- a hand-made matrix.  Fall
          // back to the positional rule, but skip a known non-score name
          // wherever it appears rather than only while it is still a prefix.
          ci = col_iter;
          col_idx = 0;
          while (1) {
            char* ce = ci;
            while (*ce != '\t' && *ce != '\n' && *ce != '\r' && *ce != '\0') ++ce;
            const uint32_t slen = S_CAST(uint32_t, ce - ci);
            if (!SscoreIsNonScoreName(ci, slen)) {
              if (first_score_col == total_header_cols) first_score_col = col_idx;
              SetBit(col_idx, score_col_bits);
            }
            if (*ce == '\0' || *ce == '\n' || *ce == '\r') break;
            ci = ce + 1;
            ++col_idx;
          }
        }
      }
      if (first_score_col == total_header_cols) {
        logerrprintfww("Error: --sscore: no score columns found in %s. Use --score-col-nums to specify them explicitly.\n", sscore_fname);
        free(score_col_bits);
        reterr = kPglRetMalformedInput;
        goto BestPgsFromSscore_ret_1;
      }
    }

    // Count score columns and collect their names AND token indices.
    // We do both here so we don't need to re-read the header later.
    score_col_ct = PopcountWords(score_col_bits, col_ctl);
    if (score_col_ct == 0) {
      logerrputs("Error: --sscore: no score columns selected.\n");
      free(score_col_bits);
      reterr = kPglRetMalformedInput;
      goto BestPgsFromSscore_ret_1;
    }

    score_col_names = S_CAST(char**, malloc(score_col_ct * sizeof(char*)));
    if (!score_col_names) {
      free(score_col_bits);
      reterr = kPglRetNomem;
      goto BestPgsFromSscore_ret_1;
    }
    memset(score_col_names, 0, score_col_ct * sizeof(char*));

    // score_tok_idxs: 0-based token index in each data row for
    // each score column.  Stored here (within the first-header scope) and
    // referenced later in the reading loop.
    uint32_t* score_tok_idxs = S_CAST(uint32_t*, malloc(score_col_ct * sizeof(uint32_t)));
    if (!score_tok_idxs) {
      free(score_col_bits);
      reterr = kPglRetNomem;
      goto BestPgsFromSscore_ret_1;
    }

    {
      char* ci = col_iter;
      uint32_t col_idx = 0;
      uint32_t sc_idx = 0;
      while (1) {
        char* ce = ci;
        while (*ce != '\t' && *ce != '\n' && *ce != '\r' && *ce != '\0') ++ce;
        const uint32_t slen = S_CAST(uint32_t, ce - ci);
        if (IsSet(score_col_bits, col_idx)) {
          score_tok_idxs[sc_idx] = col_idx;  // token index = col_idx
          // Store the full column name; PgsIdsFromScoreColNames shortens it to
          // the PGS ID below, so an .sscore written by ScoreReport reads back
          // under the same IDs it was scored under.
          char* nm = S_CAST(char*, malloc(slen + 1));
          if (!nm) {
            free(score_col_bits);
            free(score_tok_idxs);
            reterr = kPglRetNomem;
            goto BestPgsFromSscore_ret_1;
          }
          memcpy(nm, ci, slen);
          nm[slen] = '\0';
          score_col_names[sc_idx++] = nm;
        }
        if (*ce == '\0' || *ce == '\n' || *ce == '\r') break;
        ci = ce + 1;
        ++col_idx;
      }
    }
    free(score_col_bits);
    score_col_bits = nullptr;

    reterr = PgsIdsFromScoreColNames(score_col_names, score_col_ct);
    if (unlikely(reterr)) {
      free(score_tok_idxs);
      goto BestPgsFromSscore_ret_1;
    }

    // -------------------------------------------------------------------------
    // 2. Read the .sscore data rows (IID + score values).
    // -------------------------------------------------------------------------
    // First pass to count rows; then allocate; then read.
    // Actually, use a grow-as-you-go approach with realloc.
    uint32_t sscore_cap = 128;
    sscore_iids = S_CAST(char**, malloc(sscore_cap * sizeof(char*)));
    if (!sscore_iids) {
      reterr = kPglRetNomem;
      goto BestPgsFromSscore_ret_1;
    }

    // score_tok_idxs was filled above alongside the name collection.
    // IID is always token 0.  Score columns are at their original col_idx positions.
    // No need to re-read the header.

    // score_data_rowmajor: row-major [sscore_cap x score_col_ct] during reading.
    // We use row-major so that realloc (which only extends the end) doesn't
    // corrupt per-column data when the cap doubles.  After reading we transpose
    // to the column-major layout that the rest of the code expects.
    score_data_rowmajor = S_CAST(double*, malloc(
        S_CAST(uint64_t, sscore_cap) * score_col_ct * sizeof(double)));
    if (!score_data_rowmajor) {
      free(score_tok_idxs);
      reterr = kPglRetNomem;
      goto BestPgsFromSscore_ret_1;
    }

    // Read data rows.
    for (;;) {
      char* line = TextGet(&txs);
      if (!line) break;

      // Skip blank lines.
      if (*line == '\n' || *line == '\r' || *line == '\0') continue;

      // Expand arrays if needed.
      if (sscore_sample_ct == sscore_cap) {
        sscore_cap *= 2;
        char** new_iids = S_CAST(char**, realloc(sscore_iids, sscore_cap * sizeof(char*)));
        // Row-major: row goes first, so realloc just adds rows at the end.
        double* new_scores = S_CAST(double*, realloc(score_data_rowmajor,
            S_CAST(uint64_t, sscore_cap) * score_col_ct * sizeof(double)));
        if (!new_iids || !new_scores) {
          if (new_iids) sscore_iids = new_iids;
          if (new_scores) score_data_rowmajor = new_scores;
          free(score_tok_idxs);
          reterr = kPglRetNomem;
          goto BestPgsFromSscore_ret_1;
        }
        sscore_iids = new_iids;
        score_data_rowmajor = new_scores;
      }

      // Parse all tokens up to the maximum column index we need.
      // Maximum token index needed = max(0, max(score_tok_idxs)).
      uint32_t max_tok_needed = 0;
      for (uint32_t sc = 0; sc != score_col_ct; ++sc) {
        if (score_tok_idxs[sc] > max_tok_needed) {
          max_tok_needed = score_tok_idxs[sc];
        }
      }

      // Tokenize up to max_tok_needed + 1 tokens.
      const char* tok_ptrs[1024];  // generous upper bound
      if (max_tok_needed >= 1024) {
        logerrputs("Error: --sscore: too many columns (max 1023).\n");
        free(score_tok_idxs);
        reterr = kPglRetMalformedInput;
        goto BestPgsFromSscore_ret_1;
      }
      const char* tp = line;
      for (uint32_t ti = 0; ti <= max_tok_needed; ++ti) {
        tok_ptrs[ti] = tp;
        // Advance past the token.
        while (*tp != '\t' && *tp != '\n' && *tp != '\r' && *tp != '\0') ++tp;
        if (ti < max_tok_needed) {
          if (*tp != '\t') {
            logerrprintfww("Error: --sscore: too few columns on data row %u.\n",
                           sscore_sample_ct + 1);
            free(score_tok_idxs);
            reterr = kPglRetMalformedInput;
            goto BestPgsFromSscore_ret_1;
          }
          ++tp;
        }
      }

      // Extract IID.  tok_ptrs[sscore_id_skip], not [0]: with a leading FID
      // column the IID is the second token.
      {
        const char* iid_start = tok_ptrs[sscore_id_skip];
        const char* iid_end = iid_start;
        while (*iid_end != '\t' && *iid_end != '\n' && *iid_end != '\r' && *iid_end != '\0') ++iid_end;
        const uint32_t iid_slen = S_CAST(uint32_t, iid_end - iid_start);
        char* iid_copy = S_CAST(char*, malloc(iid_slen + 1));
        if (!iid_copy) {
          free(score_tok_idxs);
          reterr = kPglRetNomem;
          goto BestPgsFromSscore_ret_1;
        }
        memcpy(iid_copy, iid_start, iid_slen);
        iid_copy[iid_slen] = '\0';
        sscore_iids[sscore_sample_ct] = iid_copy;
      }

      // Extract score values (row-major: row=sscore_sample_ct, col=sc).
      for (uint32_t sc = 0; sc != score_col_ct; ++sc) {
        const char* sv = tok_ptrs[score_tok_idxs[sc]];
        double val;
        const char* sv_end = ScanadvDouble(sv, &val);
        if (!sv_end) {
          logerrprintfww("Error: --sscore: invalid score value in column '%s' on data row %u.\n",
                         score_col_names[sc], sscore_sample_ct + 1);
          free(score_tok_idxs);
          reterr = kPglRetMalformedInput;
          goto BestPgsFromSscore_ret_1;
        }
        score_data_rowmajor[S_CAST(uint64_t, sscore_sample_ct) * score_col_ct + sc] = val;
      }

      ++sscore_sample_ct;
    }
    free(score_tok_idxs);

    if (TextStreamErrcode2(&txs, &reterr)) {
      TextStreamErrPrint(sscore_fname, &txs);
      goto BestPgsFromSscore_ret_1;
    }
    CleanupTextStream2(sscore_fname, &txs, &reterr);
    PreinitTextStream(&txs);

    if (sscore_sample_ct == 0) {
      logerrputs("Error: --sscore: no data rows in the .sscore file.\n");
      reterr = kPglRetMalformedInput;
      goto BestPgsFromSscore_ret_1;
    }

    // Transpose from row-major (sscore_sample_ct x score_col_ct) to
    // column-major (score_col_ct x sscore_sample_ct) for efficient column access.
    score_data = S_CAST(double*, malloc(
        S_CAST(uint64_t, score_col_ct) * sscore_sample_ct * sizeof(double)));
    if (!score_data) {
      free(score_data_rowmajor);
      reterr = kPglRetNomem;
      goto BestPgsFromSscore_ret_1;
    }
    for (uint32_t sc = 0; sc != score_col_ct; ++sc) {
      for (uint32_t si = 0; si != sscore_sample_ct; ++si) {
        score_data[S_CAST(uint64_t, sc) * sscore_sample_ct + si] =
            score_data_rowmajor[S_CAST(uint64_t, si) * score_col_ct + sc];
      }
    }
    free(score_data_rowmajor);
    score_data_rowmajor = nullptr;

    // Check for duplicate IIDs in .sscore using a hash table (O(n)).
    {
      const uint32_t htsize = GetHtableFastSize(sscore_sample_ct);
      uint32_t* htable = S_CAST(uint32_t*, malloc(htsize * sizeof(uint32_t)));
      if (!htable) {
        reterr = kPglRetNomem;
        goto BestPgsFromSscore_ret_1;
      }
      SetAllU32Arr(htsize, htable);
      for (uint32_t ii = 0; ii != sscore_sample_ct; ++ii) {
        const char* iid = sscore_iids[ii];
        const uint32_t slen = S_CAST(uint32_t, strlen(iid));
        const uint32_t existing = IdHtableFind(iid, S_CAST(const char* const*, sscore_iids),
                                               htable, slen, htsize);
        if (existing != UINT32_MAX) {
          logerrprintfww("Error: --sscore: duplicate IID '%s' in %s.\n",
                         iid, sscore_fname);
          free(htable);
          reterr = kPglRetMalformedInput;
          goto BestPgsFromSscore_ret_1;
        }
        HtableAddNondup(iid, slen, htsize, ii, htable);
      }
      free(htable);
    }
  }

  // -------------------------------------------------------------------------
  // 2b. Correlation matrix, if asked for.  Everything it needs is in hand
  //     already, and nothing below this point is: with --give-corr-matrix as
  //     the only request there is no phenotype file to read, which is why
  //     pheno_fname is allowed to be null on this route.
  // -------------------------------------------------------------------------
  if (corr_matrix) {
    reterr = PgsCorrReport(score_data, score_col_ct, sscore_sample_ct,
                           S_CAST(const char* const*, score_col_names),
                           outname, outname_end);
    if (unlikely(reterr)) {
      goto BestPgsFromSscore_ret_1;
    }
  }
  if (!pheno_fname) {
    goto BestPgsFromSscore_ret_1;  // reterr is kPglRetSuccess
  }

  // -------------------------------------------------------------------------
  // 3. Read the phenotype file.
  //    Plain tab-delimited; first column = IID; remaining = phenotype values.
  //    Select phenotype columns by pheno_range_list (name or number).
  //    Missing = "NA" or "-9".  Classify as QT unless all non-missing values
  //    are in {0,1} or {1,2} (plink2 convention: 1=control, 2=case).
  // -------------------------------------------------------------------------
  {
    reterr = SizeAndInitTextStream(pheno_fname, bigstack_left() / 8, 1, &txs);
    if (unlikely(reterr)) {
      logerrprintfww("Error: --sscore: Failed to open phenotype file %s .\n", pheno_fname);
      goto BestPgsFromSscore_ret_1;
    }

    char* header_line = TextGet(&txs);
    if (!header_line) {
      if (!TextStreamErrcode2(&txs, &reterr)) {
        logerrprintfww("Error: --sscore: phenotype file %s is empty.\n", pheno_fname);
        reterr = kPglRetMalformedInput;
      } else {
        TextStreamErrPrint(pheno_fname, &txs);
      }
      goto BestPgsFromSscore_ret_1;
    }

    // Parse header: IID, then the phenotype names.  A leading FID column
    // is skipped here so every count below is relative to IID.
    char* ci = header_line;
    if (*ci == '#') ++ci;
    const uint32_t pheno_id_skip = SscoreIdSkip(ci);
    SscoreSkipFidCol(&ci, pheno_id_skip);

    // Count total columns.
    uint32_t total_pheno_cols = 0;
    {
      char* tmp = ci;
      while (1) {
        while (*tmp != '\t' && *tmp != '\n' && *tmp != '\r' && *tmp != '\0') ++tmp;
        ++total_pheno_cols;
        if (*tmp == '\0' || *tmp == '\n' || *tmp == '\r') break;
        ++tmp;
      }
    }
    if (total_pheno_cols < 2) {
      logerrprintfww("Error: --sscore: phenotype file %s has no phenotype columns (only IID).\n", pheno_fname);
      reterr = kPglRetMalformedInput;
      goto BestPgsFromSscore_ret_1;
    }
    const uint32_t avail_pheno_col_ct = total_pheno_cols - 1;  // exclude IID

    // Determine which phenotype columns to use.
    uintptr_t* pheno_sel_bits = S_CAST(uintptr_t*, calloc(BitCtToWordCt(avail_pheno_col_ct), sizeof(uintptr_t)));
    if (!pheno_sel_bits) {
      reterr = kPglRetNomem;
      goto BestPgsFromSscore_ret_1;
    }

    if (pheno_range_list && pheno_range_list->name_ct) {
      // --pheno-name: select by name.
      // Re-scan the header to build a sorted name list for matching.
      // For simplicity, match by direct scan.
      char* tmp_ci = ci;
      uint32_t col_i = 0;
      // Skip IID col first.
      while (*tmp_ci != '\t' && *tmp_ci != '\n' && *tmp_ci != '\r' && *tmp_ci != '\0') ++tmp_ci;
      if (*tmp_ci == '\t') ++tmp_ci;
      while (1) {
        char* ce = tmp_ci;
        while (*ce != '\t' && *ce != '\n' && *ce != '\r' && *ce != '\0') ++ce;
        const uint32_t slen = S_CAST(uint32_t, ce - tmp_ci);
        // Check if this name matches any name in pheno_range_list.
        for (uint32_t ni = 0; ni != pheno_range_list->name_ct; ++ni) {
          const char* want = &pheno_range_list->names[ni * pheno_range_list->name_max_blen];
          if (!strcmp(want, "")) continue;  // skip empty
          const uint32_t want_len = S_CAST(uint32_t, strlen(want));
          if (want_len == slen && memequal(tmp_ci, want, slen)) {
            SetBit(col_i, pheno_sel_bits);
          }
        }
        if (*ce == '\0' || *ce == '\n' || *ce == '\r') break;
        tmp_ci = ce + 1;
        ++col_i;
      }
    } else {
      // No selection: use all phenotype columns.
      for (uint32_t ci2 = 0; ci2 != avail_pheno_col_ct; ++ci2) {
        SetBit(ci2, pheno_sel_bits);
      }
    }

    pheno_col_ct = PopcountWords(pheno_sel_bits, BitCtToWordCt(avail_pheno_col_ct));
    if (pheno_col_ct == 0) {
      logerrprintfww("Error: --sscore: no phenotype columns selected from %s.\n", pheno_fname);
      free(pheno_sel_bits);
      reterr = kPglRetInconsistentInput;
      goto BestPgsFromSscore_ret_1;
    }

    // Collect phenotype column names.
    pheno_col_names = S_CAST(char**, malloc(pheno_col_ct * sizeof(char*)));
    if (!pheno_col_names) {
      free(pheno_sel_bits);
      reterr = kPglRetNomem;
      goto BestPgsFromSscore_ret_1;
    }
    memset(pheno_col_names, 0, pheno_col_ct * sizeof(char*));

    // Indices (1-based from the data row perspective, 0-based among pheno cols)
    // of the selected phenotype columns.
    uint32_t* pheno_sel_token_idxs = S_CAST(uint32_t*, malloc(pheno_col_ct * sizeof(uint32_t)));
    if (!pheno_sel_token_idxs) {
      free(pheno_sel_bits);
      reterr = kPglRetNomem;
      goto BestPgsFromSscore_ret_1;
    }

    {
      char* tmp_ci = ci;
      uint32_t col_i = 0;
      uint32_t ph_idx = 0;
      // Skip IID col first.
      while (*tmp_ci != '\t' && *tmp_ci != '\n' && *tmp_ci != '\r' && *tmp_ci != '\0') ++tmp_ci;
      if (*tmp_ci == '\t') ++tmp_ci;
      while (1) {
        char* ce = tmp_ci;
        while (*ce != '\t' && *ce != '\n' && *ce != '\r' && *ce != '\0') ++ce;
        const uint32_t slen = S_CAST(uint32_t, ce - tmp_ci);
        if (IsSet(pheno_sel_bits, col_i)) {
          char* nm = S_CAST(char*, malloc(slen + 1));
          if (!nm) {
            free(pheno_sel_token_idxs);
            free(pheno_sel_bits);
            reterr = kPglRetNomem;
            goto BestPgsFromSscore_ret_1;
          }
          memcpy(nm, tmp_ci, slen);
          nm[slen] = '\0';
          pheno_col_names[ph_idx] = nm;
          // +1 for the IID token, +pheno_id_skip for a leading FID: the
          // data ROW still holds the FID even though `ci` skipped past it.
          pheno_sel_token_idxs[ph_idx] = col_i + 1 + pheno_id_skip;
          ++ph_idx;
        }
        if (*ce == '\0' || *ce == '\n' || *ce == '\r') break;
        tmp_ci = ce + 1;
        ++col_i;
      }
    }
    free(pheno_sel_bits);

    // Read phenotype data rows.
    // pheno_data: row-major [pheno_cap x pheno_col_ct] during reading.
    uint32_t pheno_cap = 128;
    pheno_iids = S_CAST(char**, malloc(pheno_cap * sizeof(char*)));
    pheno_data = S_CAST(double*, malloc(S_CAST(uint64_t, pheno_cap) * pheno_col_ct * sizeof(double)));
    if (!pheno_iids || !pheno_data) {
      free(pheno_sel_token_idxs);
      reterr = kPglRetNomem;
      goto BestPgsFromSscore_ret_1;
    }

    for (;;) {
      char* line = TextGet(&txs);
      if (!line) break;
      if (*line == '\n' || *line == '\r' || *line == '\0') continue;

      if (pheno_sample_ct == pheno_cap) {
        pheno_cap *= 2;
        char** new_iids = S_CAST(char**, realloc(pheno_iids, pheno_cap * sizeof(char*)));
        // Row-major: just extend at the end.
        double* new_data = S_CAST(double*, realloc(pheno_data,
            S_CAST(uint64_t, pheno_cap) * pheno_col_ct * sizeof(double)));
        if (!new_iids || !new_data) {
          if (new_iids) pheno_iids = new_iids;
          if (new_data) pheno_data = new_data;
          free(pheno_sel_token_idxs);
          reterr = kPglRetNomem;
          goto BestPgsFromSscore_ret_1;
        }
        pheno_iids = new_iids;
        pheno_data = new_data;
      }

      // Find max token index needed.
      uint32_t max_tok = 0;
      for (uint32_t ph = 0; ph != pheno_col_ct; ++ph) {
        if (pheno_sel_token_idxs[ph] > max_tok) {
          max_tok = pheno_sel_token_idxs[ph];
        }
      }

      // Tokenize up to max_tok+1.
      const char* tok_ptrs[1024];
      if (max_tok >= 1024) {
        free(pheno_sel_token_idxs);
        reterr = kPglRetMalformedInput;
        goto BestPgsFromSscore_ret_1;
      }
      const char* tp = line;
      for (uint32_t ti = 0; ti <= max_tok; ++ti) {
        tok_ptrs[ti] = tp;
        while (*tp != '\t' && *tp != '\n' && *tp != '\r' && *tp != '\0') ++tp;
        if (ti < max_tok) {
          if (*tp != '\t') {
            logerrprintfww("Error: --sscore: too few columns on phenotype data row %u.\n",
                           pheno_sample_ct + 1);
            free(pheno_sel_token_idxs);
            reterr = kPglRetMalformedInput;
            goto BestPgsFromSscore_ret_1;
          }
          ++tp;
        }
      }

      // Extract IID.
      {
        const char* iid_start = tok_ptrs[pheno_id_skip];
        const char* iid_end = iid_start;
        while (*iid_end != '\t' && *iid_end != '\n' && *iid_end != '\r' && *iid_end != '\0') ++iid_end;
        const uint32_t iid_slen = S_CAST(uint32_t, iid_end - iid_start);
        char* iid_copy = S_CAST(char*, malloc(iid_slen + 1));
        if (!iid_copy) {
          free(pheno_sel_token_idxs);
          reterr = kPglRetNomem;
          goto BestPgsFromSscore_ret_1;
        }
        memcpy(iid_copy, iid_start, iid_slen);
        iid_copy[iid_slen] = '\0';
        pheno_iids[pheno_sample_ct] = iid_copy;
      }

      // Extract phenotype values.
      for (uint32_t ph = 0; ph != pheno_col_ct; ++ph) {
        const char* sv = tok_ptrs[pheno_sel_token_idxs[ph]];
        const char* sv_end = sv;
        while (*sv_end != '\t' && *sv_end != '\n' && *sv_end != '\r' && *sv_end != '\0') ++sv_end;
        const uint32_t sv_slen = S_CAST(uint32_t, sv_end - sv);

        double val;
        const char* scan_end = ScanadvDouble(sv, &val);
        if (!scan_end || (scan_end != sv_end && *scan_end != '\t' && *scan_end != '\n' && *scan_end != '\r' && *scan_end != '\0')) {
          // Not a number: check for NA / -9 (missing).
          // strequal_k handles "NA"; also accept "-9".
          if ((sv_slen == 2 && memequal(sv, "NA", 2)) ||
              (sv_slen == 2 && memequal(sv, "-9", 2))) {
            val = 0.0 / 0.0;  // NaN = missing
          } else {
            val = 0.0 / 0.0;  // treat unknown strings as missing
          }
        }
        pheno_data[S_CAST(uint64_t, pheno_sample_ct) * pheno_col_ct + ph] = val;
      }

      ++pheno_sample_ct;
    }
    free(pheno_sel_token_idxs);

    if (TextStreamErrcode2(&txs, &reterr)) {
      TextStreamErrPrint(pheno_fname, &txs);
      goto BestPgsFromSscore_ret_1;
    }
    CleanupTextStream2(pheno_fname, &txs, &reterr);
    PreinitTextStream(&txs);

    // Check for duplicate IIDs in phenotype file using a hash table (O(n)).
    {
      const uint32_t htsize = GetHtableFastSize(pheno_sample_ct);
      uint32_t* htable = S_CAST(uint32_t*, malloc(htsize * sizeof(uint32_t)));
      if (!htable) {
        reterr = kPglRetNomem;
        goto BestPgsFromSscore_ret_1;
      }
      SetAllU32Arr(htsize, htable);
      for (uint32_t ii = 0; ii != pheno_sample_ct; ++ii) {
        const char* iid = pheno_iids[ii];
        const uint32_t slen = S_CAST(uint32_t, strlen(iid));
        const uint32_t existing = IdHtableFind(iid, S_CAST(const char* const*, pheno_iids),
                                               htable, slen, htsize);
        if (existing != UINT32_MAX) {
          logerrprintfww("Error: --sscore: duplicate IID '%s' in phenotype file %s.\n",
                         iid, pheno_fname);
          free(htable);
          reterr = kPglRetMalformedInput;
          goto BestPgsFromSscore_ret_1;
        }
        HtableAddNondup(iid, slen, htsize, ii, htable);
      }
      free(htable);
    }

    // Transpose pheno_data from row-major to column-major.
    {
      double* pheno_col_major = S_CAST(double*, malloc(
          S_CAST(uint64_t, pheno_col_ct) * pheno_sample_ct * sizeof(double)));
      if (!pheno_col_major) {
        reterr = kPglRetNomem;
        goto BestPgsFromSscore_ret_1;
      }
      for (uint32_t ph = 0; ph != pheno_col_ct; ++ph) {
        for (uint32_t si = 0; si != pheno_sample_ct; ++si) {
          pheno_col_major[S_CAST(uint64_t, ph) * pheno_sample_ct + si] =
              pheno_data[S_CAST(uint64_t, si) * pheno_col_ct + ph];
        }
      }
      free(pheno_data);
      pheno_data = pheno_col_major;
    }
  }

  // -------------------------------------------------------------------------
  // 4. Read the covariate file (optional).
  // -------------------------------------------------------------------------
  if (covar_fname) {
    reterr = SizeAndInitTextStream(covar_fname, bigstack_left() / 8, 1, &txs);
    if (unlikely(reterr)) {
      logerrprintfww("Error: --sscore: Failed to open covariate file %s .\n", covar_fname);
      goto BestPgsFromSscore_ret_1;
    }

    char* header_line = TextGet(&txs);
    if (!header_line) {
      if (!TextStreamErrcode2(&txs, &reterr)) {
        logerrprintfww("Error: --sscore: covariate file %s is empty.\n", covar_fname);
        reterr = kPglRetMalformedInput;
      } else {
        TextStreamErrPrint(covar_fname, &txs);
      }
      goto BestPgsFromSscore_ret_1;
    }

    char* ci = header_line;
    if (*ci == '#') ++ci;

    const uint32_t covar_id_skip = SscoreIdSkip(ci);
    SscoreSkipFidCol(&ci, covar_id_skip);

    // Count covariate columns (all columns after IID).
    uint32_t total_cv_cols = 0;
    {
      char* tmp = ci;
      while (1) {
        while (*tmp != '\t' && *tmp != '\n' && *tmp != '\r' && *tmp != '\0') ++tmp;
        ++total_cv_cols;
        if (*tmp == '\0' || *tmp == '\n' || *tmp == '\r') break;
        ++tmp;
      }
    }
    if (total_cv_cols < 2) {
      logerrprintfww("Error: --sscore: covariate file %s has no covariate columns.\n", covar_fname);
      reterr = kPglRetMalformedInput;
      goto BestPgsFromSscore_ret_1;
    }
    covar_col_ct = total_cv_cols - 1;

    // Collect covariate column names.
    covar_col_names = S_CAST(char**, malloc(covar_col_ct * sizeof(char*)));
    if (!covar_col_names) {
      reterr = kPglRetNomem;
      goto BestPgsFromSscore_ret_1;
    }
    memset(covar_col_names, 0, covar_col_ct * sizeof(char*));

    {
      char* tmp_ci = ci;
      // Skip IID col.
      while (*tmp_ci != '\t' && *tmp_ci != '\n' && *tmp_ci != '\r' && *tmp_ci != '\0') ++tmp_ci;
      if (*tmp_ci == '\t') ++tmp_ci;
      for (uint32_t cv = 0; cv != covar_col_ct; ++cv) {
        char* ce = tmp_ci;
        while (*ce != '\t' && *ce != '\n' && *ce != '\r' && *ce != '\0') ++ce;
        const uint32_t slen = S_CAST(uint32_t, ce - tmp_ci);
        char* nm = S_CAST(char*, malloc(slen + 1));
        if (!nm) {
          reterr = kPglRetNomem;
          goto BestPgsFromSscore_ret_1;
        }
        memcpy(nm, tmp_ci, slen);
        nm[slen] = '\0';
        covar_col_names[cv] = nm;
        if (*ce == '\0' || *ce == '\n' || *ce == '\r') break;
        tmp_ci = ce + 1;
      }
    }

    // Read covariate data rows.
    // covar_data: row-major [covar_cap x covar_col_ct] during reading.
    // Transposed to column-major [covar_col_ct x covar_sample_ct] after reading.
    uint32_t covar_cap = 128;
    covar_iids = S_CAST(char**, malloc(covar_cap * sizeof(char*)));
    covar_data = S_CAST(double*, malloc(S_CAST(uint64_t, covar_cap) * covar_col_ct * sizeof(double)));
    if (!covar_iids || !covar_data) {
      reterr = kPglRetNomem;
      goto BestPgsFromSscore_ret_1;
    }

    for (;;) {
      char* line = TextGet(&txs);
      if (!line) break;
      if (*line == '\n' || *line == '\r' || *line == '\0') continue;

      if (covar_sample_ct == covar_cap) {
        covar_cap *= 2;
        char** new_iids = S_CAST(char**, realloc(covar_iids, covar_cap * sizeof(char*)));
        // Row-major: rows at the end, so realloc just adds capacity.
        double* new_data = S_CAST(double*, realloc(covar_data,
            S_CAST(uint64_t, covar_cap) * covar_col_ct * sizeof(double)));
        if (!new_iids || !new_data) {
          if (new_iids) covar_iids = new_iids;
          if (new_data) covar_data = new_data;
          reterr = kPglRetNomem;
          goto BestPgsFromSscore_ret_1;
        }
        covar_iids = new_iids;
        covar_data = new_data;
      }

      // Parse up to covar_col_ct + 1 tokens.
      const char* tok_ptrs[1024];
      const uint32_t need_ct = covar_col_ct + 1 + covar_id_skip;
      if (need_ct > 1024) {
        reterr = kPglRetMalformedInput;
        goto BestPgsFromSscore_ret_1;
      }
      const char* tp = line;
      for (uint32_t ti = 0; ti < need_ct; ++ti) {
        tok_ptrs[ti] = tp;
        while (*tp != '\t' && *tp != '\n' && *tp != '\r' && *tp != '\0') ++tp;
        if (ti < need_ct - 1) {
          if (*tp != '\t') {
            logerrprintfww("Error: --sscore: too few columns on covariate data row %u.\n",
                           covar_sample_ct + 1);
            reterr = kPglRetMalformedInput;
            goto BestPgsFromSscore_ret_1;
          }
          ++tp;
        }
      }

      // IID.
      {
        const char* iid_end = tok_ptrs[covar_id_skip];
        while (*iid_end != '\t' && *iid_end != '\n' && *iid_end != '\r' && *iid_end != '\0') ++iid_end;
        const uint32_t iid_slen = S_CAST(uint32_t, iid_end - tok_ptrs[covar_id_skip]);
        char* iid_copy = S_CAST(char*, malloc(iid_slen + 1));
        if (!iid_copy) {
          reterr = kPglRetNomem;
          goto BestPgsFromSscore_ret_1;
        }
        memcpy(iid_copy, tok_ptrs[covar_id_skip], iid_slen);
        iid_copy[iid_slen] = '\0';
        covar_iids[covar_sample_ct] = iid_copy;
      }

      // Covariate values (row-major: row=covar_sample_ct, col=cv).
      for (uint32_t cv = 0; cv != covar_col_ct; ++cv) {
        const char* sv = tok_ptrs[cv + 1 + covar_id_skip];
        double val;
        const char* scan_end = ScanadvDouble(sv, &val);
        if (!scan_end) {
          // Missing / non-numeric covariate treated as NaN.
          val = 0.0 / 0.0;
        }
        covar_data[S_CAST(uint64_t, covar_sample_ct) * covar_col_ct + cv] = val;
      }

      ++covar_sample_ct;
    }

    if (TextStreamErrcode2(&txs, &reterr)) {
      TextStreamErrPrint(covar_fname, &txs);
      goto BestPgsFromSscore_ret_1;
    }
    CleanupTextStream2(covar_fname, &txs, &reterr);
    PreinitTextStream(&txs);

    // Check for duplicate IIDs in covariate file using a hash table (O(n)).
    {
      const uint32_t htsize = GetHtableFastSize(covar_sample_ct);
      uint32_t* htable = S_CAST(uint32_t*, malloc(htsize * sizeof(uint32_t)));
      if (!htable) {
        reterr = kPglRetNomem;
        goto BestPgsFromSscore_ret_1;
      }
      SetAllU32Arr(htsize, htable);
      for (uint32_t ii = 0; ii != covar_sample_ct; ++ii) {
        const char* iid = covar_iids[ii];
        const uint32_t slen = S_CAST(uint32_t, strlen(iid));
        const uint32_t existing = IdHtableFind(iid, S_CAST(const char* const*, covar_iids),
                                               htable, slen, htsize);
        if (existing != UINT32_MAX) {
          logerrprintfww("Error: --sscore: duplicate IID '%s' in covariate file %s.\n",
                         iid, covar_fname);
          free(htable);
          reterr = kPglRetMalformedInput;
          goto BestPgsFromSscore_ret_1;
        }
        HtableAddNondup(iid, slen, htsize, ii, htable);
      }
      free(htable);
    }

    // Transpose covar_data from row-major to column-major.
    {
      double* covar_col_major = S_CAST(double*, malloc(
          S_CAST(uint64_t, covar_col_ct) * covar_sample_ct * sizeof(double)));
      if (!covar_col_major) {
        reterr = kPglRetNomem;
        goto BestPgsFromSscore_ret_1;
      }
      for (uint32_t cv = 0; cv != covar_col_ct; ++cv) {
        for (uint32_t si = 0; si != covar_sample_ct; ++si) {
          covar_col_major[S_CAST(uint64_t, cv) * covar_sample_ct + si] =
              covar_data[S_CAST(uint64_t, si) * covar_col_ct + cv];
        }
      }
      free(covar_data);
      covar_data = covar_col_major;
    }
  }

  // -------------------------------------------------------------------------
  // 5. Inner join on IID across all three inputs using hash-table lookups.
  //    Build hash tables over pheno IIDs (and covar IIDs if present), then
  //    probe with each sscore IID in order. One pass, O(n) total.
  // -------------------------------------------------------------------------
  {
    // Build hash table over pheno IIDs.
    uint32_t pheno_htsize = GetHtableFastSize(pheno_sample_ct);
    uint32_t* pheno_htable = S_CAST(uint32_t*, malloc(pheno_htsize * sizeof(uint32_t)));
    if (!pheno_htable) {
      reterr = kPglRetNomem;
      goto BestPgsFromSscore_ret_1;
    }
    SetAllU32Arr(pheno_htsize, pheno_htable);
    for (uint32_t ii = 0; ii != pheno_sample_ct; ++ii) {
      const char* iid = pheno_iids[ii];
      HtableAddNondup(iid, S_CAST(uint32_t, strlen(iid)), pheno_htsize, ii, pheno_htable);
    }

    // Build hash table over covar IIDs (if a covariate file was given).
    uint32_t covar_htsize = 0;
    uint32_t* covar_htable = nullptr;
    if (covar_fname) {
      covar_htsize = GetHtableFastSize(covar_sample_ct);
      covar_htable = S_CAST(uint32_t*, malloc(covar_htsize * sizeof(uint32_t)));
      if (!covar_htable) {
        free(pheno_htable);
        reterr = kPglRetNomem;
        goto BestPgsFromSscore_ret_1;
      }
      SetAllU32Arr(covar_htsize, covar_htable);
      for (uint32_t ii = 0; ii != covar_sample_ct; ++ii) {
        const char* iid = covar_iids[ii];
        HtableAddNondup(iid, S_CAST(uint32_t, strlen(iid)), covar_htsize, ii, covar_htable);
      }
    }

    // First pass: count the intersection size.
    // Probe sscore IIDs against pheno (and covar) hash tables.
    for (uint32_t ii = 0; ii != sscore_sample_ct; ++ii) {
      const char* iid = sscore_iids[ii];
      const uint32_t slen = S_CAST(uint32_t, strlen(iid));
      const uint32_t pi = IdHtableFind(iid, S_CAST(const char* const*, pheno_iids),
                                       pheno_htable, slen, pheno_htsize);
      if (pi == UINT32_MAX) continue;
      if (covar_fname) {
        const uint32_t ci = IdHtableFind(iid, S_CAST(const char* const*, covar_iids),
                                         covar_htable, slen, covar_htsize);
        if (ci == UINT32_MAX) continue;
      }
      ++joined_ct;
    }

    // Warn about dropped samples.
    if (joined_ct < sscore_sample_ct || joined_ct < pheno_sample_ct) {
      uint32_t dropped = sscore_sample_ct + pheno_sample_ct - 2 * joined_ct;
      if (covar_fname) {
        logerrprintfww("Warning: --sscore: %u samples in the .sscore, %u in the phenotype file, %u in the covariate file, %u after joining on IID (%u dropped).\n",
                       sscore_sample_ct, pheno_sample_ct, covar_sample_ct, joined_ct,
                       sscore_sample_ct - joined_ct + pheno_sample_ct - joined_ct);
      } else {
        logerrprintfww("Warning: --sscore: %u samples in the .sscore, %u in the phenotype file, %u after joining on IID (%u dropped).\n",
                       sscore_sample_ct, pheno_sample_ct, joined_ct, dropped);
      }
    }

    if (joined_ct == 0) {
      logerrputs("Error: --sscore: no samples remain after joining on IID. Check whether IID conventions differ between files (e.g. FID_IID vs IID-only).\n");
      free(covar_htable);
      free(pheno_htable);
      reterr = kPglRetInconsistentInput;
      goto BestPgsFromSscore_ret_1;
    }

    // Second pass: fill the index arrays in sscore order (preserves output order).
    uint32_t* join_pheno_idx = S_CAST(uint32_t*, malloc(joined_ct * sizeof(uint32_t)));
    uint32_t* join_covar_idx = covar_fname ? S_CAST(uint32_t*, malloc(joined_ct * sizeof(uint32_t))) : nullptr;
    uint32_t* join_score_idx = S_CAST(uint32_t*, malloc(joined_ct * sizeof(uint32_t)));
    if (!join_pheno_idx || !join_score_idx || (covar_fname && !join_covar_idx)) {
      free(covar_htable);
      free(pheno_htable);
      free(join_pheno_idx);
      free(join_covar_idx);
      free(join_score_idx);
      reterr = kPglRetNomem;
      goto BestPgsFromSscore_ret_1;
    }

    {
      uint32_t out_idx = 0;
      for (uint32_t ii = 0; ii != sscore_sample_ct; ++ii) {
        const char* iid = sscore_iids[ii];
        const uint32_t slen = S_CAST(uint32_t, strlen(iid));
        const uint32_t pi = IdHtableFind(iid, S_CAST(const char* const*, pheno_iids),
                                         pheno_htable, slen, pheno_htsize);
        if (pi == UINT32_MAX) continue;
        if (covar_fname) {
          const uint32_t ci = IdHtableFind(iid, S_CAST(const char* const*, covar_iids),
                                           covar_htable, slen, covar_htsize);
          if (ci == UINT32_MAX) continue;
          join_covar_idx[out_idx] = ci;
        }
        join_score_idx[out_idx] = ii;
        join_pheno_idx[out_idx] = pi;
        ++out_idx;
      }
    }

    free(covar_htable);
    free(pheno_htable);

    // -------------------------------------------------------------------------
    // 6. Classify phenotypes and validate.
    // -------------------------------------------------------------------------
    pheno_types = S_CAST(PhenoDtype*, malloc(pheno_col_ct * sizeof(PhenoDtype)));
    keep_idx = S_CAST(uint32_t*, malloc(joined_ct * sizeof(uint32_t)));
    if (!pheno_types || !keep_idx) {
      free(join_pheno_idx);
      free(join_covar_idx);
      free(join_score_idx);
      reterr = kPglRetNomem;
      goto BestPgsFromSscore_ret_1;
    }

    // The shared guards read y, so this route materialises it here too.  Sized
    // for the largest possible complete-case set, which is joined_ct; section 7
    // refills the same buffer alongside X.
    y_complete = S_CAST(double*, malloc(joined_ct * sizeof(double)));
    if (!y_complete) {
      free(join_pheno_idx);
      free(join_covar_idx);
      free(join_score_idx);
      reterr = kPglRetNomem;
      goto BestPgsFromSscore_ret_1;
    }

    // Every skip below is a warning, not an error, matching the --pfile route:
    // one unusable phenotype should not take the others down with it.  Only
    // "none of them are usable" is fatal.  Both routes now reach those skips
    // through ValidatePhenoBlock, so neither can drift out from under the other.
    uint32_t usable_pheno_ct = 0;
    for (uint32_t ph = 0; ph != pheno_col_ct; ++ph) {
      pheno_types[ph] = kPhenoDtypeCat;  // sentinel for "unusable"

      const uint32_t complete_ct =
          SscoreCompleteCases(pheno_data, pheno_sample_ct, join_pheno_idx, ph,
                              covar_data, covar_sample_ct, join_covar_idx,
                              covar_col_ct, joined_ct, keep_idx);
      // This route can say how many rows it dropped, because it knows what it
      // joined; the --pfile route has no comparable "before" count.  It has to
      // precede the skips below, which say nothing about the arithmetic.
      if (complete_ct < joined_ct) {
        logerrprintfww("Warning: --sscore: phenotype '%s': %u of the %u joined samples have a missing phenotype or covariate value and are excluded; %u remain.\n",
                       pheno_col_names[ph], joined_ct - complete_ct, joined_ct,
                       complete_ct);
      }

      // Classify: QT unless every value is in {0,1} or {1,2}.  This matches
      // LoadPhenos: initially assume CC; if any value is not in {0, 1, 2},
      // it's QT.
      //
      // plink2 convention: 2=case, anything else control -- so a 0/1-coded
      // column reads as all-controls here and falls through to QT, where an OLS
      // fit on the 0/1 outcome gives the same R2 family the CC path would have
      // taken from its own OLS companion.  That is also why the shared block's
      // case/control guard cannot fire on this route: a phenotype with no cases
      // or no controls is never labelled CC to begin with.  Vacuous is the
      // right state for it -- it stops being vacuous the day --sscore is told a
      // type instead of guessing one.
      uint32_t all_01_or_12 = 1;
      uint32_t n_case = 0;
      for (uint32_t ii = 0; ii != complete_ct; ++ii) {
        const double val = pheno_data[S_CAST(uint64_t, ph) * pheno_sample_ct + join_pheno_idx[keep_idx[ii]]];
        if (val != 0.0 && val != 1.0 && val != 2.0) {
          all_01_or_12 = 0;
          break;
        }
        if (val > 1.5) {
          ++n_case;
        }
      }
      const PhenoDtype dtype =
          (all_01_or_12 && n_case && (n_case != complete_ct))? kPhenoDtypeCc
                                                             : kPhenoDtypeQt;

      // y as the guards want it: the raw value for QT, 0/1 for case/control.
      for (uint32_t ii = 0; ii != complete_ct; ++ii) {
        const double val = pheno_data[S_CAST(uint64_t, ph) * pheno_sample_ct + join_pheno_idx[keep_idx[ii]]];
        y_complete[ii] = (dtype == kPhenoDtypeCc)? ((val > 1.5)? 1.0 : 0.0)
                                                 : val;
      }

      PgsPhenoBlock pb;
      pb.y = y_complete;
      pb.n = complete_ct;
      pb.covar_ct = covar_col_ct;
      pb.type = dtype;
      pb.name = pheno_col_names[ph];
      pb.route = "--sscore";

      uint32_t usable;
      reterr = ValidatePhenoBlock(&pb, fold_ct, &usable);
      if (unlikely(reterr)) {
        free(join_pheno_idx);
        free(join_covar_idx);
        free(join_score_idx);
        goto BestPgsFromSscore_ret_1;
      }
      if (!usable) {
        continue;
      }

      pheno_types[ph] = dtype;
      ++usable_pheno_ct;
    }

    if (usable_pheno_ct == 0) {
      logerrputs("Error: --sscore: no usable phenotype found.\n");
      free(join_pheno_idx);
      free(join_covar_idx);
      free(join_score_idx);
      reterr = kPglRetInconsistentInput;
      goto BestPgsFromSscore_ret_1;
    }

    // -------------------------------------------------------------------------
    // 7. Allocate working arrays and call FitWritePhenoBlock per phenotype.
    // -------------------------------------------------------------------------
    // Sized for the largest possible complete-case set, which is joined_ct.
    // y_complete was allocated above section 6, whose guards already needed it.
    all_rows = S_CAST(BestPgsRow*, malloc(S_CAST(uint64_t, usable_pheno_ct) * score_col_ct * sizeof(BestPgsRow)));
    X_complete = S_CAST(double*, malloc(S_CAST(uint64_t, joined_ct) * (1 + covar_col_ct) * sizeof(double)));
    score_block = S_CAST(double*, malloc(S_CAST(uint64_t, score_col_ct) * joined_ct * sizeof(double)));
    if (!all_rows || !X_complete || !score_block) {
      free(join_pheno_idx);
      free(join_covar_idx);
      free(join_score_idx);
      reterr = kPglRetNomem;
      goto BestPgsFromSscore_ret_1;
    }

    // Pre-build covar_names_arr.
    const char** covar_names_arr = nullptr;
    if (covar_col_ct) {
      covar_names_arr = S_CAST(const char**, malloc(covar_col_ct * sizeof(char*)));
      if (!covar_names_arr) {
        free(join_pheno_idx);
        free(join_covar_idx);
        free(join_score_idx);
        reterr = kPglRetNomem;
        goto BestPgsFromSscore_ret_1;
      }
      for (uint32_t cv = 0; cv != covar_col_ct; ++cv) {
        covar_names_arr[cv] = covar_col_names[cv];
      }
    }

    // --pgs-pheno-pairs: both name sets are known here.  Validated against
    // every phenotype COLUMN, not just the usable ones -- a phenotype that
    // exists but was skipped is not a typo and must not be reported as one.
    if (pairs_fname) {
      reterr = LoadPgsPhenoPairs(pairs_fname, &pairs);
      if (!reterr) {
        pair_keep_idxs = S_CAST(uint32_t*, malloc(score_col_ct * sizeof(uint32_t)));
        if (!pair_keep_idxs) {
          reterr = kPglRetNomem;
        }
      }
      if (!reterr) {
        reterr = PgsPairsValidate(&pairs,
                                  S_CAST(const char* const*, pheno_col_names), pheno_col_ct,
                                  S_CAST(const char* const*, score_col_names), score_col_ct);
      }
      if (unlikely(reterr)) {
        free(covar_names_arr);
        free(join_pheno_idx);
        free(join_covar_idx);
        free(join_score_idx);
        goto BestPgsFromSscore_ret_1;
      }
    }

    uint32_t out_pheno_pos = 0;
    for (uint32_t ph = 0; ph != pheno_col_ct; ++ph) {
      if (pheno_types[ph] == kPhenoDtypeCat) continue;  // skipped

      const PhenoDtype type_code = pheno_types[ph];
      const char* pheno_name = pheno_col_names[ph];

      // The same complete-case set section 6 classified this phenotype on.
      // Recomputed rather than kept: one pass over the joined rows, against
      // pheno_col_ct index arrays of joined_ct entries to hold them all.
      const uint32_t complete_ct =
          SscoreCompleteCases(pheno_data, pheno_sample_ct, join_pheno_idx, ph,
                              covar_data, covar_sample_ct, join_covar_idx,
                              covar_col_ct, joined_ct, keep_idx);

      // Build y_complete (phenotype values; remap 2/other -> 1/0 for CC).
      for (uint32_t ii = 0; ii != complete_ct; ++ii) {
        double val = pheno_data[S_CAST(uint64_t, ph) * pheno_sample_ct + join_pheno_idx[keep_idx[ii]]];
        if (type_code == kPhenoDtypeCc) {
          val = (val > 1.5) ? 1.0 : 0.0;
        }
        y_complete[ii] = val;
      }

      // Build X_complete: intercept + covariates.
      for (uint32_t ii = 0; ii != complete_ct; ++ii) {
        X_complete[ii] = 1.0;  // intercept
      }
      for (uint32_t cv = 0; cv != covar_col_ct; ++cv) {
        for (uint32_t ii = 0; ii != complete_ct; ++ii) {
          X_complete[S_CAST(uint64_t, 1 + cv) * complete_ct + ii] =
              covar_data[S_CAST(uint64_t, cv) * covar_sample_ct + join_covar_idx[keep_idx[ii]]];
        }
      }

      // The score columns have no missingness of their own, but the rows they
      // are taken over are this phenotype's, so score_block is per-phenotype
      // too.
      for (uint32_t sc = 0; sc != score_col_ct; ++sc) {
        for (uint32_t ii = 0; ii != complete_ct; ++ii) {
          score_block[S_CAST(uint64_t, sc) * complete_ct + ii] =
              score_data[S_CAST(uint64_t, sc) * sscore_sample_ct + join_score_idx[keep_idx[ii]]];
        }
      }

      uint32_t keep_ct = 0;
      if (pairs_fname) {
        keep_ct = PgsPairsKeepIdxs(&pairs, pheno_name,
                                   S_CAST(const char* const*, score_col_names),
                                   score_col_ct, pair_keep_idxs);
      }
      BestPgsRow* pheno_rows = &all_rows[S_CAST(uint64_t, out_pheno_pos) * score_col_ct];
      reterr = FitWritePhenoBlock(complete_ct, covar_col_ct, y_complete, X_complete,
                                  score_block,
                                  score_col_names, score_col_ct,
                                  covar_names_arr, pheno_name, type_code,
                                  prevalence, fold_ct, fold_seed,
                                  pair_keep_idxs, keep_ct, &outfile,
                                  outname, outname_end, pheno_rows);
      if (unlikely(reterr)) {
        free(covar_names_arr);
        free(join_pheno_idx);
        free(join_covar_idx);
        free(join_score_idx);
        goto BestPgsFromSscore_ret_1;
      }
      ++out_pheno_pos;
    }

    free(covar_names_arr);
    free(join_pheno_idx);
    free(join_covar_idx);
    free(join_score_idx);
  }

  // See the matching guard on the --pfile route: under a filter every named
  // phenotype can turn out to be one that was skipped, and a header-only
  // .best_pgs would read as a complete report of nothing.
  if (unlikely(!outfile)) {
    logerrputs("Error: --pgs-pheno-pairs: no pair could be evaluated; every phenotype named in\nthe pairs file was skipped.\n");
    reterr = kPglRetInconsistentInput;
    goto BestPgsFromSscore_ret_1;
  }
  if (unlikely(fclose_null(&outfile))) {
    reterr = kPglRetWriteFail;
    goto BestPgsFromSscore_ret_1;
  }

  *outname_end = '\0';
  logprintfww("--give-best-pgs (--sscore): Results written to %s.best_pgs .\n", outname);

 BestPgsFromSscore_ret_1:
  PgsPhenoPairsFree(&pairs);
  free(pair_keep_idxs);
  if (outfile) fclose(outfile);
  CleanupTextStream2("--sscore input", &txs, &reterr);

  free(score_block);
  free(X_complete);
  free(y_complete);
  free(all_rows);
  free(keep_idx);
  free(pheno_types);
  free(score_data_rowmajor);

  free(covar_data);
  if (covar_iids) {
    for (uint32_t si = 0; si != covar_sample_ct; ++si) {
      free(covar_iids[si]);
    }
    free(covar_iids);
  }
  if (covar_col_names) {
    for (uint32_t cv = 0; cv != covar_col_ct; ++cv) {
      free(covar_col_names[cv]);
    }
    free(covar_col_names);
  }

  free(pheno_data);
  if (pheno_iids) {
    for (uint32_t si = 0; si != pheno_sample_ct; ++si) {
      free(pheno_iids[si]);
    }
    free(pheno_iids);
  }
  if (pheno_col_names) {
    for (uint32_t ph = 0; ph != pheno_col_ct; ++ph) {
      free(pheno_col_names[ph]);
    }
    free(pheno_col_names);
  }

  free(score_data);
  if (sscore_iids) {
    for (uint32_t si = 0; si != sscore_sample_ct; ++si) {
      free(sscore_iids[si]);
    }
    free(sscore_iids);
  }
  if (score_col_names) {
    for (uint32_t sc = 0; sc != score_col_ct; ++sc) {
      free(score_col_names[sc]);
    }
    free(score_col_names);
  }

  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
