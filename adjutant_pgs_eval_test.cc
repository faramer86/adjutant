// Standalone unit tests for adjutant_pgs_eval.
//
// The project otherwise tests by golden-diff against pristine plink2 output,
// which cannot reach a pure function that never touches a .pgen.  Regression
// code needs its own check: "it compiled and converged" is not evidence that
// the coefficients are right.
//
// Build and run:  make pgs_eval_test && ./bin/pgs_eval_test

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "adjutant_pgs_eval.h"
#include "adjutant_pgs_match.h"

static int g_failures = 0;
static int g_checks = 0;

static void CheckClose(const char* what, double got, double want, double tol) {
  ++g_checks;
  // NaN first: `NaN > tol` is false, so without this every tolerance passes
  // silently on a NaN -- which is exactly what an unpopulated fixture key or a
  // degenerate fit produces.
  if (isnan(got) || isnan(want)) {
    printf("  FAIL %-38s got %.12g  want %.12g  (NaN)\n", what, got, want);
    ++g_failures;
    return;
  }
  const double diff = fabs(got - want);
  // Relative tolerance when want is non-zero, absolute only for exact zero.
  const double scale = (want != 0.0) ? fabs(want) : 1.0;
  if (diff / scale > tol) {
    printf("  FAIL %-38s got %.12g  want %.12g  (diff %.3g)\n",
           what, got, want, diff);
    ++g_failures;
  }
}

// A confidence interval that excludes its own point estimate is the failure
// mode the straddle-zero handling exists to prevent, so every fixture row
// asserts containment rather than only the row built to provoke it.  The
// negated form also fails on NaN.
static void CheckCiContains(const char* what, double lo, double point, double hi) {
  ++g_checks;
  if (!(lo <= point && point <= hi)) {
    printf("  FAIL %-38s %.12g <= %.12g <= %.12g violated\n",
           what, lo, point, hi);
    ++g_failures;
  }
}

// Overlap of two intervals as a fraction of the wider one.  The bootstrap is a
// check on the method, not a tolerance target: it pins down that the analytic
// interval sits in the right place and is the right width, which is all a
// closed form can be asked to agree with a resampling estimate about.
static double CiOverlapFrac(double lo1, double hi1, double lo2, double hi2) {
  const double lo = (lo1 > lo2) ? lo1 : lo2;
  const double hi = (hi1 < hi2) ? hi1 : hi2;
  const double overlap = (hi > lo) ? (hi - lo) : 0.0;
  const double w1 = hi1 - lo1;
  const double w2 = hi2 - lo2;
  return overlap / ((w1 > w2) ? w1 : w2);
}

// The near-null PGS behind the straddle-zero case.  Integer arithmetic, so the
// values are bit-identical to the ones generate_ci.R computed the reference
// interval from.
static double NullPgs(uint32_t ii) {
  return (double)((ii * 37u) % 101u);
}

// Least-squares fit of yv on the columns of `design` (column-major,
// row_ct x col_ct), restricted to the rows where use[] is nonzero, by forming
// the normal equations and solving them with Gauss-Jordan elimination.
//
// Deliberately a different route from anything in adjutant_pgs_eval.cc: no
// Frisch-Waugh, no LAPACK, no subsetting shortcut.  It exists so the
// cross-validation path's held-out predictions can be checked against a fit
// that shares none of its algebra.  Returns 1 if the system is singular.
static int SolveSubsetOls(const double* design, const double* yv,
                          uint32_t row_ct, uint32_t col_ct,
                          const unsigned char* use, double* beta_out) {
  const uint32_t kMaxCol = 8;
  if (col_ct > kMaxCol) { return 1; }
  double aa[kMaxCol * (kMaxCol + 1)];  // augmented [X'X | X'y], row-major
  for (uint32_t r = 0; r != col_ct; ++r) {
    for (uint32_t c = 0; c != col_ct; ++c) {
      double acc = 0.0;
      for (uint32_t i = 0; i != row_ct; ++i) {
        if (use[i]) {
          acc += design[(size_t)r * row_ct + i] * design[(size_t)c * row_ct + i];
        }
      }
      aa[r * (col_ct + 1) + c] = acc;
    }
    double acc = 0.0;
    for (uint32_t i = 0; i != row_ct; ++i) {
      if (use[i]) { acc += design[(size_t)r * row_ct + i] * yv[i]; }
    }
    aa[r * (col_ct + 1) + col_ct] = acc;
  }
  for (uint32_t p = 0; p != col_ct; ++p) {
    uint32_t piv = p;
    for (uint32_t r = p + 1; r != col_ct; ++r) {
      if (fabs(aa[r * (col_ct + 1) + p]) > fabs(aa[piv * (col_ct + 1) + p])) { piv = r; }
    }
    if (fabs(aa[piv * (col_ct + 1) + p]) < 1e-300) { return 1; }
    if (piv != p) {
      for (uint32_t c = 0; c <= col_ct; ++c) {
        const double t = aa[p * (col_ct + 1) + c];
        aa[p * (col_ct + 1) + c] = aa[piv * (col_ct + 1) + c];
        aa[piv * (col_ct + 1) + c] = t;
      }
    }
    const double d = aa[p * (col_ct + 1) + p];
    for (uint32_t c = 0; c <= col_ct; ++c) { aa[p * (col_ct + 1) + c] /= d; }
    for (uint32_t r = 0; r != col_ct; ++r) {
      if (r == p) { continue; }
      const double f = aa[r * (col_ct + 1) + p];
      if (f == 0.0) { continue; }
      for (uint32_t c = 0; c <= col_ct; ++c) {
        aa[r * (col_ct + 1) + c] -= f * aa[p * (col_ct + 1) + c];
      }
    }
  }
  for (uint32_t r = 0; r != col_ct; ++r) { beta_out[r] = aa[r * (col_ct + 1) + col_ct]; }
  return 0;
}

// Consume the rest of a comment line; returns 1 for the "# data:" marker that
// ends a fixture's key/value block.  Keying the end of the block off this
// marker rather than off whichever key happens to come last is what lets
// reference values be added to a fixture in any order.
static int FixtureCommentIsData(FILE* f) {
  char rest[16];
  if (!fgets(rest, sizeof(rest), f)) { return 0; }
  const int is_data = !strncmp(rest, " data", 5);
  if (!strchr(rest, '\n')) {  // comment longer than the buffer; drain it
    while (!feof(f) && fgetc(f) != '\n') {}
  }
  return is_data;
}

// ---------------------------------------------------------------------------
// Logistic fixture loader (mirrors LoadOlsFixture).
// ---------------------------------------------------------------------------
// Reference values from Tests/pgs_eval_fixtures/logistic_reference.txt.
static const uint32_t kMaxNullBeta = 16;
struct LogisticFixture {
  uint32_t n;
  uint32_t n_covar;
  double beta, se;
  // Coefficients of the covariates-only Logit fit (statsmodels null.params),
  // intercept first.  Repeated `null_beta` lines in the fixture, in order.
  double null_beta[kMaxNullBeta];
  uint32_t null_beta_ct;
  // Observed-scale R² family, from statsmodels OLS on the same 0/1 outcome.
  double r2_null, r2_full, r2_incr;
  // 95% CI on r2_incr: analytic and bootstrap, both from generate_ci.R.  The
  // null_* trio is the same statistic for a PGS with essentially no signal.
  double r2_incr_l95, r2_incr_u95;
  double r2_incr_boot_l95, r2_incr_boot_u95;
  double null_r2_incr, null_r2_incr_l95, null_r2_incr_u95;
  // Liability scale at `prevalence`, again from generate_ci.R.
  double prevalence, r2l_incr_l95, r2l_incr_u95;
  double* y;
  double* covars;  // column-major, n * n_covar
  double* pgs;
};

static int LoadLogisticFixture(const char* path, LogisticFixture* fx) {
  FILE* f = fopen(path, "r");
  if (!f) { printf("  FAIL cannot open %s\n", path); return 1; }
  char key[64];
  double val;
  const double nan_val = 0.0 / 0.0;
  fx->n = 0; fx->n_covar = 0; fx->null_beta_ct = 0;
  fx->y = nullptr; fx->covars = nullptr; fx->pgs = nullptr;
  // NaN, not zero: a key the fixture no longer carries must fail loudly rather
  // than be checked against a plausible-looking 0.
  fx->r2_incr_l95 = nan_val; fx->r2_incr_u95 = nan_val;
  fx->r2_incr_boot_l95 = nan_val; fx->r2_incr_boot_u95 = nan_val;
  fx->null_r2_incr = nan_val;
  fx->null_r2_incr_l95 = nan_val; fx->null_r2_incr_u95 = nan_val;
  fx->prevalence = nan_val;
  fx->r2l_incr_l95 = nan_val; fx->r2l_incr_u95 = nan_val;
  while (fscanf(f, "%63s", key) == 1) {
    if (key[0] == '#') { if (FixtureCommentIsData(f)) { break; } continue; }
    if (fscanf(f, "%lf", &val) != 1) { fclose(f); return 1; }
    if (!strcmp(key, "n"))        { fx->n      = (uint32_t)val; }
    else if (!strcmp(key, "n_covar"))  { fx->n_covar = (uint32_t)val; }
    else if (!strcmp(key, "beta"))     { fx->beta   = val; }
    else if (!strcmp(key, "se"))       { fx->se     = val; }
    else if (!strcmp(key, "null_beta")) {
      // Repeated, one per null-model coefficient, in order.
      if (fx->null_beta_ct == kMaxNullBeta) { fclose(f); return 1; }
      fx->null_beta[fx->null_beta_ct++] = val;
    }
    else if (!strcmp(key, "r2_null"))  { fx->r2_null = val; }
    else if (!strcmp(key, "r2_full"))  { fx->r2_full = val; }
    else if (!strcmp(key, "r2_incr"))  { fx->r2_incr = val; }
    else if (!strcmp(key, "r2_incr_l95")) { fx->r2_incr_l95 = val; }
    else if (!strcmp(key, "r2_incr_u95")) { fx->r2_incr_u95 = val; }
    else if (!strcmp(key, "r2_incr_boot_l95")) { fx->r2_incr_boot_l95 = val; }
    else if (!strcmp(key, "r2_incr_boot_u95")) { fx->r2_incr_boot_u95 = val; }
    else if (!strcmp(key, "null_r2_incr"))     { fx->null_r2_incr = val; }
    else if (!strcmp(key, "null_r2_incr_l95")) { fx->null_r2_incr_l95 = val; }
    else if (!strcmp(key, "null_r2_incr_u95")) { fx->null_r2_incr_u95 = val; }
    else if (!strcmp(key, "prevalence"))       { fx->prevalence = val; }
    else if (!strcmp(key, "r2l_incr_l95"))     { fx->r2l_incr_l95 = val; }
    else if (!strcmp(key, "r2l_incr_u95"))     { fx->r2l_incr_u95 = val; }
  }
  fx->y      = (double*)malloc(fx->n * sizeof(double));
  fx->covars = (double*)malloc((size_t)fx->n * fx->n_covar * sizeof(double));
  fx->pgs    = (double*)malloc(fx->n * sizeof(double));
  if (!fx->y || !fx->covars || !fx->pgs) { fclose(f); return 1; }
  for (uint32_t i = 0; i != fx->n; ++i) {
    if (fscanf(f, "%lf", &fx->y[i]) != 1) { fclose(f); return 1; }
    for (uint32_t j = 0; j != fx->n_covar; ++j) {
      // store column-major
      if (fscanf(f, "%lf", &fx->covars[j * fx->n + i]) != 1) { fclose(f); return 1; }
    }
    if (fscanf(f, "%lf", &fx->pgs[i]) != 1) { fclose(f); return 1; }
  }
  fclose(f);
  return 0;
}

// Reference values and data from Tests/pgs_eval_fixtures/ols_reference.txt.
struct OlsFixture {
  uint32_t n;
  uint32_t n_covar;
  double beta, se, pval, r2_null, r2_full, r2_incr;
  // 95% CI on r2_incr: analytic and bootstrap, both from generate_ci.R.  The
  // null_* trio is the same statistic for a PGS with essentially no signal.
  double r2_incr_l95, r2_incr_u95;
  double r2_incr_boot_l95, r2_incr_boot_u95;
  double null_r2_incr, null_r2_incr_l95, null_r2_incr_u95;
  double* y;
  double* covars;  // column-major, n * n_covar
  double* pgs;
};

static int LoadOlsFixture(const char* path, OlsFixture* fx) {
  FILE* f = fopen(path, "r");
  if (!f) { printf("  FAIL cannot open %s\n", path); return 1; }
  char key[64];
  double val;
  const double nan_val = 0.0 / 0.0;
  fx->n = 0; fx->n_covar = 0;
  // NaN, not zero: see LoadLogisticFixture.
  fx->r2_incr_l95 = nan_val; fx->r2_incr_u95 = nan_val;
  fx->r2_incr_boot_l95 = nan_val; fx->r2_incr_boot_u95 = nan_val;
  fx->null_r2_incr = nan_val;
  fx->null_r2_incr_l95 = nan_val; fx->null_r2_incr_u95 = nan_val;
  while (fscanf(f, "%63s", key) == 1) {
    if (key[0] == '#') { if (FixtureCommentIsData(f)) { break; } continue; }
    if (fscanf(f, "%lf", &val) != 1) { fclose(f); return 1; }
    if (!strcmp(key, "n")) { fx->n = (uint32_t)val; }
    else if (!strcmp(key, "n_covar")) { fx->n_covar = (uint32_t)val; }
    else if (!strcmp(key, "beta")) { fx->beta = val; }
    else if (!strcmp(key, "se")) { fx->se = val; }
    else if (!strcmp(key, "pval")) { fx->pval = val; }
    else if (!strcmp(key, "r2_null")) { fx->r2_null = val; }
    else if (!strcmp(key, "r2_full")) { fx->r2_full = val; }
    else if (!strcmp(key, "r2_incr")) { fx->r2_incr = val; }
    else if (!strcmp(key, "r2_incr_l95")) { fx->r2_incr_l95 = val; }
    else if (!strcmp(key, "r2_incr_u95")) { fx->r2_incr_u95 = val; }
    else if (!strcmp(key, "r2_incr_boot_l95")) { fx->r2_incr_boot_l95 = val; }
    else if (!strcmp(key, "r2_incr_boot_u95")) { fx->r2_incr_boot_u95 = val; }
    else if (!strcmp(key, "null_r2_incr"))     { fx->null_r2_incr = val; }
    else if (!strcmp(key, "null_r2_incr_l95")) { fx->null_r2_incr_l95 = val; }
    else if (!strcmp(key, "null_r2_incr_u95")) { fx->null_r2_incr_u95 = val; }
  }
  fx->y = (double*)malloc(fx->n * sizeof(double));
  fx->covars = (double*)malloc((size_t)fx->n * fx->n_covar * sizeof(double));
  fx->pgs = (double*)malloc(fx->n * sizeof(double));
  for (uint32_t i = 0; i != fx->n; ++i) {
    if (fscanf(f, "%lf", &fx->y[i]) != 1) { fclose(f); return 1; }
    for (uint32_t j = 0; j != fx->n_covar; ++j) {
      // store column-major
      if (fscanf(f, "%lf", &fx->covars[j * fx->n + i]) != 1) { fclose(f); return 1; }
    }
    if (fscanf(f, "%lf", &fx->pgs[i]) != 1) { fclose(f); return 1; }
  }
  fclose(f);
  return 0;
}

int main() {
  printf("pgs_eval_test\n");
  // Two-sided Student-t tail probabilities, checked against
  // scipy.stats.t.sf(x, df) * 2.
  CheckClose("t_pval(2.0, 10)",   plink2::PgsEvalTPval(2.0, 10.0),
             0.07338803477074039, 1e-9);
  CheckClose("t_pval(-2.0, 10)",  plink2::PgsEvalTPval(-2.0, 10.0),
             0.07338803477074039, 1e-9);
  CheckClose("t_pval(0.0, 10)",   plink2::PgsEvalTPval(0.0, 10.0),
             1.0, 1e-12);
  CheckClose("t_pval(4.5, 197)",  plink2::PgsEvalTPval(4.5, 197.0),
             1.1605958615030959e-05, 1e-7);
  // Extreme tail: PGS ranking routinely produces very small p-values; this
  // exercises the deepest continued-fraction iterations.
  // scipy: 2 * stats.t.sf(10.0, 1000.0) == 1.6670702958600065e-22
  CheckClose("t_pval(10.0, 1000)", plink2::PgsEvalTPval(10.0, 1000.0),
             1.6670702958600065e-22, 1e-7);

  // predictor_ct == 1 (intercept only, no covariates).
  // With no covariates, FWL reduces to simple linear regression of e_y on e_g,
  // where e_y = y - y_mean and e_g = g - g_mean (normalised). So r2_null == 0
  // and r2_incr == r2_full == squared Pearson correlation of y and g.
  {
    const uint32_t n = 40;
    double yy[40], gg[40];
    double ys = 0.0, gs = 0.0;
    for (uint32_t i = 0; i != n; ++i) {
      yy[i] = (double)(int32_t)(i % 7) - 3.0;
      gg[i] = (double)(int32_t)(i % 5) - 2.0;
      ys += yy[i]; gs += gg[i];
    }
    const double ym = ys / n, gm = gs / n;
    double sxy = 0.0, sxx = 0.0, syy = 0.0;
    for (uint32_t i = 0; i != n; ++i) {
      sxy += (yy[i] - ym) * (gg[i] - gm);
      sxx += (gg[i] - gm) * (gg[i] - gm);
      syy += (yy[i] - ym) * (yy[i] - ym);
    }
    const double r2_pearson = (sxy * sxy) / (sxx * syy);
    // Design matrix: single intercept column.
    double xx1[40];
    for (uint32_t i = 0; i != n; ++i) { xx1[i] = 1.0; }
    plink2::PgsEvalOlsCtx ctx1;
    if (plink2::PgsEvalOlsInit(yy, xx1, n, 1, &ctx1)) {
      printf("  FAIL predictor_ct==1 PgsEvalOlsInit\n"); ++g_failures;
    } else {
      CheckClose("pc1 r2_null", ctx1.r2_null, 0.0, 1e-15);
      plink2::PgsEvalResult res1;
      plink2::PgsEvalOlsOne(&ctx1, gg, &res1);
      if (!res1.converged) {
        printf("  FAIL predictor_ct==1 did not converge\n"); ++g_failures;
      } else {
        CheckClose("pc1 r2_full == r2_pearson", res1.r2_full, r2_pearson, 1e-12);
        CheckClose("pc1 r2_incr == r2_full",    res1.r2_incr, res1.r2_full, 1e-15);
        CheckClose("pc1 r2_null == 0",          res1.r2_null, 0.0, 1e-15);
        CheckCiContains("pc1 CI contains r2_incr",
                        res1.r2_incr_l95, res1.r2_incr, res1.r2_incr_u95);
      }
      plink2::PgsEvalOlsCleanup(&ctx1);
    }
  }

  // predictor_ct == 2 (intercept + one covariate).
  {
    const uint32_t n = 40;
    double yy[40], cc[40], gg[40];
    for (uint32_t i = 0; i != n; ++i) {
      yy[i] = (double)(int32_t)(i % 7) - 3.0;
      cc[i] = (double)(int32_t)(i % 11) - 5.0;
      gg[i] = (double)(int32_t)(i % 5) - 2.0;
    }
    // Design matrix [1, c], column-major.
    double xx2[80];
    for (uint32_t i = 0; i != n; ++i) { xx2[i] = 1.0; xx2[n + i] = cc[i]; }
    plink2::PgsEvalOlsCtx ctx2;
    if (plink2::PgsEvalOlsInit(yy, xx2, n, 2, &ctx2)) {
      printf("  FAIL predictor_ct==2 PgsEvalOlsInit\n"); ++g_failures;
    } else {
      plink2::PgsEvalResult res2;
      plink2::PgsEvalOlsOne(&ctx2, gg, &res2);
      if (!res2.converged) {
        printf("  FAIL predictor_ct==2 did not converge\n"); ++g_failures;
      } else {
        // r2_full must lie in [r2_null, 1].
        CheckClose("pc2 r2_full >= r2_null",
                   res2.r2_full >= res2.r2_null ? res2.r2_full : res2.r2_null,
                   res2.r2_full, 1e-15);
        CheckClose("pc2 r2_incr == r2_full - r2_null",
                   res2.r2_incr, res2.r2_full - res2.r2_null, 1e-15);
        CheckCiContains("pc2 CI contains r2_incr",
                        res2.r2_incr_l95, res2.r2_incr, res2.r2_incr_u95);
        // Cross-check: full normal equations on [1, c, g_std] must give same r2_full.
        {
          double xs[120];  // n * 3, column-major
          for (uint32_t i = 0; i != n; ++i) {
            xs[i] = 1.0; xs[n + i] = cc[i];
          }
          double gs2 = 0.0;
          for (uint32_t i = 0; i != n; ++i) { gs2 += gg[i]; }
          const double gm2 = gs2 / n;
          double gv2 = 0.0;
          for (uint32_t i = 0; i != n; ++i) { gv2 += (gg[i]-gm2)*(gg[i]-gm2); }
          const double gsd2 = sqrt(gv2 / (n - 1));
          for (uint32_t i = 0; i != n; ++i) { xs[2*n + i] = (gg[i]-gm2)/gsd2; }
          plink2::PgsEvalOlsCtx cf2;
          if (plink2::PgsEvalOlsInit(yy, xs, n, 3, &cf2)) {
            printf("  FAIL predictor_ct==2 cross-check init\n"); ++g_failures;
          } else {
            CheckClose("pc2 crosscheck r2_full", cf2.r2_null, res2.r2_full, 1e-12);
            plink2::PgsEvalOlsCleanup(&cf2);
          }
        }
      }
      plink2::PgsEvalOlsCleanup(&ctx2);
    }
  }

  OlsFixture fx;
  if (LoadOlsFixture("Tests/pgs_eval_fixtures/ols_reference.txt", &fx)) {
    printf("  FAIL fixture load\n");
    ++g_failures;
  } else {
    const uint32_t predictor_ct = fx.n_covar + 1;
    // Design matrix [1, covars], column-major.
    double* xx = (double*)malloc((size_t)fx.n * predictor_ct * sizeof(double));
    for (uint32_t i = 0; i != fx.n; ++i) { xx[i] = 1.0; }
    memcpy(&(xx[fx.n]), fx.covars,
           (size_t)fx.n * fx.n_covar * sizeof(double));

    plink2::PgsEvalOlsCtx ctx;
    if (plink2::PgsEvalOlsInit(fx.y, xx, fx.n, predictor_ct, &ctx)) {
      printf("  FAIL PgsEvalOlsInit\n");
      ++g_failures;
    } else {
      CheckClose("r2_null", ctx.r2_null, fx.r2_null, 1e-10);
      plink2::PgsEvalResult res;
      plink2::PgsEvalOlsOne(&ctx, fx.pgs, &res);
      CheckClose("beta",    res.beta,    fx.beta,    1e-9);
      CheckClose("se",      res.se,      fx.se,      1e-9);
      CheckClose("pval",    res.pval,    fx.pval,    1e-7);
      CheckClose("r2_full", res.r2_full, fx.r2_full, 1e-10);
      CheckClose("r2_incr", res.r2_incr, fx.r2_incr, 1e-10);
      // Analytic 95% CI, against R's independent implementation of the same
      // formula: its residuals come from a Householder QR, not from the
      // Frisch-Waugh-Lovell path under test here.
      CheckClose("r2_incr_l95", res.r2_incr_l95, fx.r2_incr_l95, 1e-9);
      CheckClose("r2_incr_u95", res.r2_incr_u95, fx.r2_incr_u95, 1e-9);
      CheckCiContains("CI contains r2_incr",
                      res.r2_incr_l95, res.r2_incr, res.r2_incr_u95);
      // Cross-check against a 2000-replicate percentile bootstrap.  Overlap,
      // deliberately not equality: the bootstrap validates the method, and
      // asserting agreement to a tight tolerance against a resampling estimate
      // would only be asserting that the seed has not changed.
      {
        const double frac = CiOverlapFrac(res.r2_incr_l95, res.r2_incr_u95,
                                          fx.r2_incr_boot_l95, fx.r2_incr_boot_u95);
        ++g_checks;
        if (!(frac > 0.8)) {
          printf("  FAIL analytic/bootstrap CI overlap %.3f: analytic [%.6g, %.6g] "
                 "bootstrap [%.6g, %.6g]\n", frac, res.r2_incr_l95, res.r2_incr_u95,
                 fx.r2_incr_boot_l95, fx.r2_incr_boot_u95);
          ++g_failures;
        }
      }
      // Straddle-zero: the same y and covariates against a PGS with essentially
      // no relation to them, so the interval on the partial correlation
      // contains zero.  Squaring that interval must give a lower bound of
      // exactly 0 -- take lo_r^2 instead and the interval excludes its own
      // point estimate, in precisely the near-null case a CI is consulted for.
      {
        double* g_null = (double*)malloc(fx.n * sizeof(double));
        for (uint32_t i = 0; i != fx.n; ++i) { g_null[i] = NullPgs(i); }
        plink2::PgsEvalResult nres;
        plink2::PgsEvalOlsOne(&ctx, g_null, &nres);
        ++g_checks;
        if (!nres.converged) {
          printf("  FAIL straddle-zero PGS did not converge\n"); ++g_failures;
        } else {
          CheckClose("straddle r2_incr",     nres.r2_incr,     fx.null_r2_incr,     1e-9);
          CheckClose("straddle r2_incr_u95", nres.r2_incr_u95, fx.null_r2_incr_u95, 1e-9);
          ++g_checks;
          if (nres.r2_incr_l95 != 0.0) {
            printf("  FAIL straddle r2_incr_l95 must be exactly 0, got %.17g\n",
                   nres.r2_incr_l95);
            ++g_failures;
          }
          CheckCiContains("straddle CI contains r2_incr",
                          nres.r2_incr_l95, nres.r2_incr, nres.r2_incr_u95);
          // The point estimate really is near null, so the case tests what its
          // name says: a positive dR2 whose interval still reaches zero.
          ++g_checks;
          if (!(nres.r2_incr > 0.0 && nres.r2_incr < 0.01)) {
            printf("  FAIL straddle case is not near-null: r2_incr %.6g\n",
                   nres.r2_incr);
            ++g_failures;
          }
        }
        free(g_null);
      }
      // Cross-check: same fit via explicit normal equations on the full
      // design [1, covars, pgs], independent of the FWL path above.
      {
        const uint32_t pc_full = predictor_ct + 1;
        double* xf = (double*)malloc((size_t)fx.n * pc_full * sizeof(double));
        memcpy(xf, xx, (size_t)fx.n * predictor_ct * sizeof(double));
        double gs = 0.0;
        for (uint32_t i = 0; i != fx.n; ++i) { gs += fx.pgs[i]; }
        const double gm = gs / fx.n;
        double gv = 0.0;
        for (uint32_t i = 0; i != fx.n; ++i) { gv += (fx.pgs[i]-gm)*(fx.pgs[i]-gm); }
        const double gsd2 = sqrt(gv / (fx.n - 1));
        for (uint32_t i = 0; i != fx.n; ++i) {
          xf[(size_t)predictor_ct * fx.n + i] = (fx.pgs[i] - gm) / gsd2;
        }
        plink2::PgsEvalOlsCtx cf;
        if (plink2::PgsEvalOlsInit(fx.y, xf, fx.n, pc_full, &cf)) {
          printf("  FAIL cross-check init\n"); ++g_failures;
        } else {
          // R2 of the full model must match the FWL r2_full exactly.
          CheckClose("crosscheck r2_full", cf.r2_null, res.r2_full, 1e-12);
          plink2::PgsEvalOlsCleanup(&cf);
        }
        free(xf);
      }
      // -----------------------------------------------------------------------
      // The cross-validated interval reaches the same Fisher z core as the
      // in-sample one, but derives its partial correlation from the reported
      // R^2 instead of from the residual dot products:
      //     r = sqrt(r2_incr / (1 - r2_null))
      // In-sample that is an identity, not an approximation -- rearranged it is
      // r^2 * (1 - R2_null) = r2_incr, the line the containment guarantee rests
      // on.  So at fold_ct == 1 the two routes must return the same interval,
      // which is the whole claim that there is one formula and not two.  It is
      // checked here rather than asserted in a comment.
      //
      // On this fixture the two routes agree bit for bit -- measured, both
      // assertions pass at a tolerance of exactly 0.  The tolerance is not 0
      // because that agreement is luck, not structure: the CV route reaches r
      // through r2_incr, and r2_incr is r2_full - r2_null, which cancels
      // catastrophically once the effect is small against strong covariates.
      // The disagreement then grows without bound as r2_incr / (1 - r2_null)
      // shrinks -- 4.4e-10 relative over a 20,000-dataset sweep across eight
      // decades of effect size, and 1.8e-7 at the worst configuration measured
      // (r2_null 0.99992, r2_incr 3.5e-13).  1e-6 clears all of that and is
      // still nowhere near the O(1) gap a genuinely different formula would
      // open, which is the only thing this assertion can catch.
      // -----------------------------------------------------------------------
      {
        double lo = 0.0, hi = 0.0;
        plink2::IncrR2CiFromRTest(sqrt(res.r2_incr / (1.0 - res.r2_null)),
                                  res.r2_null, fx.n, predictor_ct, &lo, &hi);
        CheckClose("CV CI formula == in-sample L95", lo, res.r2_incr_l95, 1e-6);
        CheckClose("CV CI formula == in-sample U95", hi, res.r2_incr_u95, 1e-6);
      }
      plink2::PgsEvalOlsCleanup(&ctx);
    }
    free(xx);
  }

  // ---------------------------------------------------------------------------
  // k-fold cross-validation.
  //
  // The one piece of new algebra is that a held-out prediction of the full
  // model can be assembled from the null model's held-out prediction plus
  // beta_g times the PGS residualised on the covariates *of the training set*,
  // evaluated on the held-out rows.  That is Frisch-Waugh applied inside each
  // training set, and it is what lets cross-validation cost k small
  // factorisations rather than k per model.  Here it is checked against
  // explicit refits of both models on each training set, through a solver that
  // shares no code with it.
  // ---------------------------------------------------------------------------
  {
    const uint32_t n = 60;
    const uint32_t predictor_ct = 2;  // intercept + one covariate
    const uint32_t fold_ct = 3;
    double yy[60], cc[60], gg[60];
    for (uint32_t i = 0; i != n; ++i) {
      cc[i] = (double)(int32_t)((i * 13u) % 17u) - 8.0;
      gg[i] = (double)(int32_t)((i * 29u) % 23u) - 11.0;
      // Real signal in g, so r2_incr comes out clearly positive out of fold.
      yy[i] = 0.4 * cc[i] + 0.6 * gg[i] + (double)(int32_t)((i * 7u) % 11u) - 5.0;
    }
    double xx[120];
    for (uint32_t i = 0; i != n; ++i) { xx[i] = 1.0; xx[n + i] = cc[i]; }

    plink2::PgsEvalOlsCtx ctx;
    if (plink2::PgsEvalOlsInit(yy, xx, n, predictor_ct, &ctx)) {
      printf("  FAIL cv PgsEvalOlsInit\n"); ++g_failures;
    } else if (plink2::PgsEvalOlsSetFolds(&ctx, fold_ct, 0, 20260802u)) {
      printf("  FAIL cv PgsEvalOlsSetFolds\n"); ++g_failures;
      plink2::PgsEvalOlsCleanup(&ctx);
    } else {
      plink2::PgsEvalResult res;
      plink2::PgsEvalOlsOne(&ctx, gg, &res);
      ++g_checks;
      if (!res.converged) {
        printf("  FAIL cv fit did not converge\n"); ++g_failures;
      } else {
        // Every fold gets n / fold_ct or that plus one, and every sample lands
        // in exactly one -- the property the pooled set being all of y depends
        // on.
        uint32_t fold_size[8] = {0};
        for (uint32_t i = 0; i != n; ++i) { ++fold_size[ctx.fold_assign[i]]; }
        for (uint32_t f = 0; f != fold_ct; ++f) {
          ++g_checks;
          if (fold_size[f] != n / fold_ct) {
            printf("  FAIL cv fold %u has %u samples, want %u\n",
                   f, fold_size[f], n / fold_ct);
            ++g_failures;
          }
        }

        // Independent out-of-fold prediction: refit both models outright on
        // each training set and predict the held-out rows.
        const double* g_used = gg;
        double gsum = 0.0;
        for (uint32_t i = 0; i != n; ++i) { gsum += g_used[i]; }
        const double gm = gsum / n;
        double gv = 0.0;
        for (uint32_t i = 0; i != n; ++i) { gv += (g_used[i]-gm)*(g_used[i]-gm); }
        const double gsd = sqrt(gv / (n - 1));
        double xg[180];  // [1, c, g_std], column-major
        for (uint32_t i = 0; i != n; ++i) {
          xg[i] = 1.0; xg[n + i] = cc[i]; xg[2*n + i] = (g_used[i] - gm) / gsd;
        }
        unsigned char use[60];
        double bnull[8], bfull[8];
        double ss_null = 0.0, ss_full = 0.0;
        int solve_fail = 0;
        for (uint32_t f = 0; f != fold_ct && !solve_fail; ++f) {
          for (uint32_t i = 0; i != n; ++i) {
            use[i] = (unsigned char)(ctx.fold_assign[i] != f);
          }
          solve_fail |= SolveSubsetOls(xx, yy, n, predictor_ct, use, bnull);
          solve_fail |= SolveSubsetOls(xg, yy, n, predictor_ct + 1, use, bfull);
          if (solve_fail) { break; }
          for (uint32_t i = 0; i != n; ++i) {
            if (use[i]) { continue; }
            double pn = 0.0, pf = 0.0;
            for (uint32_t k = 0; k != predictor_ct; ++k) {
              pn += xx[(size_t)k * n + i] * bnull[k];
            }
            for (uint32_t k = 0; k != predictor_ct + 1; ++k) {
              pf += xg[(size_t)k * n + i] * bfull[k];
            }
            ss_null += (yy[i] - pn) * (yy[i] - pn);
            ss_full += (yy[i] - pf) * (yy[i] - pf);
          }
        }
        ++g_checks;
        if (solve_fail) {
          printf("  FAIL cv independent refit is singular\n"); ++g_failures;
        } else {
          CheckClose("cv r2_null == refit out-of-fold", res.r2_null,
                     1.0 - ss_null / ctx.tss, 1e-12);
          CheckClose("cv r2_full == refit out-of-fold", res.r2_full,
                     1.0 - ss_full / ctx.tss, 1e-12);
        }
        // The identity pgs_eval_check.sh Case 21b asserts on the printed file.
        CheckClose("cv r2_incr == r2_full - r2_null", res.r2_incr,
                   res.r2_full - res.r2_null, 1e-15);
        CheckCiContains("cv CI contains r2_incr",
                        res.r2_incr_l95, res.r2_incr, res.r2_incr_u95);
        // beta is the whole-sample estimate, unaffected by the folds.
        {
          plink2::PgsEvalOlsCtx plain;
          if (plink2::PgsEvalOlsInit(yy, xx, n, predictor_ct, &plain)) {
            printf("  FAIL cv beta cross-check init\n"); ++g_failures;
          } else {
            plink2::PgsEvalResult pres;
            plink2::PgsEvalOlsOne(&plain, gg, &pres);
            CheckClose("cv beta stays in-sample", res.beta, pres.beta, 0.0);
            CheckClose("cv se stays in-sample",   res.se,   pres.se,   0.0);
            CheckClose("cv pval stays in-sample", res.pval, pres.pval, 0.0);
            // Both fits explain less of the held-out data than of the data they
            // were fitted on.  That is the whole point, and it is the assertion
            // that holds per-model: r2_incr is a difference of two quantities
            // that both shrink, so it is not itself monotone in k -- the
            // shrinkage of the *winner of a ranking* is asserted end to end in
            // pgs_eval_check.sh, on a sample large enough for it to mean
            // something.
            ++g_checks;
            if (!(res.r2_full < pres.r2_full && res.r2_null < pres.r2_null)) {
              printf("  FAIL cv out-of-fold R2 not below in-sample: null %.6g/%.6g full %.6g/%.6g\n",
                     res.r2_null, pres.r2_null, res.r2_full, pres.r2_full);
              ++g_failures;
            }
            ++g_checks;
            if (!(res.r2_incr > 0.0)) {
              printf("  FAIL cv r2_incr %.6g should stay positive for a real effect\n",
                     res.r2_incr);
              ++g_failures;
            }
            plink2::PgsEvalOlsCleanup(&plain);
          }
        }
        // A PGS with no relation to y: out of fold it predicts worse than
        // nothing, r2_incr goes negative, and there is no partial correlation
        // for the interval to be about -- so both bounds are NA rather than an
        // interval that excludes its own point estimate.
        {
          double g_noise[60];
          for (uint32_t i = 0; i != n; ++i) { g_noise[i] = NullPgs(i); }
          plink2::PgsEvalResult nres;
          plink2::PgsEvalOlsOne(&ctx, g_noise, &nres);
          ++g_checks;
          if (!nres.converged) {
            printf("  FAIL cv noise PGS did not converge\n"); ++g_failures;
          } else {
            ++g_checks;
            if (!(nres.r2_incr < 0.0)) {
              printf("  FAIL cv noise r2_incr %.6g should be negative out of fold\n",
                     nres.r2_incr);
              ++g_failures;
            }
            ++g_checks;
            if (!(isnan(nres.r2_incr_l95) && isnan(nres.r2_incr_u95))) {
              printf("  FAIL cv negative r2_incr must have NA bounds, got [%.6g, %.6g]\n",
                     nres.r2_incr_l95, nres.r2_incr_u95);
              ++g_failures;
            }
          }
        }
      }
      plink2::PgsEvalOlsCleanup(&ctx);
    }
  }

  // ---------------------------------------------------------------------------
  // Fold sizes when n is NOT divisible by k -- the case the block above cannot
  // reach with 60 samples over 3 folds.  Shuffling and then dealing round-robin
  // is meant to give every fold floor(n/k) or ceil(n/k) samples, with every
  // sample in exactly one fold; an off-by-one in either half shows up as a fold
  // one too large, an unassigned sample, or a fold index out of range, none of
  // which an exactly-divisible n can see.
  //
  // Membership in {floor, ceil} plus a total of n pins the split exactly: the
  // two sizes differ by one, so if `a` folds are large then a*(lo+1) +
  // (k-a)*lo == n forces a == n % k.  No separate count assertion is needed.
  //
  // Every pair here leaves the smallest training set above predictor_ct + 1,
  // which is what PgsEvalOlsSetFolds requires; {7, 2} does not -- holding out
  // its 4-sample fold leaves 3 training samples for 2 predictors -- which is
  // why the smallest pair is {9, 2}.
  // ---------------------------------------------------------------------------
  {
    static const uint32_t kFoldSizeCases[][2] = {
      {9, 2}, {10, 3}, {23, 5}, {31, 7}, {50, 8}, {59, 10}
    };
    const uint32_t predictor_ct = 2;
    for (uint32_t t = 0; t != sizeof(kFoldSizeCases) / sizeof(kFoldSizeCases[0]); ++t) {
      const uint32_t n = kFoldSizeCases[t][0];
      const uint32_t fold_ct = kFoldSizeCases[t][1];
      double* yy = (double*)malloc(n * sizeof(double));
      double* xx = (double*)malloc((size_t)n * predictor_ct * sizeof(double));
      for (uint32_t i = 0; i != n; ++i) {
        xx[i] = 1.0;
        xx[n + i] = (double)(int32_t)((i * 19u) % 31u) - 15.0;
        yy[i] = 0.5 * xx[n + i] + (double)(int32_t)((i * 7u) % 11u) - 5.0;
      }
      plink2::PgsEvalOlsCtx ctx;
      if (plink2::PgsEvalOlsInit(yy, xx, n, predictor_ct, &ctx)) {
        printf("  FAIL fold-size n=%u k=%u init\n", n, fold_ct); ++g_failures;
      } else if (plink2::PgsEvalOlsSetFolds(&ctx, fold_ct, 0, 20260802u + t)) {
        printf("  FAIL fold-size n=%u k=%u SetFolds\n", n, fold_ct); ++g_failures;
        plink2::PgsEvalOlsCleanup(&ctx);
      } else {
        uint32_t fold_size[16] = {0};
        uint32_t out_of_range = 0;
        for (uint32_t i = 0; i != n; ++i) {
          if (ctx.fold_assign[i] >= fold_ct) { ++out_of_range; continue; }
          ++fold_size[ctx.fold_assign[i]];
        }
        const uint32_t lo_size = n / fold_ct;
        const uint32_t hi_size = (n + fold_ct - 1) / fold_ct;
        uint32_t total = 0;
        uint32_t bad = 0;
        for (uint32_t f = 0; f != fold_ct; ++f) {
          total += fold_size[f];
          if (fold_size[f] != lo_size && fold_size[f] != hi_size) { ++bad; }
        }
        ++g_checks;
        if (out_of_range || bad || total != n) {
          printf("  FAIL fold sizes n=%u k=%u: %u out of range, %u not in {%u, %u}, "
                 "%u assigned of %u\n", n, fold_ct, out_of_range, bad, lo_size,
                 hi_size, total, n);
          ++g_failures;
        }
        plink2::PgsEvalOlsCleanup(&ctx);
      }
      free(yy);
      free(xx);
    }
  }

  // ---------------------------------------------------------------------------
  // Containment in the near-perfect-fit corner, at full double precision.
  //
  // IncrR2CiFromR clamps |r| to 0.9999999 before the Fisher z transform, which
  // caps the upper bound at
  //   tanh(atanh(0.9999999) + 1.96 / sqrt(n - p - 2))^2 * (1 - R2_NULL).
  // Once r goes past that the bound stops tracking the point estimate and ends
  // up *below* it -- measured at 9.2e-8 relative here, on the in-sample and the
  // cross-validated path alike.  A score column equal to the phenotype plus
  // 1e-5 noise gets there, which is what handing the tool the phenotype itself
  // as a model by mistake looks like.  ClampCiToPoint is what closes it.
  //
  // This has to be asserted here and not in pgs_eval_check.sh: dtoa_g prints
  // about six significant digits, so a 1e-7 breach is invisible in the file.
  // ---------------------------------------------------------------------------
  {
    const uint32_t n = 30;
    const uint32_t predictor_ct = 2;
    double yy[30], xx[60], gg[30];
    for (uint32_t i = 0; i != n; ++i) {
      xx[i] = 1.0;
      xx[n + i] = (double)(int32_t)((i * 13u) % 17u) - 8.0;
      yy[i] = 0.4 * xx[n + i] + (double)(int32_t)((i * 7u) % 11u) - 5.0;
      // Almost exactly y, so the full model's residual sum of squares lands at
      // 4.0e-10 of the null's -- well past the clamp's reach, which starts
      // biting below about 9.3e-8 at these n and p.
      gg[i] = yy[i] + 1e-5 * ((double)(int32_t)((i * 29u) % 23u) - 11.0);
    }
    static const uint32_t kNearPerfectFoldCts[] = {1, 5};
    for (uint32_t t = 0; t != sizeof(kNearPerfectFoldCts) / sizeof(kNearPerfectFoldCts[0]); ++t) {
      const uint32_t fold_ct = kNearPerfectFoldCts[t];
      plink2::PgsEvalOlsCtx ctx;
      if (plink2::PgsEvalOlsInit(yy, xx, n, predictor_ct, &ctx)) {
        printf("  FAIL near-perfect k=%u init\n", fold_ct); ++g_failures;
        continue;
      }
      if (plink2::PgsEvalOlsSetFolds(&ctx, fold_ct, 0, 20260802u)) {
        printf("  FAIL near-perfect k=%u SetFolds\n", fold_ct); ++g_failures;
        plink2::PgsEvalOlsCleanup(&ctx);
        continue;
      }
      plink2::PgsEvalResult res;
      plink2::PgsEvalOlsOne(&ctx, gg, &res);
      ++g_checks;
      if (!res.converged) {
        printf("  FAIL near-perfect k=%u did not converge\n", fold_ct);
        ++g_failures;
      } else {
        // The case only tests what its name says if the fit really is in the
        // clamped corner, so that is asserted rather than assumed: the raw
        // Fisher z bound, the one ClampCiToPoint corrects, must come out below
        // the point estimate.  A weaker gg would pass containment for the
        // boring reason and this case would be testing nothing.
        double raw_lo = 0.0, raw_hi = 0.0;
        plink2::IncrR2CiFromRTest(sqrt(res.r2_incr / (1.0 - res.r2_null)),
                                  res.r2_null, n, predictor_ct,
                                  &raw_lo, &raw_hi);
        ++g_checks;
        if (!(raw_hi < res.r2_incr)) {
          printf("  FAIL near-perfect k=%u is not in the clamped corner: "
                 "unclamped U95 %.17g, r2_incr %.17g\n",
                 fold_ct, raw_hi, res.r2_incr);
          ++g_failures;
        }
        char what[64];
        snprintf(what, sizeof(what), "near-perfect k=%u CI contains r2_incr",
                 fold_ct);
        CheckCiContains(what, res.r2_incr_l95, res.r2_incr, res.r2_incr_u95);
      }
      plink2::PgsEvalOlsCleanup(&ctx);
    }
  }

  // Stratified folds: an imbalanced binary phenotype must still put cases in
  // every fold.  30 cases over 10 folds is three per fold when stratified;
  // shuffling the 300 samples without regard to y leaves an empty fold with
  // probability around 4%, and an empty fold is a fit with nothing to score.
  {
    const uint32_t n = 300;
    const uint32_t predictor_ct = 2;
    const uint32_t fold_ct = 10;
    double* yy = (double*)malloc(n * sizeof(double));
    double* xx = (double*)malloc((size_t)n * predictor_ct * sizeof(double));
    for (uint32_t i = 0; i != n; ++i) {
      yy[i] = ((i % 10u) == 3u) ? 1.0 : 0.0;  // 30 cases
      xx[i] = 1.0;
      xx[n + i] = (double)(int32_t)((i * 19u) % 31u) - 15.0;
    }
    plink2::PgsEvalOlsCtx ctx;
    if (plink2::PgsEvalOlsInit(yy, xx, n, predictor_ct, &ctx)) {
      printf("  FAIL stratified cv init\n"); ++g_failures;
    } else if (plink2::PgsEvalOlsSetFolds(&ctx, fold_ct, 1, 20260802u)) {
      printf("  FAIL stratified cv SetFolds\n"); ++g_failures;
      plink2::PgsEvalOlsCleanup(&ctx);
    } else {
      uint32_t case_in_fold[16] = {0};
      uint32_t ctrl_in_fold[16] = {0};
      for (uint32_t i = 0; i != n; ++i) {
        if (yy[i] > 0.5) { ++case_in_fold[ctx.fold_assign[i]]; }
        else { ++ctrl_in_fold[ctx.fold_assign[i]]; }
      }
      for (uint32_t f = 0; f != fold_ct; ++f) {
        ++g_checks;
        if (case_in_fold[f] != 3 || ctrl_in_fold[f] != 27) {
          printf("  FAIL stratified fold %u has %u cases / %u controls, want 3 / 27\n",
                 f, case_in_fold[f], ctrl_in_fold[f]);
          ++g_failures;
        }
      }
      plink2::PgsEvalOlsCleanup(&ctx);
    }
    free(yy);
    free(xx);
  }

  // ---------------------------------------------------------------------------
  // Binary path: assert against statsmodels at 1e-6.
  //
  // Two fits are checked through the single PgsEvalLogisticOne call, because
  // that is what the function now does: Logit supplies beta/se, and OLS on the
  // same 0/1 outcome supplies the observed-scale R² family that R2_NULL /
  // R2_FULL / R2_INCR report.
  //
  // Tolerance is 1e-6, not 1e-9 as in the OLS path, because IRLS converges
  // iteratively (stopping when max|delta| < 1e-9 per coefficient), while the
  // OLS path is closed-form.  Do not tighten this to 1e-9: IRLS cannot
  // guarantee coefficient-level accuracy tighter than the iteration stop
  // criterion, and floating-point accumulation across 300 samples adds further
  // spread; 1e-6 is the correct bound and a mystery failure here almost
  // certainly means someone tightened it without understanding that.  The R²
  // assertions are held at 1e-10 -- they come from the closed-form OLS path,
  // and nothing about the logistic fit feeds into them.
  // ---------------------------------------------------------------------------
  LogisticFixture lfx;
  if (LoadLogisticFixture("Tests/pgs_eval_fixtures/logistic_reference.txt", &lfx)) {
    printf("  FAIL logistic fixture load\n");
    ++g_failures;
  } else {
    const uint32_t predictor_ct = lfx.n_covar + 1;  // intercept + covariates

    // Design matrix [1, covars], column-major.
    double* lxx = (double*)malloc((size_t)lfx.n * predictor_ct * sizeof(double));
    if (!lxx) {
      printf("  FAIL logistic test alloc\n"); ++g_failures;
    } else {
      for (uint32_t i = 0; i != lfx.n; ++i) { lxx[i] = 1.0; }
      memcpy(&(lxx[lfx.n]), lfx.covars, (size_t)lfx.n * lfx.n_covar * sizeof(double));

      plink2::PgsEvalLogisticCtx lctx;
      if (plink2::PgsEvalLogisticInit(lfx.y, lxx, lfx.n, predictor_ct, &lctx)) {
        printf("  FAIL logistic PgsEvalLogisticInit\n"); ++g_failures;
      } else {
        // The null model's own coefficients, which PgsEvalLogisticInit leaves in
        // lctx.beta and PgsEvalLogisticOne never touches.  Every coefficient is
        // pinned to statsmodels' null.params -- the null fit is otherwise only
        // visible through auc_null, and a rank statistic cannot pin down a fit.
        ++g_checks;
        if (lfx.null_beta_ct != predictor_ct) {
          printf("  FAIL logistic null_beta count: fixture has %u, need %u\n",
                 lfx.null_beta_ct, predictor_ct);
          ++g_failures;
        } else {
          for (uint32_t k = 0; k != predictor_ct; ++k) {
            char what[48];
            snprintf(what, sizeof(what), "logistic null beta[%u]", k);
            CheckClose(what, lctx.beta[k], lfx.null_beta[k], 1e-6);
          }
        }

        plink2::PgsEvalResult lres;
        plink2::PgsEvalLogisticOne(&lctx, lfx.pgs, -1.0, &lres);

        if (!lres.converged) {
          printf("  FAIL logistic did not converge\n"); ++g_failures;
        } else {
          CheckClose("logistic beta",     lres.beta,     lfx.beta,     1e-6);
          CheckClose("logistic se",       lres.se,       lfx.se,       1e-6);
          // Observed-scale R²: these must equal statsmodels' OLS R² on the
          // 0/1 outcome, not any pseudo-R² of the logistic fit.  A Nagelkerke
          // on this fixture gives 0.0582 / 0.1627 / 0.1045 against the observed
          // 0.0424 / 0.1200 / 0.0776 -- roughly 35% larger, and it would fail
          // loudly here.
          CheckClose("binary r2_null", lres.r2_null, lfx.r2_null, 1e-10);
          CheckClose("binary r2_full", lres.r2_full, lfx.r2_full, 1e-10);
          CheckClose("binary r2_incr", lres.r2_incr, lfx.r2_incr, 1e-10);
          // The CI reaches the binary path with no code of its own:
          // PgsEvalLogisticOne passes the caller's struct straight into
          // PgsEvalOlsOne, so a field added there arrives here.  These
          // assertions are what proves that hand-off is still structural.
          CheckClose("binary r2_incr_l95", lres.r2_incr_l95, lfx.r2_incr_l95, 1e-9);
          CheckClose("binary r2_incr_u95", lres.r2_incr_u95, lfx.r2_incr_u95, 1e-9);
          CheckCiContains("binary CI contains r2_incr",
                          lres.r2_incr_l95, lres.r2_incr, lres.r2_incr_u95);
          {
            const double frac = CiOverlapFrac(lres.r2_incr_l95, lres.r2_incr_u95,
                                              lfx.r2_incr_boot_l95,
                                              lfx.r2_incr_boot_u95);
            ++g_checks;
            if (!(frac > 0.8)) {
              printf("  FAIL binary analytic/bootstrap CI overlap %.3f: analytic "
                     "[%.6g, %.6g] bootstrap [%.6g, %.6g]\n", frac,
                     lres.r2_incr_l95, lres.r2_incr_u95,
                     lfx.r2_incr_boot_l95, lfx.r2_incr_boot_u95);
              ++g_failures;
            }
          }
          // Without --prevalence the whole liability family is absent.
          ++g_checks;
          if (!isnan(lres.r2_liability) || !isnan(lres.r2_liability_l95) ||
              !isnan(lres.r2_liability_u95)) {
            printf("  FAIL liability family must be NaN without prevalence\n");
            ++g_failures;
          }
        }

        // Liability scale: the same interval times the Lee constant, checked
        // against R rather than against LiabilityScaleR2 applied to our own
        // endpoints, which would only be restating the C code.
        plink2::PgsEvalResult pres;
        plink2::PgsEvalLogisticOne(&lctx, lfx.pgs, lfx.prevalence, &pres);
        if (!pres.converged) {
          printf("  FAIL logistic with prevalence did not converge\n"); ++g_failures;
        } else {
          CheckClose("binary r2l_incr_l95", pres.r2_liability_l95,
                     lfx.r2l_incr_l95, 1e-9);
          CheckClose("binary r2l_incr_u95", pres.r2_liability_u95,
                     lfx.r2l_incr_u95, 1e-9);
          CheckCiContains("binary liability CI contains r2l_incr",
                          pres.r2_liability_l95, pres.r2_liability,
                          pres.r2_liability_u95);
          CheckCiContains("binary CI contains r2_incr (prevalence)",
                          pres.r2_incr_l95, pres.r2_incr, pres.r2_incr_u95);
        }

        // Straddle-zero on the binary path: same near-null PGS, reaching the
        // OLS companion through PgsEvalLogisticOne.
        {
          double* g_null = (double*)malloc(lfx.n * sizeof(double));
          for (uint32_t i = 0; i != lfx.n; ++i) { g_null[i] = NullPgs(i); }
          plink2::PgsEvalResult nres;
          plink2::PgsEvalLogisticOne(&lctx, g_null, lfx.prevalence, &nres);
          ++g_checks;
          if (!nres.converged) {
            printf("  FAIL binary straddle-zero PGS did not converge\n"); ++g_failures;
          } else {
            CheckClose("binary straddle r2_incr", nres.r2_incr,
                       lfx.null_r2_incr, 1e-9);
            CheckClose("binary straddle r2_incr_u95", nres.r2_incr_u95,
                       lfx.null_r2_incr_u95, 1e-9);
            ++g_checks;
            if (nres.r2_incr_l95 != 0.0) {
              printf("  FAIL binary straddle r2_incr_l95 must be exactly 0, got %.17g\n",
                     nres.r2_incr_l95);
              ++g_failures;
            }
            // The Lee factor is positive, so a zero lower bound stays zero.
            ++g_checks;
            if (nres.r2_liability_l95 != 0.0) {
              printf("  FAIL binary straddle r2l_incr_l95 must be exactly 0, got %.17g\n",
                     nres.r2_liability_l95);
              ++g_failures;
            }
            CheckCiContains("binary straddle CI contains r2_incr",
                            nres.r2_incr_l95, nres.r2_incr, nres.r2_incr_u95);
            CheckCiContains("binary straddle liability CI contains r2l",
                            nres.r2_liability_l95, nres.r2_liability,
                            nres.r2_liability_u95);
          }
          free(g_null);
        }
        plink2::PgsEvalLogisticCleanup(&lctx);
      }
      free(lxx);
    }
    free(lfx.y); free(lfx.covars); free(lfx.pgs);
  }

  // ---------------------------------------------------------------------------
  // Separation test: complete separation, 8 rows.
  //
  // y = 0,0,0,0,1,1,1,1  x = 1,2,3,4,5,6,7,8
  // y is 0 for x <= 4 and 1 for x >= 5, so the data is completely separated.
  // Plain logistic IRLS diverges: R's glm(family=binomial) reports
  // beta ~ 45.8 with a "fitted probabilities numerically 0 or 1" warning.
  // Reporting converged = 0 (NA) is the more honest outcome.
  //
  // This test proves plain GLM *detects* separation and does not return a
  // garbage coefficient.  (Firth, which would return a finite value, has been
  // removed from this path; it is recoverable from git eeae1e2.)
  //
  // The OLS companion fit converges here (y on [1, x] has R² = 0.762, nowhere
  // near degenerate), so this also pins down the converged semantics: a usable
  // R² does not rescue a row whose logistic fit failed.
  //
  // The design matrix [intercept, x] with predictor_ct = 2 also exercises the
  // MAXV(predictor_ct, 3) sizing for dbl_2d_buf when predictor_ct = 2 < 3,
  // which is the live heap-overflow case from the OLS path review.
  // ---------------------------------------------------------------------------
  {
    const double sep_y[8] = {0,0,0,0,1,1,1,1};
    const double sep_x1[8] = {1,1,1,1,1,1,1,1};  // intercept only
    const double sep_pgs[8] = {1,2,3,4,5,6,7,8};
    const uint32_t sep_n = 8;

    plink2::PgsEvalLogisticCtx nctx;
    if (plink2::PgsEvalLogisticInit(sep_y, sep_x1, sep_n, 1, &nctx)) {
      // Init may fail on the null model itself in rare cases.
      ++g_checks;  // count as a pass: non-convergence in Init is acceptable
    } else {
      // Verify the claim above: the OLS companion really does converge on this
      // data.  Without this the separation assertion could start passing
      // because OLS bailed rather than because the logistic fit did, and the
      // test would no longer be testing what its name says.  nctx.ols is the
      // very context PgsEvalLogisticOne uses, so there is nothing to keep in
      // step here.
      plink2::PgsEvalResult ores;
      plink2::PgsEvalOlsOne(&nctx.ols, sep_pgs, &ores);
      ++g_checks;
      if (!ores.converged) {
        printf("  FAIL separation: OLS on the 0/1 outcome should converge\n");
        ++g_failures;
      }

      plink2::PgsEvalResult nres;
      plink2::PgsEvalLogisticOne(&nctx, sep_pgs, -1.0, &nres);
      ++g_checks;
      if (nres.converged) {
        printf("  FAIL separation: plain logistic should NOT converge on "
               "completely separated data\n");
        ++g_failures;
      }
      // Correct: converged == 0.
      plink2::PgsEvalLogisticCleanup(&nctx);
    }
  }

  // ---------------------------------------------------------------------------
  // AUC: two closed-form cases, no external reference needed.
  //
  // 1. Perfectly separating predictor: cases have higher probability than
  //    all controls, so AUC = 1.0 exactly.
  //
  // 2. Constant predictor: all predicted probabilities are equal, so all
  //    tied, and average-rank tie handling must give AUC = 0.5 exactly.
  //    If we skip average ranks and use raw ranks instead this comes out wrong.
  // ---------------------------------------------------------------------------
  {
    // Perfect separation: cases are samples 0,1,2; controls are 3,4,5.
    // pp_sep[0..2] > pp_sep[3..5].
    const double pp_sep[6]   = {0.9, 0.85, 0.8, 0.2, 0.15, 0.1};
    const double yy_sep[6]   = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};
    uint32_t idx_buf_sep[6];
    double auc_sep = plink2::AucMannWhitneyTest(pp_sep, yy_sep, 6, 3, idx_buf_sep);
    CheckClose("AUC perfect separation", auc_sep, 1.0, 1e-15);
  }
  {
    // Constant predictor: all pp equal, so AUC must be 0.5.
    const double pp_const[6] = {0.4, 0.4, 0.4, 0.4, 0.4, 0.4};
    const double yy_const[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};
    uint32_t idx_buf_const[6];
    double auc_const = plink2::AucMannWhitneyTest(pp_const, yy_const, 6, 3, idx_buf_const);
    CheckClose("AUC constant predictor", auc_const, 0.5, 1e-15);
  }

  // ---------------------------------------------------------------------------
  // Liability-scale R²: three closed-form tests verified in R during design.
  //
  // Let R2_obs = r2_incr and K = prevalence, P = sample case proportion.
  // The transform is: R2_liability = R2_obs * (K*(1-K))^2 / (z^2 * P*(1-P))
  // where z = dnorm(qnorm(1-K)).
  //
  // Test 1: at P == K the factor reduces exactly to K*(1-K)/z^2.
  //
  // Test 2: K=0.01, P=0.50 (50/50 ascertained): multiplier = 0.5519 < 1
  //   (shrinks because the sample over-represents cases).
  //
  // Test 3: K=0.01, P=0.01 (population sample): multiplier = 13.9371 > 1
  //   (inflates because cases are rare and the observed R² is compressed).
  //
  // Getting the formula upside-down swaps tests 2 and 3 in direction, which
  // is the exact check that distinguishes the correct formula from its inverse.
  // ---------------------------------------------------------------------------
  {
    // Test 1: at P == K the factor reduces exactly to K*(1-K)/z^2, regardless
    // of R2_obs.  We verify with R2_obs = 0.1 and K = 0.05.
    // In R: K=0.05; z=dnorm(qnorm(0.95)); K*(1-K)/z^2
    // R gives: 4.465561 (not 0.2040285; that was a transcription error during
    // plan authoring -- 0.2040285 is the factor for K=0.20 at P=K, not K=0.05).
    // The structural test (exact algebraic reduction to K*(1-K)/z^2) is what
    // matters; the reference value is computed from the same qnorm to avoid
    // circular verification.
    const double kk1 = 0.05;
    const double r2_obs1 = 0.1;
    const double r2_lia1 = plink2::LiabilityScaleR2Test(r2_obs1, kk1, kk1);
    // Expected: r2_obs1 * K*(1-K)/z^2
    const double t1 = plink2::QnormTest(1.0 - kk1);
    const double z1 = exp(-0.5 * t1 * t1) / sqrt(2.0 * M_PI);
    const double expected1 = r2_obs1 * kk1 * (1.0 - kk1) / (z1 * z1);
    CheckClose("liability P==K reduction", r2_lia1, expected1, 1e-9);
    // Verify against R's output (R: K=0.05; z=dnorm(qnorm(0.95)); K*(1-K)/z^2):
    CheckClose("liability factor K=0.05 P=K vs R", r2_lia1 / r2_obs1, 4.465561, 1e-4);
  }
  {
    // Test 2: K=0.01, P=0.50 (50/50 case-control): multiplier should be ~0.5519
    // R: K=0.01; z=dnorm(qnorm(0.99)); K^2*(1-K)^2/(z^2*0.5*0.5)
    // R gives: 0.5519...
    const double kk2 = 0.01;
    const double pp2 = 0.50;
    const double r2_obs2 = 0.1;
    const double factor2 = plink2::LiabilityScaleR2Test(r2_obs2, kk2, pp2) / r2_obs2;
    CheckClose("liability K=0.01 P=0.50 shrinks", factor2, 0.5519, 1e-3);
    ++g_checks;
    if (!(factor2 < 1.0)) {
      printf("  FAIL liability K=0.01 P=0.50 should shrink (factor=%.4f)\n", factor2);
      ++g_failures;
    }
  }
  {
    // Test 3: K=0.01, P=0.01 (population sample): multiplier should be ~13.9371
    // R: K=0.01; z=dnorm(qnorm(0.99)); K^2*(1-K)^2/(z^2*K*(1-K))
    //    = K*(1-K)/z^2 ~ 13.9371
    const double kk3 = 0.01;
    const double pp3 = 0.01;
    const double r2_obs3 = 0.1;
    const double factor3 = plink2::LiabilityScaleR2Test(r2_obs3, kk3, pp3) / r2_obs3;
    CheckClose("liability K=0.01 P=0.01 inflates", factor3, 13.9371, 1e-3);
    ++g_checks;
    if (!(factor3 > 1.0)) {
      printf("  FAIL liability K=0.01 P=0.01 should inflate (factor=%.4f)\n", factor3);
      ++g_failures;
    }
  }
  // qnorm spot checks (Wichura AS241 verification).
  CheckClose("qnorm(0.5)",   plink2::QnormTest(0.5),   0.0,                   1e-12);
  CheckClose("qnorm(0.975)", plink2::QnormTest(0.975),  1.959963984540054,     1e-9);
  CheckClose("qnorm(0.99)",  plink2::QnormTest(0.99),   2.3263478740408408,    1e-9);

  // ---------------------------------------------------------------------------
  // PgsMatchOne: validate against the committed bigsnpr fixture.
  //
  // The fixture has two sections:
  //   hand: 9 hand-designed rows covering every status.
  //   rand: 350 randomised rows (set.seed(20260802)).
  //
  // For each row we call PgsMatchOne with allow_ambiguous=0 and assert:
  //   - status matches the fixture status
  //   - out_beta matches the fixture out_beta (caller negates on reversed)
  //
  // out_beta for unmatched rows is NA in the fixture; we skip the beta check
  // for those (no output beta exists).
  //
  // The four signed betas 0.1, -0.2, 0.3, -0.4 are the critical ones: any
  // inversion of the sign convention (negating on flip, negating twice, etc.)
  // shows up immediately.
  // ---------------------------------------------------------------------------
  {
    // Fixture path relative to the project root, matching the OLS fixture.
    const char* fix_path = "Tests/pgs_match_fixtures/match_reference.txt";
    FILE* fix_f = fopen(fix_path, "r");
    if (!fix_f) {
      printf("  FAIL cannot open %s\n", fix_path);
      ++g_failures;
    } else {
      // Status string -> PgsMatchStatus mapping.
      auto ParseStatus = [](const char* s) -> plink2::PgsMatchStatus {
        if (!strcmp(s, "EXACT"))             return plink2::kPgsMatchExact;
        if (!strcmp(s, "REVERSED"))          return plink2::kPgsMatchReversed;
        if (!strcmp(s, "FLIPPED"))           return plink2::kPgsMatchFlipped;
        if (!strcmp(s, "FLIPPED_REVERSED"))  return plink2::kPgsMatchFlippedReversed;
        if (!strcmp(s, "AMBIGUOUS"))         return plink2::kPgsMatchAmbiguous;
        if (!strcmp(s, "NO_MATCH_POS"))      return plink2::kPgsMatchNoPos;
        if (!strcmp(s, "NO_MATCH_ALLELE"))   return plink2::kPgsMatchNoAllele;
        if (!strcmp(s, "DUPLICATE"))         return plink2::kPgsMatchDuplicate;
        return plink2::kPgsMatchNoAllele;  // unreachable
      };

      auto StatusName = [](plink2::PgsMatchStatus s) -> const char* {
        switch (s) {
          case plink2::kPgsMatchExact:            return "EXACT";
          case plink2::kPgsMatchReversed:         return "REVERSED";
          case plink2::kPgsMatchFlipped:          return "FLIPPED";
          case plink2::kPgsMatchFlippedReversed:  return "FLIPPED_REVERSED";
          case plink2::kPgsMatchAmbiguous:        return "AMBIGUOUS";
          case plink2::kPgsMatchNoPos:            return "NO_MATCH_POS";
          case plink2::kPgsMatchNoAllele:         return "NO_MATCH_ALLELE";
          case plink2::kPgsMatchDuplicate:        return "DUPLICATE";
          default:                                return "UNKNOWN";
        }
      };

      // We test against the target implied by the hand-designed case:
      // pos=100..700 step 100 with the documented alleles, plus pos=999.
      // For the randomised section we must reconstruct target alleles from
      // the fixture's out_a0/out_a1 (which equal the target alleles for
      // matched rows).  For unmatched rows we use dummy targets that won't
      // affect the decision (PgsMatchOne is called with the correct target
      // for matched rows, but for unmatched rows the status depends only on
      // the source alleles and the ambiguous flag; the target-allele pair
      // can be anything that provokes the same outcome).
      //
      // Strategy: for each row, derive the expected target alleles from the
      // fixture status + src alleles + out alleles, then call PgsMatchOne
      // and assert.
      //
      // For NO_MATCH_POS rows PgsMatchOne is never called by the real code
      // (the position is absent from the target); we still call it with a
      // dummy target to verify it returns NO_MATCH_ALLELE (or AMBIGUOUS),
      // then accept any status other than a false positive match.
      // Actually: NO_MATCH_POS is decided by the caller, not PgsMatchOne.
      // For testing purposes we reconstruct the situation PgsMatchOne would
      // see if called -- but since the position is absent the caller would
      // never call it.  We simply skip NO_MATCH_POS rows in the assertion
      // (the status is not a PgsMatchOne output).

      char line[512];
      int section_n = 0;   // 0=none, 1=hand, 2=rand
      int n_expected = 0;  // n_variants for current section
      int n_read = 0;

      // Hand-designed target table (same as generate_match.R):
      // index is 0-based position order.
      // tgt[i] = {a0, a1} for pos = 100*(i+1), i in [0..6], plus pos=999 at i=7.
      struct TgtRow { const char* a0; const char* a1; };
      static const TgtRow hand_tgt[8] = {
        {"A","G"}, {"A","G"}, {"A","G"}, {"A","G"},
        {"A","T"}, {"AT","A"}, {"C","G"}, {"C","T"}
      };
      auto HandTgtForPos = [&](int pos, const char** a0, const char** a1) -> int {
        const int positions[8] = {100,200,300,400,500,600,700,999};
        for (int k = 0; k < 8; ++k) {
          if (positions[k] == pos) { *a0 = hand_tgt[k].a0; *a1 = hand_tgt[k].a1; return 1; }
        }
        return 0;
      };

      while (fgets(line, (int)sizeof(line), fix_f)) {
        if (line[0] == '#') {
          // Check for section header or n_variants.
          if (strstr(line, "# section hand"))       { section_n = 1; n_read = 0; }
          else if (strstr(line, "# section rand"))  { section_n = 2; n_read = 0; }
          else if (strstr(line, "# n_variants "))   { sscanf(line, "# n_variants %d", &n_expected); }
          continue;
        }
        if (section_n == 0) continue;

        // Parse: chr pos src_a0 src_a1 src_beta status out_a0 out_a1 out_beta
        int   chr, pos;
        char  src_a0[64], src_a1[64];
        double src_beta;
        char  status_str[32];
        char  out_a0_str[64], out_a1_str[64], out_beta_str[64];
        int   nf = sscanf(line, "%d %d %63s %63s %lf %31s %63s %63s %63s",
                          &chr, &pos, src_a0, src_a1, &src_beta,
                          status_str, out_a0_str, out_a1_str, out_beta_str);
        if (nf < 9) continue;
        ++n_read;

        plink2::PgsMatchStatus exp_status = ParseStatus(status_str);

        // NO_MATCH_POS: PgsMatchOne is not called in production; skip.
        if (exp_status == plink2::kPgsMatchNoPos) continue;

        // Determine target alleles to supply to PgsMatchOne.
        const char* tgt_a0_p = nullptr;
        const char* tgt_a1_p = nullptr;

        if (section_n == 1) {
          // Hand section: look up from the fixed target table.
          const char* ta0 = nullptr; const char* ta1 = nullptr;
          if (!HandTgtForPos(pos, &ta0, &ta1)) {
            // pos=800 is NO_MATCH_POS, already skipped above.
            continue;
          }
          tgt_a0_p = ta0; tgt_a1_p = ta1;
        } else {
          // Rand section: for matched rows, out_a0/out_a1 ARE the target
          // alleles (snp_match outputs the target's allele coding).
          // For AMBIGUOUS and NO_MATCH_ALLELE, we use the source alleles as
          // a dummy target so that PgsMatchOne sees no match (it will return
          // AMBIGUOUS or NO_MATCH_ALLELE as expected).
          if (exp_status == plink2::kPgsMatchAmbiguous ||
              exp_status == plink2::kPgsMatchNoAllele) {
            // Any dummy that won't match the source works here; use "XX"
            // so none of the four candidates can equal it.
            tgt_a0_p = "XX"; tgt_a1_p = "XX";
          } else {
            tgt_a0_p = out_a0_str; tgt_a1_p = out_a1_str;
          }
        }

        // Call PgsMatchOne.
        uint32_t got_reversed = 99;
        plink2::PgsMatchStatus got_status =
            plink2::PgsMatchOne(src_a0, src_a1, tgt_a0_p, tgt_a1_p,
                                /*allow_ambiguous=*/0, /*resolve=*/nullptr,
                                &got_reversed);

        // Assert status.
        ++g_checks;
        if (got_status != exp_status) {
          printf("  FAIL pgs_match status row %d/%d (chr%d pos%d %s/%s): "
                 "got %s want %s\n",
                 n_read, n_expected, chr, pos, src_a0, src_a1,
                 StatusName(got_status), StatusName(exp_status));
          ++g_failures;
          continue;  // don't check beta on a status failure
        }

        // Assert out_beta for matched rows.
        // out_beta = src_beta if !reversed, -src_beta if reversed.
        // The fixture records the snp_match output beta, which already
        // carries the sign; our out_beta = src_beta * (reversed ? -1 : 1).
        if (exp_status == plink2::kPgsMatchExact ||
            exp_status == plink2::kPgsMatchReversed ||
            exp_status == plink2::kPgsMatchFlipped ||
            exp_status == plink2::kPgsMatchFlippedReversed) {
          double exp_beta_val;
          if (sscanf(out_beta_str, "%lf", &exp_beta_val) != 1) continue;
          double got_beta = src_beta * (got_reversed ? -1.0 : 1.0);
          ++g_checks;
          const double diff = fabs(got_beta - exp_beta_val);
          const double scale = (exp_beta_val != 0.0) ? fabs(exp_beta_val) : 1.0;
          if (diff / scale > 1e-12) {
            printf("  FAIL pgs_match beta row %d/%d (chr%d pos%d %s/%s %s): "
                   "got %.12g want %.12g\n",
                   n_read, n_expected, chr, pos, src_a0, src_a1,
                   StatusName(exp_status), got_beta, exp_beta_val);
            ++g_failures;
          }
        }
      }
      fclose(fix_f);
    }
  }

  // ---------------------------------------------------------------------------
  // PgsMatchOne with --resolve-ambig: palindromes decided by frequency.
  //
  // No file I/O -- the whole decision is a pure function of six values, so it
  // is exercised directly.
  //
  // The sign is the point.  A palindrome's exact and flipped+reversed forms
  // carry identical allele labels, so nothing but the frequency distinguishes
  // "the source effect allele is the target's a1" from "it is the target's
  // a0".  Get the comparison backwards and the beta comes out with the wrong
  // sign and no other symptom, which is why every case below asserts the beta
  // and not just the status.
  // ---------------------------------------------------------------------------
  {
    // Runs one resolution, asserts status and the reversal flag, and returns
    // the beta the CALLER would emit -- src_beta negated on reversal, exactly
    // as RunPgsModel does at its single negation site.  Returning it (rather
    // than reading a shared `reversed` out-param afterwards) keeps the beta
    // assertions tied to the call that produced them.
    const double kSrcBeta = 0.3;
    auto CheckResolve = [kSrcBeta](
        const char* what,
        const char* src_a0, const char* src_a1,
        const char* tgt_a0, const char* tgt_a1,
        double src_freq, double tgt_a1_freq, double maf_thresh,
        plink2::PgsMatchStatus want, uint32_t want_reversed) -> double {
      plink2::PgsAmbigResolve r;
      r.maf_threshold = maf_thresh;
      r.src_freq      = src_freq;
      r.tgt_a1_freq   = tgt_a1_freq;
      uint32_t got_reversed = 99;
      const plink2::PgsMatchStatus got =
          plink2::PgsMatchOne(src_a0, src_a1, tgt_a0, tgt_a1,
                              /*allow_ambiguous=*/0, &r, &got_reversed);
      ++g_checks;
      if ((got != want) || (got_reversed != want_reversed)) {
        printf("  FAIL %-46s got status %d rev %u  want status %d rev %u\n",
               what, (int)got, got_reversed, (int)want, want_reversed);
        ++g_failures;
      }
      return got_reversed ? -kSrcBeta : kSrcBeta;
    };
    auto CheckSign = [](const char* what, double beta, double want_sign) {
      ++g_checks;
      if (beta * want_sign <= 0.0) {
        printf("  FAIL %-46s beta %.12g has the wrong sign\n", what, beta);
        ++g_failures;
      }
    };
    const double kNaN = nan("");

    // The case the feature exists for: the source study reported the other
    // allele of the pair as its effect allele.  Target a1 = A at 0.85, source
    // effect allele T at 0.15 -- the same physical allele, so the source
    // effect allele is the target's a0 and the beta must flip sign.  A
    // backwards comparison yields RESOLVED_EXACT and an un-negated beta.
    CheckSign("resolve: strand-swapped beta negated",
              CheckResolve("resolve: strand-swapped palindrome -> REVERSED",
                           "A", "T", "T", "A", 0.15, 0.85, 0.4,
                           plink2::kPgsMatchResolvedReversed, 1),
              -1.0);
    // Same variant and same target frequency, source effect allele on the
    // target's a1: beta untouched.
    CheckSign("resolve: aligned beta unchanged",
              CheckResolve("resolve: aligned palindrome -> EXACT",
                           "T", "A", "T", "A", 0.85, 0.85, 0.4,
                           plink2::kPgsMatchResolvedExact, 0),
              1.0);
    // C/G behaves the same as A/T.
    CheckSign("resolve: C/G swapped beta negated",
              CheckResolve("resolve: C/G swapped -> REVERSED",
                           "G", "C", "C", "G", 0.2, 0.8, 0.4,
                           plink2::kPgsMatchResolvedReversed, 1),
              -1.0);
    CheckSign("resolve: C/G aligned beta unchanged",
              CheckResolve("resolve: C/G aligned -> EXACT",
                           "C", "G", "C", "G", 0.8, 0.8, 0.4,
                           plink2::kPgsMatchResolvedExact, 0),
              1.0);

    // Near 0.5 the two hypotheses are indistinguishable whatever the target
    // says, so both directions of target frequency must refuse.
    CheckResolve("resolve: src MAF 0.45 -> UNRESOLVED (tgt 0.45)",
                 "T", "A", "T", "A", 0.45, 0.45, 0.4,
                 plink2::kPgsMatchAmbiguousUnresolved, 0);
    CheckResolve("resolve: src MAF 0.45 -> UNRESOLVED (tgt 0.55)",
                 "T", "A", "T", "A", 0.45, 0.55, 0.4,
                 plink2::kPgsMatchAmbiguousUnresolved, 0);
    CheckResolve("resolve: src freq 0.5 -> UNRESOLVED",
                 "T", "A", "T", "A", 0.5, 0.02, 0.4,
                 plink2::kPgsMatchAmbiguousUnresolved, 0);
    // The threshold is on the MINOR allele frequency: a source frequency of
    // 0.58 has MAF 0.42, so it is refused at 0.4 and allowed at 0.5.  A
    // threshold compared against the raw frequency would let 0.58 through.
    CheckResolve("resolve: src freq 0.58 -> UNRESOLVED at 0.4",
                 "T", "A", "T", "A", 0.58, 0.05, 0.4,
                 plink2::kPgsMatchAmbiguousUnresolved, 0);
    CheckResolve("resolve: src freq 0.58 resolves at threshold 0.5",
                 "T", "A", "T", "A", 0.58, 0.05, 0.5,
                 plink2::kPgsMatchResolvedReversed, 1);

    // The SAME threshold applies to the target frequency.  The signal is
    // |2f - 1| on each side, so a target near 0.5 is exactly as uninformative
    // as a source near 0.5 -- and a target of exactly 0.5 makes the two
    // distances equal, which without this guard fell through to
    // RESOLVED_REVERSED and negated the beta on zero evidence.  The source
    // frequency here is 0.01, far below any threshold, so only the target
    // guard can produce these.
    CheckResolve("resolve: tgt freq exactly 0.5 -> UNRESOLVED",
                 "T", "A", "T", "A", 0.01, 0.5, 0.4,
                 plink2::kPgsMatchAmbiguousUnresolved, 0);
    // ...and it is the tie that is refused, not a threshold coincidence: 0.5
    // is still refused when the threshold is as small as the CLI allows.
    CheckResolve("resolve: tgt freq 0.5 -> UNRESOLVED at threshold 0.01",
                 "T", "A", "T", "A", 0.01, 0.5, 0.01,
                 plink2::kPgsMatchAmbiguousUnresolved, 0);
    // The boundary, approached from both sides of 0.5 and from inside/outside.
    CheckResolve("resolve: tgt MAF 0.4 (freq 0.6) -> UNRESOLVED at 0.4",
                 "T", "A", "T", "A", 0.01, 0.6, 0.4,
                 plink2::kPgsMatchAmbiguousUnresolved, 0);
    CheckResolve("resolve: tgt MAF 0.4 (freq 0.4) -> UNRESOLVED at 0.4",
                 "T", "A", "T", "A", 0.01, 0.4, 0.4,
                 plink2::kPgsMatchAmbiguousUnresolved, 0);
    CheckSign("resolve: tgt MAF 0.39 (freq 0.61) beta negated",
              CheckResolve("resolve: tgt MAF 0.39 (freq 0.61) -> REVERSED",
                           "T", "A", "T", "A", 0.01, 0.61, 0.4,
                           plink2::kPgsMatchResolvedReversed, 1),
              -1.0);
    CheckSign("resolve: tgt MAF 0.39 (freq 0.39) beta unchanged",
              CheckResolve("resolve: tgt MAF 0.39 (freq 0.39) -> EXACT",
                           "T", "A", "T", "A", 0.01, 0.39, 0.4,
                           plink2::kPgsMatchResolvedExact, 0),
              1.0);

    // No usable frequency on either side: refuse rather than assume.
    CheckResolve("resolve: target freq NaN -> UNRESOLVED",
                 "T", "A", "T", "A", 0.15, kNaN, 0.4,
                 plink2::kPgsMatchAmbiguousUnresolved, 0);
    CheckResolve("resolve: source freq NaN -> UNRESOLVED",
                 "T", "A", "T", "A", kNaN, 0.15, 0.4,
                 plink2::kPgsMatchAmbiguousUnresolved, 0);

    // Palindromic on both sides but not the same pair: no shared allele, so
    // no amount of frequency agreement makes this a match.
    CheckResolve("resolve: A/T source vs C/G target -> NO_MATCH_ALLELE",
                 "T", "A", "C", "G", 0.15, 0.15, 0.4,
                 plink2::kPgsMatchNoAllele, 0);
    CheckResolve("resolve: C/G source vs A/T target -> NO_MATCH_ALLELE",
                 "G", "C", "A", "T", 0.15, 0.15, 0.4,
                 plink2::kPgsMatchNoAllele, 0);

    // A non-palindromic variant never reaches the frequency logic, even when
    // handed frequencies that would have flipped it.
    CheckResolve("resolve: non-palindromic still EXACT on alleles",
                 "A", "G", "A", "G", 0.15, 0.85, 0.4,
                 plink2::kPgsMatchExact, 0);
    CheckResolve("resolve: non-palindromic still REVERSED on alleles",
                 "G", "A", "A", "G", 0.15, 0.85, 0.4,
                 plink2::kPgsMatchReversed, 1);

    // Without the resolve struct the default is unchanged: dropped.
    ++g_checks;
    {
      uint32_t r0 = 99;
      if (plink2::PgsMatchOne("A", "T", "T", "A", 0, nullptr, &r0) !=
          plink2::kPgsMatchAmbiguous) {
        printf("  FAIL resolve: palindrome not AMBIGUOUS without resolve\n");
        ++g_failures;
      }
    }
  }

  printf("%d checks, %d failures\n", g_checks, g_failures);
  return g_failures ? 1 : 0;
}
