#ifndef __ADJUTANT_PGS_EVAL_H__
#define __ADJUTANT_PGS_EVAL_H__

// Rank PGS models by incremental variance explained over a covariates-only
// null model.
//
// Deliberately takes plain double arrays: no .pgen, no bigstack, no PhenoCol.
// That boundary is what lets these numerics be checked against statsmodels and
// R with plain fixtures, which is the only real evidence that a regression
// implementation is correct.

#include "include/plink2_base.h"
#include "adjutant_cmdline.h"
#include "adjutant_common.h"
#include "adjutant_matrix.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// PgsEvalOlsOne writes every field on every path, converged included, so this
// struct has no partially-initialised state.  Keep it that way: the logistic
// path passes the caller's struct straight into PgsEvalOlsOne rather than
// copying selected fields out of a local, so a field added here is filled in on
// both paths without a copy-out to keep in sync.
typedef struct PgsEvalResultStruct {
  double beta;      // per standard deviation of the PGS
  double se;
  double pval;      // two-sided
  // Observed-scale R^2 of the covariates-only model, of the model with the PGS
  // added, and their difference.  For a binary phenotype these come from an OLS
  // fit on the 0/1 outcome (see PgsEvalLogisticOne), not from the logistic fit,
  // so r2_incr is the same statistic for both trait types.
  double r2_null;
  double r2_full;
  double r2_incr;
  // Analytic 95% CI on r2_incr.  Whether the winner of a ranking is separable
  // from the runner-up is not answerable from a point estimate, and the tool's
  // whole output is a winner.  NaN when sample_ct <= predictor_ct + 2, which
  // leaves no degrees of freedom for the Fisher z standard error, and under
  // cross-validation when r2_incr < 0 (see PgsEvalOlsSetFolds).
  double r2_incr_l95;
  double r2_incr_u95;
  // Binary-trait fields.  PgsEvalOlsOne sets them to NaN; only the logistic
  // path fills them in, and only the CC output columns read them.
  double r2_liability;  // liability-scale R^2 (Lee 2012); NaN without prevalence
  double r2_liability_l95;
  double r2_liability_u95;
  double auc_null;      // AUC of the covariates-only model
  double auc_full;      // AUC with the PGS added
  uint32_t converged;   // 0 if the column was degenerate and nothing was fit
} PgsEvalResult;

// Everything that depends only on y and the covariates, computed once and
// reused for every PGS column.
typedef struct PgsEvalOlsCtxStruct {
  const double* xx;   // borrowed: column-major, sample_ct * predictor_ct,
                      // column 0 is the intercept
  const double* yy;   // borrowed, like xx; only the fold_ct > 1 path reads it
  double* xtx_inv;    // predictor_ct * predictor_ct
  double* e_y;        // residuals of y on xx
  double* work;       // scratch, >= sample_ct + 2 * predictor_ct
                      // Layout: [0, predictor_ct) = xtv, [predictor_ct, 2*predictor_ct) = coef,
                      // [2*predictor_ct, 2*predictor_ct+sample_ct) = e_g.  e_g MUST sit after
                      // both xtv and coef: Residualize's first loop writes xtv while still
                      // reading vv, which is what makes in-place residualisation safe.
  double tss;
  double rss_null;
  double r2_null;
  uint32_t sample_ct;
  uint32_t predictor_ct;
  // --- k-fold cross-validation state.  fold_ct == 1 and every pointer null
  // unless PgsEvalOlsSetFolds has been called; at fold_ct == 1 PgsEvalOlsOne
  // does not read any of it, which is what makes --folds 1 bit-identical to
  // the no-flag output. ---
  uint32_t fold_ct;
  uint32_t* fold_assign;   // sample_ct: which fold each sample is held out in
  double* fold_xtx_inv;    // fold_ct * predictor_ct^2, one per training set
  double* fold_beta;       // fold_ct * predictor_ct, null-model coefficients
  double* cv_work;         // 3 * sample_ct + 2 * predictor_ct
                           // Layout: [0, n) = out-of-fold null residuals
                           // y - yhat_null (fixed by PgsEvalOlsSetFolds),
                           // [n, 2n) = per-column incremental prediction,
                           // [2n, 3n) = standardised PGS column,
                           // [3n, 3n+p) = X'g, [3n+p, 3n+2p) = g-on-X coefs.
  double r2_null_oof;      // 1 - sum(cv_work[0, n)^2) / tss
} PgsEvalOlsCtx;

// Two-sided p-value for Student's t on df degrees of freedom.
double PgsEvalTPval(double tt, double df);

// Factor X'X, residualise y, and record TSS / RSS_null / R2_null.  `xx` and
// `yy` are borrowed and must outlive the context.  Returns 1 on allocation
// failure or if X is rank-deficient -- which is what catches a constant
// covariate colliding with the intercept column.
BoolErr PgsEvalOlsInit(const double* yy, const double* xx, uint32_t sample_ct,
                       uint32_t predictor_ct, PgsEvalOlsCtx* ctx);

void PgsEvalOlsCleanup(PgsEvalOlsCtx* ctx);

// Switch this context to k-fold cross-validated R^2.  Call after
// PgsEvalOlsInit and before the first PgsEvalOlsOne; fold_ct == 1 is a no-op,
// so callers can call it unconditionally.
//
// Assigns folds, factors X'X on each of the fold_ct training sets, and records
// the pooled out-of-fold null residuals.  All of that depends only on y and X,
// so it is paid once rather than once per PGS column.
//
// From then on PgsEvalOlsOne reports r2_null / r2_full / r2_incr over the
// pooled out-of-fold predictions -- each sample predicted by the fit that did
// not see it -- while beta / se / pval stay the in-sample estimates from all
// the data.  That mixture is deliberate: the coefficient answers "how big is
// the effect", for which the whole sample is the most efficient estimate,
// while the cross-validated R^2 answers "how well does this generalise", which
// an in-sample R^2 overstates by exactly the amount the ranking exploits.
//
// stratify_binary: split cases and controls separately, so an imbalanced
// phenotype cannot put nearly all its cases in one fold.  Pass 1 for a 0/1 y.
// seed: fold assignment is deterministic given (seed, y, sample order).
//
// Returns 1 on allocation failure, if a training set has too few samples for
// the predictors, or if X'X is singular on one -- the last is what catches a
// covariate that is constant within a fold.
BoolErr PgsEvalOlsSetFolds(PgsEvalOlsCtx* ctx, uint32_t fold_ct,
                           uint32_t stratify_binary, uint64_t seed);

// Fit y ~ 1 + covars + gg by Frisch-Waugh-Lovell.  Algebraically identical to
// refitting the full model, not an approximation, and needs only this one PGS
// column -- so peak memory does not grow with the number of models.
//
// `gg` is standardised internally, so beta is per standard deviation.
// Sets result->converged = 0 and leaves the estimates as NaN on early return.
// Early-return conditions (all leave converged = 0):
//   - zero-variance PGS (all weights matched the same genotype);
//   - collinear PGS (gg_dot <= 1e-10 * (sample_ct - 1) after residualisation);
//   - degenerate fit (rss_full <= 0 due to floating-point cancellation when
//     the PGS nearly perfectly explains all remaining variance);
//   - insufficient degrees of freedom (sample_ct <= predictor_ct + 1);
//   - with folds set, a PGS collinear with the covariates inside one training
//     set even though the whole-sample fit was fine.
//
// THREADING: ctx->work and ctx->cv_work are scratch and are overwritten on
// every call despite the const pointer.  PgsEvalOlsOne is NOT reentrant; use
// one ctx per thread.
//
// PRECONDITION: all elements of yy and xx passed to PgsEvalOlsInit must be
// finite.  A NaN in either passes PgsEvalOlsInit undetected (dpocon returns
// NaN and NaN < threshold is false) and produces NaN estimates with
// converged = 1.
void PgsEvalOlsOne(const PgsEvalOlsCtx* ctx, const double* gg,
                   PgsEvalResult* result);

// ---------------------------------------------------------------------------
// Logistic regression path
// ---------------------------------------------------------------------------
// Logistic fit of a binary phenotype.  Frisch-Waugh-Lovell does not apply, so
// unlike the OLS path each model needs a genuine refit; only the null model is
// shared.  `yy` must be 0/1.
typedef struct PgsEvalLogisticCtxStruct {
  const double* xx;      // borrowed: column-major, sample_ct * predictor_ct
  const double* yy;      // borrowed: 0/1
  double* beta;          // predictor_ct null-model coefficients.
                         // PgsEvalLogisticOne never touches them (it allocates
                         // its own full-model beta) and nothing else in the
                         // main binary reads them either: they are kept so
                         // adjutant_pgs_eval_test.cc can pin the null fit
                         // against statsmodels, which is the only direct check
                         // the shared null model gets.
  double* work;          // see PgsEvalLogisticInit for the layout
  // The OLS companion fit, on the same yy and xx, which supplies the
  // observed-scale R² family.  Owned: built by PgsEvalLogisticInit and released
  // by PgsEvalLogisticCleanup, so it cannot be built on a different design
  // matrix than the logistic fit.
  PgsEvalOlsCtx ols;
  double auc_null;       // AUC of the covariates-only model
  uint32_t sample_ct;
  uint32_t predictor_ct;
  uint32_t case_ct;
} PgsEvalLogisticCtx;

BoolErr PgsEvalLogisticInit(const double* yy, const double* xx,
                            uint32_t sample_ct, uint32_t predictor_ct,
                            PgsEvalLogisticCtx* ctx);
void PgsEvalLogisticCleanup(PgsEvalLogisticCtx* ctx);

// Fit y ~ 1 + covars + gg by Newton-Raphson IRLS -- a standard binomial GLM,
// the same estimator as R's glm(family=binomial) and statsmodels Logit.
//
// Two models are fitted per call.  The logistic fit supplies beta, se, pval and
// auc_full.  The R² family (r2_null, r2_full, r2_incr, and hence r2_liability)
// comes from PgsEvalOlsOne on the same 0/1 outcome, called internally through
// ctx->ols.  Lee et al. (2012) derive the liability-scale transform from an
// observed-scale R² produced by a linear model, so that is the R² the transform
// is applied to; reporting an odds ratio from the logistic fit alongside it is
// standard practice.
//
// Sets result->converged = 0 if gg is degenerate, if the logistic fit does not
// converge (which is what complete separation produces), or if the OLS fit
// takes an early return.  The last case is deliberate: r2_incr is the sort key,
// so a row without one is useless even if the logistic coefficients came out
// fine.  With converged = 0 the estimates are not to be read: some are NaN, and
// once the OLS fit has run beta / se / pval / the R² family hold its values
// rather than the logistic fit's.  Both writers emit NA for the whole row.
//
// prevalence: population case proportion K for liability-scale R² (Lee 2012).
// Pass prevalence < 0 (e.g. -1.0) when not provided; result->r2_liability will
// then be NaN.  Prevalence is passed here rather than stored on the context
// because it is user-supplied per-analysis and does not depend on y, whereas
// the context carries y-dependent state fit once and reused per PGS column.
void PgsEvalLogisticOne(PgsEvalLogisticCtx* ctx, const double* gg,
                        double prevalence, PgsEvalResult* result);

// ---------------------------------------------------------------------------
// Test-visible wrappers for internal metrics functions.
// These forward to the static implementations in adjutant_pgs_eval.cc.
// They are only compiled into pgs_eval_test; the main binary does not
// reference them and the linker will discard them if they are unreferenced.
// ---------------------------------------------------------------------------
double AucMannWhitneyTest(const double* pp, const double* yy,
                          uint32_t sample_ct, uint32_t case_ct,
                          uint32_t* idx_buf);
double LiabilityScaleR2Test(double r2_observed, double prevalence,
                             double case_proportion);
double QnormTest(double pp);
// The Fisher z core the cross-validated path reaches with its own partial
// correlation.  Exposed so the test can feed it the in-sample r and check that
// it returns the in-sample interval -- i.e. that there really is one formula.
void IncrR2CiFromRTest(double rr, double r2_null, uint32_t sample_ct,
                       uint32_t predictor_ct, double* lo, double* hi);

// Everything checkable without the scores.  Plink2Core calls this before
// ScoreReport so a misconfigured run fails in a second rather than after the
// scoring pass.  BestPgsReport still calls it, so the unit stays correct if
// used standalone.
PglErr BestPgsValidate(const uintptr_t* sample_include, uint32_t sample_ct,
                       const PhenoCol* pheno_cols, const char* pheno_names,
                       uintptr_t max_pheno_name_blen, uint32_t pheno_ct,
                       const PhenoCol* covar_cols, uint32_t covar_ct,
                       uint32_t fold_ct,
                       uint32_t* usable_pheno_idxs, uint32_t* usable_pheno_ct_ptr);

// --pgs-pheno-pairs: the (phenotype, PGS) combinations to evaluate.
//
// A filter, not a selector.  Every phenotype in the --pheno file is still
// loaded; these pairs restrict which combinations get fitted, at the one place
// both input routes share.  Consequently the flag is mutually exclusive with
// --pheno-name / --pheno-col-nums: the pairs file already says which
// phenotypes matter, and a second flag saying so too would be two sources of
// truth for one question.
typedef struct PgsPhenoPairsStruct {
  char** pheno_names;   // pair_ct strings, owned
  char** pgs_names;     // pair_ct strings, owned
  uint32_t pair_ct;
} PgsPhenoPairs;

void PgsPhenoPairsPreinit(PgsPhenoPairs* pp);
void PgsPhenoPairsFree(PgsPhenoPairs* pp);

// Two whitespace-separated columns per line; blank lines and lines whose first
// non-space character is '#' are skipped.  A line with one column, or a file
// with no pairs at all, is an error.
PglErr LoadPgsPhenoPairs(const char* fname, PgsPhenoPairs* pp);

// Every name in the file must exist, or the run stops.  A silently dropped
// typo would leave a report that looks complete and is not.
PglErr PgsPairsValidate(const PgsPhenoPairs* pp,
                        const char* const* pheno_names, uint32_t pheno_ct,
                        const char* const* score_col_names, uint32_t score_col_ct);

// Score-column indices paired with `pheno_name`, in score-column order (not
// pairs-file order, so the row list does not depend on how the file was
// written).  Writes at most score_col_ct entries; returns how many.
uint32_t PgsPairsKeepIdxs(const PgsPhenoPairs* pp, const char* pheno_name,
                          const char* const* score_col_names,
                          uint32_t score_col_ct, uint32_t* keep_idxs);

// --give-corr-matrix: Pearson correlation matrix across the PGS columns of a
// column-major [model_ct x sample_ct] score matrix, written to
// <outname><outname_end>.pgs_corr as a full square matrix.
//
// Both input routes call this one helper with their own matrix, which is what
// stops "the correlation between two PGSs" acquiring two definitions.
//
// Computed over EVERY sample in the matrix — no phenotype, no covariate, no
// complete-case restriction.  A column that is constant or holds a non-finite
// value gets NA across its whole row and column, not 0.
// Errors if sample_ct < 2.
PglErr PgsCorrReport(const double* score_data, uint32_t model_ct,
                     uint32_t sample_ct, const char* const* model_names,
                     char* outname, char* outname_end);

// Read back the .sscore file that ScoreReport already wrote, rank the PGS
// models by incremental R² over the covariates-only null, and write
// <outname><outname_end>.best_pgs.
//
// sscore_fname is the full path (possibly .zst) produced by ScoreReport.
// sample_ct is the count of set bits in sample_include — i.e. the number of
// data rows in the .sscore file.
// pheno_names and covar_names are flat char arrays with stride
// max_pheno_name_blen / max_covar_name_blen respectively.
// outname / outname_end: same convention as ScoreReport; on return outname
// contains the path of the last file written.
//
// give_best_pgs / corr_matrix: which of the two reports to produce; at least
// one must be set.  With corr_matrix alone the phenotype validation and the
// whole fitting loop are skipped, so pheno_cols may be empty.
PglErr BestPgsReport(const char* sscore_fname,
                     const uintptr_t* sample_include, uint32_t sample_ct,
                     const PhenoCol* pheno_cols, const char* pheno_names,
                     uintptr_t max_pheno_name_blen, uint32_t pheno_ct,
                     const PhenoCol* covar_cols, const char* covar_names,
                     uintptr_t max_covar_name_blen, uint32_t covar_ct,
                     double prevalence, uint32_t fold_ct, uint64_t fold_seed,
                     uint32_t give_best_pgs, uint32_t corr_matrix,
                     uint32_t no_na, const char* pairs_fname,
                     char* outname, char* outname_end);

// --sscore standalone input mode: read a user-supplied .sscore, phenotype
// file, and optional covariate file; join on IID; rank and write .best_pgs.
// No .pgen/.psam involved.
//
// score_col_range: from --score-col-nums (pass nullptr or a range_list with
//   name_ct=0 to use the auto-detection heuristic).
// pheno_fname: nullptr means "scores only" — legal only with corr_matrix set,
//   in which case the .pgs_corr is written and nothing is fitted.
// pheno_range_list: from --pheno-name or --pheno-col-nums (nullptr or
//   name_ct=0 means select all phenotype columns).
// pairs_fname: from --pgs-pheno-pairs; nullptr means the full cross-product.
PglErr BestPgsFromSscore(const char* sscore_fname,
                         const RangeList* score_col_range,
                         const char* pheno_fname,
                         const RangeList* pheno_range_list,
                         const char* covar_fname,
                         double prevalence, uint32_t fold_ct, uint64_t fold_seed,
                         uint32_t corr_matrix, uint32_t no_na,
                         const char* pairs_fname,
                         char* outname, char* outname_end);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __ADJUTANT_PGS_EVAL_H__
