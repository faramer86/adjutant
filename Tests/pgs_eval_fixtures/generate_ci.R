#!/usr/bin/env Rscript
# Confidence-interval reference values for adjutant_pgs_eval, in R.
#
#     Rscript Tests/pgs_eval_fixtures/generate_ci.R
#
# Run once; the output is committed.  R is not needed to build or run adjutant,
# or to run pgs_eval_check.sh.
#
# Deliberately a *second* implementation, in a second language, of the interval
# adjutant computes -- the formula is written out here from the definition
# rather than ported from adjutant_pgs_eval.cc, and the residuals come from R's
# lm() (Householder QR) rather than from the Frisch-Waugh-Lovell path under
# test.  A port would only prove the C compiles.
#
# It reads the data already committed in ols_reference.txt and
# logistic_reference.txt and rewrites those files with the CI keys inserted
# before the '# data:' marker; it never touches the data block or the values
# generate.py / generate_logistic.py produced.  Re-running is idempotent.

# --- The interval -----------------------------------------------------------
# Fisher z transform of the partial correlation between y and g, both
# residualised on the covariates:
#
#   r  = cor(resid(y ~ X), resid(g ~ X))
#   z  = atanh(r),  SE(z) = 1 / sqrt(n - p - 2)
#   CI(r)    = tanh(z +/- qnorm(0.975) * SE)
#   CI(dR2)  = CI(r)^2 * (1 - R2_null)
#
# Squaring an interval that straddles zero has minimum 0, not lo^2.
#
# X carries its own intercept column, so R² is computed from residuals here
# rather than taken from summary.lm(): given a no-intercept formula summary.lm
# reports the *uncentered* R², 1 - RSS/sum(y^2), which is a different statistic
# from the one adjutant reports and silently inflates every interval below.
Resid <- function(X, v) qr.resid(qr(X), v)
Tss <- function(y) sum((y - mean(y))^2)
R2Null <- function(y, X) 1 - sum(Resid(X, y)^2) / Tss(y)

AnalyticCi <- function(y, X, g) {
  n <- length(y)
  p <- ncol(X)                      # X includes the intercept column
  e_y <- Resid(X, y)
  e_g <- Resid(X, g)
  r <- sum(e_y * e_g) / sqrt(sum(e_y * e_y) * sum(e_g * e_g))
  se <- 1 / sqrt(n - p - 2)
  lo_r <- tanh(atanh(r) - qnorm(0.975) * se)
  hi_r <- tanh(atanh(r) + qnorm(0.975) * se)
  lo_sq <- if (lo_r <= 0 && hi_r >= 0) 0 else min(lo_r^2, hi_r^2)
  hi_sq <- max(lo_r^2, hi_r^2)
  scale <- 1 - R2Null(y, X)
  c(lo = lo_sq * scale, hi = hi_sq * scale)
}

# --- The cross-check --------------------------------------------------------
# Percentile bootstrap on the same statistic.  This is a check on the *method*,
# not a tolerance target: it says the analytic interval is in the right place
# and the right width, which is all a closed form can be asked to agree with a
# resampling estimate about.
BootstrapCi <- function(y, X, g, B = 2000, seed = 20260802) {
  set.seed(seed)
  n <- length(y)
  Xg <- cbind(X, g)
  d <- numeric(B)
  for (b in seq_len(B)) {
    i <- sample.int(n, n, replace = TRUE)
    yb <- y[i]
    rss_null <- sum(Resid(X[i, , drop = FALSE], yb)^2)
    rss_full <- sum(Resid(Xg[i, , drop = FALSE], yb)^2)
    d[b] <- (rss_null - rss_full) / Tss(yb)
  }
  q <- unname(quantile(d, c(0.025, 0.975)))
  c(lo = q[1], hi = q[2])
}

# Liability-scale transform (Lee et al. 2012).  A positive constant given K and
# P, so it maps interval endpoints straight across.
LiabilityFactor <- function(K, P) {
  t <- qnorm(1 - K)
  z <- dnorm(t)
  (K * (1 - K))^2 / (z^2 * P * (1 - P))
}

# The near-null PGS for the straddle-zero case.  Integer arithmetic so R and C
# produce bit-identical values from the same one-line formula.
NullPgs <- function(n) ((seq_len(n) - 1) * 37) %% 101

# --- Fixture I/O ------------------------------------------------------------
ReadFixture <- function(path) {
  lines <- readLines(path)
  data_at <- grep("^# data", lines)[1]
  hdr <- lines[seq_len(data_at)]
  n <- as.integer(sub("^n ", "", grep("^n ", hdr, value = TRUE)[1]))
  n_covar <- as.integer(sub("^n_covar ", "", grep("^n_covar ", hdr, value = TRUE)[1]))
  m <- matrix(scan(text = lines[-seq_len(data_at)], quiet = TRUE),
              ncol = n_covar + 2, byrow = TRUE)
  stopifnot(nrow(m) == n)
  list(lines = lines, data_at = data_at, n = n, n_covar = n_covar,
       y = m[, 1], covars = m[, 2:(n_covar + 1), drop = FALSE],
       pgs = m[, n_covar + 2])
}

# Insert `block` immediately before the '# data:' marker, dropping any earlier
# run of the same keys so re-running does not accumulate duplicates.  Comment
# lines are matched whole, not by first token: every comment's first token is
# "#", so treating that as a key would delete the fixture's own provenance
# header.
WriteFixture <- function(fx, path, block) {
  keys <- setdiff(sub(" .*", "", block), "#")
  first <- sub(" .*", "", fx$lines)
  keep <- !(first != "#" & first %in% keys) & !(fx$lines %in% block)
  keep[fx$data_at:length(keep)] <- TRUE   # never filter the data block
  lines <- fx$lines[keep]
  at <- grep("^# data", lines)[1]
  writeLines(c(lines[seq_len(at - 1)], block, lines[at:length(lines)]), path)
}

Fmt <- function(key, v) sprintf("%s %.17g", key, v)

Emit <- function(path, extra_prevalence = NA) {
  fx <- ReadFixture(path)
  X <- cbind(1, fx$covars)
  ana <- AnalyticCi(fx$y, X, fx$pgs)
  boot <- BootstrapCi(fx$y, X, fx$pgs)
  g0 <- NullPgs(fx$n)
  ana0 <- AnalyticCi(fx$y, X, g0)
  e_y <- Resid(X, fx$y)
  e_g0 <- Resid(X, g0)
  r2_incr0 <- (sum(e_y * e_g0)^2 / sum(e_g0 * e_g0)) / Tss(fx$y)

  block <- c(
    sprintf("# CI reference added by generate_ci.R -- %s", R.version.string),
    Fmt("r2_incr_l95", ana[["lo"]]),
    Fmt("r2_incr_u95", ana[["hi"]]),
    "# percentile bootstrap, B=2000, seed 20260802: a check on the method, not",
    "# a tolerance target.  The test asserts overlap, never equality.",
    Fmt("r2_incr_boot_l95", boot[["lo"]]),
    Fmt("r2_incr_boot_u95", boot[["hi"]]),
    "# straddle-zero case: pgs_null[i] = ((i * 37) %% 101), i = 0..n-1, which",
    "# has almost no relation to y, so the interval on r contains zero and the",
    "# lower bound on dR2 must be exactly 0.",
    Fmt("null_r2_incr", r2_incr0),
    Fmt("null_r2_incr_l95", ana0[["lo"]]),
    Fmt("null_r2_incr_u95", ana0[["hi"]]))

  if (!is.na(extra_prevalence)) {
    P <- mean(fx$y)
    f <- LiabilityFactor(extra_prevalence, P)
    block <- c(block,
      "# liability scale (Lee 2012) at this prevalence: the observed-scale",
      "# interval times a positive constant.",
      Fmt("prevalence", extra_prevalence),
      Fmt("r2l_incr_l95", f * ana[["lo"]]),
      Fmt("r2l_incr_u95", f * ana[["hi"]]))
  }

  WriteFixture(fx, path, block)
  cat(sprintf("%s\n  analytic  [%.6f, %.6f]\n  bootstrap [%.6f, %.6f]\n",
              basename(path), ana[["lo"]], ana[["hi"]],
              boot[["lo"]], boot[["hi"]]))
  cat(sprintf("  straddle  r2_incr %.3e  CI [%.9g, %.6f]\n",
              r2_incr0, ana0[["lo"]], ana0[["hi"]]))
}

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
Emit(file.path(here, "ols_reference.txt"))
Emit(file.path(here, "logistic_reference.txt"), extra_prevalence = 0.05)
