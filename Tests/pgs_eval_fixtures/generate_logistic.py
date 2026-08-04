#!/usr/bin/env python3
"""statsmodels reference values for the binary path.  Run once; output is
committed.

    python3 Tests/pgs_eval_fixtures/generate_logistic.py

This overwrites the file, so the confidence-interval keys generate_ci.R adds go
with it.  Re-run `Rscript Tests/pgs_eval_fixtures/generate_ci.R` afterwards.

Two fits on the same data, matching what PgsEvalLogisticOne does:
  - Logit, for the coefficient and its standard error, plus the covariates-only
    fit's own coefficients (null.params), which pin down the null model that
    PgsEvalLogisticInit leaves in ctx->beta;
  - OLS on the same 0/1 outcome, for the observed-scale R2 family that the
    binary path now reports (Nagelkerke was removed -- Lee et al. (2012) derive
    the liability-scale transform from a linear-model R2).

Note: Firth (logistf) reference generation was removed when Firth was dropped
from the binary path (see git eeae1e2 for the Firth implementation, and the
removal commit for the rationale).
"""
import numpy as np, statsmodels.api as sm, os

OUT = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20260801)
n, n_covar = 300, 2
covars = rng.normal(size=(n, n_covar))
pgs = rng.normal(size=n)
eta = -0.7 + 0.5 * covars[:, 0] + 0.8 * pgs
y = rng.binomial(1, 1.0 / (1.0 + np.exp(-eta)))

X_null = sm.add_constant(covars)
X_full = np.column_stack([X_null, (pgs - pgs.mean()) / pgs.std(ddof=1)])
null = sm.Logit(y, X_null).fit(disp=0)
full = sm.Logit(y, X_full).fit(disp=0)

# Observed-scale R2: plain OLS of the 0/1 outcome on the same two designs.
ols_null = sm.OLS(y, X_null).fit()
ols_full = sm.OLS(y, X_full).fit()

with open(f"{OUT}/logistic_reference.txt", "w") as f:
    f.write("# statsmodels Logit + OLS reference -- generated, do not edit\n")
    # The versions the committed values came from.  A regeneration on a
    # different stack has already moved beta/se by 1 ULP; without this line
    # there is no way to tell a numerics change from a dependency change.
    f.write(f"# numpy {np.__version__}  statsmodels {sm.__version__}\n")
    f.write(f"n {n}\nn_covar {n_covar}\n")
    f.write(f"beta {full.params[-1]:.17g}\n")
    f.write(f"se {full.bse[-1]:.17g}\n")
    # Null-model coefficients, intercept first: one line each, in order.
    for b in null.params:
        f.write(f"null_beta {b:.17g}\n")
    f.write(f"r2_null {ols_null.rsquared:.17g}\n")
    f.write(f"r2_full {ols_full.rsquared:.17g}\n")
    f.write(f"r2_incr {ols_full.rsquared - ols_null.rsquared:.17g}\n")
    f.write("# data: y then covars then pgs\n")
    for i in range(n):
        f.write(" ".join(f"{v:.17g}" for v in [y[i], *covars[i], pgs[i]]) + "\n")

print(f"Written logistic_reference.txt ({n} samples, {n_covar} covariates)")
