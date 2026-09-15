# Using CERES with errors on errors

A ready-to-run fit and dedicated validation are available in
[`errors-on-errors-CERES`](../errors-on-errors-CERES/README.md).

Replace the example's `Minimizer: MINUIT` and `MINUIT` command block with:

```yaml
Minimizer: CERES
CERES:
  offset: 10
  tolerance: 1e-7
  strategy: 0
  covariance: 1
  doErrors: Hesse
```

Keep the source-level `:N`/`:E` and `@eps=` annotations in `steering.txt`.
External nuisance parameters are registered automatically; do not also add
these sources to the YAML `Parameters` section. The example's MINUIT commands
set their initial values explicitly; those commands do not apply to CERES,
which starts the automatically registered external sources at zero.

CERES performs a minimization, whereas the supplied MINUIT regression only
evaluates a fixed parameter point. Thus the resulting minima need not match
the example's reference chi-square. Check convergence, and increase the
profiling iteration count as needed (the paper reports eight profiling iterations in its
production configuration).

The exported nuisance residuals contain the logarithmic EoE penalties.
Bartlett preparation uses the active minimizer's free parameter names and
uncertainties. It excludes external nuisances from the PDF-parameter count,
computes theory derivatives, and restores the fitted parameter values.

`doErrors: Hesse` requires `covariance: 1`. The symmetric PDF bands use the
CERES covariance and the common Bartlett CI scale. `parsout_0` uses the same
corrected-uncertainty convention as MINUIT: the CI scale for PDF parameters,
and the conditional quadratic nuisance uncertainty for external sources.
CERES's `ceres.out.txt` covariance and `pars.yaml` errors retain the raw CERES
covariance convention. Bartlett factors remain local quadratic, order-epsilon-
squared approximations.
