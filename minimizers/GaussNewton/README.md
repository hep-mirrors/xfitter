# GaussNewton minimizer

A Levenberg–Marquardt (damped Gauss–Newton) minimizer for xFitter. It minimises the
official xFitter χ² with the nuisance shifts profiled, using a Jacobian of the theory
predictions. It is selected with

```yaml
Minimizer: GaussNewton
```

and nothing in this directory runs otherwise.

**The minimizer never decides what the χ² is.** Every step is accepted or rejected with
the χ² returned by `fcn` (plus the Gaussian parameter priors, which `fcn` does not include
for non-MINUIT minimizers). The module only computes a search direction and the
covariance. A bug in the derivative model can slow the convergence or bias the
covariance; it cannot move the χ² minimum, and `checkDerivatives: 1` tests the model
against finite differences of `fcn`.

## 1. Steering

```yaml
GaussNewton:
  maxIterations: 50
  tolerance: 1.e-4            # stop when the Newton decrement, Eq. (11), is below this (chi2 units)
  lambda0: 1.e-3              # initial damping
  jacobian: FiniteDifference  # the only one in this version
  fdScheme: forward           # forward | central
  fdStep: 1.e-3               # h_j = fdStep * max(|p_j|, YAML step_j)
  threads: 0                  # forked processes for the Jacobian columns and the observed Hessian (0: none)
  covariance: 1               # C = 2 H^-1, Eq. (10)
  hessian: observed           # observed (second differences of the chi2) | gaussnewton
  hessianScheme: central      # central: O(h^2), P(P+1) chi2 calls | forward: O(h), P(P+3)/2 calls, as HESSE
  hessianDelta: 1.e-3         # chi2 change targeted by the second-difference steps
  doErrors: Hesse             # Hesse | None: symmetric PDF error bands from C
  checkDerivatives: 0         # 1: compare the model with finite differences of fcn at the start
```

All keys are optional. Unknown keys give a warning.

Outputs, in the xFitter output directory: `parsout_0`, `pars.yaml` and
`gaussnewton.out.txt` (summary, covariance, correlations), as for CERES, plus the usual
files of the final `fcn` call. The iteration log goes to stdout.

## 2. Code layout

| Directory | Job | xFitter headers |
|---|---|---|
| `core/` | numerics: LM loop, point model, Schur complement | none (unit-tested in `test/`) |
| `bridge/` | the only reads of xFitter state; `fcn` and theory-only evaluations | yes |
| `jacobian/` | J_t = ∂t/∂p | through `bridge/` |
| `util/` | fork pool, timer, output files | yes |
| top level | `GaussNewtonMinimizer` (glue), `ProfiledChi2` (the fit as a `core::Problem`), `DerivativeCheck` | yes |

`bridge/xfitter_gn_bridge.f` reads the commons through the real include files, so the
common-block layouts cannot drift. It writes nothing.

One iteration:

```
p --fcn--> official chi2; snapshot of d, t, sigma, Gamma, b, ...           (bridge)
  --FiniteDifferenceJacobian--> J_t                                         (jacobian)
  --PointModel--> r, dr/dt, dr/db, log term                                 (core)
  --ProfiledNuisances--> g_eff, H_eff                                       (core)
  --LevenbergMarquardt--> trial step, judged with the official chi2 (fcn)   (core)
```

The snapshot is taken right after `fcn`, because the Jacobian evaluations overwrite the
theory.

Unit tests of `core/`: configure with `-DXFITTER_GN_TESTS=ON` and run
`build/minimizers/GaussNewton/test/gn_core_tests`.

## 3. Mathematics

Notation: p are the minimised parameters, t(p) the theory of the fitted points i
(FitSample, diagonal errors), k the nuisance-parameter sources, b their shifts, π_k the
prior scale (`SysPriorScale`). χ² is in χ² units (not NLL = χ²/2).

**(1) Objective**

χ²(p) = min_b [ Σ_i r_i(t, b)² + Σ_k π_k b_k² + Σ_i L_i(t) ] + Σ_j ((p_j − μ_j)/s_j)²

The last term is the Gaussian parameter priors from `parameters.yaml`.

**(2) Systematic scaling.** Γ_ki = β_ki m_k(t_i) with m = d_i (NoRescale, or forced
additive point), t_i (Linear) or √(t_i d_i) (Poisson). So ∂Γ_ki/∂t_i = c_k Γ_ki / t_i with
c_k = 0, 1, ½.

**(3) Errors.** σ_i² = A_i |d_i t_i| + B_i t_i² + const, with A_i = e_stat,poisson² +
e_uncor,poisson² and B_i = e_uncor,mult² (relative errors). ∂σ_i²/∂t_i = A_i |d_i| + 2 B_i t_i,
and 0 where xFitter replaces t by d (t ≤ 0). The value of σ_i is taken from xFitter.

**(4) Residual.** r_i = (d_i − t_i + Σ_k Γ_ki b_k) / σ_i.

**(5) Derivatives.**
a_i ≡ ∂r_i/∂t_i = (−1 + Σ_k c_k Γ_ki b_k / t_i)/σ_i − r_i (∂σ_i²/∂t_i)/(2σ_i²),
∂r_i/∂b_k = Γ_ki/σ_i.

**(6) Log term** (`PoissonCorr`). L_i = log(σ_i²/α_i²). L′_i = (∂σ_i²/∂t_i)/σ_i², and
L″_i = (∂²σ_i²/∂t_i²)/σ_i² − (L′_i)². Only max(L″_i, 0) enters H (design decision D2).

**(7) Jacobian on (p, b).** J_p = diag(a) J_t, J_b = (Γ/σ)ᵀ, prior rows √π_k for b.

**(8) Normal equations.**
g_p = 2 J_pᵀ r + J_tᵀ L′,  g_b = 2 (Γ/σ) r + 2 π b,
H_pp = 2 J_pᵀ J_p + J_tᵀ diag(L″₊) J_t,  H_pb = 2 J_pᵀ (Γ/σ)ᵀ,  H_bb = 2 ((Γ/σ)(Γ/σ)ᵀ + diag π).

**(9) Profiling (Schur complement).**
H_eff = H_pp − H_pb H_bb⁻¹ H_bp,  g_eff = g_p − H_pb H_bb⁻¹ g_b.
xFitter's b is the exact profile for the current t, so g_b ≈ 0 and g_eff is the exact
gradient of the profiled χ² (envelope theorem). The log reports max |H_bb⁻¹ g_b|.
The parameter priors add 2(p_j − μ_j)/s_j² to g and 2/s_j² to H.

**(10) Covariance.** C = 2 H⁻¹ (Δχ² = 1). By default H is the observed Hessian of the
objective (1) at the minimum, from second differences of the official χ²: central
diagonal and off-diagonal terms, f(p ± h_j e_j) and f(p ± (h_j e_j + h_k e_k)), with
steps for Δχ² ≈ 10⁻³ set from H_eff and refined once; P(P+1) χ² calls. The forward
off-diagonal scheme of MINUIT HESSE (P(P+3)/2 calls) is available. With `hessian: gaussnewton`, H = H_eff
of Eq. (9) at no extra cost (design decision D10). If H is not positive definite, the
pseudo-inverse over the positive eigenvalues is used, with a warning.

**(11) LM step** (hylu `fitModel`).
Solve (H + λ D) dx = −g with D = diag(max(|H_ii|, 10⁻⁸ max_j |H_jj|)). Predicted decrease pred = −(g·dx + ½ dxᵀ H dx).
ρ = (χ²(p) − χ²(p + dx)) / pred with both χ² from `fcn`. Accept if ρ > 10⁻⁴; then λ/3 if
ρ > 0.75 and 2λ if ρ < 0.25. On rejection 4λ; at most 20 trials per iteration.
Stop when the Newton decrement ½ gᵀ H⁻¹ g < `tolerance`.

**(12) Bounds.** Projected steps: x + dx is clipped to the box, and the predicted decrease
uses the clipped step. A parameter at a bound whose gradient points out of the box is held
fixed (active set). There is no variable transformation, so C keeps its meaning. A
parameter left at a bound is reported.

**(13) Supported configurations.**

| Feature | |
|---|---|
| Diagonal points, Nuisance sources, NoRescale/Linear/Poisson scaling, PoissonCorr, parameter priors, bounds, `ControlFitSplit`, `CovarToNuisance` | supported |
| Covariance-matrix points (correlation files, `:C` sources, `CorChi2Type: Matrix`) | stops: use `CovarToNuisance` |
| External sources (fit parameters entering only through the shifts) | stops |
| Offset sources, `ExtraSystRescale`, `FirstIterationRescale`, asymmetric errors | warning: the derivative model is approximate; the χ² stays exact |
| Polarisation shifts | stops |

## 4. Design decisions

- **D1. The χ² is never recomputed by the module for step acceptance.** Only `fcn` decides.
  The module χ² from the snapshot is compared with `fcn` at every model evaluation
  (`checks:` line of the log), which tests the bridge.
- **D2. Log-term curvature.** L″ can be negative; a negative contribution could make H
  indefinite. Only the non-negative part is kept, as hylu does. The gradient is exact.
- **D3. Finite-difference step.** h_j = `fdStep`·max(|p_j|, YAML step), one-sided away from
  a bound. Forward differences need one theory evaluation per parameter, central two.
- **D4. Theory-only Jacobian.** The columns use `update_theory_iteration` (no χ², no
  nuisance solve). The derivatives of σ, Γ and the profiling come from Eqs. (2)–(9).
- **D5. Parameter priors.** `fcn` does not include them for non-MINUIT minimizers (the
  MINUIT extra-parameter list is reset at start, as in CERES); the module adds them.
- **D6. Stopping rule.** The Newton decrement does not depend on the damping λ, unlike the
  predicted decrease of a damped step.
- **D7. Damping floor.** Marquardt damping λ diag|H| does not restrain a parameter with
  almost zero curvature; the floor 10⁻⁸ max|H_jj| does.
- **D8. Non-finite derivatives.** If a perturbed point of the Jacobian gives a non-finite
  theory (e.g. an invalid PDF), those entries are set to zero with a warning, and a
  non-finite step counts as a rejected trial. `fcn` returns 10¹⁰ for a NaN χ², so such a
  step is never accepted.
- **D9. History check.** At the start the theory is recomputed at the same p; a difference
  gives a warning (e.g. an `"=expr"` parameter that uses a sum-rule normalisation reads the
  previous iteration's value).
- **D10. Observed Hessian for the covariance.** The Gauss–Newton Hessian drops the term
  2 Σ r ∂²r. On `examples/CERES-fit` this term dominates along the least constrained
  directions (nearly degenerate negative-gluon parameters): along the three weakest
  eigenvectors of H_eff the second difference of the χ² is 4–48 times the model one, and
  the Gauss–Newton errors of Agp, Bg, Bgp are 6–8 times the observed ones. The steps are
  still judged with the official χ², so the minimum is not affected; only the covariance
  needs the observed Hessian. Central off-diagonal differences give the same errors to
  4 digits for steps targeting Δχ² = 10⁻², 10⁻³ and 10⁻⁴; forward differences (as HESSE)
  drift by a few percent with the step, and three HESSE runs of the same MINUIT job
  differ by up to 12% for these parameters. The χ² calls are independent. They are made in
  three batches (diagonal, refined diagonal steps, off-diagonal), which run in `threads`
  forked processes (`fcn` writes no file for flag 2). On `examples/GaussNewton-fit` with
  `threads: 2`, the fit takes 410 s instead of 776 s. The errors differ by at most
  4·10⁻⁵ relative: in a child, `fcn` starts from a different, history-dependent state at
  round-off level. With `threads: 0` the output is byte-identical to before.

## 5. Validation

Configuration: `examples/CERES-fit` (HERA I+II inclusive, NNLO, RT, QCDNUM, 14 parameters,
169 nuisance sources, PoissonCorr). All three minimisers start from the same point, displaced
from the minimum along the error ellipse (χ² = 1397.32, Δχ² = 34). Laptop, one fit at a
time, about 3 s per χ² evaluation.

| | iterations | fcn calls | theory evaluations | wall time | χ² | parameters vs MINUIT | errors / HESSE |
|---|---|---|---|---|---|---|---|
| GaussNewton: minimisation | 6 | 9 | 98 | 313 s | 1363.325 | ≤ 0.03σ | |
| GaussNewton: + observed Hessian | | + 210 | | 927 s | | | 0.983 – 1.013 |
| CERES | 10 | 245 | | 719 s | 1363.33 | ≤ 0.03σ | 0.74 – 3.9 |
| MINUIT MIGRAD | | 1150 | | 3307 s | 1363.326 | reference | |
| MINUIT MIGRAD + HESSE | | 1298 | | ≈ 3750 s | | | 1 |

- σ are the HESSE errors. Two independent MINUIT fits (from this start and from the
  minimum) differ by up to 0.028σ, which is the precision of the comparison.
- GaussNewton reaches the minimum 10× faster than MIGRAD, and gives HESSE-quality errors
  4× faster than MIGRAD + HESSE.
- With `hessian: gaussnewton` the errors would be 1.0–5.3 times HESSE (D10); CERES has the
  same problem, made worse by a factor-2 convention in its covariance.

**Built-in checks** (`checkDerivatives: 1` at the minimum): the model gradient agrees with
central finite differences of `fcn` within 3·10⁻³ (forward-difference Jacobian); the module
χ² equals `fcn` within 10⁻¹² and max |H_bb⁻¹ g_b| ≈ 10⁻¹³ at every model evaluation.

Unit tests (`gn_core_tests`): LM on linear least squares, Rosenbrock and bounded problems;
point-model derivatives, profiled gradient and Schur complement against finite differences
of a toy xFitter χ² with exact profiling; observed Hessian against analytic Hessians.

Not supported in this version (stops with a message): covariance-matrix points (use
`CovarToNuisance`). Gains in this version come from fewer iterations; the cost per
iteration (1 + P theory evaluations) is reduced in the next merge requests (Broyden
updates, factorised Jacobian).
