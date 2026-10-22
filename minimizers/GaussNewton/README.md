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
  jacobian: FiniteDifference  # FiniteDifference | Factorised (J_t = K F, Eq. 14)
  fdScheme: forward           # forward | central
  fdStep: 1.e-3               # h_j = fdStep * max(|p_j|, YAML step_j)
  threads: 0                  # forked processes for the Jacobian columns and the observed Hessian (0: none)
  spanTolerance: 1.e-3        # Factorised: K F_j used if F_j is in the explored span to this, Eq. (17)
  covariance: 1               # C = 2 H^-1, Eq. (10)
  hessian: observed           # observed (second differences of the chi2) | gaussnewton
  hessianScheme: central      # central: O(h^2), P(P+1) chi2 calls | forward: O(h), P(P+3)/2 calls, as HESSE
  hessianDelta: 1.e-3         # chi2 change targeted by the second-difference steps
  doErrors: Hesse             # Hesse: symmetric bands from C | Pumplin: iterative Hessian, Eq. (18),
                              #   and 2P asymmetric eigenvector members, as MINUIT ITERATE | None
  pumplinIterations: 10       # Pumplin: at most this many iterations
  pumplinDelta: 1.            # Pumplin: chi2 rise of the steps along the eigenvectors (MINUIT UP)
  pumplinTolerance: 0.02      # Pumplin: stop when max |H_y/2 - 1| < this, Eq. (18)
  checkDerivatives: 0         # 1: compare the model with finite differences of fcn at the start
```

All keys are optional. Unknown keys give a warning.

Outputs, in the xFitter output directory: `parsout_0`, `pars.yaml` and
`gaussnewton.out.txt` (summary, covariance, correlations), as for CERES, plus the usual
files of the final `fcn` call. The iteration log goes to stdout.

## 2. Code layout

| Directory | Job | xFitter headers |
|---|---|---|
| `core/` | numerics: LM loop, point model, Schur complement, linear response K | none (unit-tested in `test/`) |
| `bridge/` | the only reads of xFitter state; `fcn`, theory-only evaluations, input PDFs f0(p) | yes |
| `jacobian/` | J_t = ∂t/∂p | through `bridge/` |
| `util/` | fork pool, timer, output files | yes |
| top level | `GaussNewtonMinimizer` (glue), `ProfiledChi2` (the fit as a `core::Problem`), `DerivativeCheck` | yes |

`bridge/xfitter_gn_bridge.f` reads the commons through the real include files, so the
common-block layouts cannot drift. It writes nothing.

One iteration:

```
p --fcn--> official chi2; snapshot of d, t, sigma, Gamma, b, ...           (bridge)
  --FiniteDifferenceJacobian or FactorisedJacobian--> J_t                  (jacobian)
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
| `ExtraSystRescale` | supported: exact gradient, Eq. (13b); Gauss–Newton Hessian |
| Offset sources, `FirstIterationRescale`, asymmetric errors | warning: the derivative model is approximate; the χ² stays exact |
| Polarisation shifts | stops |

**(13b) `ExtraSystRescale`.** xFitter profiles the shifts twice: b¹ with the nominal errors,
then b with the Poisson stat part of σ² rescaled, S_i → S_i E_i, E_i = exp(−Σ_k β_ki b¹_k)
over the linear nuisance sources. The χ² then depends on t also through b¹(t):

dχ²/dt_i = γ_i + Σ_k w_k ∂b¹_k/∂t_i,  w_k = −Σ_j (∂Φ/∂σ_j²) S_j E_j β_kj,
∂b¹/∂t = −Φ¹_bb⁻¹ Φ¹_bt,  ∂Φ/∂σ_j² = (L_j-flag − r_j²)/σ_j²,

with Φ¹ the first-pass objective (σ¹, b¹) and γ of Eq. (8) computed with ∂σ²/∂t including E.
The module recomputes b¹ from the first-pass errors (`GetPointErrors`) and checks that the
rescaled σ equals xFitter's (`checks:` line of the log). The Hessian stays Gauss–Newton.

**(14) Factorised Jacobian** (`jacobian: Factorised`). At fixed α_s, masses, Q0 and data,
the theory of DIS data is linear in the input PDFs f0 = xf(x, Q0) (evolution and
coefficient functions are linear). Then

J_t = K F,  K = ∂t/∂f0,  F = ∂f0/∂p.

f0(p) is sampled on fixed x nodes for each decomposition and flavour (`bridge/InputPdfs`):
only the dependent parameters, parameterisations and decompositions (sum rules) are
recomputed, no evolution and no theory. F costs 2P such samples (central differences),
which is negligible against one theory evaluation. K is never formed as a matrix.

**(15) Exact pairs.** For two evaluated points a, b with the same non-PDF parameters,
w = t(b) − t(a) = K u with u = f0(b) − f0(a), exactly for any distance between the points.
Pairs come from consecutive χ² evaluations of the fit (accepted and rejected trials) and
from the finite-difference columns. They are stored as an orthonormal basis Q of the
explored directions and Z = K Q (`core/LinearResponse`): with a = Qᵀu and r = u − Q a
(Gram–Schmidt applied twice), a pair adds q = r/|r| and K q = (w − Z a)/|r| if
|r| > 10⁻⁴ |u|. The threshold keeps nearly dependent pairs (e.g. repeated trials along one
direction) out: they would divide the part of w that the nodes do not describe by a
small |r|.

**(16) Response.** K F_j ≈ Z Qᵀ F_j. The error is K (F_j − Q Qᵀ F_j).

**(17) Column rule.** Column j is Z Qᵀ F_j if |F_j − Q Qᵀ F_j| ≤ `spanTolerance` |F_j|,
otherwise a finite-difference column, which also adds the pair of its two points. Columns
of parameters that do not change f0 (α_s, masses, normalisations: F_j = 0) are always
finite differences, and Q, Z are dropped when such a parameter moves (K depends on it).
This is exact, but with a free α_s it removes the gain.
**(18) Pumplin errors** (`doErrors: Pumplin`, hep-ph/0008191, as MINUIT `ITERATE`). With
the eigen-decomposition H = Σ ε_j v_j v_jᵀ, the displacements U_j = v_j √(2/ε_j) raise the χ²
by 1 in the quadratic approximation, and in the coordinates y of x = x̂ + U y the Hessian is
2·1. An iteration measures it, H_y: as MINUIT ITERATE, the step along each eigenvector is
searched (up to 6 batches) so that the average rise of the χ² at ±step is within
[½, 2]·`pumplinDelta`; a non-finite or very large χ² (invalid PDF) shortens the step, a flat
or concave direction lengthens it, and a direction without a usable rise keeps its
curvature. The off-diagonal terms use these steps (central differences); a pair with a χ²
more than 100·`pumplinDelta` above the minimum is not used. Then H ← U⁻ᵀ H_y U⁻¹ (with the
eigenvalues of H_y kept ≥ 0.02, so the errors grow at most 10× per iteration), until
max |H_y/2 − 1| < `pumplinTolerance`; otherwise the most self-consistent iteration is used.
The start is the Gauss–Newton Hessian; the first iteration then replaces the observed
Hessian. C = 2 H⁻¹ is written as the covariance, and the 2P members x̂ ± U_j of the final H
give the same files as MINUIT (`pdfs_q2val_sNN{m,p}`, LHAPDF `hessian` set with 2P + 1
members, asymmetric theory
uncertainties). The χ² calls of an iteration are independent and run in `threads`
processes.

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

- **D11. Explored-subspace K instead of nodal probes.** The probe approach (K on a basis
  of input functions, one theory evaluation per x node and flavour) needs 7 × 199 ≈ 1400
  theory evaluations for the QCDNUM grid of `examples/CERES-fit`, 14 times a whole fit
  with finite differences (98). The fit only needs K on the
  directions F_j, which change slowly from one iteration to the next. The pairs of
  Eq. (15) are exact and come free with every χ² evaluation, so after the first
  iteration most columns need no theory evaluation. No probe decomposition, no change of
  the steering outside the `GaussNewton` block.
- **D12. Linearity check and exact confirmation.** The first time a column comes from K F
  (and again after each reset), the column with the largest σ-weighted norm is also
  computed by finite differences (one theory evaluation, two with `fdScheme: central`);
  the fit uses the finite-difference column. If they differ by more than 2% (σ-weighted
  norm), the theory is not linear in f0 (hadronic data, ratios in `TheorEval`
  expressions, ...): a warning is printed and the fit continues with finite differences
  only. A Jacobian error does not change the χ² values (D1) but moves the point where the
  fit stops, because the gradient is K F-based. Therefore, when the Newton decrement falls
  below `tolerance`, or an accepted step gains less than it, while the Jacobian contains
  K F columns, the loop takes a finite-difference Jacobian at the same point and tests
  again; from then on the fit uses finite differences. Every minimum is confirmed with
  exact derivatives, for one or two finite-difference Jacobians.
- **D13. Input-PDF nodes and metric.** 100 nodes uniform in log x on [10⁻⁶, 0.1) and 50
  uniform in x on [0.1, 0.98]; the flavours −6..6 and the photon of every decomposition.
  Each value is divided by max(|xf|, 10⁻³ max_x |xf|) of the first point, so the span test
  measures relative changes of the PDFs, as the data do. With plain xf values the test is
  dominated by the smallest x when the gluon is steep: on the N3LO HERA configuration
  (Bg ≈ −0.6, xg ≈ 4000 at x = 10⁻⁶) columns in the span to 10⁻⁵ deviated by 5–10% from
  finite differences, and the fit stopped 0.47σ from the CERES minimum. The theory is
  linear there: along parameters whose f0 path is two-dimensional, a third point is
  predicted from two pairs to 10⁻¹³. With the relative metric the span residual tracks the
  column error (ratio 1–10), and the fit reaches the CERES minimum. The nodes only measure
  whether F_j is in the explored span; they need not coincide with the evolution grid
  (denser nodes up to 1 − x = 10⁻⁴ changed nothing), because the pairs (u, w) relate whole
  functions. DIS theory at x depends on f0 only at x′ ≥ x, so nodes below the smallest x of
  the data would only cost extra finite-difference columns.

## 5. Validation

Configuration: `examples/CERES-fit` (HERA I+II inclusive, NNLO, RT, QCDNUM, 14 parameters,
169 nuisance sources, PoissonCorr). All three minimisers start from the same point, displaced
from the minimum along the error ellipse (χ² = 1397.32, Δχ² = 34). Laptop, one fit at a
time, about 3 s per χ² evaluation.

| | iterations | fcn calls | theory evaluations | wall time | χ² | parameters vs MINUIT | errors / HESSE |
|---|---|---|---|---|---|---|---|
| GaussNewton: minimisation | 6 | 9 | 98 | 313 s | 1363.325 | ≤ 0.03σ | |
| GaussNewton: + observed Hessian | | + 210 | | 927 s | | | 0.983 – 1.013 |
| GaussNewton, `jacobian: Factorised`: minimisation | 6 | 10 | 55 | 201 s | 1363.325 | ≤ 0.03σ | |
| CERES | 10 | 245 | | 719 s | 1363.33 | ≤ 0.03σ | 0.74 – 3.9 |
| MINUIT MIGRAD | | 1150 | | 3307 s | 1363.326 | reference | |
| MINUIT MIGRAD + HESSE | | 1298 | | ≈ 3750 s | | | 1 |

- σ are the HESSE errors. Two independent MINUIT fits (from this start and from the
  minimum) differ by up to 0.028σ, which is the precision of the comparison.
- GaussNewton reaches the minimum 10× faster than MIGRAD, and gives HESSE-quality errors
  4× faster than MIGRAD + HESSE.
- With `hessian: gaussnewton` the errors would be 1.0–5.3 times HESSE (D10); CERES has the
  same problem, made worse by a factor-2 convention in its covariance.
- **Factorised Jacobian** (Eqs. 14–17), with the exact confirmation of D12 (minima equal
  to the finite-difference ones within 0.002σ):

  | configuration, start | Jacobian | iterations | theory evaluations | time to minimum |
  |---|---|---|---|---|
  | `CERES-fit`, ellipse (χ² = 1397) | finite differences | 6 | 98 | 313 s |
  | | factorised | 6 | 55 (28 for the exact confirmation) | 201 s |
  | `CERES-fit`, all parameters 1σ off (χ² ≈ 3.9·10⁴) | finite differences | 8 | 126 | 398 s |
  | | factorised | 12 | 115 | 390 s |
  | N3LO HERA (APFELxx, `ExtraSystRescale`), 1σ off | finite differences | 5 | 84 | 230 s |
  | | factorised | 6 | 54 | 170 s |

  The gain is 1.0–1.6×: after the first iteration almost every column comes from K F, but
  the confirmation of the minimum costs one or two finite-difference Jacobians, and from a
  far start the last iterations crawl in the flat negative-gluon directions with either
  Jacobian. Without the confirmation the fits were 3× faster, but on the N3LO
  configuration the minimum was off by up to 0.47σ (D13).

**Built-in checks** (`checkDerivatives: 1` at the minimum): the model gradient agrees with
central finite differences of `fcn` within 3·10⁻³ (forward-difference Jacobian); the module
χ² equals `fcn` within 10⁻¹² and max |H_bb⁻¹ g_b| ≈ 10⁻¹³ at every model evaluation.

Unit tests (`gn_core_tests`): LM on linear least squares, Rosenbrock and bounded problems;
point-model derivatives, profiled gradient and Schur complement against finite differences
of a toy xFitter χ² with exact profiling; observed Hessian against analytic Hessians;
linear response: exact in the span, error bounded by the span residual, toy J = K F.

Not supported in this version (stops with a message): covariance-matrix points (use
`CovarToNuisance`). With `jacobian: FiniteDifference` an iteration costs 1 + P theory
evaluations; `jacobian: Factorised` saves up to 40% of them for DIS data (Sec. 5).
