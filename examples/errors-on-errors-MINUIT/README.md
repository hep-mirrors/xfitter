# Full MINUIT errors-on-errors fit

This example configures the HERAPDF2.0, NNLO, RTOPT HERA fit with errors on
errors (EoE) on seven selected systematic sources. It performs a full MINUIT
minimisation and calculates Pumplin PDF uncertainties. All seven selected
nuisance parameters are fitted externally.

## Inputs and settings

| Files | Uncertainty treatment |
| --- | --- |
| `parameters.yaml`, `steering.txt` | Seven-source EoE, all selected sources `:E@eps=0.6` |
| `constants.yaml` | Electroweak constants and quark masses |

The fit uses the seven combined HERA inclusive datasets with Q² ≥ 3.5 GeV²:
1145 data points and 14 free PDF parameters, giving 1131 degrees of freedom.
The setup uses QCDNUM evolution, RTOPT neutral-current coefficient functions,
Q₀² = 1.9 GeV², αs(MZ) = 0.118, mc = 1.43 GeV and mb = 4.50 GeV.

The EoE sources, in their declaration order, are:

1. `sysHZComb1053`
2. `proc_tb21`
3. `sysHZComb1064`
4. `sysHZComb1114`
5. `sysHZComb1144`
6. `sysHZComb1111`
7. `sysHZComb1113`

Each has ε = 0.6. The other 162 correlated sources retain Gaussian
constraints. MINUIT registers the seven external nuisances after the PDF
parameter definitions; the `set par 23` through `set par 29` commands supply
their starting values. These are starting guesses, not fixed parameters.
Recheck the indices if the PDF definitions or source order are changed.

The fit runs `migrad`, followed by the final function evaluation, and uses:

```yaml
doErrors: |
  Pumplin
  SET ERRDEF 0.5
  ITERATE 10
```

This block belongs under `MINUIT`. Pumplin's `ITERATE 10` controls the
uncertainty calculation. It is separate from `n_iterations` in steering,
which concerns internal nuisance profiling. All seven EoE nuisances in this
example are external; the remaining internally profiled nuisances are
Gaussian. Bartlett corrections are enabled for the EoE fit.