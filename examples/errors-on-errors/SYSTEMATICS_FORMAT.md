# Errors-on-Errors HERA7 regression

This example is a fast regression test for the errors-on-errors (EoE)
implementation and its Bartlett corrections. It evaluates the HERA inclusive
DIS fit at a known parameter point; it deliberately does not run a MINUIT
minimisation.

## Regression configuration

The seven sources are the forward-selected HERA7 set used in the EoE PDF
study. Two sources are represented as external MINUIT parameters and five are
profiled internally:

```fortran
&Systematics
  ListOfSources =
    'sysHZComb1053:E@eps=0.6',
    'proc_tb21:E@eps=0.6',
    'sysHZComb1064:N@eps=0.6',
    'sysHZComb1114:N@eps=0.6',
    'sysHZComb1144:N@eps=0.6',
    'sysHZComb1111:N@eps=0.6',
    'sysHZComb1113:N@eps=0.6'

  n_iterations = 4
  Enable_Bartlett = .true.
&End
```

All seven sources use `epsilon = 0.6`. The five `:N` sources receive exactly
four Newton updates during the FCN evaluation. Four updates are a bounded,
near-converged regression setting; production configurations with profiled EoE
sources in the study used eight updates.

The active MINUIT command stream is:

```text
set par 23 -16.40759961
set par 24  -6.47988938
call fcn 3
```

For this parameter and source ordering, indices 23 and 24 are
`sysHZComb1053` and `proc_tb21`, respectively. Recheck the MINUIT parameter
definitions if the example configuration is changed.

The PDF parameters and external systematics are initialised at the completed
N=4 solution. The second fields of the free PDF parameters are the fitted
errors from that solution. In this evaluation-only example they provide the
finite-difference scales used by the Bartlett calculation; they are not errors
estimated by this run.

The test therefore covers mixed external/profiled EoE handling, four internal
profiling updates, Bartlett factors, corrected nuisance-parameter errors and
output. It does not test convergence of the outer PDF fit.

## Running a full fit manually

For a slower physics validation, replace the evaluation-only command stream
with:

```text
set par 23 -16.40759961
set par 24  -6.47988938
call fcn 1
set str 2
migrad
call fcn 3
```

Pumplin error-band generation is intentionally not part of this regression.

## Source syntax

Each configured source uses:

```text
'source-name:form[:scaling][:type]@eps=value'
```

- `:N` selects internal nuisance-parameter profiling.
- `:E` selects an external MINUIT parameter.
- Optional scaling modifiers are `:M`, `:A` and `:P`.
- Optional source-type modifiers are `:D` and `:T`.
- `epsilon` must be non-negative.

EoE is supported for `:N` and `:E` sources. An `@eps=` value attached to a
covariance (`:C`) or offset (`:O`) source is ignored with a warning.

Setting `n_iterations` explicitly requests exactly that many profiling
updates. If it is omitted, profiling instead runs until the configured
tolerance is reached, subject to `max_iterations`; the built-in defaults are
`tolerance = 1e-8` and `max_iterations = 50`.

`Enable_Bartlett = .true.` enables the goodness-of-fit and confidence-interval
Bartlett factors reported in `Results.txt`.
