# CERES errors-on-errors fit

This example fits the seven-source HERA configuration with CERES, computes
Bartlett corrections, and produces symmetric PDF bands at Q2 = 4 GeV2.
It has 1145 data points, 14 free PDF parameters, two external nuisances,
and five iteratively profiled EoE nuisances. Each selected source has
`epsilon = 0.6`; profiling uses eight updates.

The inputs derive from `../errors-on-errors/`. Keep their data selection,
PDF parameterization, and selected systematic sources synchronized when
updating either example. The deliberate differences are the minimizer,
eight instead of four profiling updates, and PDF-band output.

## Running locally

Build and install xFitter with CERES and QCDNUM available, source the dependency
setup, and run from the installation/repository root:

```bash
source setup.sh
./tools/test.sh errors-on-errors-CERES
```

The test runs in `temp/errors-on-errors-CERES/`. Its numerical checks are
reported in `test.log` and `validation.xml`; the fit log is `xfitter.log`,
and generated tables are in `output/`. Python 3 is required for validation,
with no additional Python packages.

The example uses two derivative workers. Change `CERES.threads` to match
available resources; zero runs the same forward derivatives serially.
Unlike the MINUIT example, this performs a full fit. External nuisances are
registered automatically from steering and start at zero; do not also add
them to `Parameters` in YAML.

This example is in `omitTests` in `tools/test.sh`, so an ordinary run does
not execute it. Explicitly naming it overrides that omission.

## Validation

`validate.py` checks:

- successful CERES convergence and available Bartlett corrections;
- raw/corrected chi-square, 1131 degrees of freedom, and 14 PDF parameters;
- exactly the seven expected EoE sources and their internal/external treatment;
- finite correction factors consistent with summed source contributions;
- restored fitted parameter values and correctly scaled parameter errors;
- a finite covariance and all 16 finite symmetric PDF-band tables.

`expected.json` records the numerical reference and absolute tolerances:
raw chi-square 1312.79 +/- 0.10, corrected chi-square 1305.90 +/- 0.20,
GoF factor 0.994747 +/- 0.0001, and CI factor 1.003563 +/- 0.0005.
The tolerances allow small solver/library differences while rejecting
missing corrections or a substantially different minimum. Review these
numbers if the physics inputs change; `--copy` is intentionally unsupported.

There are no checked-in generated output tables. Validation does not require
specific eigenvector signs or exact nuisance-table snapshots. Recheck an
existing run without fitting again using:

```bash
python3 examples/errors-on-errors-CERES/validate.py temp/errors-on-errors-CERES
```

The fast minimizer/Bartlett regression can be run independently:

```bash
./tools/test-bartlett-minimizer.sh
```

See [CERES uncertainty conventions](../errors-on-errors/CERES.md) for the
difference between corrected `parsout_0` uncertainties and raw CERES output.

## GitLab CI

`job-ceres-eoe` in `.gitlab-ci.yml` runs automatically in scheduled pipelines
and is an optional manual job in merge-request pipelines. Trigger it for
CERES/EoE changes. A failed scheduled run fails the pipeline; the manual MR
job is non-blocking. The job uses the existing CERES-capable image and checks
that the CERES module was built, publishes JUnit results, and retains logs
and output tables for one week even on failure.

The fast Bartlett regression also runs in `job-install-full`, providing
routine coverage without requiring the full CERES fit. The CI configuration
does not create a GitLab schedule: it participates in schedules configured
for this project and branch.
