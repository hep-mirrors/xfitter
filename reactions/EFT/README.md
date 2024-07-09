# EFT-reaction
Fitting SM/EFT parameters in the xFitter framework

## installation

## Example: fit top quark mass with xFitter + EFT reaction

1. Let xFitter know how the theoretical prediction depends on `mt`. 

TermName = 'Kmt'
TermSource = 'EFT'
TermInfo = 'ListEFTParam=deltamt:FileName=/path/to/EFT_file.yaml'

TheorExpr  = 'theory1725/*Kmt'

Here we include the top quark mass dependence via a K factor `Kmt`.
The argument `ListEFTParam` tells xFitter that only one parameter `deltamt` is involved. 
Details of how `Kmt` depends on `deltamt` is given in another YAML file given by the `FileName` argument.

2. The YAML file `/path/to/EFT_file.yaml` provides necessary inputs such that xFitter knows how to
calculate `Kmt` for any top quark mass.
This file should be provided by the user.
An example is provided at `examples/fit_mt.yaml`

3. Provide the initial values of delta_mt in `parameters.yaml`.


