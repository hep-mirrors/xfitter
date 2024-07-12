# EFT-reaction
Fitting SM/EFT parameters in the xFitter framework (together with PDF and alphas)

## installation
Given xFitter properly installed, the EFT reaction may be added by:

1. downloading the source code of the EFT reaction
   and adding the EFT directory under /path_to_xfitter/xfitter_master/reactions ;
2. appending a new line
   add_subdirectory(reactions/EFT)
   to /path_to_xfitter/xfitter_master/CMakeLists.txt ;
3. recompiling xFitter by executing
  cd /path_to_xfitter
  source setup.sh
  cd xfitter_master
  ./make.sh install .

## Example:
Please check examples/ttbar, which fits top quark mass to inclusive
ttbar production rates.  Here we highlight details relevant for the
`EFT` reaction.

1. Let xFitter know how the theoretical prediction depends on `mt`. 

  TermName = 'SMNNLO'
  TermSource = 'EFT'
  TermInfo =
    'ListEFTParam=deltamt:FileName=fit_mt.yaml:NoCentral=False:AbsOutput=True',
  TheorExpr = 'SMNNLO'

Here the observable (inclusive ttbar x-xsection) is calculated by the
`EFT` reaction.  The argument `ListEFTParam` tells xFitter that only
one parameter `deltamt` is involved.  Details of how `SMNNLO` depends
on `deltamt` is given in another YAML file given by the `FileName`
argument (`fit_mt.yaml` in this case).

2. The YAML file `fit_mt.yaml` provides necessary inputs such that
xFitter knows how to calculate `SMNNLO` for any top quark mass.

3. Provide the initial values of `deltamt` in `parameters.yaml`.

  Parameters:
    mtp    : [173.0, 1.0] # the top quark pole mass.
    deltamt: "=mtp-172.5" # deltamt is defined by mtp, and is fitted using the EFT reaction

Here `deltamt` is not an independent parameter, but is defined as the
difference between the top quark pole mass `mtp` (which we want to
fit) and 172.5GeV.

`deltamt` is introduced since we expect

   sigma(172.5+deltamt) = sigma(172.5) + x * deltamt + y * (deltamt)^2

to be a good approximation as long as deltamt is not very large
(e.g. abs(deltamt) < 10 GeV).  On the other hand, we will NOT expect the
inclusive ttbar x-section may be expanded as

    sigma(mt) ~ sigma(mt=0) + x * mt + y * (mt)^2

for all mt in e.g. 0GeV ~ 180GeV