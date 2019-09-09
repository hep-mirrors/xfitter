One chi2 iteration for ALL data using NLO RTOPT LHAPDF=NNPDF30_nlo_as_0118
it tests also FlipCharge and FlipUD evolutions
it tests also storing PDF from LHAPDF evolution in LHAPDF6 format

This tests aims to enable all data available in xFitter datafiles.
As of 2.07.2019 only these datasets are not covered by this test:

datafiles/hera/zeus/diffractiveDis/0812.2003: not supported in xfitter-2.2, use xfitter-2.0.1
datafiles/hera/h1/jets/0904.3870: not supported in xfitter-2.1 and xfitter-2.2 (normalised jet cross sections), use xfitter-2.0.1
datafiles/hera/h1/jets/1406.4709: not supported in xfitter-2.1 and xfitter-2.2 (normalised jet cross sections), use xfitter-2.0.1
datafiles/lhc/atlas/topProduction/1407.0371: not supported in xfitter-2.1 and xfitter-2.2 (DiffTop & fastNLO), use xfitter-2.0.1
datafiles/lhc/cms/topProduction/1211.2220: not supported in xfitter-2.1 and xfitter-2.2 (DiffTop & fastNLO), use xfitter-2.0.1
LHeC: pseudodata covered in dedicated test (LHeC)
AFB: pseudodata covered in dedicated test (AFB)

Some other data sets which are superseeded are explicitly listed and commented out in steering.txt.
