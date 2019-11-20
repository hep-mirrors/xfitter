This test demonstrate how to use the Hathor reaction for charm, bottom and top quark hadroproduction,
with different heavy-quark masses, PDFs and scale settings.
In addition, it tests specifying PDF evolution per reaction and demonstrates how to combine several
reactions in one expression using KFactor.

There are two data files:

1) ttbar.dat calculates chi2 for total top quark pair production cross sections at sqrt(s)=13 TeV 
for different settings from CMS publication Eur.Phys.J. C79 (2019) no.5, 368. The idea is that 
all these different settings should give chi2 ~ 0 since they correspond to the values of top quark 
mass extracted using different PDFs/scales (see Table 6 fom the paper)

2) ccbar-bbbar.dat calculates chi2 for charm and bottom quark pair production cross sections at 
sqrt(s)=13 TeV (for dummy data points)
