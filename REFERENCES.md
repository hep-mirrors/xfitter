If you use the xFitter package in a scientific publication, please 
consider adding the following references.  The main citations list contains the papers 
which should be cited for any use of the xFitter program. In addition, some citations 
are required depending on the modules, data and theory tables used in the program.


#  Main citation

1) "HERAFitter, Open Source QCD Fit Project"
By S. Alekhin at al., DESY Report 14-188, 7 Nov 2014, Published in EPJC (2015), 75: 304,
e-Print: arXiv:1410.4412 [hep-ex]

Also add a reference to the xFitter web portal:  [xfitter.org](xfitter.org)


# Citations depending on the usage


## QCDNUM  

"Fast QCD Evolution and Convolution", M. Botje,
NIKHEF-10-002, May 2010. 74pp. 
Published in Comput.Phys.Commun.182:490-532,2011. 
e-Print: arXiv:1005.1481 [hep-ph]


## LHAPDF6 

'LHAPDF6: parton density access in the LHC precision era'
A. Buckley et. al,
Published in Eur.Phys.J. C75 (2015) 3, 132, http://arxiv.org/abs/1412.7420

## APFEL  

"APFEL: A PDF Evolution Library with QED corrections"
V. Bertone, S. Carrazza and J. Rojo,
Published in Comput. Phys. Commun.  185 (2014) 1647,
e-Print: arXiv:1310.1394 [hep-ph].

When using APFEL with NLO QED corrections (TheoryType = 'DGLAP_APFEL_QED') , please cite:
"The photon PDF from high-mass Drell Yan data at the LHC"
xFitter Developers' Team (F. Giuli et al.),
e-Print: arXiv:1701.08553


## MINUIT  

F. James, M. Roos (CERN). Jul 1975. 38 pp. 
Published in Comput.Phys.Commun. 10 (1975) 343-367 
CERN-DD-75-20 
DOI: 10.1016/0010-4655(75)90039-9

## HERA data and HERAPDFs 

"Combination of measurements of inclusive deep inelastic e±p scattering cross sections and QCD analysis of HERA data."
By H1 and ZEUS Collaborations (H. Abramowicz et al.). DESY-15-039, Jun 19, 2015. 160 pp. 
Published in EPJC (2015) 75: 12.
e-Print: arXiv:1506.06042 [hep-ex]

"Combined Measurement and QCD Analysis of the Inclusive e+- p Scattering Cross Sections at HERA."
By H1 and ZEUS Collaborations (F.D. Aaron et al.). DESY-09-158, Oct 2009. 61pp. 
e-Print: arXiv:0911.0884 [hep-ex]

"Inclusive Deep Inelastic Scattering at High Q2 with Longitudinally Polarised Lepton Beams at HERA"
By H1 Collaboration (F.D. Aaron et al.) DESY-12-107
e-Print: arXiv:1206.7007 [hep-ex] 


## RT- variable flavour number scheme (for HF_SCHEME =  'RT' or 'RT FAST')

"An Ordered analysis of heavy flavor production in deep inelastic scattering"
R.S. Thorne, R.G. Roberts. RAL-TR-97-049, Sep 1997. 48pp. 
Published in Phys.Rev.D57:6871-6898,1998. 
e-Print: hep-ph/9709442

"A Variable-flavor number scheme for NNLO"
R.S. Thorne. CAVENDISH-HEP-2006-01, Jan 2006. 17pp. 
Published in Phys.Rev.D73:054019,2006. 
e-Print: hep-ph/0601245

 RT- variable flavour number scheme (for HF_SCHEME =  'RT OPT' or 'RT OPT FAST')
----------------------------------------------------------------------------------

"The Effect of Changes of Variable Flavour Number Scheme on PDFs and Predicted Cross Sections"
R.S. Thorne.  Jan 2012.  14pp.
Published in Phys.Rev. D86 (2012) 074017
e-Print: arXiv:1201.6180 [hep-ph]

## ACOT-variable flavour number scheme ( for HF_SCHEME = 'ACOT Full' or 'ACOT Chi')  

"Leptoproduction of Heavy Quarks II -- A Unified QCD Formulation of 
Charged and Neutral Current Processes from Fixed-target to Collider Energies"
M.A.G. Aivazis, John C. Collins, Fredrick I. Olness, Wu-Ki Tung
Journal reference:      Phys.Rev.D50:3102-3118,1994
DOI:    10.1103/PhysRevD.50.3102
Cite as:        arXiv:hep-ph/9312319v2


## ABM-fixed flavour number scheme 
---------------------------------------------------------------------------
	
     http://www-zeuthen.desy.de/~alekhin/OPENQCDRAD/ 

## FONLL general-mass variable-flavor-number scheme from APFEL

"APFEL: A PDF Evolution Library with QED corrections"
V. Bertone, S. Carrazza and J. Rojo,
Published in Comput. Phys. Commun.  185 (2014) 1647,
e-Print: arXiv:1310.1394 [hep-ph].

"Heavy quarks in deep-inelastic scattering"
S. Forte, E. Laenen, P. Nason and J. Rojo,
Published in Nucl. Phys. B 834 (2010) 116,
e-Print: arXiv:1001.2312 [hep-ph].

##  Chi2 definitions 

'H12000':
 CHI2SettingsName = 'StatScale', 'UncorSysScale', 'CorSysScale', 'UncorChi2Type', 'CorChi2Type'
 Chi2Settings     = 'NoRescale'  , 'NoRescale', 'Linear'     , 'Diagonal'     , 'Hessian'


"Measurement and QCD analysis of neutral and charged current cross-sections at HERA"
By H1 Collaboration (C. Adloff et al.). DESY-03-038, Apr 2003. 54pp. 
Published in Eur.Phys.J.C30:1-32,2003. 
e-Print: hep-ex/0304003

"A ZEUS next-to-leading-order QCD analysis of data on deep inelastic scattering"
By ZEUS Collaboration (S. Chekanov et al.). DESY-02-105, Aug 2002. 50pp. 
Published in Phys.Rev.D67:012007,2003. 
e-Print: hep-ex/0208023

'HERAPDF1.0': 
   CHI2SettingsName = 'StatScale', 'UncorSysScale', 'CorSysScale', 'UncorChi2Type', 'CorChi2Type'
   Chi2Settings     = 'Poisson'  , 'Linear',        'Linear'     , 'Diagonal'     , 'Hessian'

Add two references above (for 'H12000'), in addition: 

"Measurement of the Inclusive ep Scattering Cross Section at Low Q^2 and x at HERA"
F.D. Aaron et al. DESY-08-171, 2009. 90pp. 
Published in Eur.Phys.J.C63:625-678,2009. 
e-Print: arXiv:0904.0929 [hep-ex]


'HERAPDF2.0', 'H12012':
   CHI2SettingsName = 'StatScale', 'UncorSysScale', 'CorSysScale', 'UncorChi2Type', 'CorChi2Type'
   Chi2Settings     = 'Poisson'  , 'Linear',        'Linear'     , 'Diagonal'     , 'Hessian'
   Chi2ExtraParam = 'PoissonCorr' 

Add two references above (for 'H12000'), in addition: 

"Inclusive Deep Inelastic Scattering at High Q2 with Longitudinally Polarised Lepton Beams at HERA"
By H1 Collaboration (F.D. Aaron et al.)	DESY-12-107
e-Print: arXiv:1206.7007 [hep-ex]



## PDF Uncertainties 

"Multivariate fitting and the error matrix in global analysis of data"
J. Pumplin, D.R. Stump, W.K. Tung. MSU-HEP-07100, CERN-TH-2000-249, Aug 2000. 14pp. 
Published in Phys.Rev.D65:014011,2001. 
e-Print: hep-ph/0008191

"New generation of parton distributions with uncertainties from global QCD analysis"
J. Pumplin, D.R. Stump, J. Huston, H.L. Lai, Pavel M. Nadolsky, W.K. Tung, MSU-HEP-011101, 
Jan 2002. 44pp. 
e-Print: hep-ph/0201195

## APPLGRID 

"A posteriori inclusion of parton density functions in NLO QCD final-state calculations at hadron colliders: The APPLGRID Project"
Tancredi Carli, Dan Clements, Amanda Cooper-Sarkar, Claire Gwenlan, Gavin P. Salam, 
Frank Siegert, Pavel Starovoitov,  Mark Sutton. 2010. 33pp. 
Published in Eur.Phys.J.C66:503-524,2010. 
e-Print: arXiv:0911.2985 [hep-ph]

## FastNLO 

"FastNLO: Fast pQCD calculations for PDF fits"
T. Kluge, K. Rabbertz, M. Wobisch. DESY-06-186, FERMILAB-CONF-06-352-E, Sep 2006. 8pp. 
Presented at 14th International Workshop on Deep Inelastic Scattering (DIS 2006), Tsukuba, Japan, 20-24 Apr 2006. 
Published in *Tsukuba 2006, Deep inelastic scattering* 483-486 
e-Print: hep-ph/0609285

"Theory-Data Comparisons for Jet Measurements in Hadron-Induced Processes"
M. Wobisch, D. Britzger, T. Kluge, K. Rabbertz, F. Stober.
DESY 11-150, FERMILAB-PUB-11-418-PPD, Sep 2011. arXiv:1109.1310v1 [hep-ph]

D. Britzger, K. Rabbertz, F. Stober, M. Wobisch, "New features in version 2 of the fastNLO project",
in the proceedings of the XX International Workshop on Deep Inelastic Scattering (DIS12), 26-30th March 2012 hep-ph/1208.3641.

## APFELgrid 

"APFELgrid: a high performance tool for parton density determinations"
V. Bertone, S. Carrazza and N. P. Hartland,
Published in Comput. Phys. Commun. 212 (2017) 205.
e-Print: arXiv:1605.02070 

## HATHOR 

from http://www-zeuthen.desy.de/~moch/hathor/

## PDF Reweighting (if used) 

Description of NNPDF method to create NNPDF PDF sets:	

- "A first unbiased global NLO determination of parton distributions and their uncertainties"
Richard D. Ball, Luigi Del Debbio, Stefano Forte, Alberto Guffanti, Jose I. Latorre, Juan Rojo, Maria Ubiali
Published in Nucl.Phys. B838 (2010) 136-206
e-Print: arXiv:1002.4407 [hep-ph] 

Description of reweighting method:
 	
- "Reweighting NNPDFs: the W lepton asymmetry"
NNPDF Collaboration
Published in Nucl.Phys. B849 (2011) 112-143, Erratum-ibid. B854 (2012) 926-927, Erratum-ibid. B855 (2012) 927-928
e-Print: arXiv:1012.0836 [hep-ph]
 	
- "Reweighting and Unweighting of Parton Distributions and the LHC W lepton asymmetry data"
Richard D. Ball, Valerio Bertone, Francesco Cerutti, Luigi Del Debbio, Stefano Forte, Alberto Guffanti, Nathan P. Hartland, Jose I. Latorre, Juan Rojo, Maria Ubiali
Published in Nucl.Phys. B855 (2012) 608-638
e-Print: arXiv:1108.1758 [hep-ph]

Description of reweighting approach for Hessian PDF error sets:
"Study of Monte Carlo approach to experimental uncertainty propagation with MSTW 2008 PDFs"
G. Watt, R.S. Thorne
e-Print: arXiv:1205.4024 [hep-ph]


##  PDF profiling  

"PDF reweighting in the Hessian matrix approach"
H. Paukkunen and P. Zurita, PDF reweighting in the Hessian matrix approach,
JHEP 12 (2014) 100, arXiv:1402.6623 [hep-ph].

"QCD analysis of W- and Z-boson production at Tevatron"
XFitter developers, S. Camarda et al., QCD analysis of W- and Z-boson production at Tevatron,
Eur. Phys. J. C 75 (2015) 458, arXiv:1503.05221 [hep-ph].


## MNR code 

M. Mangano, P. Nason and G. Ridolfi
Published in Nucl. Phys. B 373 (1992) 295
from http://www.ge.infn.it/~ridolfi/hvqlibx.tgz



## Citations for data tables 

Please use the citations as given in headers of the files which are used in the fit.

