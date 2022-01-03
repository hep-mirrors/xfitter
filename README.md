# xFitter   --- PDF fit program from HERA.

xFitter is an open source QCD fit framework desinged to extract PDFs 
and assess the impact of new data. The xFitter project is is a 
common initiative by the H1 and ZEUS collaborations and extended 
to the LHC collaborations to provide precision QCD analyses. 
xFitter has been used as one of the main software packages for the 
determination of the HERA proton parton densities (PDFs), HERAPDFs. 
For further details please check the xfitter.org web page.


The current package includes code to fit DIS inclusive cross section 
data, Drell-Yan, jet and ttbar processes (using APPLGRID and FastNLO
interfaces). The program is distributed under the GPL v3 license, see
LICENCE file for more details. The program uses the QCD evolution 
package QCDNUM developed by M. Botje and includes other parts of the code:
-- VFNS from R. Thorne, G. Watt (MSTW) @ LO, NLO, NNLO
-- VFNS from F. Olness (ACOT) @ LO, NLO and NNLO, NNNLO corrections for FL 
-- VFNS from APFEL (FONLL) @ LO, NLO and NNLO
-- FFNS from S. Alekhin (ABM) @ NLO, NNLO (pole and running heavy quark masses)
-- PDF error estimation from J. Pumplin
-- Bayesian reweighting tool from A. Guffanti (a la NNPDF) and based on
   EIGENVECTORS from G. Watt (a la MSTW).
-- total ttbar production cross sections via HATHOR (S. Moch et al.)
-- differential ttbar production cross sections with DiffTop (M. Guzzi, S. Moch et al.)
-- MNR calculation for heavy quark production (Mangano, Nason and Ridolfi, 
   implemented by O.Zenaiev)

If the results obtained with the program are to be included in a scientific 
publication, please use the citations as suggested by the REFERENCES file. 

For support information, please visit https://wiki-zeuthen.desy.de/xFitter/xFitter

For documentation of the code, please use https://gitlab.cern.ch/fitters/xfitter/-/wikis/home


# Installation  Instructions:
Please refer to the [INSTALLATION](INSTALLATION) file.
=====================

# Brief Introduction
=====================
  ## Steering cards
    The software behaviour is controlled by two files with steering commands.
    These files have predefined names:

      steering.txt   --  controls main "stable" (un-modified during 
                         minimisation) parameters. The file also contains
                         names of data files to be fitted to, definition 
                         of kinematic cuts                        
                        
      parameters.yaml -- new unified parameter file, to control parameters transfered to the
                         reaction interfaces. 


  ## Inclusion of data files
  
    Inclusion of the data files is controlled by &InFiles namelist in the 
    steering.txt file. For example, by default the following four HERA-I
    files are included:

&InFiles
    NInputFiles = 7

    InputFileNames(1) = 'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_920.dat'
    InputFileNames(2) = 'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_820.dat'
    InputFileNames(3) = 'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_575.dat'
    InputFileNames(4) = 'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_460.dat'
    InputFileNames(5) = 'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCem.dat'
    InputFileNames(6) = 'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_CCep.dat'
    InputFileNames(7) = 'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_CCem.dat'

&End

    To include more files:
       -- Increase NInputFiles
       -- Specify  InputFileNames()

another option would be:

    NInputFiles = 7
    InputFileNames =
 'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_920.dat'
 'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_820.dat'
 'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_575.dat'
 'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_460.dat'
 'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCem.dat'
 'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_CCep.dat'
 'datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_CCem.dat'

-> then the order does matter of the files listed.

   Inclusion of the statistical or systematic correlations of the data 
   in the fit is done via &InCorr namelist:

&InCorr
  ! Number of correlation (statistical, systematical or full) files
    NCorrFiles = 1
    CorrFileNames(1) = 'datafiles/hera/H1_NormInclJets_HighQ2_99-07___H1_NormInclJets_HighQ2_99-07.corr'
&End
   in this case the statistical correlations for H1_NormInclJets_HighQ2_99-07 data file are included
   (the method also allows to include correlations between data sets via file names, i.e:
    H1_NormInclJets_HighQ2_99-07___H1_InclJets_HighQ2_99-00.dat.corr)


   As additional option for data sets with covariance matrix, it is possible
   to convert covariance matrix to nuisance parameter representation
   (following the prescription suggested by J. Gao and P. Nadolsky  in arXiv:1401.0013):

&CovarToNuisance
   ! Global switch for using nuisance param representation for covariance mat.
  LConvertCovToNui = .true.

   ! Tolerance -- zero means exact transformation
  Tolerance = 0.0

   ! The following lines allow to adjust error scaling properties 
   ! (default: :M - multiplicative, A - additive)
  DataName     = 'CMS electon Asymmetry rapidity', 'CMS W muon asymmetry'
  DataSystType = ':A', ':A'
&End

