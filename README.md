# xFitter   --- A PDF fit program from HERA.

xFitter is an open source QCD fit framework desinged to extract PDFs and assess the impact of new data. 
The xFitter project is is a common initiative by the H1 and ZEUS collaborations and extended 
to the LHC collaborations to provide precision QCD analyses. 
xFitter has been used as one of the main software packages for the 
determination of the HERA proton parton densities (PDFs), HERAPDFs. 
For further details please check the xfitter.org web page.


The current package includes code to fit DIS inclusive cross section 
data, Drell-Yan, jet and ttbar processes (using APPLGRID and FastNLO
interfaces). The program is distributed under the GPL v3 license, see
[LICENCE](LICENCE) file for more details. The program uses the QCD evolution 
package QCDNUM developed by M. Botje and includes other parts of the code:
- VFNS from R. Thorne, G. Watt (MSTW) @ LO, NLO, NNLO
- VFNS from F. Olness (ACOT) @ LO, NLO and NNLO, NNNLO corrections for FL 
- VFNS from APFEL (FONLL) @ LO, NLO and NNLO
- FFNS from S. Alekhin (ABM) @ NLO, NNLO (pole and running heavy quark masses)
- PDF error estimation from J. Pumplin
- Bayesian reweighting tool from A. Guffanti (a la NNPDF) and based on EIGENVECTORS from G. Watt (a la MSTW).
- total ttbar production cross sections via HATHOR (S. Moch et al.)
- differential ttbar production cross sections with DiffTop (M. Guzzi, S. Moch et al.)
- MNR calculation for heavy quark production (Mangano, Nason and Ridolfi, implemented by O.Zenaiev)

If the results obtained with the program are to be included in a scientific 
publication, please use the citations as suggested by the [REFERENCES](REFERENCES) file. 

For support information, please visit https://wiki-zeuthen.desy.de/xFitter/xFitter

For documentation of the code, please use https://gitlab.cern.ch/fitters/xfitter/-/wikis/home


## Installation  Instructions:
Please refer to the [INSTALLATION](INSTALLATION) file.
