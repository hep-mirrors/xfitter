# xFitter - A PDF fit program from HERA.

xFitter is an open source QCD fit framework desinged to extract PDFs and assess the impact of new data. 
The xFitter project is is a common initiative by the H1 and ZEUS collaborations and extended 
to the LHC collaborations to provide precision QCD analyses. 
xFitter has been used as one of the main software packages for the 
determination of the HERA proton parton densities (PDFs), HERAPDFs. 
For further details and information, please check the [xfitter](xfitter.org) web page.

The current package includes code to fit DIS inclusive cross section data,
Drell-Yan, jet, ttbar and other processes (using APPLGRID and FastNLO interfaces).

## Installation  Instructions:
Please refer to the [INSTALLATION](INSTALLATION) file.

## Basic Usage
   The software behaviour is controlled by two files with steering commands.

      steering.txt   --  controls the main "stable" (un-modified during minimisation) parameters.
      		     The file also contains
                         names of data files to be fitted to, definition 
                         of kinematic cuts                        
                        
      parameters.yaml -- new unified parameter file, to control parameters transfered to the
                         reaction interfaces.
			 
For documentation of the code, please use https://gitlab.cern.ch/fitters/xfitter/-/wikis/home


## Citation Policy
If the results obtained with the program are to be included in a scientific publication,
 add a reference to the xFitter paper  together with the xFitter web portal: xfitter.org

"HERAFitter, Open Source QCD Fit Project"
By S. Alekhin at al., DESY Report 14-188, 7 Nov 2014, Published in EPJC (2015), 75: 304,
e-Print: arXiv:1410.4412 [hep-ex]

In addition to the main xFitter paper, additional references to be cited
for specific modules and software packages can be found in the [REFERENCES](REFERENCES) file. 

## Licence
xFitter is distributed under the GPL v3 license, see [LICENCE](LICENCE) for more details.
