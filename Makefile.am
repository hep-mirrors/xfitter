ACLOCAL_AMFLAGS=-I m4
AUTOMAKE_OPTIONS = foreign
SUBDIRS = minuit/src interfaces/src DY/src DIPOLE/src RT/src EW/src common common/linalg \
          common/num_utils \
          tools/draw tools/MakeLHAPDF FastNLO/src DiffDIS/src ACOT/src SACOT/src ABM/src FONLL/src Cascade/src \
          genetic/mixmax_r004 genetic/src QEDevol/src \
          include interfaces/include FastNLO/include FastNLO/include/fastnlotk DiffDIS/include \
          DY/include tools/draw/include  \
          pdf2yaml tools/process/ tools/ 

DIST_SUBDIRS = $(SUBDIRS)  Hathor/src  HVQMNR/src doc/tex/manual  ACOT SACOT datafiles QEDevol/include \
           datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042 doc/logo examples Hathor/interface \
           FONLL HVQMNR HVQMNR/include/ python/ input_steering/

if ENABLE_HATHOR
  SUBDIRS+= Hathor/src
endif

if ENABLE_HVQMNR
  SUBDIRS+= HVQMNR/src
endif

SUBDIRS+= common
if BUILD_DOCS
  DOCS=doc/tex/manual doc/tex/paper
endif
#$(info ${DOCS})
SUBDIRS+= src $(DOCS)

include $(top_srcdir)/aminclude.am



# Specify which SUBDIRS to distribute:
#DIST_SUBDIRS = minuit/src interfaces/src DY/src DIPOLE/src RT/src EW/src \
#          tools/draw tools/MakeLHAPDF FastNLO/src DiffDIS/src ACOT/src SACOT/src ABM/src FONLL/src Cascade/src \
#          genetic/mixmax_r004 genetic/src common/num_utils \
#          include interfaces/include FastNLO/include FastNLO/include/fastnlotk DiffDIS/include \
#          DY/include tools/draw/include input_steering \
#          common common/linalg Hathor/src pdf2yaml tools/process src


# DOXYGEN SUPPORT
# ensure the distribution of the doxygen configuration file
#EXTRA_DIST = doxygen.cfg

# Tell which program should run the .sh scripts.
#SH_LOG_COMPILER = $(SHELL) -ex
#TESTS_ENVIRONMENT = $(SHELL) 
#TESTS= ./tools/check.sh

EXTRA_DIST= README INSTALLATION LICENCE REFERENCES \
  ewparam.txt minuit.in.txt steering.txt testPythonInterface.py
bin_SCRIPTS = xfitter-config tools/xfitter-getdata.sh
