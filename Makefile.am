AUTOMAKE_OPTIONS = foreign
SUBDIRS = minuit/src interfaces/src DY/src DIPOLE/src RT/src EW/src common common/linalg \
          common/num_utils tools/DrawResultsPackage \
          tools/DrawPdfs tools/MakeLHAPDF FastNLO/src DiffDIS/src ACOT/src ABM/src FONLL/src Cascade/src \
          genetic/mixmax_r004 genetic/src QEDevol/src \
          include interfaces/include FastNLO/include FastNLO/include/fastnlotk DiffDIS/include \
          DY/include tools/DrawResultsPackage/include tools/DrawPdfs/include input_steering \
          NNPDF/include pdf2yaml tools/PostProcess/ 
if ENABLE_HATHOR
  SUBDIRS+= Hathor/src
endif

if ENABLE_HVQMNR
  SUBDIRS+= HVQMNR/src
endif

SUBDIRS+= NNPDF/src

SUBDIRS+= common
if BUILD_DOCS
  DOCS=doc/tex/manual doc/tex/paper
endif
#$(info ${DOCS})
SUBDIRS+= src $(DOCS)

include $(top_srcdir)/aminclude.am

# Specify which SUBDIRS to distribute:
DIST_SUBDIRS = minuit/src interfaces/src DY/src DIPOLE/src RT/src EW/src tools/DrawResultsPackage \
          tools/DrawPdfs tools/MakeLHAPDF FastNLO/src DiffDIS/src ACOT/src ABM/src FONLL/src Cascade/src \
          genetic/mixmax_r004 genetic/src common/num_utils \
          include interfaces/include FastNLO/include FastNLO/include/fastnlotk DiffDIS/include \
          DY/include tools/DrawResultsPackage/include tools/DrawPdfs/include input_steering \
          NNPDF/include NNPDF/src common common/linalg Hathor/src pdf2yaml tools/PostProcess src


# DOXYGEN SUPPORT
# ensure the distribution of the doxygen configuration file
#EXTRA_DIST = doxygen.cfg

# Tell which program should run the .sh scripts.
#SH_LOG_COMPILER = $(SHELL) -ex
#TESTS_ENVIRONMENT = $(SHELL) 
#TESTS= ./tools/check.sh

EXTRA_DIST= README README-MacOSX LICENCE REFERENCES ewparam.txt minuit.in.txt steering.txt 
bin_SCRIPTS = xfitter-config
