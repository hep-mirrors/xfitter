AUTOMAKE_OPTIONS = foreign
SUBDIRS = minuit/src interfaces/src DY/src DIPOLE/src RT/src EW/src tools/DrawResultsPackage \
	  tools/DrawPdfs tools/MakeLHAPDF FastNLO/src DiffDIS/src ACOT/src ABM/src Cascade/src \
          genetic/mixmax_r004 genetic/src \
          include interfaces/include FastNLO/include DiffDIS/include \
          DY/include tools/DrawResultsPackage/include tools/DrawPdfs/include input_steering \
          NNPDF/include common common/linalg
if ENABLE_HATHOR
  SUBDIRS+= Hathor/src
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
          tools/DrawPdfs tools/MakeLHAPDF FastNLO/src DiffDIS/src ACOT/src ABM/src Cascade/src \
          genetic/mixmax_r004 genetic/src \
          include interfaces/include FastNLO/include DiffDIS/include \
          DY/include tools/DrawResultsPackage/include tools/DrawPdfs/include input_steering datafiles/bcdms \
          NNPDF/include NNPDF/src common common/linalg Hathor/src src

# Also datafiles
DIST_SUBDIRS+= datafiles/hera datafiles/lhc/atlas/Jets2010 datafiles/lhc/atlas/WZ2010 \
          datafiles/lhc/cms datafiles/tevatron 

# DOXYGEN SUPPORT
# ensure the distribution of the doxygen configuration file
#EXTRA_DIST = doxygen.cfg

# Tell which program should run the .sh scripts.
#SH_LOG_COMPILER = $(SHELL) -ex
#TESTS_ENVIRONMENT = $(SHELL) 
#TESTS= ./tools/check.sh

EXTRA_DIST= README README-MacOSX LICENCE REFERENCES ewparam.txt minuit.in.txt steering.txt 
