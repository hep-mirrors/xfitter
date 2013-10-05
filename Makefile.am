AUTOMAKE_OPTIONS = foreign
SUBDIRS = minuit/src interfaces/src DY/src DIPOLE/src RT/src EW/src tools/DrawResultsPackage \
          tools/MakeLHAPDF FastNLO/src DiffDIS/src ACOT/src ABM/src Cascade/src \
          genetic/mixmax_r004 genetic/src \
          include interfaces/include FastNLO/include DiffDIS/include \
          DY/include tools/DrawResultsPackage/include input_steering datafiles/bcdms \
          datafiles/hera datafiles/lhc/atlas/Jets2010 datafiles/lhc/atlas/WZ2010 \
          datafiles/lhc/cms datafiles/tevatron theoryfiles theoryfiles/updf \
          theoryfiles/cms/CMS-EWK-10-010 theoryfiles/cms/EWK-11-005 \
          theoryfiles/cms/SPM_12_001 theoryfiles/atlas/Jets2010/R04 \
          theoryfiles/atlas/Jets2010/R06 theoryfiles/atlas/WZ2010 \
          NNPDF/include common common/linalg
if ENABLE_HATHOR
  SUBDIRS+= Hathor/src
endif

SUBDIRS+= NNPDF/src

SUBDIRS+= common
SUBDIRS+= src doc/tex/manual doc/tex/paper

# DOXYGEN SUPPORT
include $(top_srcdir)/aminclude.am

# ensure the distribution of the doxygen configuration file
EXTRA_DIST = doxygen.cfg

# Tell which program should run the .sh scripts.
#SH_LOG_COMPILER = $(SHELL) -ex
#TESTS_ENVIRONMENT = $(SHELL) 
#TESTS= ./tools/check.sh

