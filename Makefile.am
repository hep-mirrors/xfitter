AUTOMAKE_OPTIONS = foreign
SUBDIRS = minuit/src interfaces/src DY/src DIPOLE/src RT/src EW/src tools/DrawResultsPackage \
          tools/MakeLHAPDF FastNLO/src DiffDIS/src ACOT/src ABM/src Cascade/src
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
TESTS_ENVIRONMENT = $(SHELL) 
TESTS= ./tools/check.sh

