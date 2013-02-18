AUTOMAKE_OPTIONS = foreign
SUBDIRS = minuit/src interfaces/src DY/src DIPOLE/src RT/src HS/src tools/DrawResultsPackage \
          tools/MakeLHAPDF FastNLO/src DiffDIS/src ACOT/src ABKM/src Cascade/src
if ENABLE_HATHOR
  SUBDIRS+= Hathor/src
endif

SUBDIRS+= NNPDF/src

SUBDIRS+= src doc/tex

# DOXYGEN SUPPORT
include $(top_srcdir)/aminclude.am

# ensure the distribution of the doxygen configuration file
EXTRA_DIST = doxygen.cfg

TESTS= ./tools/check.sh

