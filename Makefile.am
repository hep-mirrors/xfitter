AUTOMAKE_OPTIONS = foreign
SUBDIRS = minuit  DY DIPOLE RT HS tools/H1FitterDrawResults \
          tools/MakeLHAPDF FastNLO/src
if ENABLE_HATHOR
  SUBDIRS+= Hathor/src
endif
SUBDIRS+= src doc/tex

# DOXYGEN SUPPORT
include $(top_srcdir)/aminclude.am

# ensure the distribution of the doxygen configuration file
EXTRA_DIST = doxygen.cfg


