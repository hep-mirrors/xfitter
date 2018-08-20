
# Created by AddMinimizer.py on 2018-08-17

if ENABLE_CERES

  AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include $(CERES_CPPFLAGS)  -I$(srcdir)/../../BaseMinimizer/include -Wall -fPIC -Wno-deprecated 

  lib_LTLIBRARIES = libCERESMinimizer_xfitter.la
  libCERESMinimizer_xfitter_la_SOURCES = CERESMinimizer.cc

  libCERESMinimizer_xfitter_la_LDFLAGS = $(CERES_LDFLAGS) -lbaseminimizer_xfitter -L$(libdir)

endif


datadir = ${prefix}/yaml/minimizers/CERES
data_DATA = ../yaml/parameters.yaml
dist_noinst_HEADERS = ../include ../yaml

