
# Created by AddEvolution.py on 2018-08-20


#ifdef LHAPDF_ENABLED

  AM_CXXFLAGS = -I$(srcdir)/../include -I$(srcdir)/../../BaseEvolution/include -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated $(LHAPDF_CPPFLAGS)

  lib_LTLIBRARIES = liblhapdf_xfitter.la
  liblhapdf_xfitter_la_SOURCES = EvolutionLHAPDF.cc

# liblhapdf_xfitter_la_LDFLAGS = place_if_needed  

#endif

datadir = ${prefix}/yaml/evolutions/LHAPDF
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
 
