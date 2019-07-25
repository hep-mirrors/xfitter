
#ifdef LHAPDF_ENABLED

  AM_CXXFLAGS = -I$(srcdir)/../include -I$(srcdir)/../../BaseEvolution/include -I$(srcdir)/../../../include -Wall -fPIC -Wno-deprecated $(LHAPDF_CPPFLAGS)

  lib_LTLIBRARIES = liblhapdf_xfitter.la
  liblhapdf_xfitter_la_SOURCES = EvolutionLHAPDF.cc

# liblhapdf_xfitter_la_LDFLAGS = place_if_needed  

#endif

datadir = ${prefix}/yaml/evolutions/LHAPDF
dist_noinst_HEADERS=../include
 
