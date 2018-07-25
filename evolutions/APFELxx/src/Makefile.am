if ENABLE_APFELXX
   AM_CXXFLAGS = $(APFEL_CPPFLAGS) \
	-I$(srcdir)/../include \
	-I$(srcdir)/../../BaseEvolution/include \
	-I$(srcdir)/../../../pdfdecompositions/BasePdfDecomposition/include \
	-I$(srcdir)/../../../include \
	-Wall -fPIC -Wno-deprecated 

   lib_LTLIBRARIES = libapfelxx_xfitter.la
   libapfelxx_xfitter_la_SOURCES = EvolutionAPFELxx.cc

   libapfelxx_xfitter_la_LDFLAGS = $(APFELXX_LDFLAGS)  #-L$(top_srcdir)/evolutions/BaseEvolution/src/.libs -lbaseevolution_xfitter
endif

datadir = ${prefix}/yaml/evolutions/APFELxx
data_DATA = ../yaml/parameters.yaml
dist_noinst_HEADERS = ../include ../yaml
