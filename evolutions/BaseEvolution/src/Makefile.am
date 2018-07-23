AM_CXXFLAGS = -I$(srcdir)/../include -I$(srcdir)/../../../pdfdecompositions/BasePdfDecomposition/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libbaseevolution_xfitter.la

libbaseevolution_xfitter_la_SOURCES = BaseEvolution.cc

datadir = ${prefix}/yaml/evolutions/BaseEvolution
data_DATA = ../yaml/parameters.yaml
dist_noinst_HEADERS = ../include ../yaml
