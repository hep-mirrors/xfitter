AM_CXXFLAGS = -I$(srcdir)/../include  -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libbasepdfdecomposition_xfitter.la

libbasepdfdecomposition_xfitter_la_SOURCES = BasePdfDecomposition.cc

datadir = ${prefix}/yaml/pdfdecompositions/BasePdfDecomposition
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
