
# Created by AddPdfDecomposition.py on 2018-08-07

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include -I$(srcdir)/../../../pdfparams/BasePdfParam/include/ -I$(srcdir)/../../BasePdfDecomposition/include -Wall -fPIC -Wno-deprecated

lib_LTLIBRARIES=libSU3_PionPdfDecomposition_xfitter.la
libSU3_PionPdfDecomposition_xfitter_la_SOURCES = SU3_PionPdfDecomposition.cc

datadir = ${prefix}/yaml/pdfdecompositions/SU3_Pion
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
