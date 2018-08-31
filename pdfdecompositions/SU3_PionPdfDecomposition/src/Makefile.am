
# Created by AddPdfDecomposition.py on 2018-08-07

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include -I$(srcdir)/../../../pdfparams/BasePdfParam/include/   -I$(srcdir)/../../BasePdfDecomposition/include -Wall -fPIC -Wno-deprecated \
              -I$(srcdir)/../../../pdfparams/PolySqrtPdfParam/include/ #Hacks - remove later

lib_LTLIBRARIES=libSU3_PionPdfDecomposition_xfitter.la
libSU3_PionPdfDecomposition_xfitter_la_SOURCES = SU3_PionPdfDecomposition.cc

datadir = ${prefix}/yaml/pdfdecompositions/SU3_Pion
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml

# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx These hacks are temporary 
libSU3_PionPdfDecomposition_xfitter_la_LDFLAGS=-lPolySqrtPdfParam_xfitter -lBasePdfParam_xfitter -L$(libdir)
