
# Created by AddPdfDecomposition.py on 2018-08-07

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include -I$(srcdir)/../../../pdfparams/BasePdfParam/include/   -I$(srcdir)/../../BasePdfDecomposition/include -Wall -fPIC -Wno-deprecated \
              -I$(srcdir)/../../../pdfparams/HERAPDF_PdfParam/include/ #Hacks - remove later

lib_LTLIBRARIES=libGRV_PionPdfDecomposition_xfitter.la
libGRV_PionPdfDecomposition_xfitter_la_SOURCES = GRV_PionPdfDecomposition.cc

datadir = ${prefix}/yaml/pdfdecompositions/GRV_Pion
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml

# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx  HERAPDF is temporary 
libGRV_PionPdfDecomposition_xfitter_la_LDFLAGS=-lHERAPDF_PdfParam_xfitter -lBasePdfParam_xfitter -L$(libdir)
