
# Created by AddPdfDecomposition.py on 2018-07-12
if ENABLE_LHAPDF
   AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include $(LHAPDF_CPPFLAGS) \
                 -Wall -fPIC -Wno-deprecated 

   lib_LTLIBRARIES = liblhapdfdecomposition_xfitter.la
   liblhapdfdecomposition_xfitter_la_SOURCES = LHAPDFDecomposition.cc
endif

datadir = ${prefix}/yaml/pdfdecompositions/LHAPDFDecomposition
data_DATA = ../yaml/parameters.yaml
dist_noinst_HEADERS = ../include ../yaml
