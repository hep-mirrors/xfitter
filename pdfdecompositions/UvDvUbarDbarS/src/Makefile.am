
# Created by AddPdfDecomposition.py on 2018-07-11

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include -I$(srcdir)/../../../pdfparams/BasePdfParam/include/ -I$(srcdir)/../../BasePdfDecomposition/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libuvdvubardbars_xfitter.la
libuvdvubardbars_xfitter_la_SOURCES = UvDvUbarDbarS.cc

datadir = ${prefix}/yaml/pdfdecompositions/UvDvUbarDbarS
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
