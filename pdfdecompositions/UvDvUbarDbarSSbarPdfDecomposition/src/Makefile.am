
# Created by AddPdfDecomposition.py on 2019-02-25

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../pdfparams/BasePdfParam/include/ -Wall -fPIC -Wno-deprecated

#lib_LTLIBRARIES = libUvDvUbarDbarSSbarPdfDecomposition_xfitter.la
lib_LTLIBRARIES = libuvdvubardbarssbarpdfdecomposition_xfitter.la
libuvdvubardbarssbarpdfdecomposition_xfitter_la_SOURCES = UvDvUbarDbarSSbarPdfDecomposition.cc

datadir = ${prefix}/yaml/pdfdecompositions/UvDvUbarDbarSSbar
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
