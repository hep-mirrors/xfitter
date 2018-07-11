
# Created by AddPdfParam.py on 2018-07-11

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../BasePdfParam/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libHERAPDF_PdfParam_xfitter.la
libHERAPDF_PdfParam_xfitter_la_SOURCES = HERAPDF_PdfParam.cc

datadir = ${prefix}/yaml/pdfparams/HERAPDF_
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
