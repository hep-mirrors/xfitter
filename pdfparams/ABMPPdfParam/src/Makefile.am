
# Created by AddPdfParam.py on 2019-02-25

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../BasePdfParam/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libABMPPdfParam_xfitter.la
libABMPPdfParam_xfitter_la_SOURCES = ABMPPdfParam.cc

datadir = ${prefix}/yaml/pdfparams/ABMP
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
libABMPPdfParam_xfitter_la_LDFLAGS=-lBasePdfParam_xfitter -L$(libdir)
