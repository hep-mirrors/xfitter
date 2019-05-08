
# Created by AddPdfParam.py on 2019-02-25

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../BasePdfParam/include -Wall -fPIC -Wno-deprecated

lib_LTLIBRARIES = libABMPgluonPdfParam_xfitter.la
libABMPgluonPdfParam_xfitter_la_SOURCES = ABMPgluonPdfParam.cc

datadir = ${prefix}/yaml/pdfparams/ABMPgluon
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
libABMPgluonPdfParam_xfitter_la_LDFLAGS=-lBasePdfParam_xfitter -L$(libdir)
