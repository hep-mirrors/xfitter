
# Created by AddPdfParam.py on 2019-02-25

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../BasePdfParam/include -Wall -fPIC -Wno-deprecated

lib_LTLIBRARIES = libABMPvalencePdfParam_xfitter.la
libABMPvalencePdfParam_xfitter_la_SOURCES = ABMPvalencePdfParam.cc

datadir = ${prefix}/yaml/pdfparams/ABMPvalence
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
libABMPvalencePdfParam_xfitter_la_LDFLAGS=-lBasePdfParam_xfitter -L$(libdir)
