
# Created by AddPdfParam.py on 2019-02-25

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../BasePdfParam/include -Wall -fPIC -Wno-deprecated

lib_LTLIBRARIES = libABMPseaPdfParam_xfitter.la
libABMPseaPdfParam_xfitter_la_SOURCES = ABMPseaPdfParam.cc

datadir = ${prefix}/yaml/pdfparams/ABMPsea
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
libABMPseaPdfParam_xfitter_la_LDFLAGS=-lBasePdfParam_xfitter -L$(libdir)
