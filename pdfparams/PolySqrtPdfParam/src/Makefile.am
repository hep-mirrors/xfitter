
# Created by AddPdfParam.py on 2018-08-16

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../BasePdfParam/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libPolySqrtPdfParam_xfitter.la
libPolySqrtPdfParam_xfitter_la_SOURCES = PolySqrtPdfParam.cc

dist_noinst_HEADERS = ../include
libPolySqrtPdfParam_xfitter_la_LDFLAGS=-lBasePdfParam_xfitter -L$(libdir)
