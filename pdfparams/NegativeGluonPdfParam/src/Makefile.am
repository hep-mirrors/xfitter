
# Created by AddPdfParam.py on 2018-09-26

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include -Wall -fPIC -Wno-deprecated

lib_LTLIBRARIES = libNegativeGluonPdfParam_xfitter.la
libNegativeGluonPdfParam_xfitter_la_SOURCES = NegativeGluonPdfParam.cc

datadir = ${prefix}/yaml/pdfparams/NegativeGluonPdfParam
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml

libNegativeGluonPdfParam_xfitter_la_LDFLAGS=-lBasePdfParam_xfitter -L$(libdir)
