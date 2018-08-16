AM_CXXFLAGS =-O0 -g -I$(srcdir)/../include -I$(srcdir)/../../../include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libBasePdfParam_xfitter.la

libBasePdfParam_xfitter_la_SOURCES = BasePdfParam.cc

datadir = ${prefix}/yaml/pdfparams/BasePdfParam
data_DATA = ../yaml/parameters.yaml


dist_noinst_HEADERS = ../include ../yaml
