AM_CXXFLAGS = -I$(srcdir)/../../../include -Wall -fPIC -Wno-deprecated

lib_LTLIBRARIES = libbaseminimizer_xfitter.la

libbaseminimizer_xfitter_la_SOURCES = BaseMinimizer.cc

datadir = ${prefix}/yaml/minimizers/BaseMinimizer
data_DATA = ../yaml/parameters.yaml
dist_noinst_HEADERS = ../yaml
