
# Created by AddMinimizer.py on 2018-08-17

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include   -I$(srcdir)/../../../common  -I$(srcdir)/../../BaseMinimizer/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libMINUITMinimizer_xfitter.la
libMINUITMinimizer_xfitter_la_SOURCES = MINUITMinimizer.cc

datadir = ${prefix}/yaml/minimizers/MINUIT
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml

libMINUITMinimizer_xfitter_la_LDFLAGS = -lbaseminimizer_xfitter -L$(libdir)
