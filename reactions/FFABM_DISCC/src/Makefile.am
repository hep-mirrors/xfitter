
# Created by AddReaction.py on 2017-10-09

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -I$(srcdir)/../../../reactions/BaseDISCC/include -Wall -fPIC -Wno-deprecated

lib_LTLIBRARIES = libffabm_discc_xfitter.la
libffabm_discc_xfitter_la_SOURCES = ReactionFFABM_DISCC.cc

libffabm_discc_xfitter_la_LDFLAGS = -lmyabm -lbasediscc_xfitter -L$(top_srcdir)/reactions/BaseDISCC/src/.libs -L$(top_srcdir)/ABM/src/.libs

datadir = ${prefix}/yaml/reactions/FFABM_DISCC
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml

