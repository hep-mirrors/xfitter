
# Created by AddReaction.py on 2018-07-16

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libafb_xfitter.la
libafb_xfitter_la_SOURCES = ReactionAFB.cc

libafb_xfitter_la_LDFLAGS =  -lgsl -lgslcblas -lm

datadir = ${prefix}/yaml/reactions/AFB
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
 
