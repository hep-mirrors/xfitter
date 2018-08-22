
# Created by AddReaction.py on 2018-08-03

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libkmatrix_xfitter.la
libkmatrix_xfitter_la_SOURCES = ReactionKMatrix.cc

# libkmatrix_xfitter_la_LDFLAGS = place_if_needed  

datadir = ${prefix}/yaml/reactions/KMatrix
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
 