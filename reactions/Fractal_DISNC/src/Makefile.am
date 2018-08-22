
# Created by AddReaction.py on 2017-05-15

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libfractal_disnc_xfitter.la
libfractal_disnc_xfitter_la_SOURCES = ReactionFractal_DISNC.cc

# libfractal_disnc_xfitter_la_LDFLAGS = place_if_needed  

datadir = ${prefix}/yaml/reactions/Fractal_DISNC
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
