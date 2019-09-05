
# Created by AddReaction.py on 2019-09-05

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated

lib_LTLIBRARIES = libubardbar_xfitter.la
libubardbar_xfitter_la_SOURCES = ReactionUbarDbar.cc

# libubardbar_xfitter_la_LDFLAGS = place_if_needed

datadir = ${prefix}/yaml/reactions/UbarDbar
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
 