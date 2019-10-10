
# Created by AddReaction.py on 2019-10-07

if ENABLE_DYTURBO
  AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated
  AM_CXXFLAGS += $(DYTURBO_CPPFLAGS) $(LHAPDF_CPPFLAGS) # $(CUBA_CPPFLAGS) 
  lib_LTLIBRARIES = libdyturbo_xfitter.la
  libdyturbo_xfitter_la_SOURCES = ReactionDYTurbo.cc

  libdyturbo_xfitter_la_LDFLAGS = $(DYTURBO_LDFLAGS) $(CUBA_LDFLAGS)
endif

datadir = ${prefix}/yaml/reactions/DYTurbo
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
 
