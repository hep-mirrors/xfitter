
# Created by AddReaction.py on 2017-10-05

AM_CXXFLAGS = $(QCDNUM_CPPFLAGS)\
  -I$(srcdir)/../include\
  -I$(srcdir)/../../../include\
  -I$(srcdir)/../../../interfaces/include\
  -I$(srcdir)/../../../evolutions/BaseEvolution/include\
  -Wall -fPIC -Wno-deprecated

lib_LTLIBRARIES = libbasediscc_xfitter.la
libbasediscc_xfitter_la_SOURCES = ReactionBaseDISCC.cc

# libbasediscc_xfitter_la_LDFLAGS = place_if_needed  

dist_noinst_HEADERS = ../include ../yaml
