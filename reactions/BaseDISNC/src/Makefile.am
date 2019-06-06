
# Created by AddReaction.py on 2017-04-08

AM_CXXFLAGS = $(QCDNUM_CPPFLAGS)\
  -I$(srcdir)/../include\
  -I$(srcdir)/../../../include\
  -I$(srcdir)/../../../interfaces/include\
  -I$(srcdir)/../../../evolutions/BaseEvolution/include\
  -I$(srcdir)/../../../evolutions/QCDNUM/include\
  -Wall -fPIC -Wno-deprecated

lib_LTLIBRARIES = libbasedisnc_xfitter.la
libbasedisnc_xfitter_la_SOURCES = ReactionBaseDISNC.cc

# libbasedisnc_xfitter_la_LDFLAGS = place_if_needed  

dist_noinst_HEADERS = ../include ../yaml
