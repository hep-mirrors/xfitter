
# Created by AddReaction.py on 2017-10-05

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libbasediscc_xfitter.la
libbasediscc_xfitter_la_SOURCES = ReactionBaseDISCC.cc

# libbasediscc_xfitter_la_LDFLAGS = place_if_needed  
# link to main xFitter library (for IntegrateDIS)
libbasediscc_xfitter_la_LDFLAGS = -lxfmain -L$(srcdir)/../../../src/.libs

dist_noinst_HEADERS = ../include ../yaml
