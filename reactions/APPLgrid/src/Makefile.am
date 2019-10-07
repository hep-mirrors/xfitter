if ENABLE_APPLGRID
  AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated 

  lib_LTLIBRARIES = libapplgrid_xfitter.la
  libapplgrid_xfitter_la_SOURCES = ReactionAPPLgrid.cc

  libapplgrid_xfitter_la_CPPFLAGS = $(APPLGRID_CPPFLAGS) $(ROOT_CFLAGS)

  libapplgrid_xfitter_la_LDFLAGS = $(APPLGRID_LDFLAGS)
endif

datadir = ${prefix}/yaml/reactions/APPLgrid
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
