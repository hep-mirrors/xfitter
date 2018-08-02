# Created by AddReaction.py on 2017-03-28

if ENABLE_APPLGRID
  AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated 

  lib_LTLIBRARIES = libapplgrid_xfitter.la
  libapplgrid_xfitter_la_SOURCES = ReactionAPPLgrid.cc

  libapplgrid_xfitter_la_LDFLAGS:=
  libapplgrid_xfitter_la_CPPFLAGS = $(APPLGRID_CPPFLAGS) $(ROOT_CFLAGS)
if ENABLE_LHAPDF
  AM_CXXFLAGS+=$(LHAPDF_CPPFLAGS)
  libapplgrid_xfitter_la_CPPFLAGS+=$(LHAPDF_CPPFLAGS)
  libapplgrid_xfitter_la_LDFLAGS+=$(LHAPDF_LDFLAGS)
endif
endif

datadir = ${prefix}/yaml/reactions/APPLgrid
data_DATA = ../yaml/parameters.yaml
