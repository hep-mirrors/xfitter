# Created by AddReaction.py on 2017-03-28

if ENABLE_APPLGRID
  AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated 

  lib_LTLIBRARIES = libapplgrid_xfitter.la
  libapplgrid_xfitter_la_SOURCES = ReactionAPPLgrid.cc

# libapplgrid_xfitter_la_LDFLAGS = place_if_needed  
  libapplgrid_xfitter_la_CPPFLAGS = $(APPLGRID_CPPFLAGS) $(ROOT_CFLAGS)
endif

 
