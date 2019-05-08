
# Created by AddReaction.py on 2017-01-02

if HAVE_ROOT

AM_CXXFLAGS = $(ROOT_CFLAGS) -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated

lib_LTLIBRARIES = libbasehvqmnr_xfitter.la
libbasehvqmnr_xfitter_la_SOURCES = ReactionBaseHVQMNR.cc MNR.cc MNRFrag.cc MNRGrid.cc hvqcrsx.f

# libbasehvqmnr_xfitter_la_LDFLAGS = place_if_needed
libbasehvqmnr_xfitter_la_LDFLAGS = $(ROOT_LIBS)

endif

dist_noinst_HEADERS = ../include ../yaml
