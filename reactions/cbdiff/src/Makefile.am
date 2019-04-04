
# Created by AddReaction.py on 2019-02-01

if HAVE_ROOT

AM_CXXFLAGS = $(ROOT_CFLAGS) -I$(srcdir)/../../../reactions/BaseHVQMNR/include -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated

lib_LTLIBRARIES = libcbdiff_xfitter.la
libcbdiff_xfitter_la_SOURCES = Reactioncbdiff.cc

# libcbdiff_xfitter_la_LDFLAGS = place_if_needed
libcbdiff_xfitter_la_LDFLAGS = -lbasehvqmnr_xfitter -L$(srcdir)/../../../reactions/BaseHVQMNR/src/.libs

datadir = ${prefix}/yaml/reactions/cbdiff
data_DATA = ../yaml/parameters.yaml

endif

dist_noinst_HEADERS = ../include ../yaml

