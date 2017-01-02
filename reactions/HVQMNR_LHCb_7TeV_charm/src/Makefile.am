
# Created by AddReaction.py on 2017-01-02

AM_CXXFLAGS = $(ROOT_CFLAGS) -I$(srcdir)/../../../reactions/BaseHVQMNR/include -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libhvqmnr_lhcb_7tev_charm_xfitter.la
libhvqmnr_lhcb_7tev_charm_xfitter_la_SOURCES = ReactionHVQMNR_LHCb_7TeV_charm.cc

libhvqmnr_lhcb_7tev_charm_xfitter_la_LDFLAGS = -lbasehvqmnr_xfitter -L$(srcdir)/../../../reactions/BaseHVQMNR/src/.libs

 
