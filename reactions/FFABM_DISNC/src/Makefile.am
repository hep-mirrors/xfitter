
# Created by AddReaction.py on 2017-09-29

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -I$(srcdir)/../../../reactions/BaseDISNC/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libffabm_disnc_xfitter.la
libffabm_disnc_xfitter_la_SOURCES = ReactionFFABM_DISNC.cc

libffabm_disnc_xfitter_la_LDFLAGS =  -lbasedisnc_xfitter -L$(top_srcdir)/reactions/BaseDISNC/src/.libs   

dist_noinst_HEADERS = ../include ../yaml

