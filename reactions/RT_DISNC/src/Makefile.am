
# Created by AddReaction.py on 2017-04-10

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -I$(srcdir)/../../../reactions/BaseDISNC/include  -I$(top_srcdir)/RT/src/ -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = librt_disnc_xfitter.la
librt_disnc_xfitter_la_SOURCES = ReactionRT_DISNC.cc

librt_disnc_xfitter_la_LDFLAGS =  -lbasedisnc_xfitter -lmyrt -L$(top_srcdir)/reactions/BaseDISNC/src/.libs    -L$(top_srcdir)/RT/src/.libs/

datadir = ${prefix}/yaml/reactions/RT_DISNC
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
