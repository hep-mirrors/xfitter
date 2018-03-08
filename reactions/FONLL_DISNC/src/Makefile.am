
# Created by AddReaction.py on 2017-11-29

if ENABLE_APFEL
   AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include $(APFEL_CPPFLAGS) -I$(srcdir)/../../../reactions/BaseDISNC/include -Wall -fPIC -Wno-deprecated 

   lib_LTLIBRARIES = libfonll_disnc_xfitter.la
   libfonll_disnc_xfitter_la_SOURCES = ReactionFONLL_DISNC.cc

   libfonll_disnc_xfitter_la_LDFLAGS =  -lbasedisnc_xfitter -L$(top_srcdir)/reactions/BaseDISNC/src/.libs
endif

datadir = ${prefix}/yaml/reactions/FONLL_DISNC
data_DATA = ../yaml/parameters.yaml
dist_noinst_HEADERS = ../include ../yaml
