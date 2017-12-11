
# Created by AddReaction.py on 2017-11-29
if ENABLE_APFEL
   AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include $(APFEL_CPPFLAGS) -I$(srcdir)/../../../reactions/BaseDISCC/include -Wall -fPIC -Wno-deprecated 

   lib_LTLIBRARIES = libfonll_discc_xfitter.la
   libfonll_discc_xfitter_la_SOURCES = ReactionFONLL_DISCC.cc

   libfonll_discc_xfitter_la_LDFLAGS =  -lbasediscc_xfitter -L$(top_srcdir)/reactions/BaseDISCC/src/.libs
endif

datadir = ${prefix}/yaml/reactions/FONLL_DISCC
data_DATA = ../yaml/parameters.yaml
