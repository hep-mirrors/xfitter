
# Created by AddReaction.py on 2017-10-28

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libkfactor_xfitter.la
libkfactor_xfitter_la_SOURCES = ReactionKFactor.cc

# libkfactor_xfitter_la_LDFLAGS = place_if_needed  

datadir = ${prefix}/yaml/reactions/KFactor
data_DATA = ../yaml/parameters.yaml

 