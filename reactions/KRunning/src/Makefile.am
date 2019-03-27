
# Created by AddReaction.py on 2019-01-16

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libkrunning_xfitter.la
libkrunning_xfitter_la_SOURCES = ReactionKRunning.cc

# libkrunning_xfitter_la_LDFLAGS = place_if_needed  

datadir = ${prefix}/yaml/reactions/KRunning
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
 