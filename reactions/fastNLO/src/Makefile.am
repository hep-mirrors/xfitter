
# Created by AddReaction.py on 2016-12-06

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include   -Wall -I$(top_builddir)/FastNLO/include -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libfastnlo_xfitter.la
libfastnlo_xfitter_la_SOURCES = ReactionfastNLO.cc speaker.cc fastNLOTools.cc fastNLOTable.cc fastNLOCoefficients.cc fastNLOCoeffBase.cc fastNLOCoeffData.cc fastNLOCoeffMult.cc fastNLOCoeffAddBase.cc fastNLOCoeffAddFix.cc fastNLOCoeffAddFlex.cc fastNLOReader.cc fastNLOPDFLinearCombinations.cc FastNLOInterface.cc FastNLOxFitter.cc
libfastnlo_xfitter_la_LDFLAGS = #-lxffnlo

# libfastnlo_xfitter_la_LDFLAGS = place_if_needed  

dist_noinst_HEADERS = ../include ../yaml
