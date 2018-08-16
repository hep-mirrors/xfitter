
# Created by AddEvolution.py on 2018-08-14

AM_CXXFLAGS = $(QCDNUM_CPPFLAGS) \
        -I$(srcdir)/../include  -I$(srcdir)/../../../include \
        -I$(srcdir)/../../BaseEvolution/include \
        -I$(srcdir)/../../../pdfdecompositions/BasePdfDecomposition/include \
        -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = libqcdnum_xfitter.la
libqcdnum_xfitter_la_SOURCES = EvolutionQCDNUM.cc

# libqcdnum_xfitter_la_LDFLAGS = place_if_needed  

datadir = ${prefix}/yaml/evolutions/QCDNUM
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
 
