##########################################################################
#
#   File:         Makefile.am
#
#                 Copyright (C) M.Sutton (sutt@cern.ch) 
#
#   Description:  generate the Makefile for the applwrap
#                 
#
#   Created:    M.Sutton (sutt@cern.ch) Sun  Sun  9 Dec 2018 14:52:45 CET 
#   
#   $Id: makefile Sun 17 Aug 2025 15:32:39 BST
#                   
#########################################################################


APPLCFLAGS = $(shell applgrid-config --cxxflags)
APPLLIBS   = $(shell applgrid-config --ldflags)


# NB For some *stupid* reason, although hoppet requests these two functions, 
#    the linker doesn't find them on lxplus - becuase the STUPID GFORTRAN
#    installation is broken 
# FRTLIB   = $(shell gfortran -print-file-name=libgfortran.a)
# FRTLLIB  = -u _gfortran_runtime_error_at  -u _gfortran_os_error -L$(subst /libgfortran.a, ,$(FRTLIB) ) -lgfortran 
# FRTLLIB  = -L$(subst /libgfortran.a, ,$(FRTLIB) ) -lgfortran 
# gfortran not needed for the library
FRTLIB  = 
FRTLLIB = 



AM_CXXFLAGS = -O2 -Wall -fPIC  -I. -I.. -I$(srcdir) $(APPLCFLAGS) $(EIGEN_CXXFLAGS)
AM_CFLAGS   = -O2 -Wall -fPIC  -I. -I.. -I$(srcdir) $(APPLCFLAGS) $(EIGEN_CXXFLAGS)

AM_LDFLAGS = $(ROOTARCH) -O2  $(FCLIBS) $(APPLLIBS) $(LHAPDFLIB) $(FRTLLIB) 

AM_FCFLAGS = -c
AM_FCLIBS  = $(FRTLLIB)
AM_SOFLAGS = -shared



lib_LTLIBRARIES= libapplwrap.la libfapplwrap.la 
libapplwrap_la_SOURCES = \
	applwrap.cxx smooth.cxx \
	applwrap.h smooth.h solve.h fit.h 

libfapplwrap_la_SOURCES = \
	fapplwrap.cxx \
	applwrap.cxx smooth.cxx \
	applwrap.h smooth.h solve.h fit.h 


library_includedir=$(includedir)/

library_include_HEADERS=applwrap.h 


clean-local:
	rm -rf *.o *.lo 

.cxx.o : 
	$(CXX) $(AM_CXXFLAGS) -c $<
.c.o : 
	$(CC) $(AM_CFLAGS) -c $<


