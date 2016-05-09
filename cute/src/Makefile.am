AUTOMAKE_OPTIONS = foreign

# The name of the library                                                                                                                                                         
noinst_LIBRARIES = libCuteInterface.a

# Next files will be processed by g77 (and automatically gcc for *.c) compiler                                                                                                     
libCuteInterface_a_SOURCES = CuteInterface.C

AM_CXXFLAGS = -I$(srcdir)/../CuTe-1.1.0/src -I$(srcdir)/../include -I$(srcdir)/../../include -std=c++11

AM_CPPFLAGS=$(LHAPDF_CPPFLAGS)
#AM_FFLAGS+=$(LHAPDF_CPPFLAGS)
