AUTOMAKE_OPTIONS = foreign

# The name of the library
noinst_LTLIBRARIES = libHathor.la

# Next files will be processed by g77 (and automatically gcc for *.c) compiler
libHathor_la_SOURCES = xFitterPdf.cc HathorInterface.cc

AM_CXXFLAGS = -I$(srcdir)/../include -I$(HATHOR_ROOT) -I$(srcdir)/../../interfaces/include -Wall -fPIC -Wno-deprecated $(ROOT_CFLAGS)
