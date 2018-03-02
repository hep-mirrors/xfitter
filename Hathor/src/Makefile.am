AUTOMAKE_OPTIONS = foreign

# The name of the library
noinst_LTLIBRARIES = libHathor.la

# Next files will be processed by g77 (and automatically gcc for *.c) compiler
libHathor_la_SOURCES = xFitterPdf.cc HathorInterface.cc

# HATHOR_ROOT points to the directory containing libHathor.a:
#   (1) in Hathor-1.5 this is the same common directory which contains all headers and sources,
#   (2) in Hathor-2.0 the library is Hathor-2.0/lib while the headers are in Hathor-2.0/include,
# therefore add both header places below
AM_CXXFLAGS = -I$(srcdir)/../include -I$(HATHOR_ROOT) -I$(HATHOR_ROOT)/../include $(LHAPDF_CPPFLAGS) -I$(srcdir)/../../interfaces/include -Wall -fPIC -Wno-deprecated $(ROOT_CFLAGS)
