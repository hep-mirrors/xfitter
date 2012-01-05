AUTOMAKE_OPTIONS = foreign

# The name of the library
noinst_LIBRARIES = libHathor.a

# Next files will be processed by g77 (and automatically gcc for *.c) compiler
libHathor_a_SOURCES = H1FitterPdf.cc HathorInterface.cc

AM_CXXFLAGS = -I$(srcdir)/../include -I$(HATHOR_ROOT) -Wall -fPIC -Wno-deprecated $(ROOT_CFLAGS)
