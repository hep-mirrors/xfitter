AUTOMAKE_OPTIONS = foreign

# The name of the library
noinst_LIBRARIES = libFastNLO.a

# Next files will be processed by g77 (and automatically gcc for *.c) compiler
libFastNLO_a_SOURCES = FastNLOBlockB.cc FastNLOReader.cc FastNLOInterface.cc FastNLOHeraFitter.cc speaker.cc
#LHAPDF.cc Alphas.cc

AM_CXXFLAGS = -I$(srcdir)/../include -I$(srcdir)/../../interfaces/include -Wall -fPIC -Wno-deprecated $(ROOT_CFLAGS)
