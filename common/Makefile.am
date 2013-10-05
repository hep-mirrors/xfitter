AUTOMAKE_OPTIONS = foreign
pkginclude_HEADERS = FTNFitPars.h FitPars_base.h TextReader.h Xstring.h dbgtools.h fiasco.h
# The name of the library
noinst_LIBRARIES = libHFcommon.a

libHFcommon_a_SOURCES = FitPars_base.cc FTNFitPars.cc Xstring.cc
# libHFcommon_a_SOURCES += offset/g_offset.cc

AM_CXXFLAGS = -I$(srcdir)/linalg -fPIC -Wno-deprecated $(ROOT_CFLAGS) -w

