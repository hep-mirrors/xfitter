AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libmyrt.la
libmyrt_la_SOURCES = mstw2008.f mstw2008_wrap.f ext_pdfs.cc ext_pdfs_qcdnum.cc

AM_FFLAGS = -I$(srcdir)/../../include -ffixed-line-length-132
AM_CXXFLAGS = -I$(srcdir)/../../include -I$(srcdir)/../../interfaces/include

if ENABLE_OPENMP
AM_FFLAGS+= -fopenmp
endif

dist_noinst_HEADERS = ext_pdfs.h
