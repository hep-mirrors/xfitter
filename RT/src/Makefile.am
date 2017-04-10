AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libmyrt.la

libmyrt_la_SOURCES = mstw2008.f mstw2008_wrap.f 
AM_FFLAGS = -I$(srcdir)/../../include -ffixed-line-length-132

if ENABLE_OPENMP
AM_FFLAGS+= -fopenmp
endif

