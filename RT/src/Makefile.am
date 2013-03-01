AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libmyrt.a

libmyrt_a_SOURCES = mstw2008.f mstw2008_wrap.f 
AM_FFLAGS = -I$(srcdir)/../../include -ffixed-line-length-132
#AM_FFLAGS = -I$(srcdir)/../../include -fno-automatic -finit-local-zero -ffixed-line-length-132
if ENABLE_OPENMP
AM_FFLAGS+= -fopenmp
endif

