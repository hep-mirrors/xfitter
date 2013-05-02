AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libinterfaces.a

# explicit dependence:
hf_pdf_calls.o: ../../include/fcn.inc ../../include/ntot.inc ../../include/scales.inc ../../include/steering.inc

libinterfaces_a_SOURCES = hf_pdf_calls.f Cache_pdf_calls.cxx
AM_FFLAGS = -I$(srcdir)/../../include -ffixed-line-length-132
AM_CXXFLAGS = -I$(srcdir)/../include 
if ENABLE_OPENMP
AM_FFLAGS+= -fopenmp 
AM_CXXFLAGS+= -fopenmp
endif
