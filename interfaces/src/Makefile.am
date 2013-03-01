AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libinterfaces.a

libinterfaces_a_SOURCES = hf_pdf_calls.f Cache_pdf_calls.cxx
AM_FFLAGS = -I$(srcdir)/../../include -ffixed-line-length-132
AM_CXXFLAGS = -I$(srcdir)/../include 
if ENABLE_OPENMP
AM_FFLAGS+= -fopenmp -frecursive
AM_CXXFLAGS+= -fopenmp
endif
