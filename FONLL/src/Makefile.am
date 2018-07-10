AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libmyfonll.la

# Force recompilation if include files change
*.o: ../../include/couplings.inc

libmyfonll_la_SOURCES = ext_pdfs.cc

AM_CXXFLAGS = -I$(srcdir)/../../include 

dist_noinst_HEADERS = ext_pdfs.h
