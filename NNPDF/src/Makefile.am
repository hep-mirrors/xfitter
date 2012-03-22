AUTOMAKE_OPTIONS = foreign

# The name of the library
noinst_LIBRARIES = libNNPDF.a

# Next files will be processed by g77 (and automatically gcc for *.c) compiler
libNNPDF_a_SOURCES = pdfs.cc utils.cc nnpdfrw.cc 

AM_CXXFLAGS = -I$(srcdir)/../include -Wall -fPIC -Wno-deprecated $(ROOT_CFLAGS)
