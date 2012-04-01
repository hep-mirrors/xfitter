AUTOMAKE_OPTIONS = foreign

# The name of the library
noinst_LIBRARIES = libNNPDF.a

if ENABLE_NNPDFWEIGHT

# Next files will be processed by g77 (and automatically gcc for *.c) compiler
  libNNPDF_a_SOURCES = pdfs.cc utils.cc nnpdfrw.cc 

   AM_CXXFLAGS = -I$(srcdir)/../include -Wall -fPIC -Wno-deprecated $(LHAPDF_CPPFLAGS) $(ROOT_CFLAGS)
else
  libNNPDF_a_SOURCES = nnpdf_dum.f
endif

