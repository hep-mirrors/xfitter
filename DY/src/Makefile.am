AUTOMAKE_OPTIONS = foreign

# The name of the library
noinst_LIBRARIES = libmydy.a

# Next files will be processed by g77 (and automatically gcc for *.c) compiler
libmydy_a_SOURCES = PhysPar.cc \
                    finterface.cc \
                    BinMatrix.cc IntSteps.cc PDFconv.cc DYcalc.cc
AM_FFLAGS = -I$(srcdir)/../include $(NOAUTOFCFLAG) -finit-local-zero -ffixed-line-length-132

AM_CXXFLAGS = -I$(srcdir)/../include -I$(srcdir)/../../interfaces/include  -O  -fPIC -Wno-deprecated -w
