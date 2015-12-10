AUTOMAKE_OPTIONS = foreign

# The name of the library
noinst_LIBRARIES = libFastNLO.a

# Next files will be processed by g77 (and automatically gcc for *.c) compiler
libFastNLO_a_SOURCES =  speaker.cc fastNLOTools.cc fastNLOBase.cc fastNLOTable.cc fastNLOCoefficients.cc fastNLOCoeffBase.cc fastNLOCoeffData.cc fastNLOCoeffMult.cc fastNLOCoeffAddBase.cc fastNLOCoeffAddFix.cc fastNLOCoeffAddFlex.cc fastNLOReader.cc fastNLOPDFLinearCombinations.cc FastNLOInterface.cc FastNLOxFitter.cc 
#LHAPDF.cc Alphas.cc

AM_CXXFLAGS = -I$(srcdir)/../include -I$(srcdir)/../include/fastnlotk -I$(srcdir)/../../interfaces/include -Wall -fPIC -Wno-deprecated $(ROOT_CFLAGS)
