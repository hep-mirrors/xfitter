AUTOMAKE_OPTIONS = foreign

# The name of the library
noinst_LIBRARIES = libFastNLO.a

# Next files will be processed by g77 (and automatically gcc for *.c) compiler
libFastNLO_a_SOURCES =  speaker.cc fastNLOTools.cc fastNLOInterpolBase.cc fastNLOInterpolCatmullRom.cc fastNLOInterpolLagrange.cc fastNLOInterpolLinear.cc fastNLOInterpolOneNode.cc fastNLOBase.cc fastNLOTable.cc fastNLOCoefficients.cc fastNLOCreate.cc read_steer.cc fastNLOCoeffBase.cc fastNLOCoeffData.cc fastNLOCoeffMult.cc fastNLOCoeffAddBase.cc fastNLOCoeffAddFix.cc fastNLOCoeffAddFlex.cc fastNLOReader.cc fastNLOPDFLinearCombinations.cc Alphas.cc CRunDec.cc FastNLOInterface.cc FastNLOHeraFitter.cc 
#LHAPDF.cc Alphas.cc

AM_CXXFLAGS = -I$(srcdir)/../include -I$(srcdir)/../../interfaces/include -Wall -fPIC -Wno-deprecated $(ROOT_CFLAGS)
