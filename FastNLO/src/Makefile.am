AUTOMAKE_OPTIONS = foreign

# The name of the library
#noinst_LIBRARIES = #libFastNLO.a

# Next files will be processed by g77 (and automatically gcc for *.c) compiler
#libFastNLO_a_SOURCES= # speaker.cc fastNLOTools.cc fastNLOBase.cc fastNLOTable.cc fastNLOCoefficients.cc fastNLOCoeffBase.cc fastNLOCoeffData.cc fastNLOCoeffMult.cc fastNLOCoeffAddBase.cc fastNLOCoeffAddFix.cc fastNLOCoeffAddFlex.cc fastNLOReader.cc fastNLOPDFLinearCombinations.cc FastNLOInterface.cc FastNLOxFitter.cc 
#LHAPDF.cc Alphas.cc

#AM_CXXFLAGS = -I$(srcdir)/../include -I$(srcdir)/../include/fastnlotk -I$(srcdir)/../../interfaces/include -Wall -fPIC -Wno-deprecated $(ROOT_CFLAGS)

# shared library for 'theory reacton interface'
lib_LTLIBRARIES = libxffnlo.la libxfrtfnlo.la
libxffnlo_la_SOURCES= speaker.cc fastNLOTools.cc fastNLOBase.cc fastNLOTable.cc fastNLOCoefficients.cc fastNLOCoeffBase.cc fastNLOCoeffData.cc fastNLOCoeffMult.cc fastNLOCoeffAddBase.cc fastNLOCoeffAddFix.cc fastNLOCoeffAddFlex.cc fastNLOReader.cc fastNLOPDFLinearCombinations.cc FastNLOInterface.cc FastNLOxFitter.cc #LHAPDF.cc Alphas.cc $(libFastNLO_a_SOURCES)
libxffnlo_la_CFLAGS=
libxffnlo_la_CPPFLAGS=
libxffnlo_la_CXXFLAGS=-I$(srcdir)/../include -I$(srcdir)/../include/fastnlotk -I$(srcdir)/../../interfaces/include -Wall #$(AM_CXXFLAGS)
libxffnlo_la_LDFLAGS= #-shared

# shared library for 'reaction theory'
# 'implementation' of ReactionFastNLO
libxfrtfnlo_la_SOURCES= ReactionFastNLO.cc
libxfrtfnlo_la_CFLAGS=
libxfrtfnlo_la_CPPFLAGS=
libxfrtfnlo_la_CXXFLAGS=-I$(srcdir)/../include -I$(srcdir)/../include/fastnlotk -I$(top_builddir)/include
libxfrtfnlo_la_LDFLAGS= -lxffnlo

