#bin_PROGRAMS = CuTe
#CuTe_SOURCES = main.cpp ConfigFile.cpp apkernels.h polynom.h fixed.h lhafunc.h resumed.h ConfigFile.h gPDF.h vecmath.h Hardfunc.h param.h crosssection.h  konst.h plot.h

AUTOMAKE_OPTIONS = foreign

# The name of the library                                                                                                                                                         
noinst_LIBRARIES = libcute.a

libcute_a_SOURCES = ConfigFile.cpp plot.C lhafunc.C param.C

#AM_CXXFLAGS = -I$(srcdir)/../CuTe-1.1.0/src -I$(srcdir)/../include
AM_CXXFLAGS = -I$(srcdir)/../../../include

AM_CPPFLAGS=$(LHAPDF_CPPFLAGS) -std=c++11
#AM_FFLAGS+=$(LHAPDF_CPPFLAGS)


#CuTe_CXXFLAGS = $(INC_LHAPDF) $(INC_GSL)
#CuTe_LDADD = $(LIB_LHAPDF) $(LIB_GSL)
#CuTe_LDFLAGS = $(LDD_LHAPDF) $(LDD_GSL)


libcute_a_CXXFLAGS = $(INC_LHAPDF) $(INC_GSL) -I$(srcdir)/../../../include
