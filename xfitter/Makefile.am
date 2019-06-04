AUTOMAKE_OPTIONS = foreign
#pkginclude_HEADERS = #offset.h
# Explicit depence for all files:
*.o: ../include/*.inc 

# the main target, i.e. the executable to build
bin_PROGRAMS = xfitter

xfitter_SOURCES = cmain.cxx

# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#xfitter_LDFLAGS = -L${srcdir}/../lib -lxfmain -lxffnlo -lmyrt -lHFcommon -lmyew ../minuit/src/.libs/minuit.o -lmyminuit -lgenetic -lQEDevol -lnum_utils
xfitter_LDFLAGS = -L${srcdir}/../lib -lxfmain -lHFcommon -lmyew ../minuit/src/.libs/minuit.o -lmyminuit -lgenetic -lQEDevol -lnum_utils ../tools/draw/src/pdferrors.o
#TODO move stat functions from pdferrors.cc to some place where they make sense

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../interfaces/include -I$(srcdir)/../FastNLO/include
AM_CXXFLAGS += -I$(srcdir)/../common -I$(srcdir)/../common/linalg
#AM_CXXFLAGS += -I$(srcdir)/../tools/draw/include

# FORTRAN compiler flags to be used for the main target
AM_FFLAGS = $(DEFS) -cpp -I$(srcdir)/../include -fno-automatic -I$(srcdir)/../minuit/src

# 
AM_CPPFLAGS = -I$(srcdir)/../include   -DPREFIX="\"$(prefix)\""

if ENABLE_APPLGRID 
  AM_CPPFLAGS+=$(APPLGRID_CPPFLAGS) $(ROOT_CFLAGS)
endif

if ENABLE_APFEL 
  AM_CPPFLAGS+=$(APFEL_CPPFLAGS)
  AM_FFLAGS+=$(APFEL_CPPFLAGS)
endif


if ENABLE_APFELGRID
  AM_CPPFLAGS+=$(APFELGRID_CPPFLAGS)
endif

# xFitter main library:
#noinst_LIBRARIES = libxfitter.a
#libxfitter_a_SOURCES = 
LDADD = #libxfitter.a
###################LDADD += libxfitter.a

# libraries to be linked to the main target
LDADD += $(QCDNUMLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) -lyaml-cpp


if ENABLE_LHAPDF 
  LDADD+=$(LHAPDF_LDFLAGS)
endif

if HAVE_ROOT
  LDADD+=$(ROOT_LIBS)
endif

# libraries to be linked to the main target
LDADD+= $(LIBS) $(FLIBS)

if ENABLE_OPENMP
  LDADD+= -lgomp
  AM_FFLAGS+= -fopenmp
endif

if ENABLE_LHAPDF 
  AM_CPPFLAGS+=$(LHAPDF_CPPFLAGS)
  AM_FFLAGS+=$(DEFS) $(LHAPDF_CPPFLAGS)
endif

# LDADD += ../Cascade/src/libcasutil.a
# LDADD += ../common/libHFcommon.a
# F2PYLIBS += ../Cascade/src/libcasutil.a
# F2PYLIBS += ../common/libHFcommon.a

if ENABLE_UPDF
   LDADD += $(CASCADELIBS)  $(PYTHIALIBS) 
endif

if ENABLE_HATHOR
# this is obsolete: code is in reactions/Hathor now
#  LDADD+= ../Hathor/src/libHathor.a \
  LDADD+= $(HATHOR_LIBS)
endif

# Link with HVQMNR library, otherwise with dummy stub
if ENABLE_HVQMNR
  LDADD += ../HVQMNR/src/libHVQMNR.a $(ROOT_LIBS)
#  libxfitter_a_SOURCES += hvqmnr_dum.f	need more
else
#  libxfitter_a_SOURCES += hvqmnr_dum.f
#  libxfmain_la_SOURCES += hvqmnr_dum.f	
endif

XFLIBS=-lstdc++

LDADD+=$(XFLIBS)

GFRTLIB  = -lgfortran

APPLGRID_1 = $(subst -m64,,$(APPLGRID_LDFLAGS))
APPLGRID_2 = $(subst -pthread,,$(APPLGRID_1))
APPLGRID_3 = $(subst -rdynamic,,$(APPLGRID_2))

if ENABLE_APPLGRID 
  LDADD+=$(APPLGRID_LDFLAGS) $(GFRTLIB)
#  libagtheory_la_LDFLAGS= $(APPLGRID_LDFLAGS) $(GFRTLIB)
endif

if ENABLE_APFELGRID 
  LDADD+=$(APFELGRID_LDFLAGS)
endif

if ENABLE_APFEL
  LDADD+=$(APFEL_LDFLAGS) -lmyfonll
endif

if ENABLE_MELA
  LDADD+=$(MELA_LDFLAGS)
endif

if ENABLE_GENETIC
  LDADD+= ../genetic/mixmax_r004/libmixmax.a
endif

