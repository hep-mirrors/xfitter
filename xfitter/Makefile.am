AUTOMAKE_OPTIONS = foreign
#pkginclude_HEADERS = #offset.h
# Explicit depence for all files:
*.o: ../include/*.inc 

# the main target, i.e. the executable to build
bin_PROGRAMS = xfitter

xfitter_SOURCES = cmain.cxx

xfitter_LDFLAGS = -L${pkglibdir} -L${srcdir}/../src -L${top_builddir}/src -L${top_srcdir}/src -L${top_srcdir}/FastNLO/src  -L${top_srcdir}/RT/src -lxfmain -lxffnlo -lmyrt -lHFcommon -lmyew -lmyminuit -lgenetic -lQEDevol -lnum_utils -linterfaces
#xfitter_LDADD = 
# all source files which go into the main target

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../interfaces/include -I$(srcdir)/../FastNLO/include
AM_CXXFLAGS += -I$(srcdir)/../common -I$(srcdir)/../common/linalg
AM_CXXFLAGS += -I$(srcdir)/../tools/draw/include

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
  LDADD+= ../Hathor/src/libHathor.a \
  $(HATHOR_LIBS)
endif

# Link with HVQMNR library, otherwise with dummy stub
if ENABLE_HVQMNR
  LDADD += ../HVQMNR/src/libHVQMNR.a $(ROOT_LIBS)
#  libxfitter_a_SOURCES += hvqmnr_dum.f	need more
else
#  libxfitter_a_SOURCES += hvqmnr_dum.f
#  libxfmain_la_SOURCES += hvqmnr_dum.f	
endif

XFLIBS =  -lstdc++  ../tools/draw/src/pdferrors.o

#  ../common/libHFcommon.a \
#  ../EW/src/libmyew.a \
#  ../interfaces/src/libinterfaces.a\
#  ../minuit/src/libmyminuit.a $(CERNLIBS) @FLIBS@ \
# ../genetic/src/libgenetic.a \
# ../QEDevol/src/libQEDevol.a \
# ../common/num_utils/libnum_utils.a 
# ../FastNLO/src/libFastNLO.a 


LDADD += $(XFLIBS)

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
  LDADD+=$(APFEL_LDFLAGS)
endif

if ENABLE_MELA
  LDADD+=$(MELA_LDFLAGS)
endif

if ENABLE_GENETIC
  LDADD+= ../genetic/mixmax_r004/libmixmax.a
endif

