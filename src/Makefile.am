AUTOMAKE_OPTIONS = foreign
dist_noinst_HEADERS = offset.h
# Explicit depence for all files:
*.o: ../include/*.inc

# the main shared library
lib_LTLIBRARIES = libxfmain.la

libxfmain_la_SOURCES =  hfbanner.cc init_theory.f fcn.f GetChisquare.f minuit_ini.f read_data.f store_output.f\
read_steer.f error_bands_pumplin.f systematics.f\
mc_errors.f\
dataset_tools.f theory_dispatcher.f error_logging.f\
GetDummyXsection.f main.f\
prep_corr.f GetPDFUncType.f rediagonalize.f\
ExtraParConstr.cc TheorEval.cc ftheor_eval.cc  xfitter_cpp_base.cc xfitter_pars.cc xfitter_steer.cc Profiler.cc\
tinyexpr.c expression_utils.cc dependent_pars.cc\
fortran_interface.cc\
offset_fns.f g_offset.cc\
c_interface.f\
lhapdf6_output.cc\
dis_sigma_ht.f\
fortran_io.f lhapdferrors.cc\
chi2scan.cc\
ReactionTheory.cc\
TermData.cc\
QCDNUM_Manager.cc\
IntegrateDIS.cc
# optional sources
# these should be moved to indiviual shared libraries in the future
##APFEL and MELA commented out because they need to be updated for the new evolutions interface
#if ENABLE_APFEL
#  libxfmain_la_SOURCES+= apfel.f
#else
#  libxfmain_la_SOURCES+= apfel_dum.f
#endif
#if ENABLE_MELA
#  libxfmain_la_SOURCES+= mela.f
#else
#  libxfmain_la_SOURCES+= mela_dum.f
#endif

#Broken during "TermData" rewrite
#difftop_fastnlo.f
#store_h1qcdfunc.f

#
# flags for libfxmain
#
libxfmain_la_CXXFLAGS=\
  -I$(top_srcdir)/interfaces/include\
  -I$(top_srcdir)/FastNLO/include\
  -I$(top_srcdir)/common\
  -I$(top_srcdir)/common/linalg\
  -I$(top_srcdir)/tools/draw/include
libxfmain_la_CFLAGS=
libxfmain_la_CPPFLAGS=$(QCDNUM_CPPFLAGS) -I$(top_srcdir)/include -DPREFIX="\"$(prefix)\""
libxfmain_la_FFLAGS=$(DEFS) -cpp -I$(top_srcdir)/include  $(NOAUTOFCFLAG) -I$(top_srcdir)/minuit/src
libxfmain_la_LDFLAGS=-ldl

# xFitter main library:
LDADD = libxfitter.a
# libraries to be linked to the main target
LDADD += $(QCDNUMLIBS) $(LAPACK_LIBS) $(BLAS_LIBS)

if ENABLE_LHAPDF
  LDADD+=$(LHAPDF_LDFLAGS)
  libxfmain_la_CPPFLAGS+=$(LHAPDF_CPPFLAGS)
  libxfmain_la_FFLAGS+=$(DEFS) $(LHAPDF_CPPFLAGS)
  libxfmain_la_SOURCES+=CheckForPDF.cxx
else
  libxfmain_la_SOURCES+= lhapdf_dum.f
endif

if HAVE_ROOT
  LDADD+=$(ROOT_LIBS)
  libxfmain_la_CXXFLAGS+=$(ROOT_CFLAGS)
endif

# libraries to be linked to the main target
LDADD+=$(LIBS) $(FLIBS)

if ENABLE_OPENMP
  LDADD+= -lgomp
  libxfmain_la_FFLAGS+= -fopenmp
endif

XFLIBS=../EW/src/libmyew.a\
 ../interfaces/src/libinterfaces.a\
 ../minuit/src/libmyminuit.a $(CERNLIBS) @FLIBS@ \
 ../genetic/src/libgenetic.a\
 ../QEDevol/src/libQEDevol.a\
 ../common/num_utils/libnum_utils.a\
 ../tools/draw/src/pdferrors.o

LDADD+=$(XFLIBS)

GFRTLIB=-lgfortran

APPLGRID_1 = $(subst -m64,,$(APPLGRID_LDFLAGS))
APPLGRID_2 = $(subst -pthread,,$(APPLGRID_1))
APPLGRID_3 = $(subst -rdynamic,,$(APPLGRID_2))

#if ENABLE_MELA
#  LDADD+=$(MELA_LDFLAGS)
#endif

if ENABLE_GENETIC
  LDADD+= ../genetic/mixmax_r004/libmixmax.a
endif
