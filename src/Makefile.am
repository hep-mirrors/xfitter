AUTOMAKE_OPTIONS = foreign
pkginclude_HEADERS = offset.h
# Explicit depence for all files:
*.o: ../include/*.inc 

# the main target, i.e. the executable to build
bin_PROGRAMS = xfitter

# all source files which go into the main target
xfitter_SOURCES = main.f hfbanner.cc init_theory.f fcn.f GetChisquare.f minuit_ini.f read_data.f store_output.f\
sumrules.f read_steer.f error_bands_pumplin.f pdf_param.f systematics.f\
evolution.f mc_errors.f ep_jets_fastnlo.f difftop_fastnlo.f dipole.f store_h1qcdfunc.f\
dis_sigma.f  dataset_tools.f theory_dispatcher.f dy_cc_sigma.f pp_jets_applgrid.f error_logging.f \
getncxskt.f Getgridkt.f GetDummyXsection.f \
GetDiffDisXsection.f FixModelParams.f  prep_corr.f GetPDFUncType.f rediagonalize.f qedevol.f \
ExtraParConstr.cc

# WS 2013-01-19 for the Offset method
xfitter_SOURCES += offset_fns.f g_offset.cc
# LHAPDF6 output
xfitter_SOURCES += lhapdf6_output.c c_interface.f

#LHAPDFErrors
xfitter_SOURCES += fortran_io.f lhapdferrors.cc

#Chi2scan
xfitter_SOURCES += chi2scan.cc

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../interfaces/include -I$(srcdir)/../FastNLO/include
AM_CXXFLAGS += -I$(srcdir)/../common -I$(srcdir)/../common/linalg
AM_CXXFLAGS += -I$(srcdir)/../tools/draw/include

if ENABLE_APFEL
  xfitter_SOURCES+= apfel.f
else 
  xfitter_SOURCES+= apfel_dum.f
endif

if ENABLE_MELA
  xfitter_SOURCES+= mela.f
else 
  xfitter_SOURCES+= mela_dum.f
endif

if ENABLE_APPLGRID 
  xfitter_SOURCES+= fappl_grid.cxx applgrids.f TheorEval.cc ftheor_eval.cc CommonGrid.cc
else 
  xfitter_SOURCES+= applgrids_dum.f ftheor_eval_dum.cc
endif

if ENABLE_LHAPDF 
else
  xfitter_SOURCES+= lhapdf_dum.f
endif

if ENABLE_HATHOR
  xfitter_SOURCES+= ttbar_hathor.f
else
  xfitter_SOURCES+= ttbar_hathor_dum.f	
endif

# FORTRAN compiler flags to be used for the main target
AM_FFLAGS = $(DEFS) -cpp -I$(srcdir)/../include -fno-automatic -I$(srcdir)/../minuit/src

# 
AM_CPPFLAGS = -I$(srcdir)/../include 

if ENABLE_APPLGRID 
  AM_CPPFLAGS+=$(APPLGRID_CPPFLAGS) $(ROOT_CFLAGS)
endif

if ENABLE_APFEL 
  AM_CPPFLAGS+=$(APFEL_CPPFLAGS)
  AM_FFLAGS+=$(APFEL_CPPFLAGS)
endif

# libraries to be linked to the main target
LDADD = $(QCDNUMLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) 

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

LDADD += ../Cascade/src/libcasutil.a
LDADD += ../common/libHFcommon.a

if ENABLE_UPDF
   LDADD += $(CASCADELIBS)  $(PYTHIALIBS) 
endif

if ENABLE_HATHOR
  LDADD+= ../Hathor/src/libHathor.a \
  $(HATHOR_ROOT)/libHathor.a
endif

# Link with HVQMNR library, otherwise with dummy stub
if ENABLE_HVQMNR
  LDADD += ../HVQMNR/src/libHVQMNR.a $(ROOT_LIBS)
else
  xfitter_SOURCES += hvqmnr_dum.f	
endif

LDADD+=  ../DY/src/libmydy.a\
 ../RT/src/libmyrt.a \
 ../ACOT/src/libmyacot.a \
 ../FONLL/src/libmyfonll.a \
 ../EW/src/libmyew.a \
 ../ABM/src/libmyabm.a \
 ../FastNLO/src/libFastNLO.a \
 ../DIPOLE/src/libmydipole.a \
 ../DiffDIS/src/libDiffDIS.a\
 ../interfaces/src/libinterfaces.a\
 ../minuit/src/libmyminuit.a $(CERNLIBS) @FLIBS@ \
  -lstdc++ \
 ../genetic/src/libgenetic.a \
 ../QEDevol/src/libQEDevol.a \
 ../common/num_utils/libnum_utils.a \
 ../tools/draw/src/libpdferrors.la


#GFRTLLIB = $(shell gfortran  -print-file-name=libgfortran.a )
#GFRTLIB  = -L$(subst /libgfortran.a, ,$(GFRTLLIB)) -lgfortran
GFRTLIB  = -lgfortran

if ENABLE_APPLGRID 
  LDADD+=$(APPLGRID_LDFLAGS) $(GFRTLIB)
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
