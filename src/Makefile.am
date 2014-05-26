AUTOMAKE_OPTIONS = foreign
pkginclude_HEADERS = offset.h
# Explicit depence for all files:
# *.o: ../include/*.inc 

# the main target, i.e. the executable to build
bin_PROGRAMS = FitPDF

# all source files which go into the main target
FitPDF_SOURCES = main.f init_theory.f fcn.f GetChisquare.f minuit_ini.f read_data.f store_output.f\
sumrules.f read_steer.f error_bands_pumplin.f pdf_param.f systematics.f\
evolution.f mc_errors.f ep_jets_fastnlo.f dipole.f store_h1qcdfunc.f\
dis_sigma.f  dataset_tools.f theory_dispatcher.f dy_cc_sigma.f pp_jets_applgrid.f error_logging.f \
getncxskt.f Getgridkt.f GetDummyXsection.f \
GetDiffDisXsection.f FixModelParams.f reweighting.f prep_corr.f GetPDFUncType.f 

# WS 2013-01-19 for the Offset method
FitPDF_SOURCES += offset_fns.f g_offset.cc
# LHAPDF6 output
FitPDF_SOURCES += lhapdf6_output.c c_interface.f

#LHAPDFErrors
FitPDF_SOURCES += fortran_io.f lhapdferrors.cc

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../interfaces/include 
AM_CXXFLAGS += -I$(srcdir)/../common -I$(srcdir)/../common/linalg

if ENABLE_APPLGRID 
  FitPDF_SOURCES+= fappl_grid.cxx applgrids.f TheorEval.cc ftheor_eval.cc
else 
  FitPDF_SOURCES+= applgrids_dum.f ftheor_eval_dum.cc
endif

if ENABLE_LHAPDF 
else
  FitPDF_SOURCES+= lhapdf_dum.f
endif

if ENABLE_HATHOR
  FitPDF_SOURCES+= ttbar_hathor.f
else
  FitPDF_SOURCES+= ttbar_hathor_dum.f	
endif

# FORTRAN compiler flags to be used for the main target
AM_FFLAGS = -I$(srcdir)/../include -fno-automatic -I$(srcdir)/../minuit/src

# 
AM_CPPFLAGS = -I$(srcdir)/../include 

if ENABLE_APPLGRID 
  AM_CPPFLAGS+=$(APPLGRID_CPPFLAGS) $(ROOT_CFLAGS)
endif


# libraries to be linked to the main target
LDADD = $(QCDNUMLIBS)
if ENABLE_OPENMP
  LDADD+= -lgomp
  AM_FFLAGS+= -fopenmp
endif

if ENABLE_LHAPDF 
  AM_CPPFLAGS+=$(LHAPDF_CPPFLAGS)
  AM_FFLAGS+=$(LHAPDF_CPPFLAGS)
endif

LDADD += ../Cascade/src/libcasutil.a
LDADD += ../NNPDF/src/libNNPDF.a
LDADD += ../common/libHFcommon.a

if ENABLE_UPDF
   LDADD += $(CASCADELIBS)  $(PYTHIALIBS) 
endif

if ENABLE_HATHOR
  LDADD+= ../Hathor/src/libHathor.a \
  $(HATHOR_ROOT)/libHathor.a
endif

LDADD+=  ../DY/src/libmydy.a\
 ../RT/src/libmyrt.a \
 ../ACOT/src/libmyacot.a \
 ../EW/src/libmyew.a \
 ../ABM/src/libmyabm.a \
 ../FastNLO/src/libFastNLO.a \
 ../DIPOLE/src/libmydipole.a \
 ../DiffDIS/src/libDiffDIS.a\
 ../interfaces/src/libinterfaces.a\
 ../minuit/src/libmyminuit.a $(CERNLIBS) @FLIBS@ \
  -lstdc++ \
 ../genetic/src/libgenetic.a

GFRTLLIB = $(shell gfortran  -print-file-name=libgfortran.a )
GFRTLIB  = -L $(subst /libgfortran.a, ,$(GFRTLLIB)) -lgfortran

if ENABLE_APPLGRID 
  LDADD+=$(APPLGRID_LDFLAGS) $(GFRTLIB) $(ROOT_LIBS)
endif

if ENABLE_GENETIC
  LDADD+= ../genetic/mixmax_r004/libmixmax.a
endif

if ENABLE_LHAPDF 
  LDADD+=$(LHAPDF_LDFLAGS)
endif

if ENABLE_NNPDFWEIGHT
  LDADD+=$(ROOT_LIBS)
endif
