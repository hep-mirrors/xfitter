AUTOMAKE_OPTIONS = foreign

# the main target, i.e. the executable to build
bin_PROGRAMS = FitPDF

# all source files which go into the main target
FitPDF_SOURCES = main.f init_theory.f fcn.f GetChisquare.f minuit_ini.f read_data.f store_output.f\
sumrules.f read_steer.f error_bands_pumplin.f pdf_param.f systematics.f\
evolution.f mc_errors.f ep_jets_fastnlo.f dipole.f store_h1qcdfunc.f\
dis_sigma.f  dataset_tools.f theory_dispatcher.f dy_cc_sigma.f pp_jets_applgrid.f error_logging.f \
getncxskt.f Getgridkt.f \
GetDiffDisXsection.f FixModelParams.f nnpdfreweighting.f GetPointScaledErrors.f GetCovChisquare.f prep_corr.f

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../interfaces/include

# FTN_sock_DDIS.cc
  
if ENABLE_APPLGRID 
  FitPDF_SOURCES+= fappl_grid.cxx applgrids.f
else 
  FitPDF_SOURCES+= applgrids_dum.f 
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
AM_FFLAGS = -I$(srcdir)/../include -fno-automatic

# 
AM_CPPFLAGS = -I$(srcdir)/../include 
#	      -I../DiffDIS/include

if ENABLE_APPLGRID 
  AM_CPPFLAGS+=$(APPLGRID_CPPFLAGS) $(ROOT_CFLAGS)
endif


# libraries to be linked to the main target
LDADD = $(QCDNUMLIBS)

if ENABLE_LHAPDF 
  AM_CPPFLAGS+=$(LHAPDF_CPPFLAGS)
  AM_FFLAGS+=$(LHAPDF_CPPFLAGS)
endif

LDADD += ../Cascade/src/libcasutil.a
LDADD += ../NNPDF/src/libNNPDF.a

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
 ../HS/src/libmyhs.a \
 ../ABKM/src/libmyabkm.a \
 ../FastNLO/src/libFastNLO.a \
 ../DIPOLE/src/libmydipole.a \
 ../DiffDIS/src/libDiffDIS.a\
 ../interfaces/src/libinterfaces.a\
 ../minuit/src/libmyminuit.a $(CERNLIBS) @FLIBS@ \
  -lstdc++ 

GFRTLLIB = $(shell gfortran  -print-file-name=libgfortran.a )
GFRTLIB  = -L $(subst /libgfortran.a, ,$(GFRTLLIB)) -lgfortran

if ENABLE_APPLGRID 
  LDADD+=$(APPLGRID_LDFLAGS) $(GFRTLIB) $(ROOT_LIBS)
endif

if ENABLE_LHAPDF 
  LDADD+=$(LHAPDF_LDFLAGS)
endif

if ENABLE_NNPDFWEIGHT
  LDADD+=$(ROOT_LIBS)
endif
