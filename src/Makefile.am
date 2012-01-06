AUTOMAKE_OPTIONS = foreign

# the main target, i.e. the executable to build
bin_PROGRAMS = FitPDF

# all source files which go into the main target
FitPDF_SOURCES = main.f init_theory.f fcn.f minuit_ini.f read_data.f store_output.f\
sumrules.f read_steer.f error_bands_pumplin.f pdf_param.f systematics.f\
evolution.f mc_errors.f ep_jets_fastnlo.f dipole.f\
  dis_nc_sigma.f dis_f2.f dis_cc_sigma.f dataset_tools.f theory_dispatcher.f dy_cc_sigma.f pp_jets_applgrid.f
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
if ENABLE_APPLGRID 
  AM_CPPFLAGS+=$(APPLGRID_CPPFLAGS) $(ROOT_CFLAGS)
endif

if ENABLE_LHAPDF 
  AM_CPPFLAGS+=$(LHAPDF_CPPFLAGS)
  AM_FFLAGS+=$(LHAPDF_CPPFLAGS)
endif

# libraries to be linked to the main target
LDADD = $(QCDNUMLIBS) \
 ../DY/src/libmydy.a\
 ../RT/src/libmyrt.a \
 ../ACOT/src/libmyacot.a \
../HS/src/libmyhs.a \
 ../FastNLO/src/libFastNLO.a \
 ../DIPOLE/src/libmydipole.a \
 ../minuit/src/libmyminuit.a $(CERNLIBS) @FLIBS@ \
-lgsl -lgslcblas -lstdc++

if ENABLE_HATHOR
  LDADD+= ../Hathor/src/libHathor.a \
  $(HATHOR_ROOT)/libHathor.a
endif

GFRTLLIB = $(shell gfortran  -print-file-name=libgfortran.a )
GFRTLIB  = -L $(subst /libgfortran.a, ,$(GFRTLLIB)) -lgfortran

#if ENABLE_APPLGRID 
  LDADD+=$(APPLGRID_LDFLAGS) $(GFRTLIB)
#endif
#if ENABLE_LHAPDF 
  LDADD+=$(LHAPDF_LDFLAGS)
#endif


LDADD+=$(ROOT_LIBS)
