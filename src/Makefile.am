AUTOMAKE_OPTIONS = foreign

# the main target, i.e. the executable to build
bin_PROGRAMS = FitPDF

# all source files which go into the main target
FitPDF_SOURCES = main.f init_theory.f fcn.f minuit_ini.f qcdnum_ini.f covarmatrix.f fail_cuts.f read_covarmatrix.f read_data.f store_pdfs.f\
sumrules.f read_steer.f error_bands_pumplin.f formff.f pdf_param.f systematics.f\
get_partons.f ssdint.f writefittedpoints.f\
alnorm.f logshift.f zeroth.f evolution.f func1.f func22.f func24.f pdflength.f\
  dis_nc_sigma.f dis_cc_sigma.f dataset_tools.f theory_dispatcher.f dy_cc_sigma.f  theory_init.f pp_jets_applgrid.f
if ENABLE_APPLGRID 
  FitPDF_SOURCES+= fappl_grid.cxx applgrids.f
else 
  FitPDF_SOURCES+= applgrids_dum.f 
endif

# FORTRAN compiler flags to be used for the main target
AM_FFLAGS = -I$(srcdir)/../include -fno-automatic

# 
AM_CPPFLAGS = -I$(srcdir)/../include 
if ENABLE_APPLGRID 
  AM_CPPFLAGS+=$(APPLGRID_CPPFLAGS) $(ROOT_CFLAGS)
endif

# libraries to be linked to the main target
LDADD = $(QCDNUMLIBS) \
 ../DY/src/libmydy.a\
 ../RT/src/libmyrt.a \
 ../minuit/src/libmyminuit.a $(CERNLIBS) @FLIBS@ \
-lgsl -lgslcblas -lstdc++

GFRTLLIB = $(shell gfortran  -print-file-name=libgfortran.a )
GFRTLIB  = -L $(subst /libgfortran.a, ,$(GFRTLLIB)) -lgfortran

#if ENABLE_APPLGRID 
  LDADD+=$(APPLGRID_LDFLAGS) $(GFRTLIB)
#endif


