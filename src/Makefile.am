AUTOMAKE_OPTIONS = foreign

# the main target, i.e. the executable to build
bin_PROGRAMS = FitPDF

# source code needed for the fastNLO interface, to fit jets
#FASTNLOSOURCES = fn-interface.f fnh1001.f H1fnh1001.f fnh1003.f H1fnh1003.f fnh1004.f H1fnh1004.f fnh2003.f fnh2003_DIS.f fnh2003_LQ2H1.f H1fnh2003.f

# all source files which go into the main target
FitPDF_SOURCES = main.f fcn.f minuit_ini.f qcdnum_ini.f covarmatrix.f fail_cuts.f read_covarmatrix.f read_data.f store_pdfs.f\
sumrules.f formff.f read_steer.f error_bands_pumplin.f pdf_param.f systematics.f\
get_partons.f ssdint.f writefittedpoints.f\
alnorm.f logshift.f zeroth.f evolution.f func1.f func22.f func24.f pdflength.f\
  dis_nc_sigma.f dis_cc_sigma.f dataset_tools.f theory_dispatcher.f dy_cc_sigma.f unresolved.f theory_init.f pp_jets_applgrid.f
if ENABLE_APPLGRID 
  FitPDF_SOURCES+= fappl_grid.cxx applgrids.f
else 
  FitPDF_SOURCES+= applgrids_dum.f 
endif
#$(FASTNLOSOURCES)

# FORTRAN compiler flags to be used for the main target
AM_FFLAGS = -I$(srcdir)/../include -fno-automatic
#AM_FFLAGS = -I$(srcdir)/../include -fno-automatic -finit-local-zero -ffixed-line-length-132

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
#GFRTLIB  = -L/usr/lib/gcc/x86_64-redhat-linux6E/4.4.0 -lgfortran

#if ENABLE_APPLGRID 
  LDADD+=$(APPLGRID_LDFLAGS) $(GFRTLIB)
#  LDADD+=$(APPLGRID_LDFLAGS) -L/usr/lib/gcc/i386-redhat-linux5E/4.1.2 -lgfortran
#endif


