AUTOMAKE_OPTIONS = foreign
dist_noinst_HEADERS = offset.h
# Explicit depence for all files:
*.o: ../include/*.inc 




# the main shared library
lib_LTLIBRARIES = libxfmain.la

libxfmain_la_SOURCES =  hfbanner.cc init_theory.f fcn.f GetChisquare.f minuit_ini.f read_data.f store_output.f\
sumrules.f read_steer.f error_bands_pumplin.f pdf_param.f systematics.f\
evolution.f mc_errors.f ep_jets_fastnlo.f difftop_fastnlo.f dipole.f store_h1qcdfunc.f\
dis_sigma.f  dataset_tools.f theory_dispatcher.f dy_cc_sigma.f pp_jets_applgrid.f error_logging.f \
getncxskt.f Getgridkt.f GetDummyXsection.f main.f \
GetDiffDisXsection.f FixModelParams.f  prep_corr.f GetPDFUncType.f rediagonalize.f qedevol.f \
ExtraParConstr.cc CommonGrid.cc TheorEval.cc ftheor_eval.cc  xfitter_cpp_base.cc xfitter_pars.cc

# WS 2013-01-19 for the Offset method
libxfmain_la_SOURCES += offset_fns.f g_offset.cc
# LHAPDF6 output
libxfmain_la_SOURCES += lhapdf6_output.c c_interface.f
# WS 2016-01-13 DIS higher twist contributions
libxfmain_la_SOURCES += dis_sigma_ht.f
#LHAPDFErrors
libxfmain_la_SOURCES += fortran_io.f lhapdferrors.cc
#Chi2scan
libxfmain_la_SOURCES += chi2scan.cc
# ReactionTheory
libxfmain_la_SOURCES += ReactionTheory.cc
# optional sources
# these should be moved to indiviual shared libraries in the future
if ENABLE_APFEL
  libxfmain_la_SOURCES+= apfel.f
else
  libxfmain_la_SOURCES+= apfel_dum.f
endif
if ENABLE_MELA
  libxfmain_la_SOURCES+= mela.f
else
  libxfmain_la_SOURCES+= mela_dum.f
endif
if ENABLE_APPLGRID
  libxfmain_la_SOURCES+= fappl_grid.cxx applgrids.f
else
  libxfmain_la_SOURCES+= applgrids_dum.f
endif


if ENABLE_APFELGRID
  libxfmain_la_SOURCES+= APFELgridGeneration.cc
endif

if ENABLE_LHAPDF 
  libxfmain_la_SOURCES+= CheckForPDF.cxx
else
  libxfmain_la_SOURCES+= lhapdf_dum.f 
endif
if ENABLE_HATHOR
  libxfmain_la_SOURCES+= ttbar_hathor.f
else
  libxfmain_la_SOURCES+= ttbar_hathor_dum.f
endif

#
# flags for libfxmain
#
libxfmain_la_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../interfaces/include -I$(srcdir)/../FastNLO/include 
libxfmain_la_CXXFLAGS += -I$(srcdir)/../common -I$(srcdir)/../common/linalg
libxfmain_la_CXXFLAGS += -I$(srcdir)/../tools/draw/include
libxfmain_la_CFLAGS=
libxfmain_la_CPPFLAGS= -I$(srcdir)/../include -DPREFIX="\"$(prefix)\""
libxfmain_la_FFLAGS= $(DEFS) -cpp -I$(srcdir)/../include  $(NOAUTOFCFLAG)  -I$(srcdir)/../minuit/src
if ENABLE_LHAPDF
  libxfmain_la_CPPFLAGS+=$(LHAPDF_CPPFLAGS)
  libxfmain_la_FFLAGS+=$(DEFS) $(LHAPDF_CPPFLAGS)
endif
# libraries, where this library depends upon
libxfmain_la_LDFLAGS = -ldl
# libxfmain_la_LDFLAGS = -ldl  -L./ -lfortran_thunking -L/usr/local/cuda/lib64 -lcudart -lcublas

# FORTRAN compiler flags to be used for the main target
AM_FFLAGS = $(DEFS) -cpp -I$(srcdir)/../include $(NOAUTOFCFLAG) -I$(srcdir)/../minuit/src

# 
AM_CPPFLAGS = -I$(srcdir)/../include 

if ENABLE_APPLGRID 
  libxfmain_la_CPPFLAGS+=$(APPLGRID_CPPFLAGS) $(ROOT_CFLAGS)
endif

if ENABLE_APFEL 
  libxfmain_la_CPPFLAGS+=$(APFEL_CPPFLAGS)
  libxfmain_la_FFLAGS+=$(APFEL_CPPFLAGS)
endif


if ENABLE_APFELGRID
  libxfmain_la_CPPFLAGS+=$(APFELGRID_CPPFLAGS)
endif

# xFitter main library:
LDADD = libxfmain.la

# libraries to be linked to the main target
LDADD += $(QCDNUMLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) 

# separate list for python: 
F2PYLIBS = libxfmain.la 
F2PYLIBS +=  $(QCDNUMLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(LIBS) $(FLIBS)

if ENABLE_LHAPDF 
  LDADD+=$(LHAPDF_LDFLAGS)
  F2PYLIBS +=  $(LHAPDF_LDFLAGS)
endif

if HAVE_ROOT
  LDADD+=$(ROOT_LIBS)
endif

# libraries to be linked to the main target
LDADD+= $(LIBS) $(FLIBS)

if ENABLE_OPENMP
  LDADD+= -lgomp
  libxfmain_la_FFLAGS+= -fopenmp
endif

if ENABLE_LHAPDF 
  libxfmain_la_CPPFLAGS+=$(LHAPDF_CPPFLAGS)
  libxfmain_la_FFLAGS+=$(DEFS) $(LHAPDF_CPPFLAGS)
endif

# LDADD += ../Cascade/src/libcasutil.la
# LDADD += ../common/libHFcommon.la
# F2PYLIBS += ../Cascade/src/libcasutil.la
# F2PYLIBS += ../common/libHFcommon.la

if ENABLE_UPDF
   LDADD += $(CASCADELIBS)  $(PYTHIALIBS) 
   F2PYLIBS +=  $(CASCADELIBS)  $(PYTHIALIBS) 
endif

if ENABLE_HATHOR
  LDADD+= ../Hathor/src/libHathor.la \
  $(HATHOR_LIBS)
endif

# Link with HVQMNR library, otherwise with dummy stub
if ENABLE_HVQMNR
  LDADD += ../HVQMNR/src/libHVQMNR.la $(ROOT_LIBS)
else
  libxfmain_la_SOURCES += hvqmnr_dum.f	
endif

XFLIBS =  ../DY/src/libmydy.la\
 ../Cascade/src/libcasutil.la \
 ../common/libHFcommon.la \
 ../RT/src/libmyrt.la \
 ../ACOT/src/libmyacot.la \
 ../SACOT/src/libmysacot.la \
 ../FONLL/src/libmyfonll.la \
 ../EW/src/libmyew.la \
 ../ABM/src/libmyabm.la \
 ../FastNLO/src/libxffnlo.la \
 ../DIPOLE/src/libmydipole.la \
 ../DiffDIS/src/libDiffDIS.a\
 ../interfaces/src/libinterfaces.la\
 ../minuit/src/libmyminuit.la $(CERNLIBS) @FLIBS@ \
  -lstdc++ \
 ../genetic/src/libgenetic.la \
 ../QEDevol/src/libQEDevol.la \
 ../common/num_utils/libnum_utils.la \
 ../tools/draw/src/libpdferrors.la


LDADD += $(XFLIBS)
F2PYLIBS += $(XFLIBS)

GFRTLIB  = -lgfortran

APPLGRID_1 = $(subst -m64,,$(APPLGRID_LDFLAGS))
APPLGRID_2 = $(subst -pthread,,$(APPLGRID_1))
APPLGRID_3 = $(subst -rdynamic,,$(APPLGRID_2))

if ENABLE_APPLGRID 
  LDADD+=$(APPLGRID_LDFLAGS) $(GFRTLIB)
#  libagtheory_la_LDFLAGS= $(APPLGRID_LDFLAGS) $(GFRTLIB)
  F2PYLIBS +=  $(APPLGRID_3)
endif

if ENABLE_APFELGRID 
  LDADD+=$(APFELGRID_LDFLAGS)
  F2PYLIBS+=$(APFELGRID_LDFLAGS) 
endif

if ENABLE_APFEL
  LDADD+=$(APFEL_LDFLAGS)
  F2PYLIBS +=  $(APFEL_LDFLAGS)
endif

if ENABLE_MELA
  LDADD+=$(MELA_LDFLAGS)
  F2PYLIBS +=  $(MELA_LDFLAGS)
endif

if ENABLE_GENETIC
  LDADD+= ../genetic/mixmax_r004/libmixmax.la
endif


if HAVE_BOOST
   lib_LTLIBRARIES += libxfitter_fit.la
   libxfitter_fit_la_SOURCES = xFitter_fit.cc
   libxfitter_fit_la_CPPFLAGS = $(BOOST_PYTHON_CPPFLAGS) 
   libxfitter_fit_la_LDFLAGS = $(BOOST_PYTHON_LIBS)  $(F2PYLIBS) 
endif

