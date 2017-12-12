AUTOMAKE_OPTIONS = foreign
dist_noinst_HEADERS = offset.h
# Explicit depence for all files:
*.o: ../include/*.inc 

# # the main target, i.e. the executable to build
# bin_PROGRAMS = xfitter
# noinst_LIBRARIES = libxfitter.a
# 
# xfitter_SOURCES = main.f
# xfitter_LDFLAGS = -L${pkglibdir} -L${srcdir}/../src -L${top_builddir}/src -L${top_srcdir}/src -L${top_srcdir}/FastNLO/src -lxfmain -lxffnlo
# #xfitter_LDADD = 
# # all source files which go into the main target
# libxfitter_a_SOURCES = 
# 


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
  libxfitter_a_SOURCES+= CheckForPDF.cxx
else
  libxfitter_a_SOURCES+= lhapdf_dum.f 
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
LDADD = libxfitter.a
###################LDADD += libxfitter.a

# libraries to be linked to the main target
LDADD += $(QCDNUMLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) 

# separate list for python: 
F2PYLIBS = 
###########F2PYLIBS += libxfitter.a 
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

# LDADD += ../Cascade/src/libcasutil.a
# LDADD += ../common/libHFcommon.a
# F2PYLIBS += ../Cascade/src/libcasutil.a
# F2PYLIBS += ../common/libHFcommon.a

if ENABLE_UPDF
   LDADD += $(CASCADELIBS)  $(PYTHIALIBS) 
   F2PYLIBS +=  $(CASCADELIBS)  $(PYTHIALIBS) 
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
  libxfmain_la_SOURCES += hvqmnr_dum.f	
endif

XFLIBS =  ../DY/src/libmydy.a\
 ../Cascade/src/libcasutil.a \
 ../common/libHFcommon.a \
 ../RT/src/libmyrt.a \
 ../ACOT/src/libmyacot.a \
 ../SACOT/src/libmysacot.a \
 ../FONLL/src/libmyfonll.a \
 ../EW/src/libmyew.a \
 ../ABM/src/libmyabm.a \
 ../FastNLO/src/libFastNLO.a \
 ../DIPOLE/src/libmydipole.a \
 ../DiffDIS/src/libDiffDIS.a\
 ../interfaces/src/libinterfaces.a\
 ../minuit/src/libmyminuit.a $(CERNLIBS) @FLIBS@ \
 ../genetic/src/libgenetic.a \
 ../QEDevol/src/libQEDevol.a \
 ../common/num_utils/libnum_utils.a \
 ../tools/draw/src/pdferrors.o


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
  LDADD+= ../genetic/mixmax_r004/libmixmax.a
endif

EXTRA_DIST = xfitter_py.F
python:
	f2py -I../include  --f77flags="--fixed-line-length-132 -cpp" --quiet -c xfitter_py.F libxfitter.a  -m xfitter $(F2PYLIBS)  --fcompiler=$(F77)

