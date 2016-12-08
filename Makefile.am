ACLOCAL_AMFLAGS=-I m4
AUTOMAKE_OPTIONS = foreign
SUBDIRS = minuit/src interfaces/src DY/src DIPOLE/src RT/src EW/src common common/linalg \
          common/num_utils \
          tools/draw tools/MakeLHAPDF FastNLO/src DiffDIS/src ACOT/src SACOT/src ABM/src FONLL/src Cascade/src \
          genetic/mixmax_r004 genetic/src QEDevol/src \
          include interfaces/include FastNLO/include FastNLO/include/fastnlotk DiffDIS/include \
          DY/include tools/draw/include input_steering  \
          pdf2yaml tools/process reactions/testZMVFNS/src reactions/fastNLO/src/ 
if ENABLE_HATHOR
  SUBDIRS+= Hathor/src
endif

if ENABLE_HVQMNR
  SUBDIRS+= HVQMNR/src
endif

SUBDIRS+= common
if BUILD_DOCS
  DOCS=doc/tex/manual doc/tex/paper
endif
#$(info ${DOCS})
SUBDIRS+= src xfitter $(DOCS) 

include $(top_srcdir)/aminclude.am

# Specify which SUBDIRS to distribute:
DIST_SUBDIRS = minuit/src interfaces/src DY/src DIPOLE/src RT/src EW/src \
          tools/draw tools/MakeLHAPDF FastNLO/src DiffDIS/src ACOT/src SACOT/src ABM/src FONLL/src Cascade/src \
          genetic/mixmax_r004 genetic/src common/num_utils \
          include interfaces/include FastNLO/include FastNLO/include/fastnlotk DiffDIS/include \
          DY/include tools/draw/include input_steering \
          common common/linalg Hathor/src pdf2yaml tools/process reactions/testZMVFNS/src reactions/fastNLO/src src xfitter

# Tell which program should run the .sh scripts.
#SH_LOG_COMPILER = $(SHELL) -ex
#TESTS_ENVIRONMENT = $(SHELL) 
#TESTS= ./tools/check.sh

EXTRA_DIST= README README-MacOSX LICENCE REFERENCES ewparam.txt minuit.in.txt steering.txt 
bin_SCRIPTS = xfitter-config

# the main target, i.e. the executable to build
# bin_PROGRAMS = xfitter
# xfitter_SOURCES = ./src/main.f
# xfitter_LDFLAGS = -L${pkglibdir} -L${srcdir}/src -L${top_builddir}/src -L${top_srcdir}/src -L${top_srcdir}/FastNLO/src -lxfmain -lxffnlo
# xfitter_LDADD =  
# noinst_LIBRARIES = #./src/libxfitter.a 
# AM_CXXFLAGS = -I$(srcdir)/include  -I$(srcdir)/interfaces/include -I$(srcdir)/FastNLO/include
# AM_CXXFLAGS += -I$(srcdir)/common -I$(srcdir)/common/linalg
# AM_CXXFLAGS += -I$(srcdir)/tools/draw/include
# AM_FFLAGS = $(DEFS) -cpp -I$(srcdir)/include -fno-automatic -I$(srcdir)/minuit/src
# AM_CPPFLAGS = -I$(srcdir)/include   -DPREFIX="\"$(prefix)\""
# # xFitter main library:
# LDADD =
# XFLIBS =  ./DY/src/libmydy.a\
#  ./Cascade/src/libcasutil.a \
#  ./common/libHFcommon.a \
#  ./RT/src/libmyrt.a \
#  ./ACOT/src/libmyacot.a \
#  ./SACOT/src/libmysacot.a \
#  ./FONLL/src/libmyfonll.a \
#  ./EW/src/libmyew.a \
#  ./ABM/src/libmyabm.a \
#  ./FastNLO/src/libFastNLO.a \
#  ./DIPOLE/src/libmydipole.a \
#  ./DiffDIS/src/libDiffDIS.a\
#  ./interfaces/src/libinterfaces.a\
#  ./minuit/src/libmyminuit.a $(CERNLIBS) @FLIBS@ \
#   -lstdc++ \
#  ./genetic/src/libgenetic.a \
#  ./QEDevol/src/libQEDevol.a \
#  ./common/num_utils/libnum_utils.a \
#  ./tools/draw/src/pdferrors.o
# 
# 
# LDADD += $(XFLIBS)
# #F2PYLIBS += $(XFLIBS)
# if ENABLE_GENETIC
#   LDADD+= ../genetic/mixmax_r004/libmixmax.a
# endif
# xfitter_LDADD += $(LDADD)
