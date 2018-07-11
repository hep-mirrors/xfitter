ACLOCAL_AMFLAGS=-I m4
AUTOMAKE_OPTIONS = foreign
SUBDIRS = minuit/src interfaces/src DY/src DIPOLE/src RT/src EW/src common common/linalg \
	common/num_utils \
	tools/draw tools/MakeLHAPDF FastNLO/src DiffDIS/src ACOT/src SACOT/src ABM/src FONLL/src Cascade/src \
	genetic/mixmax_r004 genetic/src QEDevol/src \
	include interfaces/include FastNLO/include FastNLO/include/fastnlotk DiffDIS/include \
	DY/include tools/draw/include \
	pdf2yaml tools/process \
	reactions/KFactor/src reactions/Fractal_DISNC/src reactions/BaseDISCC/src reactions/Hathor/src reactions/BaseDISNC/src \
	reactions/RT_DISNC/src reactions/FFABM_DISNC/src reactions/FFABM_DISCC/src reactions/APPLgrid/src reactions/BaseHVQMNR/src \
	reactions/HVQMNR_LHCb_7TeV_beauty/src reactions/HVQMNR_LHCb_7TeV_charm/src reactions/testZMVFNS/src reactions/fastNLO/src/ \
	reactions/FONLL_DISCC/src reactions/FONLL_DISNC/src pdfparams/BasePdfParam/src \
        pdfdecompositions/BasePdfDecomposition/src evolutions/BaseEvolution/src

DIST_SUBDIRS = $(SUBDIRS) Hathor/src HVQMNR/src ACOT SACOT datafiles QEDevol/include \
           datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042 doc/logo examples Hathor/interface \
           FONLL HVQMNR HVQMNR/include/ python/ input_steering/ tools/

if ENABLE_HATHOR
  SUBDIRS+= Hathor/src
endif

if ENABLE_HVQMNR
  SUBDIRS+= HVQMNR/src
endif

if BUILD_DOCS
  DOCS=doc/tex/manual doc/tex/paper
endif
#$(info ${DOCS})
SUBDIRS+= src xfitter $(DOCS) 

include $(top_srcdir)/aminclude.am

reactiondir = $(prefix)/lib/
reaction_DATA = Reactions.txt

# Tell which program should run the .sh scripts.
#SH_LOG_COMPILER = $(SHELL) -ex
#TESTS_ENVIRONMENT = $(SHELL) 
#TESTS= ./tools/check.sh

EXTRA_DIST= README INSTALLATION LICENCE REFERENCES \
   minuit.in.txt steering.txt testPythonInterface.py parameters.yaml Reactions.txt
bin_SCRIPTS = xfitter-config tools/xfitter-getdata.sh
