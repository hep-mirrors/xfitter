ACLOCAL_AMFLAGS=-I m4
AUTOMAKE_OPTIONS = foreign
SUBDIRS=src\
  xfitter\
  minuit/src\
  DY/src\
  DIPOLE/src\
  RT/src\
  EW/src\
  common\
  common/linalg\
  common/num_utils\
  pdf2yaml\
  FastNLO/src\
  DiffDIS/src\
  ACOT/src\
  SACOT/src\
  ABM/src\
  FONLL/src\
  Cascade/src\
  genetic/mixmax_r004\
  genetic/src\
  QEDevol/src\
  include\
  FastNLO/include\
  FastNLO/include/fastnlotk\
  DiffDIS/include\
  DY/include\
  tools/draw/include\
  tools/process\
  tools/draw\
  tools/MakeLHAPDF\
  pdfparams/BasePdfParam/src\
  pdfparams/Expression/src\
  pdfparams/HERAPDF_PdfParam/src\
  pdfparams/PolySqrtPdfParam/src\
  pdfparams/NegativeGluonPdfParam/src\
  pdfdecompositions/BasePdfDecomposition/src\
  pdfdecompositions/UvDvUbarDbarS/src\
  pdfdecompositions/SU3_PionPdfDecomposition/src\
  evolutions/BaseEvolution/src\
  evolutions/QCDNUM/src\
  evolutions/APFELxx/src\
  evolutions/FlipUD/src\
  evolutions/FlipCharge/src\
  minimizers/BaseMinimizer/src\
  minimizers/MINUITMinimizer/src\
  reactions/APPLgrid/src\
  reactions/KFactor/src\
  reactions/Fractal_DISNC/src\
  reactions/BaseDISCC/src\
  reactions/Hathor/src\
  reactions/BaseDISNC/src\
  reactions/RT_DISNC/src\
  reactions/FFABM_DISNC/src\
  reactions/FFABM_DISCC/src\
  reactions/BaseHVQMNR/src\
  reactions/HVQMNR_LHCb_7TeV_beauty/src\
  reactions/HVQMNR_LHCb_7TeV_charm/src\
  reactions/testZMVFNS/src\
  reactions/fastNLO/src/\
  reactions/FONLL_DISCC/src\
  reactions/FONLL_DISNC/src\
  reactions/AFB/src\
  reactions/KMatrix/src

if ENABLE_CERES
  SUBDIRS+=minimizers/CERESMinimizer/src
endif

DIST_SUBDIRS = $(SUBDIRS) Hathor/src HVQMNR/src ACOT SACOT datafiles QEDevol/include \
           datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042 doc/logo examples Hathor/interface \
           FONLL HVQMNR HVQMNR/include/ python/ input_steering/ tools/

if ENABLE_LHAPDF
  SUBDIRS+= evolutions/LHAPDF/src\
            pdfdecompositions/LHAPDFDecomposition/src
endif

if ENABLE_HATHOR
  SUBDIRS+= Hathor/src
endif

if ENABLE_HVQMNR
  SUBDIRS+= HVQMNR/src
endif

if BUILD_DOCS
  DOCS=doc/tex/manual doc/tex/paper
endif
SUBDIRS+=$(DOCS)

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
