ACLOCAL_AMFLAGS=-I m4
AUTOMAKE_OPTIONS = foreign
SUBDIRS=src\
  include\
  minuit/src\
  EW/src\
  common\
  common/linalg\
  common/num_utils\
  pdf2yaml\
  DiffDIS/src\
  ABM/src\
  RT/src\
  genetic/mixmax_r004\
  genetic/src\
  QEDevol/src\
  DiffDIS/include\
  tools/draw/include\
  tools/process reactions/DYTurbo/src reactions/HathorSingleTop/src\
  tools/draw\
  tools/MakeLHAPDF\
  pdfparams/BasePdfParam/src\
  pdfparams/Factor/src\
  pdfparams/Expression/src\
  pdfparams/HERAPDF_PdfParam/src\
  pdfparams/PolySqrtPdfParam/src\
  pdfparams/NegativeGluonPdfParam/src\
  pdfparams/ABMPgluonPdfParam/src\
  pdfparams/ABMPseaPdfParam/src\
  pdfparams/ABMPvalencePdfParam/src \
  pdfdecompositions/BasePdfDecomposition/src\
  pdfdecompositions/UvDvUbarDbarS/src\
  pdfdecompositions/UvDvUbarDbarSSbarPdfDecomposition/src\
  pdfdecompositions/SU3_PionPdfDecomposition/src\
  evolutions/BaseEvolution/src\
  evolutions/APFEL/src\
  evolutions/QCDNUM/src\
  evolutions/APFELxx/src\
  evolutions/FlipUD/src\
  evolutions/FlipCharge/src\
  minimizers/BaseMinimizer/src\
  minimizers/MINUITMinimizer/src\
  reactions/APPLgrid/src\
  reactions/BaseDISCC/src\
  reactions/BaseDISNC/src\
  reactions/Fractal_DISNC/src\
  reactions/KMatrix/src\
  reactions/KFactor/src\
  reactions/RT_DISNC/src\
  reactions/FONLL_DISNC/src\
  reactions/FONLL_DISCC/src\
  reactions/FFABM_DISNC/src\
  reactions/FFABM_DISCC/src\
  reactions/BaseHVQMNR/src\
  reactions/cbdiff/src\
  reactions/fastNLO/src\
  reactions/BaseHVQMNR/src\
  reactions/HVQMNR_LHCb_7TeV_charm/src\
  reactions/HVQMNR_LHCb_7TeV_beauty/src\
  reactions/Hathor/src\
  reactions/HathorSingleTop/src\
  reactions/AFB/src
#

SUBDIRS+=xfitter
if ENABLE_CERES
  SUBDIRS+=minimizers/CERESMinimizer/src
endif

DIST_SUBDIRS = $(SUBDIRS)  datafiles QEDevol/include \
           datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042 doc/logo  \
           python/  tools/

if ENABLE_LHAPDF
  SUBDIRS+= evolutions/LHAPDF/src\
            pdfdecompositions/LHAPDFDecomposition/src

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
    steering.txt  parameters.yaml Reactions.txt constants.yaml
bin_SCRIPTS = xfitter-config
