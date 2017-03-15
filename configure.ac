# Process this file with autoconf to produce a configure script.

AC_INIT([XFitter],
m4_esyscmd(git describe --always 2>&1 | awk '/Not a git/  {print "2.0.0\033@<:@1;32m FrozenFrog\033@<:@0m"}; !/Not a git/ {print $0}'  | tr -d '\n'),
[xfitter-help@desy.de],[xfitter],[http://xfitter.org])
m4_ifndef([AC_PACKAGE_URL], AC_DEFINE([PACKAGE_URL],["http://xfitter.org"]))

#Suppress verbose output when compiling (use make V=99 for verbose output)
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

#Make sure long file names are saved by make dist:
AM_INIT_AUTOMAKE([1.9 tar-ustar])
AC_PREFIX_DEFAULT(${PWD})

# Checks for programs.

# Try to respect users' Fortran compiler variables
if test "x$FC" == "x"; then
  if test "x$F77" == "x"; then
    FC="$GFORTRAN"
  else
    FC="$F77"
  fi
fi

# Reset F77:
F77=$FC

AC_PROG_FC
AC_PROG_F77

AC_PROG_CPP
AC_PROG_CXX
AC_PROG_LIBTOOL


AM_PROG_CC_C_O
# Fortran configures
AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS
AC_F77_WRAPPERS
#AC_FC_FREEFORM

#Use libtoolize to build libraries
LT_INIT
AC_CONFIG_MACRO_DIR([m4])

# lapack macros
m4_include([m4/ax_blas.m4])
m4_include([m4/ax_lapack.m4])
# compare versions
m4_include([m4/ax_compare_version.m4])
# openmp suport
m4_include([m4/ax_openmp.m4])

#
# A lot of code written with 132 char. strings.
#
if test "uname -m" = "64"; then
   FCFLAGS="$FCFLAGS  -cpp  -fPIC -mcmodel=medium"
   FFLAGS="$FFLAGS    -cpp  -fPIC -mcmodel=medium"
else
   echo "Disabling mcmodel=medium option"
   FCFLAGS="$FCFLAGS   -fPIC -cpp"
   FFLAGS="$FFLAGS     -fPIC -cpp"
fi    


# fortran detection:
testtest=`echo $FC  | awk '/ifort/ {print "yes"}'`

if test "$testtest" = "yes"; then
# ifortran:
  AC_MSG_WARN([Intel compiler detected, setting options])
  FCFLAGS="$FCFLAGS -132 -Difort=1"
  FFLAGS="$FFLAGS -132 -Difort=1"
  NOAUTOFCFLAG="-noauto"
else
# gfortran most likely:
  AC_MSG_NOTICE([gfortran compiler options])
  FCFLAGS="$FCFLAGS  -ffixed-line-length-132 -Difort=0"
  FFLAGS="$FFLAGS  -ffixed-line-length-132 -Difort=0"
  NOAUTOFCFLAG="-fno-automatic"
fi

AC_SUBST(NOAUTOFCFLAG)

CFLAGS="$CFLAGS -fPIC"
CPPFLAGS="$CPPFLAGS -fPIC"

# check for lapack and blas
AX_BLAS
if test "$ax_blas_ok" = "no"; then
   AC_MSG_ERROR([Cannot find BLAS libraries])
fi

AX_LAPACK
if test "$ax_lapack_ok" = "no"; then
   AC_MSG_ERROR([Cannot find LAPACK libraries])
fi

AC_ARG_ENABLE([process],
                 [AS_HELP_STRING([--enable-process],[enable xfitter-process tool])])
: ${enable_process=yes} 

AS_IF([test "$enable_process" = "yes"], [
        AC_MSG_CHECKING([for libyaml])
        AC_LANG_PUSH(C)
        SAVED_LIBS=$LIBS
        LIBS="-lyaml $LIBS"
                AC_LINK_IFELSE(
                               [AC_LANG_PROGRAM([#include<yaml.h>],[yaml_document_t doc;])],
                                [AC_MSG_RESULT([yes])]
                                [AM_CONDITIONAL([HAVE_YAML],[true])],
                                [AC_MSG_RESULT([no])]
                                [AM_CONDITIONAL([HAVE_YAML],[false])]
                                [AC_MSG_WARN([libyaml is not installed, xfitter-process tool is not available.])]
                                )
                LIBS=$SAVED_LIBS
        AC_LANG_POP(C)],
[AM_CONDITIONAL([HAVE_YAML], [false])]
)

AC_ARG_ENABLE([openmp],
                 [AS_HELP_STRING([--enable-openmp],[enable openmp support])])
: ${enable_openmp=no} 

AS_IF([test "$enable_openmp" = "yes"], [
    AX_OPENMP([AM_CONDITIONAL(ENABLE_OPENMP,true)],
    [AM_CONDITIONAL(ENABLE_OPENMP,false)])
],[
AM_CONDITIONAL(ENABLE_OPENMP, test true = false)
])
AC_SUBST(OPENMP_CFLAGS)
dnl OpenMP checker only defines for C when compiling both C and C++
OPENMP_CXXFLAGS=$OPENMP_CFLAGS
AC_SUBST(OPENMP_CXXFLAGS)

# Check for trapFPE flag
AC_ARG_ENABLE([trapFPE],
  [AC_HELP_STRING(--enable-trapFPE, [Stop of floating point errors (default=no)])])
if test x$enable_trapFPE == xyes; then
  echo "Add floating point trap"
  FFLAGS="$FFLAGS -ffpe-trap='invalid,overflow,zero' "
fi

# Check for check bounds flag
AC_ARG_ENABLE([checkBounds],
  [AC_HELP_STRING(--enable-checkBounds, [add -fbounds-check flag for compilation (default=no)])])
if test x$enable_checkBounds == xyes; then
  echo "Add check bounds flag"
  FFLAGS="$FFLAGS -fbounds-check"
fi


# Checks for sources.
AC_CONFIG_SRCDIR([src/main.f])

## Checks for libraries.
#AC_MSG_CHECKING([for CERN libraries])
#if test "x$CERN_ROOT" = "x"; then
#  AC_MSG_ERROR([CERN_ROOT environment variable is not set!
#  CERN_ROOT/lib must be the path to -lkernlib,   -lmathlib, -lpacklib])
#else
#  AC_LANG(Fortran)
#  LDFLAGS+="-L$CERN_ROOT/lib64"
#  AC_SEARCH_LIBS([dinv],[kernlib],[AC_MSG_RESULT([Using -L$CERN_ROOT/lib -lkernlib])],
#  [AC_MSG_ERROR([-lkernlib not found in CERN_ROOT])])
#  AC_SEARCH_LIBS([dgammf],[mathlib],[AC_MSG_RESULT([Using -L$CERN_ROOT/lib -lmathlib])],
#  [AC_MSG_ERROR([-lmathlib not found in CERN_ROOT])])
#  AC_SEARCH_LIBS([dsyevd],[lapack],[AC_MSG_RESULT([Using -L$CERN_ROOT/lib -llapack])],
#  [AC_MSG_ERROR([-llapack not found in CERN_ROOT])])
##  AC_SUBST([CERNLIBS], ["-L$CERN_ROOT/lib -lmathlib  -lkernlib -lpacklib  -llapack -lblas"])
##  AC_MSG_RESULT([Using -L$CERN_ROOT/lib -lmathlib   -lkernlib -lpacklib  -llapack -lblas])
#fi

#AC_MSG_CHECKING([for QCDNUM libraries])
#if test "x$QCDNUM_ROOT" = "x"; then
#   AC_MSG_ERROR([QCDNUM enviroment not set ! QCDNUM_ROOT/lib must point to libhqstf.a  libmbutil.a  libqcdnum.a  libzmstf.a])
#else
#  AC_SUBST([QCDNUMLIBS], ["-L$QCDNUM_ROOT/lib -lhqstf -lzmstf -lqcdnum -lmbutil"])
#  AC_MSG_RESULT([Using -L$QCDNUM_ROOT/lib -lhqstf -lzmstf -lqcdnum -lmbutil])
#fi

# check for QCDNUM installation
AC_MSG_CHECKING([for QCDNUM installation])
qcdnum_config=`which qcdnum-config`
if test x$qcdnum_config == x; then
  AC_MSG_ERROR([Unable to find qcdnum-config.])
else 
  QCDNUMVERS=`qcdnum-config --version`
  QCDNUMLIBS=`qcdnum-config --ldflags`
  AC_MSG_RESULT([Using $qcdnum_config version $QCDNUMVERS])
  AC_SUBST(QCDNUMLIBS)
  AX_COMPARE_VERSION([17-01-13],[gt],[$QCDNUMVERS],[AC_MSG_ERROR([QCDNUM version outdated (17-01-13 or later needed).])],[])
fi

# check for lhapdf installation
AC_ARG_ENABLE([lhapdf],
  [AC_HELP_STRING(--enable-lhapdf, [use lhapdf (default=no)])])
if test $enable_lhapdf; then
  AC_MSG_CHECKING([for lhapdf installation])
  lhapdf_config=`which lhapdf-config`
  if test x$lhapdf_config == x; then
    AC_MSG_ERROR([Unable to find lhapdf-config.])
  else 
    LHAPDFVERS=`lhapdf-config --version`
    LHAPDF_CPPFLAGS=`lhapdf-config --cppflags`
    LHAPDF_LDFLAGS=`lhapdf-config --ldflags`
    AC_MSG_RESULT([Using $lhapdf_config version $LHAPDFVERS])
    AC_SUBST(LHAPDF_FAMILY)
    AC_SUBST(LHAPDF_CPPFLAGS)
    AC_SUBST(LHAPDF_LDFLAGS)
    AC_DEFINE([LHAPDF_ENABLED],[1],[Define if LHAPDF is enabled])

    lhapdf_family=`lhapdf-config --version | sed 's/\..*$//'`
    if test "x$lhapdf_family" != "x5"; then 
      AC_DEFINE([LHAPDF_FAMILY],[6],[Define value of LHAPDF family])
    fi
  fi
fi
AM_CONDITIONAL(ENABLE_LHAPDF, [test $enable_lhapdf])

# check for apfel installation
AC_ARG_ENABLE([apfel],
  [AC_HELP_STRING(--enable-apfel, [use apfel (default=no)])])
if test $enable_apfel; then
  AC_MSG_CHECKING([for apfel installation])
  apfel_config=`which apfel-config`
  if test x$apfel_config == x; then
    AC_MSG_ERROR([Unable to find apfel-config.])
  else 
    APFELVERS=`apfel-config --version`
    APFEL_CPPFLAGS=`apfel-config --cppflags`
    APFEL_LDFLAGS=`apfel-config --ldflags`
    AC_MSG_RESULT([Using $apfel_config version $APFELVERS])
    AC_SUBST(APFEL_CPPFLAGS)
    AC_SUBST(APFEL_LDFLAGS)
    AC_DEFINE([APFEL_ENABLED],[1],[Define if APFEL is enabled])
    AX_COMPARE_VERSION([2.6.1],[gt],[$APFELVERS],[AC_MSG_ERROR([APFEL version outdated (2.6.1 or later needed).])],[])
  fi
fi
AM_CONDITIONAL(ENABLE_APFEL, [test $enable_apfel])

# check for mela installation
AC_ARG_ENABLE([mela],
  [AC_HELP_STRING(--enable-mela, [use mela (default=no)])])
if test $enable_mela; then
  AC_MSG_CHECKING([for mela installation])
  mela_config=`which mela-config`
  if test x$mela_config == x; then
    AC_MSG_ERROR([Unable to find mela-config.])
  else 
    MELAVERS=`mela-config --version`
    MELA_LDFLAGS=`mela-config --ldflags`
    AC_MSG_RESULT([Using $mela_config version $MELAVERS])
    AC_SUBST(MELA_LDFLAGS)
    AC_DEFINE([MELA_ENABLED],[1],[Define if MELA is enabled])
    AX_COMPARE_VERSION([2.0.0],[gt],[$MELAVERS],[AC_MSG_ERROR([MELA version outdated (2.0.0 or later needed).])],[])
  fi
fi
AM_CONDITIONAL(ENABLE_MELA, [test $enable_mela])

# check for applgrid installation
AC_ARG_ENABLE([applgrid],
  [AC_HELP_STRING(--enable-applgrid, [use applgrid for fast pdf convolutions (default=no)])])
if test x$enable_applgrid == xyes; then
  AC_MSG_CHECKING([for applgrid installation])
  applgrid_config=`which applgrid-config`
  if test x$applgrid_config == x; then
    AC_MSG_ERROR([Unable to find applgrid-config.])
  else 
    APPLGRIDVERS=`applgrid-config --version`
    APPLGRID_CPPFLAGS=`applgrid-config --cxxflags`
    APPLGRID_LDFLAGS=`applgrid-config --ldcflags`
    AC_MSG_RESULT([Using $applgrid_config version $APPLGRIDVERS])
    AC_SUBST(APPLGRID_CPPFLAGS)
    AC_SUBST(APPLGRID_LDFLAGS)
    AC_DEFINE([APPLGRID_ENABLED],[1],[Define if applgrid is enabled])
  fi
fi
AM_CONDITIONAL(ENABLE_APPLGRID, [test x$enable_applgrid == xyes])

# check if applgrid support CKM replacement
ag_ckm=no
ag_version=`applgrid-config --version`
if test x$enable_applgrid == xyes; then
  AC_MSG_CHECKING([if APPLgrid supports CKM replacement])
  AX_COMPARE_VERSION([1.4.33],[gt],[$ag_version],[ag_ckm=no],[ag_ckm=yes])
  if test x$ag_ckm == xyes; then
    AC_MSG_RESULT([yes])
    AC_DEFINE([APPLGRID_CKM],[1],[Define if APPLGrid supports CKM definition])
  else
    AC_MSG_WARN([The APPLgrid version does not allow to set the CKM matrix ( needs v. >= 1.4.33 )])
  fi
fi
AM_CONDITIONAL(APPLGRID_CKM, [test x$ag_ckm == xyes])

# check if applgrid support dynamic scale emulation
ag_dynscale=no
ag_version=`applgrid-config --version`
if test x$enable_applgrid == xyes; then
  AC_MSG_CHECKING([if APPLgrid supports dynamic scale emulation])
  AX_COMPARE_VERSION([1.4.43],[gt],[$ag_version],[ag_dynscale=no],[ag_dynscale=yes])
  if test x$ag_dynscale == xyes; then
    AC_MSG_RESULT([yes])
    AC_DEFINE([APPLGRID_DYNSCALE],[1],[Define if APPLGrid supports dynamic scale emulation])
  else
    AC_MSG_WARN([The APPLgrid version does not allow dynamic scale emulation ( needs v. >= 1.4.43 )])
  fi
fi
AM_CONDITIONAL(APPLGRID_DYNSCALE, [test x$ag_dynscale == xyes])

# check for APFELgrid installation
AC_ARG_ENABLE([apfelgrid],
  [AC_HELP_STRING(--enable-apfelgrid, [use apfelgrid for fast pdf convolutions (default=no)])])
if test x$enable_apfelgrid == xyes; then
  AC_MSG_CHECKING([for apfelgrid installation])
  apfelgrid_config=`which apfelgrid-config`
  if test x$apfelgrid_config == x; then
    AC_MSG_ERROR([Unable to find apfelgrid-config.])
  else
    APFELGRIDVERS=`apfelgrid-config --version`
    APFELGRID_CPPFLAGS=`apfelgrid-config --cxxflags`
    APFELGRID_LDFLAGS=`apfelgrid-config --ldflags`
    AC_MSG_RESULT([Using $apfelgrid_config version $APFELGRIDVERS])
    AC_SUBST(APFELGRID_CPPFLAGS)
    AC_SUBST(APFELGRID_LDFLAGS)
    AC_DEFINE([APFELGRID_ENABLED],[1],[Define if apfelgrid is enabled])
  fi
fi
AM_CONDITIONAL(ENABLE_APFELGRID, [test x$enable_apfelgrid == xyes])

# check for genetic installation
AC_ARG_ENABLE([genetic],
    [AC_HELP_STRING(--enable-genetic, [use genetic for general minimia search (defaults=no)])])
if test x$enable_genetic == xyes; then
  AC_MSG_CHECKING([for genetic installation])
# FIXME: add a condition (version(gcc) > xxx) here
fi
AM_CONDITIONAL(ENABLE_GENETIC, [test x$enable_genetic == xyes])

# check for hathor installation
AC_ARG_ENABLE([hathor],
  [AC_HELP_STRING(--enable-hathor, [use hathor for ttbar cross section predictions (default=no)])])
if test x$enable_hathor == xyes; then
  AC_MSG_CHECKING([for hathor installation])
  if test "x$HATHOR_ROOT" = "x"; then
    AC_MSG_ERROR([HATHOR_ROOT environment variable is not set!.])
  else 
    AC_MSG_RESULT([Using $HATHOR_ROOT])
  fi
fi
AM_CONDITIONAL(ENABLE_HATHOR, [test x$enable_hathor == xyes])

# check for HVQMNR installation
AC_ARG_ENABLE([hvqmnr],
  [AC_HELP_STRING(--enable-hvqmnr, [use hvqmnr for heavy-flavour pp cross section predictions (default=no)])])
if test x$enable_hvqmnr == xyes; then
  AC_MSG_CHECKING([for HVQMNR installation])
fi
AM_CONDITIONAL(ENABLE_HVQMNR, [test x$enable_hvqmnr == xyes])

# check for ROOT libraries
AC_ARG_ENABLE([root],[AS_HELP_STRING([--enable-root],[use ROOT libraries])]): ${enable_root=yes}

AS_IF([test "$enable_root" = "yes"], [
	    AC_MSG_CHECKING([for ROOT installation])
	    root_config=`which root-config`
	    if test x$root_config == x; then
	      AC_MSG_ERROR([Unable to find root-config, install ROOT or configure with --disable-root])
	    else 
	      AC_MSG_RESULT([Using $root_config])
	      root_ok=1
	      ROOT_CFLAGS=`root-config --cflags`
	      ROOT_LDFLAGS=`root-config --ldflags`
	      ROOT_LIBS=`root-config --libs`
	      AC_SUBST(ROOT_CFLAGS)
	      AC_SUBST(ROOT_LDFLAGS)
	      AC_SUBST(ROOT_LIBS)
	      AC_DEFINE([ROOT_ENABLED],[1],[Define if ROOT is enabled])
	    fi
	  ],
	  [
	    AC_MSG_WARN([ROOT libraries are disabled, xfitter-draw not available])
	    if test x$enable_applgrid == xyes; then
	       AC_MSG_ERROR([Root is required for APPLGRID])
	    fi
	    ])
AM_CONDITIONAL([HAVE_ROOT],test $root_ok)


#AC_MSG_CHECKING([for root installation])
#root_config=`which root-config`
#if test x$root_config == x; then
#  AC_MSG_WARN([Unable to find root-config.])
#  if test x$enable_applgrid == xyes; then
#     AC_MSG_ERROR([Root is required for APPLGRID])
#  fi
#  if test x$enable_nnpdfWeight == xyes; then
#     AC_MSG_ERROR([Root is required for NNPDF])
#  fi
#else 
#  AC_MSG_RESULT([Using $root_config])
#  root_ok=1
#  ROOT_CFLAGS=`root-config --cflags`
#  ROOT_LDFLAGS=`root-config --ldflags`
#  ROOT_LIBS=`root-config --libs`
#  AC_SUBST(ROOT_CFLAGS)
#  AC_SUBST(ROOT_LDFLAGS)
#  AC_SUBST(ROOT_LIBS)
#  AC_DEFINE([ROOT_ENABLED],[1],[Define if ROOT is enabled])
#fi
#AM_CONDITIONAL([HAVE_ROOT],test $root_ok)


# UPDFs

AC_ARG_ENABLE([updf],
  [AC_HELP_STRING(--enable-updf, [use uPDF evolution (default=no)])])
AM_CONDITIONAL(ENABLE_UPDF, [test x$enable_updf == xyes])

if test x$enable_updf == xyes; then

   AC_MSG_CHECKING([for CASCADE libraries])
   if test "x$CASCADE_ROOT" = "x"; then
      AC_MSG_ERROR([CASCADE enviroment not set ! CASCADE_ROOT/lib must point to libbases.als.  Check Cascade/README for more details.])
   else
      AC_SUBST([CASCADE_SOURCE], ["$CASCADE_ROOT"])
      AC_SUBST([CASCADELIBS], ["-L$CASCADE_ROOT/lib -lbases"])
      AC_MSG_RESULT([Using -L$CASCADE_ROOT/lib  -lbases])
   fi

   AC_MSG_CHECKING([for PYTHIA libraries])
   if test "x$PYTHIA_ROOT" = "x"; then
      AC_MSG_ERROR([PYTHIA enviroment not set ! PYTHIA_ROOT/lib must point to ....])
   else
      AC_SUBST([PYTHIA_SOURCE], ["$PYTHIA_ROOT"])
      AC_SUBST([PYTHIALIBS], ["-L$PYTHIA_ROOT/lib -lpythia6 -lpythia6_dummy $PYTHIA_ROOT/lib/pydata.o "])
      AC_MSG_RESULT([Using -L$PYTHIA_ROOT/lib -lpythia6 -lpythia6_dummy $PYTHIA_ROOT/lib/pydata.o ])
   fi
      AC_SUBST([CASCADE_SOURCE], ["$CASCADE_ROOT"])
      AC_SUBST([CASCADELIBS], ["-L$CASCADE_ROOT/lib -lbases"])
      AC_SUBST([PYTHIA_SOURCE], ["$PYTHIA_ROOT"])
      AC_SUBST([PYTHIALIBS], ["-L$PYTHIA_ROOT/lib -lpythia6 -lpythia6_dummy $PYTHIA_ROOT/lib/pydata.o "])
fi

#AC_F77_LIBRARY_LDFLAGS

# Check if we need to build documentation
AC_ARG_ENABLE([doc],
  [AC_HELP_STRING(--enable-doc, [Build documentation (default=no)])])
if test x$enable_doc = xyes; then
  # Check for presence of pdfLaTeX
  AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex)
  if test -z "$PDFLATEX"; then
    AC_MSG_ERROR([Pdflatex is required to create the user manual.])
  fi
  
  # Check for presence of pdfLaTeX
  AC_CHECK_PROG(BIBTEX, bibtex, bibtex)
  if test -z "$BIBTEX"; then
    AC_MSG_ERROR([Bibtex is required to create the user manual.])
  fi
fi
AM_CONDITIONAL([BUILD_DOCS], [test x$enable_doc = xyes])

# DOXYGEN SUPPORT

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN([$PACKAGE_NAME],[doxygen.cfg])

# Output
AC_CONFIG_FILES([include/Makefile 
                 src/Makefile
                 common/Makefile
                 common/linalg/Makefile
                 common/num_utils/Makefile
                 interfaces/src/Makefile
                 interfaces/include/Makefile
                 FastNLO/src/Makefile
                 FastNLO/include/Makefile
                 FastNLO/include/fastnlotk/Makefile
                 minuit/src/Makefile
                 EW/src/Makefile
                 RT/src/Makefile
                 DiffDIS/src/Makefile
                 DiffDIS/include/Makefile
                 DIPOLE/src/Makefile
                 DY/src/Makefile
                 DY/include/Makefile
                 Hathor/src/Makefile
                 Hathor/interface/Makefile
                 HVQMNR/Makefile
                 HVQMNR/src/Makefile
                 HVQMNR/include/Makefile
                 ACOT/src/Makefile
                 ACOT/Makefile
                 SACOT/src/Makefile
                 SACOT/Makefile 
                 ABM/src/Makefile
		 FONLL/src/Makefile
                 FONLL/Makefile
                 Cascade/src/Makefile
                 genetic/mixmax_r004/Makefile
                 genetic/src/Makefile
                 QEDevol/src/Makefile
                 QEDevol/include/Makefile
                 pdf2yaml/Makefile
                 tools/process/Makefile
                 tools/draw/Makefile
                 tools/draw/src/Makefile
                 tools/draw/include/Makefile
                 tools/MakeLHAPDF/Makefile
                 input_steering/Makefile
                 tools/Makefile
                 datafiles/Makefile
                 datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/Makefile
                 doc/logo/Makefile
                 examples/Makefile
                 python/Makefile
		 xfitter-config
                 Makefile])

if test x$enable_doc = xyes; then
AC_CONFIG_FILES([doc/tex/manual/Makefile
                 doc/tex/paper/Makefile])
fi
AC_OUTPUT
