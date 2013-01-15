# Process this file with autoconf to produce a configure script.
AC_INIT([herafitter], [0.2.0], [herafitter-help@desy.de])

AM_INIT_AUTOMAKE
AC_PREFIX_DEFAULT(${PWD})


# Checks for programs.
AC_PROG_FC([gfortran])
AC_PROG_F77([gfortran])

AC_PROG_CPP
AC_PROG_CXX
AC_PROG_RANLIB

# Fortran configures
AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS
AC_F77_WRAPPERS
#AC_FC_FREEFORM

#
# A lot of code written with 132 char. strings.
#
FCFLAGS="$FCFLAGS -ffixed-line-length-132"
FFLAGS="$FFLAGS -ffixed-line-length-132"

# Check for trapFPE flag
AC_ARG_ENABLE([trapFPE],
  [AC_HELP_STRING(--enable-trapFPE, [Stop of floating point errors (default=no)])])
if test x$enable_trapFPE == xyes; then
  echo "Add floating point trap"
  FFLAGS="$FFLAGS -ffpe-trap='invalid,overflow' "
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

# Checks for libraries.
AC_MSG_CHECKING([for CERN libraries])
if test "x$CERN_ROOT" = "x"; then
  AC_MSG_ERROR([CERN_ROOT environment variable is not set!
  CERN_ROOT/lib must be the path to -lkernlib,   -lmathlib, -lpacklib])
else
  AC_SUBST([CERNLIBS], ["-L$CERN_ROOT/lib -lmathlib  -lkernlib -lpacklib"])
  AC_MSG_RESULT([Using -L$CERN_ROOT/lib -lmathlib   -lkernlib -lpacklib ])
fi

AC_MSG_CHECKING([for QCDNUM libraries])
if test "x$QCDNUM_ROOT" = "x"; then
   AC_MSG_ERROR([QCDNUM enviroment not set ! QCDNUM_ROOT/lib must point to libhqstf.a  libmbutil.a  libqcdnum.a  libzmstf.a])
else
  AC_SUBST([QCDNUMLIBS], ["-L$QCDNUM_ROOT/lib -lhqstf -lzmstf -lqcdnum -lmbutil"])
  AC_MSG_RESULT([Using -L$QCDNUM_ROOT/lib -lhqstf -lzmstf -lqcdnum -lmbutil])
fi


# NNPDF reweighting
AC_ARG_ENABLE([nnpdfWeight],
  [AC_HELP_STRING(--enable-nnpdfWeight, [use NNPDF weighting (default=no)])])
if test x$enable_nnpdfWeight == xyes; then
  enable_lhapdf=1
  echo "Force LHAPDF installation"
fi
AM_CONDITIONAL(ENABLE_NNPDFWEIGHT, [test x$enable_nnpdfWeight == xyes])

# check for lhapdf installation
AC_ARG_ENABLE([lhapdf],
  [AC_HELP_STRING(--enable-lhapdf, [use lhapdf (default=no)])])
if test $enable_lhapdf; then
  AC_MSG_CHECKING([for lhapdf installation])
  lhapdf_config=`which lhapdf-config`
  if test x$lhapdf_config == x; then
    AC_MSG_ERROR([Unable to find lhapdf-config.])
  else 
    AC_MSG_RESULT([Using $lhapdf_config])
    LHAPDF_CPPFLAGS=`lhapdf-config --cppflags`
    LHAPDF_LDFLAGS=`lhapdf-config --ldflags`
    AC_SUBST(LHAPDF_CPPFLAGS)
    AC_SUBST(LHAPDF_LDFLAGS)
  fi
fi
AM_CONDITIONAL(ENABLE_LHAPDF, [test $enable_lhapdf])

# check for applgrid installation
AC_ARG_ENABLE([applgrid],
  [AC_HELP_STRING(--enable-applgrid, [use applgrid for fast pdf convolutions (default=no)])])
if test x$enable_applgrid == xyes; then
  AC_MSG_CHECKING([for applgrid installation])
  applgrid_config=`which applgrid-config`
  if test x$applgrid_config == x; then
    AC_MSG_ERROR([Unable to find applgrid-config.])
  else 
    AC_MSG_RESULT([Using $applgrid_config])
    APPLGRID_CPPFLAGS=`applgrid-config --cxxflags`
    APPLGRID_LDFLAGS=`applgrid-config --ldcflags`
    AC_SUBST(APPLGRID_CPPFLAGS)
    AC_SUBST(APPLGRID_LDFLAGS)
  fi
fi
AM_CONDITIONAL(ENABLE_APPLGRID, [test x$enable_applgrid == xyes])

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

AC_MSG_CHECKING([for root installation])
root_config=`which root-config`
if test x$root_config == x; then
  AC_MSG_WARN([Unable to find root-config.])
  if test x$enable_applgrid == xyes; then
     AC_MSG_ERROR([Root is required for APPLGRID])
  fi
  if test x$enable_nnpdfWeight == xyes; then
     AC_MSG_ERROR([Root is required for NNPDF])
  fi
else 
  AC_MSG_RESULT([Using $root_config])
  root_ok=1
  ROOT_CFLAGS=`root-config --cflags`
  ROOT_LDFLAGS=`root-config --ldflags`
  ROOT_LIBS=`root-config --libs`
  AC_SUBST(ROOT_CFLAGS)
  AC_SUBST(ROOT_LDFLAGS)
  AC_SUBST(ROOT_LIBS)
fi

AM_CONDITIONAL([HAVE_ROOT],test $root_ok)


# UPDFs
AC_ARG_ENABLE([uPDF],
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
fi

#AC_F77_LIBRARY_LDFLAGS

# Check for presence of pdfLaTeX
AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex)
if test -z "$PDFLATEX"; then
  AC_MSG_WARN([Unable to create PDF version of the user manual.])
fi

# Check for presence of pdfLaTeX
AC_CHECK_PROG(BIBTEX, bibtex, bibtex)
if test -z "$BIBTEX"; then
  AC_MSG_WARN([Unable to create PDF version of the user manual.])
fi

AM_CONDITIONAL([HAVE_PDFLATEX], test -n "$PDFLATEX")
AM_CONDITIONAL([HAVE_BIBTEX], test -n "$BIBTEX")

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
AC_CONFIG_FILES([src/Makefile
                 interfaces/src/Makefile
                 tools/DrawResultsPackage/Makefile
                 tools/MakeLHAPDF/Makefile
                 FastNLO/src/Makefile
                 minuit/src/Makefile
                 HS/src/Makefile
                 RT/src/Makefile
                 DiffDIS/src/Makefile
                 DIPOLE/src/Makefile
                 DY/src/Makefile
                 Hathor/src/Makefile
                 ACOT/src/Makefile
                 ABKM/src/Makefile
                 NNPDF/src/Makefile
                 Cascade/src/Makefile
                 doc/tex/Makefile
                 Makefile])
# AC_CONFIG_SUBDIRS([])
AC_OUTPUT
