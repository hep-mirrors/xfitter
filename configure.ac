# Process this file with autoconf to produce a configure script.
AC_INIT([h1fitter], [0.0.1], [glazov@mail.desy.de])
AM_INIT_AUTOMAKE
AC_PREFIX_DEFAULT(${HOME})

# Checks for programs.
AC_PROG_F77
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_RANLIB

# Checks for sources.
AC_CONFIG_SRCDIR([src/main.f])

# Checks for libraries.
AC_MSG_CHECKING([for CERN libraries])
if test "x$CERN_ROOT" = "x"; then
  AC_MSG_ERROR([CERN_ROOT environment variable is not set!
  CERN_ROOT/lib must be the path to -lkernlib, -lphtools, -llapack3, -lmathlib, -lpacklib])
else
  AC_SUBST([CERNLIBS], ["-L$CERN_ROOT/lib -lmathlib -lphtools -llapack3 -lkernlib -lpacklib"])
  AC_MSG_RESULT([Using -L$CERN_ROOT/lib -lmathlib -lphtools -llapack3 -lkernlib -lpacklib ])
fi

AC_MSG_CHECKING([for QCDNUM libraries])
if test "x$QCDNUM_ROOT" = "x"; then
   AC_MSG_ERROR([QCDNUM enviroment not set ! QCDNUM_ROOT/lib must point to libhqstf.a  libmbutil.a  libqcdnum.a  libzmstf.a])
else
  AC_SUBST([QCDNUMLIBS], ["-L$QCDNUM_ROOT/lib -lhqstf -lzmstf -lqcdnum -lmbutil"])
  AC_MSG_RESULT([Using -L$QCDNUM_ROOT/lib -lhqstf -lzmstf -lqcdnum -lmbutil])
fi

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

if test x$enable_applgrid == xyes; then
  AC_MSG_CHECKING([for root installation, needed for applgrid])
  root_config=`which root-config`
  if test x$root_config == x; then
    AC_MSG_ERROR([Unable to find root-config.])
  else 
    AC_MSG_RESULT([Using $root_config])
    ROOT_CFLAGS=`root-config --cflags`
    ROOT_LDFLAGS=`root-config --ldflags`
    AC_SUBST(ROOT_CFLAGS)
    AC_SUBST(ROOT_LDFLAGS)
  fi
fi

AC_F77_LIBRARY_LDFLAGS

# Output
AC_CONFIG_FILES([src/Makefile
                 Makefile])
AC_CONFIG_SUBDIRS([minuit DY RT])
AC_OUTPUT
