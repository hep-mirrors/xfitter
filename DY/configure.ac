# Process this file with autoconf to produce a configure script.
AC_INIT([libmydy], [1.0], [voica@mail.desy.de])
AM_INIT_AUTOMAKE([libmydy], [1.0])
AC_PREFIX_DEFAULT(${HOME})

# Checks for programs.
AC_PROG_FC([gfortran])
#AC_PROG_F77([gfortran])
AC_PROG_CPP
AC_PROG_RANLIB
AC_PROG_CXX

# Fortran configures
AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS
#AC_F77_WRAPPERS
#AC_FC_FREEFORM

# Output
AC_CONFIG_FILES([src/Makefile
                 Makefile])
AC_OUTPUT
