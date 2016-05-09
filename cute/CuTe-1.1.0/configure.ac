	# Prelude.
AC_INIT([CuTe], [1.1.0], [dwilhelm@students.uni-mainz.de])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
	# Checks for programs.
AC_PROG_CXX
	# Checks for librarys.
		# Checks for gsl
AC_MSG_NOTICE([***** Check for gsl  *****])
AC_CHECK_LIB([m],[cos])
if	test "$GSL_CONFIG" != "";
then 	AC_MSG_NOTICE([GSL_CONFIG = $GSL_CONFIG])
	AC_SUBST([LIB_GSL],`$GSL_CONFIG --libs`)
	AC_SUBST([INC_GSL],`$GSL_CONFIG --cflags`)
	AC_SUBST([LDD_GSL],["-Xlinker -R`$GSL_CONFIG --prefix`/lib"])
	AC_MSG_NOTICE([INC_GSL = $INC_GSL])
	AC_MSG_NOTICE([LIB_GSL = $LIB_GSL])
	AC_MSG_NOTICE([LDD_GSL = $LDD_GSL])
else	AC_CHECK_LIB([gslcblas],[cblas_dgemm],[],[AC_MSG_ERROR([cannot find library gslcblas])])
	AC_CHECK_LIB([gsl],[gsl_blas_dgemm],[],[AC_MSG_ERROR([cannot find library gsl])])
	AC_SUBST([LIB_GSL],[""])
	AC_SUBST([INC_GSL],[""])
	AC_SUBST([IDD_GSL],[""])	
fi
		# Checks for LHAPDF
AC_MSG_NOTICE([***** Check for LHAPDF  *****])
if	test "$LHAPDF_CONFIG" != "";
then 	AC_MSG_NOTICE([LHAPDF_CONFIG = $LHAPDF_CONFIG])
	LHAPDF_INC_DIR=`$LHAPDF_CONFIG --incdir`
	LHAPDF_LIB_DIR=`$LHAPDF_CONFIG --libdir`
	AC_SUBST([LIB_LHAPDF],["-L$LHAPDF_LIB_DIR -lLHAPDF"])
	AC_SUBST([INC_LHAPDF],["-I$LHAPDF_INC_DIR"])
	AC_SUBST([LDD_LHAPDF],["-Xlinker -R$LHAPDF_LIB_DIR"])
	AC_MSG_NOTICE([LIB_LHAPDF = $LIB_LHAPDF])
	AC_MSG_NOTICE([INC_LHAPDF = $INC_LHAPDF])
	AC_MSG_NOTICE([LDD_LHAPDF = $LDD_LHAPDF])
else	AC_LANG([Fortran])
	AC_SUBST([LHAPDF_CONFIG],[lhapdf-config])
	AC_MSG_NOTICE([LHAPDF_CONFIG = $LHAPDF_CONFIG])
	AC_CHECK_LIB([LHAPDF],[initpdfsetm],[],[AC_MSG_ERROR([cannot find library LHAPDF])])
	AC_SUBST([LIB_LHAPDF],[""])
	AC_SUBST([INC_LHAPDF],[""])
	AC_SUBST([FLAG_LHAPDF],[""])		
fi
LHAPDF_PDF_DIR=`$LHAPDF_CONFIG --pdfsets-path`
AC_DEFINE_UNQUOTED([LHAPDF_PDF_DIR],["$LHAPDF_PDF_DIR"],[Path to the directory containing the PDF set data files])
AC_MSG_NOTICE([LHAPDF_PDF_DIR = $LHAPDF_PDF_DIR])
AC_MSG_NOTICE([**********])


	# Output files.
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

