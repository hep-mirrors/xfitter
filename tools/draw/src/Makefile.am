#AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = -I$(srcdir)/../include
#	-DSHAREDIR=\"${datadir}\"

lib_LTLIBRARIES = libpdferrors.la
libpdferrors_la_SOURCES = pdferrors.cc

if HAVE_ROOT

# the main target, i.e. the executable to build
bin_PROGRAMS = xfitter-draw

# source files for the xfitter-draw executable
xfitter_draw_SOURCES = xfitter-draw.cc PdfsPainter.cc CommandParser.cc DataPainter.cc DrawLogo.cc ShiftPainter.cc FitPainter.cc ParPainter.cc Chi2scanPainter.cc Outdir.cc PdfData.cc Dataset.cc Chi2.cc Par.cc Chi2scanData.cc pdferrors.cc

# libraries to be linked to the main target
xfitter_draw_LDADD = $(ROOT_LIBS)  ../../../common/libHFcommon.a libpdferrors.la

xfitter_draw_CPPFLAGS = -I$(srcdir)/../include -I$(srcdir)/../../../common $(ROOT_CFLAGS) -w

endif
