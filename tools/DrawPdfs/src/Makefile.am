#AUTOMAKE_OPTIONS = foreign

if HAVE_ROOT

# the main target, i.e. the executable to build
bin_PROGRAMS = DrawPdfs

# source files for the DrawPdfs executable
DrawPdfs_SOURCES = DrawPdfs.cc PdfsPainter.cc CommandParser.cc DataPainter.cc DrawLogo.cc ShiftPainter.cc FitPainter.cc ParPainter.cc Outdir.cc PdfData.cc Dataset.cc Chi2.cc Par.cc

DrawPdfs_CPPFLAGS = -I$(srcdir)/../include $(ROOT_CFLAGS) -w

# libraries to be linked to the main target
DrawPdfs_LDADD = $(ROOT_LIBS)

endif
