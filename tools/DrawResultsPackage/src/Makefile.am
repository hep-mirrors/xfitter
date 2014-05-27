#AUTOMAKE_OPTIONS = foreign

if HAVE_ROOT

# the main target, i.e. the executable to build
bin_PROGRAMS = DrawResults StoreGraphs StoreGraphsMC

# all source files which go into the main target
DrawResults_SOURCES = DrawResults.cc Painter.cc \
  Output.cc DataSet.cc PdfTable.cc

DrawResults_CPPFLAGS = -I$(srcdir)/../include $(ROOT_CFLAGS) -w

# libraries to be linked to the main target
DrawResults_LDADD = $(ROOT_LIBS)

# all source files 
StoreGraphs_SOURCES = StoreGraphs.cc PdfTable.cc

StoreGraphs_CPPFLAGS = -I$(srcdir)/../include $(ROOT_CFLAGS) -w

# libraries to be linked to the main target
StoreGraphs_LDADD = $(ROOT_LIBS)

# all source files 
StoreGraphsMC_SOURCES = StoreGraphsMC.cc PdfTable.cc

StoreGraphsMC_CPPFLAGS = -I$(srcdir)/../include $(ROOT_CFLAGS) -w

# libraries to be linked to the main target
StoreGraphsMC_LDADD = $(ROOT_LIBS)

endif
