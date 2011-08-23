#AUTOMAKE_OPTIONS = foreign

# the main target, i.e. the executable to build
bin_PROGRAMS = DrawResults StoreGraphs

# all source files which go into the main target
DrawResults_SOURCES = H1FitterDrawResults.cc H1FitterPainter.cc \
  H1FitterOutput.cc DataSet.cc PdfTable.cc

DrawResults_CPPFLAGS = -I$(srcdir)/../include $(ROOT_CFLAGS)

# libraries to be linked to the main target
DrawResults_LDADD = $(ROOT_LIBS)

# all source files which go into th
StoreGraphs_SOURCES = StoreGraphs.cc PdfTable.cc

StoreGraphs_CPPFLAGS = -I$(srcdir)/../include $(ROOT_CFLAGS)

# libraries to be linked to the main target
StoreGraphs_LDADD = $(ROOT_LIBS)
