#AUTOMAKE_OPTIONS = foreign

# the main target, i.e. the executable to build
bin_PROGRAMS = DrawResults

# all source files which go into the main target
DrawResults_SOURCES = H1FitterDrawResults.cc H1FitterPainter.cc \
  H1FitterOutput.cc DataSet.cc

DrawResults_CPPFLAGS = -I$(srcdir)/../include $(ROOT_CFLAGS)

# libraries to be linked to the main target
DrawResults_LDADD = $(ROOT_LIBS)
