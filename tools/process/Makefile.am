#AUTOMAKE_OPTIONS = foreign

if HAVE_YAML
SUBDIRS = ../../pdf2yaml

# the main target, i.e. the executable to build
bin_PROGRAMS = xfitter-process

# source files for the xfitter-process executable
xfitter_process_SOURCES = xfitter-process.c utils.c profile.c rotate.c custom.c symmetrize.c reweight.c 

pkginclude_HEADERS = rotate.h  utils.h  xfitter-process.h

xfitter_process_CFLAGS = -I$(srcdir)/../../include
xfitter_process_CFLAGS += -Wall
xfitter_process_CPPFLAGS =

if ENABLE_LHAPDF
  xfitter_process_SOURCES += interpolation.cc
  xfitter_process_SOURCES += combine.c
  xfitter_process_CPPFLAGS += $(LHAPDF_CPPFLAGS)
endif


# libraries to be linked to the main target
xfitter_process_LDADD = -L../../pdf2yaml/ -lpdf2yaml -lyaml

if ENABLE_LHAPDF 
  xfitter_process_LDADD+=$(LHAPDF_LDFLAGS)
endif


if HAVE_ROOT
   xfitter_process_SOURCES += rootplot.cc 
   xfitter_process_CPPFLAGS += -I$(srcdir)/../../include  $(ROOT_CFLAGS) -w
   xfitter_process_LDADD += $(ROOT_LIBS) 
else 
   xfitter_process_SOURCES += rootplot_dummy.cc
endif


endif
