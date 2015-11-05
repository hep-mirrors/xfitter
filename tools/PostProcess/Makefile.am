#AUTOMAKE_OPTIONS = foreign

if HAVE_YAML
SUBDIRS = ../../pdf2yaml

# the main target, i.e. the executable to build
bin_PROGRAMS = postproc

# source files for the DrawPdfs executable
postproc_SOURCES = post_proc.c utils.c profile.c rotate.c custom.c symmetrize.c reweight.c 


postproc_CFLAGS = -I$(srcdir)/../../include
postproc_CFLAGS += -Wall

# libraries to be linked to the main target
postproc_LDADD = -L../../pdf2yaml/ -lpdf2yaml -lyaml


if HAVE_ROOT
   postproc_SOURCES += rootplot.cc 
   postproc_CPPFLAGS = -I$(srcdir)/../../include  $(ROOT_CFLAGS) -w
   postproc_LDADD += $(ROOT_LIBS) 
else 
   postproc_SOURCES += rootplot_dummy.cc
endif

 
endif
