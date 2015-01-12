#AUTOMAKE_OPTIONS = foreign

if HAVE_YAML
SUBDIRS = ../../pdf2yaml

# the main target, i.e. the executable to build
bin_PROGRAMS = postproc

# source files for the DrawPdfs executable
postproc_SOURCES = post_proc.c profile.c rotate.c custom.c symmetrize.c reweight.c

postproc_CFLAGS = -O0 -g -I$(srcdir)/../../include
postproc_CFLAGS += -Wall

# libraries to be linked to the main target
postproc_LDADD = -L../../pdf2yaml/ -lpdf2yaml -lyaml

endif
