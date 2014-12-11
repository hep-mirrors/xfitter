#AUTOMAKE_OPTIONS = foreign

if HAVE_YAML

# the main target, i.e. the executable to build
bin_PROGRAMS = postproc

# source files for the DrawPdfs executable
postproc_SOURCES = post_proc.c profile.c rotate.c custom.c

postproc_CFLAGS = -g -I$(srcdir)/../../include $(YAML_CFLAGS)

# libraries to be linked to the main target
postproc_LDADD = -L../../pdf2yaml/ -lpdf2yaml $(YAML_LIBS)

endif
