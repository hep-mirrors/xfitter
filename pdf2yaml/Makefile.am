AUTOMAKE_OPTIONS = foreign

if HAVE_YAML

noinst_LIBRARIES = libpdf2yaml.a
libpdf2yaml_a_SOURCES = c2yaml.c list.c pdf2yaml.c

noinst_PROGRAMS = tests
tests_SOURCES= tests.c
tests_LDADD= -L. -lpdf2yaml -lyaml
tests_DEPENDENCIES= libpdf2yaml.a

AM_CFLAGS = -fpic -I$(srcdir)/../include

#TEST= tests
endif
