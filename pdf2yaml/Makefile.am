AUTOMAKE_OPTIONS = foreign

if HAVE_YAML

noinst_LIBRARIES = libpdf2yaml.a

libpdf2yaml_a_SOURCES = c2yaml.c  \
list.c \
pdf2yaml.c

AM_CFLAGS = -I$(srcdir)/../include

endif
