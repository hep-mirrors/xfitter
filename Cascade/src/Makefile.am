AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libcasutil.a

if ENABLE_UPDF
  libcasutil_a_SOURCES = cascade-util.F cascade-misc.F
else
  libcasutil_a_SOURCES = cascade-dummy.F
endif


AM_FFLAGS = -I$(srcdir)/../../include -I$(CASCADE_ROOT)/include/cascade  -fno-automatic -ffixed-line-length-132
