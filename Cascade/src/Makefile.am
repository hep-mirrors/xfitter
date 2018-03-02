AUTOMAKE_OPTIONS = foreign

dist_noinst_HEADERS = SMallx.inc cascade-util_nbx.inc
noinst_LTLIBRARIES = libcasutil.la

if ENABLE_UPDF
  libcasutil_la_SOURCES = cascade-util.F cascade-misc.F evolve_tmd.F  
else
  libcasutil_la_SOURCES = cascade-dummy.F
endif


AM_FFLAGS = -O0 -I$(srcdir)/../../include -I$(CASCADE_ROOT)/include/cascade $(NOAUTOFCFLAG)
#AM_FFLAGS = -O0 -I$(srcdir)/../../include -I$(CASCADE_ROOT)/include/cascade -finit-local-zero -fno-automatic -fno-align-commons -ffixed-line-length-132
#AM_FFLAGS = -O0 -I$(srcdir)/../../include -I$(CASCADE_ROOT)/include/cascade -finit-local-zero -fno-automatic -fno-align-commons 
#FFLAGS = -O0 -I$(srcdir)/../../include -I$(CASCADE_ROOT)/include/cascade -finit-local-zero -fno-automatic -fno-align-commons -ffixed-line-length-132
