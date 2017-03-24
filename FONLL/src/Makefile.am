AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libmyfonll.la

# Force recompilation if include files change
*.o: ../../include/couplings.inc

libmyfonll_la_SOURCES = FONLL_wrap.f
AM_FFLAGS = -I$(srcdir)/../../include $(NOAUTOFCFLAG) -ffixed-line-length-132

