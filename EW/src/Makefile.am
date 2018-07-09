AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libmyew.la

#
wrap_ew.o: ../../include/couplings.inc ../../include/polarity.inc ../../include/steering.inc

libmyew_la_SOURCES = formff.f wrap_ew.f hadr5n12.f
AM_FFLAGS = -I$(srcdir)/../../include $(NOAUTOFCFLAG)
#AM_FFLAGS = -I$(srcdir)/../../include -fno-automatic -finit-local-zero -ffixed-line-length-132


