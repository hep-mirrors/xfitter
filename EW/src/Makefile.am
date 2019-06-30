AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libmyew.a

#
wrap_ew.o: ../../include/couplings.inc ../../include/polarity.inc ../../include/steering.inc

libmyew_a_SOURCES = formff.f wrap_ew.f hadr5n12.f
AM_FFLAGS = -I$(srcdir)/../../include $(NOAUTOFCFLAG)
#AM_FFLAGS = -I$(srcdir)/../../include -fno-automatic -finit-local-zero -ffixed-line-length-132


