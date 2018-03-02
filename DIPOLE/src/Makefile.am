AUTOMAKE_OPTIONS = foreign
dist_noinst_HEADERS = xsect_dipole_proton.h
#
# The name of the library
noinst_LTLIBRARIES = libmydipole.la

# Force recompile
fit_sf_bes.o: ../../include/steering.inc ../../include/dipole.inc

# Next files will be processed by g77 (and automatically gcc for *.c) compiler
libmydipole_la_SOURCES = fit_simps.f \
			fit_sf_bes.f \
			xsect_dipole_proton.c
AM_FFLAGS = -I$(srcdir)/../../include $(NOAUTOFCFLAG)
# -finit-local-zero -ffixed-line-length-132
