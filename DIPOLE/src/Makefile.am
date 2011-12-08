AUTOMAKE_OPTIONS = foreign

# The name of the library
noinst_LIBRARIES = libmydipole.a

# Next files will be processed by g77 (and automatically gcc for *.c) compiler
libmydipole_a_SOURCES = fit_simps.f \
			fit_sf_bes.f \
			xsect_dipole_proton.c
AM_FFLAGS = -I$(srcdir)/../../include -fno-automatic
# -finit-local-zero -ffixed-line-length-132
