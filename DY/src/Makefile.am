AUTOMAKE_OPTIONS = foreign

# The name of the library
noinst_LIBRARIES = libmydy.a

# Next files will be processed by g77 (and automatically gcc for *.c) compiler
libmydy_a_SOURCES = PhysPar.cc \
			zee_integ.cc \
			wpenu_integ.cc \
			wmenu_integ.cc \
			dy_integ_utils.cc \
			wpm_eta_bins.cc \
			z_eta_bins.cc \
			z_eta_bc.cc \
			z_y_bins.cc \
			z_y_bc.cc
AM_FFLAGS = -I$(srcdir)/../../include -fno-automatic -finit-local-zero -ffixed-line-length-132

AM_CXXFLAGS = -I$(srcdir)/../../include -O -Wall -fPIC -Wno-deprecated
