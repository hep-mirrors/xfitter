AUTOMAKE_OPTIONS = foreign

if ENABLE_HVQMNR

# ROOT is needed
if HAVE_ROOT

# The name of the library
  noinst_LIBRARIES = libHVQMNR.a

# Next files will be processed by g77 (and automatically gcc for *.c, g++ for *.cxx) compiler
  libHVQMNR_a_SOURCES = hvqcrsx.f hvqmnr.f hvqmnr_calculation.cxx hvqmnr_grid.cxx hvqmnr_mnr.cxx hvqmnr_frag.cxx

# C++ compiler flags
  AM_CXXFLAGS = -I$(srcdir)/../../include -I$(srcdir)/../include $(ROOT_CFLAGS) $(LHAPDF_CPPFLAGS)

# FORTRAN compiler flags
  AM_FFLAGS = -I$(srcdir)/../../include $(LHAPDF_CPPFLAGS) -xf77-cpp-input

else
  @echo "ROOT is needed for HVQMNR!"
  exit 1
endif

endif
