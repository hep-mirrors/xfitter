AUTOMAKE_OPTIONS = foreign
# Explicit depence for all files:
*.o: ../include/*.inc 


AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../interfaces/include -I$(srcdir)/../FastNLO/include
AM_CXXFLAGS += -I$(srcdir)/../common -I$(srcdir)/../common/linalg
AM_CXXFLAGS += -I$(srcdir)/../tools/draw/include


lib_LTLIBRARIES = libagtheory.la

if ENABLE_APPLGRID
  libagtheory_la_SOURCES= ReactionAPPLGrid.cc
else 
#  libxfitter_a_SOURCES+= applgrids_dum.f 
# REMOVE ftheor_eval_dum.cc
endif

# FORTRAN compiler flags to be used for the main target
AM_FFLAGS = $(DEFS) -cpp -I$(srcdir)/../include -fno-automatic -I$(srcdir)/../minuit/src

# 
AM_CPPFLAGS = -I$(srcdir)/../include 

if ENABLE_APPLGRID 
#  libagtheory_la_LDFLAGS= $(APPLGRID_LDFLAGS) $(GFRTLIB)
endif

