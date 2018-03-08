
# Created by AddReaction.py on 2017-08-07

if ENABLE_HATHOR

# LHAPDF is required by Hathor-2.0
if ENABLE_LHAPDF
  
  # HATHOR_ROOT points to the directory containing libHathor.a:
  #   (1) in Hathor-1.5 this is the same common directory which contains all headers and sources,
  #   (2) in Hathor-2.0 the library is Hathor-2.0/lib while the headers are in Hathor-2.0/include,
  # therefore add both header places below
  HATHOR_CXXFLAGS = -I$(HATHOR_ROOT) -I$(HATHOR_ROOT)/../include

  # Hathor-2.0 requires LHAPDF
  HATHOR_CXXFLAGS += $(LHAPDF_CPPFLAGS)

  AM_CXXFLAGS = $(HATHOR_CXXFLAGS) -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated -O0

  lib_LTLIBRARIES = libhathor_xfitter.la
  libhathor_xfitter_la_SOURCES = ReactionHathor.cc HathorPdfxFitter.cc

  # add Hathor library(ies)
  libhathor_xfitter_la_LDFLAGS = $(HATHOR_LIBS)

endif

endif

dist_noinst_HEADERS = ../include ../yaml
