AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libmyfonll.a

# Force recompilation if include files change
*.o: ../../include/couplings.inc

libmyfonll_a_SOURCES = FONLL_wrap.f
AM_FFLAGS = -I$(srcdir)/../../include -fno-automatic -ffixed-line-length-132

