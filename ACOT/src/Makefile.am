AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libmyacot.a

# Force recompilation if include files change:
*.o: ../../include/couplings.inc ../../include/ntot.inc ../../include/qcdnumhelper.inc ../../include/steering.inc

libmyacot_a_SOURCES =sf_acot_wrap.f  alphas70.f  cc116.f fgen122qcdnum.f  qcdnumPDF.f  tot130.f  unified121.f  utilities113.f  c2terms.f dinteg.f ingo.f ingoUtil.f
AM_FFLAGS = -I$(srcdir)/../../include $(NOAUTOFCFLAG)

