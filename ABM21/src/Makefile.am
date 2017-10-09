AUTOMAKE_OPTIONS = foreign
pkginclude_HEADERS = APSCOM6. CONSTCOM. PDFCOM. PRECCOM.
# Force recompilation if include files are modified:
*.o: APSCOM6. CONSTCOM. PDFCOM. PRECCOM.


noinst_LIBRARIES = libmyabm.a

libmyabm_a_SOURCES = sf_abkm_wrap.f alphasvfn.f asympcoef.f cq.f dishq.f dislt.f disqpm.f f2ccbmsn.f gauss.f hqcoef.f hqnnlocoef.f lpcoef.f ome.f spline.f grid.F initgridconst.F asy-hpcc.f hplog.f split.f s2nlo.f
AM_FFLAGS = -I$(srcdir)/../../include $(NOAUTOFCFLAG)

#AM_FFLAGS = -I$(srcdir)/../../include -fno-automatic -finit-local-zero -ffixed-line-length-132
