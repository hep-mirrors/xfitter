AUTOMAKE_OPTIONS = foreign
pkginclude_HEADERS = APSCOM6. CONSTCOM. OMECOM. PDFCOM. PRECCOM.
# Force recompilation if include files are modified:
*.o: APSCOM6. CONSTCOM. OMECOM. PDFCOM. PRECCOM.


noinst_LIBRARIES = libmyabm.a

libmyabm_a_SOURCES = sf_abkm_wrap.f alphasvfn.f asympcoef.f cq.f dishq.f dislt.f disqpm.f f2ccbmsn.f gauss.f hqcoef.f hqnnlocoef.f lpcoef.f ome.f spline.f wrapNNPDFgrid.f NNPDFwrapinit.f grid.F initgridconst.F
AM_FFLAGS = -I$(srcdir)/../../include -fno-automatic
#AM_FFLAGS = -I$(srcdir)/../../include -fno-automatic -finit-local-zero -ffixed-line-length-132
