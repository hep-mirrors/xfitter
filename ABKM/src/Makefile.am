AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libmyabkm.a

libmyabkm_a_SOURCES = sf_abkm_wrap.f alphasvfn.f asympcoef.f cq.f dishq.f dislt.f disqpm.f f2ccbmsn.f gauss.f hqcoef.f hqnnlocoef.f lpcoef.f ome.f spline.f wrapNNPDFgrid.f NNPDFwrapinit.f grid.F initgridconst.F
AM_FFLAGS = -I$(srcdir)/../../include -fno-automatic
#AM_FFLAGS = -I$(srcdir)/../../include -fno-automatic -finit-local-zero -ffixed-line-length-132

