AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libnum_utils.a
pkginclude_HEADERS = c205body.inc  d107rd.inc  fact.inc  feqn.inc  finv.inc

libnum_utils_a_SOURCES = dchpws.F ddilog.F dgammf.F ranlux.F rnorml.F rnormx.F \
                         rnpssn.F snleq64.F dbesk0.F dgquad.F d107d1.F \
			 ranmar.F \
			 dfact.F dinv.F deqinv.F dfeqn.F deqn.F dfinv.F \
			 rzero.F datime.F \
                   ucopy.F timex.F timest.F rgs56p.F radapt.F expint.F \
                   divdif.F dgs56p.F dgauss.F dgausn.F dfreq.F dadapt.F abend.F                                                            
AM_FFLAGS = -I$(srcdir)/../../include -ffixed-line-length-132
#AM_FFLAGS = -I$(srcdir)/../../include -fno-automatic -finit-local-zero -ffixed-line-length-132
if ENABLE_OPENMP
AM_FFLAGS+= -fopenmp
endif

