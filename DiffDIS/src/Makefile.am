AUTOMAKE_OPTIONS = foreign
pkginclude_HEADERS = JScif.h TRc.h trs.h 
# defparams.inc DDIS.coca xdr.c
# The name of the library
noinst_LIBRARIES = libDiffDIS.a


# Next files will be processed by g77 (and automatically gcc for *.c) compiler
libDiffDIS_a_SOURCES = Coca.cc   GrdHeader.cc  Param.cc    ddis.cc   fileut.cc  hdpdf.cc     qcdpars.cc     tblipol.cc\
	DataTable.cc  LambertW1.cc  dilog.cc  gZDPDF.cc  ifun.cc      qgrid_base.cc\
	FFstrFns.cc   MDipol.cc     aqcd.cc     emsg.cc   gauss.cc   pdf_base.cc  tblint.cc      xgrid_base.cc\
	FTN_sock_DDIS.cc\
	JScif.c  grvpi.c TRc.c\
	cheavy.f  coeff2.f  coeff3.f js_fns.f sfunws.f js_fns.f jacksmith.f \
	FASTLI2.f FASTLI3.f FASTS12.f TR_alpha.f WATE96.f WGPLG.f \
	grvpiho.f grvpilo.f GetPdf.f 	

AM_CXXFLAGS = -I$(srcdir)/../include -I$(srcdir)/../../common  -fPIC -Wno-deprecated $(ROOT_CFLAGS) -w
# AM_CPPFLAGS = -I$(srcdir)/../include -Wall -Wno-deprecated -fPIC $(ROOT_CFLAGS)
