AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libmysacot.a

# Force recompilation if include files change:
*.o: ../../include/couplings.inc ../../include/ntot.inc ../../include/qcdnumhelper.inc ../../include/steering.inc

libmysacot_a_SOURCES =  0mDisPac.F 1uti.F \
                        alfaem.F alfewk.F ewcpl2an.F ewcplsc.F filcpl.F setewk.F stewcpl2.F stewcpl.F \
                        alambd.F alamf.F alepi.F alfset.F alphem.F alpior.F alpqcd.F amass.F amhatf.F amumin.F anom.F ch.F cnvl1.F evoluf.F g.F lamcwz.F masssync.F mtmul.F namqcd.F nfl.F nfltot.F parqcd.F qcdget.F qcdin.F qcdout.F qcdset.F qwarn.F qzbrnt.F rtalf.F SetAmhatf.F setl1.F setlam.F setqcd.F sud.F trnlam.F zbrlam.F zcnvlm.F \
                        pdf.F pdfh.F \
                        Acot05.F GmJSa.F GmJS.F GotSK1Q.F KinForm.F lhXsec.F SFnGm.F SfZm05.F StrFnGm.F ZmNNLO.F \
                        sf_sacotchi_wrap.F
AM_FFLAGS = -I$(srcdir)/../../include $(NOAUTOFCFLAG) -ffixed-line-length-132

