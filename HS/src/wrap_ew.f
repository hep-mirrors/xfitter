      subroutine wrap_ew
     $     (q2_in, sweff_out, deltar_out, 
     $     cau_out,cad_out,cvu_out, cvd_out,
     $     polarity_in, charge_in)


      include 'polarity.inc'


      double precision q2_in, sweff_out, deltar_out

      double precision cau_out, cad_out, cvd_out, cvu_out
      double precision DELTAR,AGF0,DRHOT,DALPMZ,XGMT,ALPQCD,BTOP4,DRPIW2
      double precision ALPFFQ,AKAPPA,GMUFFQ,SWEFF2,sm
      double precision epMz,epsin2thw,cau,cvu,cad,cvd

      COMMON /HDELTR/ DELTAR,AGF0,DRHOT,DALPMZ,XGMT,ALPQCD,BTOP4,DRPIW2
      COMMON /FORMFF/ ALPFFQ,AKAPPA,GMUFFQ,SWEFF2
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
c      COMMON /GSW/    SW,CW,SW2,CW2
c     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
c     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      common/eprc_manu/sm
      common/couplings/epsin2thw,epMz,
     $     cau,cvu,cad,cvd

      integer i 	
      double precision charge_in,polarity_In, polari, llept, lqua
      COMMON /PARAM/  POLARI,LLEPT,LQUA      



      polari = polarity_in

      if (charge_in.gt.0.) then
         LLEPT=1.
      else
         LLEPT=-1.
      endif     




      cau=cau_ew
      cad=cad_ew
      cvu=cvu_ew
      cvd=cvd_ew


      call eprc_effective(q2_in)

      sweff_out=sweff2
      deltar_out=deltar

      cau_out=cau
      cad_out=cad
      cvd_out=cvd
      cvu_out=cvu
      
cv      print*, 'am here', cau,p(95),p(91)

      RETURN
      END
