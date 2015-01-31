

      Subroutine GetH1qcdfuncOutput(charge_in,polarity_in)
!(xsh,q2sh, charge, polarity
!     $     f2pRT, f2qcd, f2gammash,f2gam, f2gamz,f2z,
!     $     flpRT, flqcd, flgammash,flgam, flgamz,flz,
!     $     xf3qcd, xf3gamz,xf3z)
C----------------------------------------------------------------
C
C---------------------------------------------------------------
      implicit none

#include "steering.inc"
#include "couplings.inc"
#include "qcdnumhelper.inc"


 
      double precision charge_in, polarity_in
      double precision xsh,q2sh

cv shiraz 

      integer nxboo, nq2boo
      parameter (nxboo=66)
      double precision xboo(nxboo)
      data xboo/
     $     0.00001d0,
     $     0.000010588d0, 0.00001409d0,  0.00001875d0,  0.000024952d0,
     $     0.000033205d0, 0.000044187d0, 0.000058801d0, 0.000078249d0,
     $     0.00010413d0,  0.00013857d0,  0.00018440d0,  0.00024539d0,
     $     0.00032655d0,  0.00043456d0,  0.00057829d0,  0.00076955d0,
     $     0.00102410d0,  0.00136280d0,  0.00181350d0,  0.00241330d0,
     $     0.00321150d0,  0.00427370d0,  0.00568720d0,  0.00756830d0,
     $     0.00898370d0,  0.010664d0,    0.012658d0,    0.015025d0,
     $     0.017835d0,    0.021171d0,    0.025130d0,    0.029830d0,
     $     0.035409d0,    0.042031d0,    0.049891d0,    0.055932d0,
     $     0.062704d0,    0.070297d0,    0.078809d0,    0.088351d0,
     $     0.099049d0,    0.11104d0,     0.12449d0,     0.13956d0,
     $     0.15646d0,     0.17540d0,     0.19664d0,     0.22045d0,
     $     0.24715d0,     0.27707d0,     0.31062d0,     0.34823d0,
     $     0.39039d0,     0.43767d0,     0.49066d0,     0.55007d0,
     $     0.61667d0,     0.65294d0,     0.69134d0,     0.73200d0,
     $     0.77505d0,     0.82064d0,     0.86890d0,     0.92000d0,
     $     0.989999d0/

      parameter (nq2boo=37)
      double precision q2oo(nq2boo)
      data q2oo /  
     $     1.d0, 1.20d0, 1.500d0, 2.000d0, 2.500d0, 4.000d0,
     &     6.5d0, 10.00d0, 13.300d0, 17.800d0, 23.700d0,
     &     31.60d0, 42.00d0, 56.000d0, 75.000d0, 100.00d0,
     &     133.0d00,178.0d0,250.0d0,350.0d0,500.0d0,
     &     650.0d00,1000.0d0,1200.0d0,2000.0d0,3000.0d0,
     &     5000.0d0,8000.0d0,10000.d0,12000.d0,20000.d0,
     &     30000.d0,50000.d0,100000.d0,200000.d0,500000.d0,
     &     1000000.d0/

      integer ns, nshiraz, ishx,ishq2, nsh2
      parameter (ns=2442)
      double precision dxoo(ns),dq2oo(ns)



      double precision F2CCepSH(ns), F2ccemSH(ns),xF3ccemSH(ns)
      double precision F2ncepSH(ns), F2ncemSH(ns),xF3ccepSH(ns)
      double precision FLCCepSH(ns), FLccemSH(ns),xF3ncemSH(ns)
      double precision FLncepSH(ns), FLncemSH(ns),xF3ncepSH(ns)

      

cv shiraz


      double precision ve,ae,au,ad,vu,vd,A_u,A_d,B_u,B_d,pz


      double precision F2,xF3,FL
      double precision F2qcd,xF3qcd,FLqcd
      double precision FLGam, F2Gam, xF3Gam
      double precision FLGamz, F2Gamz, xF3Gamz
      double precision FLzz, F2zz, xF3zz


      double precision x_in,q2_in




C RT:
      Double precision f2RT_ou,flRT_ou,
     $     f2cRT_ou,flcRT_ou,f2bRT_ou,
     $     flbRT_ou
      logical UseKFactors

c EW param

      double precision sin2th_eff, xkappa, epsilon
      double precision deltar,sweff, sin2thw2
      double precision cau, cad, cvu, cvd


      open (74,file=TRIM(OutDirName)//'/f2nc.qcdfunc',status='unknown')
      open (84,file=TRIM(OutDirName)//'/flnc.qcdfunc',status='unknown')
      open (94,file=TRIM(OutDirName)//'/xf3nc.qcdfunc',status='unknown')
      open (96,file=TRIM(OutDirName)//'/ccem.qcdfunc', status='unknown')
      open (98,file=TRIM(OutDirName)//'/ccep.qcdfunc',status='unknown')

C prepare bins:
      nshiraz=0
      do ishq2=1,37
         q2_in=q2oo(ishq2)
         do ishx=1, 66
            x_in=xboo(ishx)
            nshiraz=nshiraz+1
            dxoo(nshiraz)=x_in
            dq2oo(nshiraz)=q2_in
         enddo
      enddo

C QCDNUM ZMVFNS, caclulate FL, F2 and xF3 for d- and u- type quarks all bins:

cv NC
      CALL ZMSTFUN(2,CNEP2F,DXoo,DQ2oo,F2NcepSH,nshiraz,0)
      CALL ZMSTFUN(2,CNEM2F,DXoo,DQ2oo,F2NCEMSH,nshiraz,0)
            
      CALL ZMSTFUN(1,CNEP2F,DXoo,DQ2oo,FLNcepSH,nshiraz,0)
      CALL ZMSTFUN(1,CNEM2F,DXoo,DQ2oo,FLNCEMSH,nshiraz,0)
      
      CALL ZMSTFUN(3,CNEP3F,DXoo,DQ2oo,XF3NcepSH,nshiraz,0)    
      CALL ZMSTFUN(3,CNEM3F,DXoo,DQ2oo,XF3NCEMSH,nshiraz,0) 
            
cv CC

c      if (charge_in.gt.0) then
      CALL ZMSTFUN(1,CCEP2F,DXoo,DQ2oo,FLCCepSH,nshiraz,0)
      CALL ZMSTFUN(2,CCEP2F,DXoo,DQ2oo,F2CCepSH,nshiraz,0)
      CALL ZMSTFUN(3,CCEP3F,DXoo,DQ2oo,XF3CCepSH,nshiraz,0)    
c      else
      CALL ZMSTFUN(1,CCEM2F,DXoo,DQ2oo,FLCCemSH,nshiraz,0)
      CALL ZMSTFUN(2,CCEM2F,DXoo,DQ2oo,F2CCemSH,nshiraz,0)      
      CALL ZMSTFUN(3,CCEM3F,DXoo,DQ2oo,XF3CCemSH,nshiraz,0) 
c      endif

            
CV call GAMMA for RT
c            
c      CALL ZMSTFUN(1,PROT,DXoo,DQ2oo,FLGAMMASH,nshiraz,0)      
c      CALL ZMSTFUN(2,PROT,DXoo,DQ2oo,F2GAMMASH,nshiraz,0)
      

      ve = -0.5d0 + 2.*sin2thw
      ae = -0.5d0         

C
C and quarks
C         
         
      au = 0.5d0
      ad = -0.5d0
                  
      vu = au - (4.d0/3.d0)*sin2thw
      vd = ad + (2.d0/3.d0)*sin2thw
      
      print*,'polarity', polarity_in



C
C xF3, F2, FL from QCDNUM:
C
      nsh2=0
      
      
      do ishq2=1,37
         q2_in=q2oo(ishq2)
         do ishx=1, 66
            x_in=xboo(ishx)
            nsh2=nsh2+1

            PZ = 4.d0 * sin2thw * cos2thw * (1.+Mz**2/Q2_in)
            PZ = 1./Pz
            if (charge_in.gt.0) then
               A_u = e2u        ! gamma
     $              + (-ve-polarity_in*ae)*PZ*2.*euq*vu !gamma-Z
     $              + (ve**2 + ae**2+2*polarity_in*ve*ae)*PZ**2*(vu**2+au**2) !Z
               
               A_d = e2d 
     $              + (-ve-polarity_in*ae)*PZ*2.*edq*vd 
     $              + (ve**2 + ae**2+2*polarity_in*ve*ae)*PZ**2*(vd**2+ad**2)
               
               B_u = (ae+polarity_in*ve)*PZ*2.*euq*au !gamma-Z
     $              + (-2.*ve*ae-polarity_in*(ve**2+ae**2))*(PZ**2)*2.*vu*au !Z
               B_d = (ae+polarity_in*ve)*PZ*2.*edq*ad 
     $              + (-2.*ve*ae-polarity_in*(ve**2+ae**2))*(PZ**2)*2.*vd*ad
            else
               A_u = e2u        ! gamma
     $              + (-ve+polarity_in*ae)*PZ*2.*euq*vu !gamma-Z
     $              + (ve**2 + ae**2-2*polarity_in*ve*ae)*PZ**2*(vu**2+au**2) !Z
               
               A_d = e2d 
     $        + (-ve+polarity_in*ae)*PZ*2.*edq*vd 
     $              + (ve**2 + ae**2-2*polarity_in*ve*ae)*PZ**2*(vd**2+ad**2)
               
               B_u = (-ae+polarity_in*ve)*PZ*2.*euq*au !gamma-Z
     $              + (2.*ve*ae-polarity_in*(ve**2+ae**2))*(PZ**2)*2.*vu*au !Z
               B_d = (-ae+polarity_in*ve)*PZ*2.*edq*ad 
     $              + (2.*ve*ae-polarity_in*(ve**2+ae**2))*(PZ**2)*2.*vd*ad
               
            endif

cv for un polarised case should reduce to:
cv         A_u = e2u - ve*PZ*2.*euq*vu +(ve**2 + ae**2)*PZ**2*(vu**2+au**2)
cv         A_d = e2d - ve*PZ*2.*edq*vd +(ve**2 + ae**2)*PZ**2*(vd**2+ad**2)
cv         B_u = -ae*PZ*2.*euq*au + 2.*ve*ae*(PZ**2)*2.*vu*au
cv         B_d = -ae*PZ*2.*edq*ad + 2.*ve*ae*(PZ**2)*2.*vd*ad



cv            mode=1
            XF3qcd = B_U*XF3NCEPSH(nsh2) + B_D*XF3NCEMSH(nsh2)
            F2qcd   = A_U*F2NCEPSH(nsh2)   + A_D*F2NCEMSH(nsh2)
            FLqcd   = A_U*FLNCEPSH(nsh2)   + A_D*FLNCEMSH(nsh2)
              
*     save QCDNUM values (for debugging purposes)
c            f2qcd = f2
c            flqcd = fl
c            xf3qcd = xf3

            F2gam=e2u*F2NCEPSH(nsh2)+e2d*F2NCEMSH(nsh2)
            F2gamz=2.d0*euq*vu*F2NCEPSH(nsh2)+2.d0*edq*vd*F2NCEMSH(nsh2)
            F2zz=(vu**2+au**2)*F2NCEPSH(nsh2)+(vd**2+ad**2)*F2NCEMSH(nsh2)
            
            Flgam=e2u*FLNCEPSH(nsh2)+e2d*FLNCEMSH(nsh2)
            Flgamz=2.d0*euq*vu*FLNCEPSH(nsh2)+2.d0*edq*vd*FLNCEMSH(nsh2)
            Flzz=(vu**2+au**2)*FLNCEPSH(nsh2)+(vd**2+ad**2)*FLNCEMSH(nsh2)
            
            XF3gam=0.d0
            xF3gamz=2.d0*euq*au*XF3NCEPSH(nsh2)+2.d0*edq*ad*XF3NCEMSH(nsh2)
            xF3zz=2.d0*vu*au*XF3NCEPSH(nsh2) + 2.d0*vd*ad*XF3NCEMSH(nsh2)



            if (F2gam.lt.0.d0)   F2gam=0.d0
            if (F2gamz.lt.0.d0)  F2gamz=0.d0
            if (F2zz.lt.0.d0)    F2zz=0.d0


            if (FLgam.lt.0.d0)   FLgam=0.d0
            if (FLgamz.lt.0.d0)  FLgamz=0.d0
            if (FLzz.lt.0.d0)    FLzz=0.d0
                  

            if (xF3gam.lt.0.d0)   xF3gam=0.d0
            if (xF3gamz.lt.0.d0)  xF3gamz=0.d0
            if (xF3zz.lt.0.d0)    xF3zz=0.d0


           call  mstwnc_wrap(
     $        x_in,q2_in,1,f2RT_ou,
     $        f2cRT_ou,f2bRT_ou,flRT_ou,flcRT_ou,flbRT_ou
           ! Input:
     $          ,1,nsh2    ! fcn flag, data point index
     $          ,F2Gam,FLGam
     $          ,.false.
     $          )
           
c           print*,'NC SH F2', x,q2, f2pRT, f2qcd, f2gam, f2gamz,f2zz
c           print*,'NC SH FL', x,q2, flpRT, flqcd, flgam, flgamz,flzz
c           print*,'NC SH F3', x,q2, xf3qcd, xf3gamz,xf3zz
           if (f2rt_ou.lt.0.d0)  f2RT_ou=0.d0
           if (flrt_ou.lt.0.d0)  flRT_ou=0.d0
           if (f2qcd.lt.0.d0)  f2qcd=0.d0
           if (flqcd.lt.0.d0)  flqcd=0.d0
           if (xf3qcd.lt.0.d0) xf3qcd=0.d0


           
           write(74,888) x_in,q2_in, f2RT_ou, f2qcd,
     $          f2gam, f2gam, f2gamz,f2zz
           write(84,888) x_in,q2_in, flRT_ou, flqcd, 
     $          flgam, flgam, flgamz,flzz
           write(94,889) x_in,q2_in, xf3qcd, xf3gamz,xf3zz


           if (f2ccemSH(nsh2).lt.0.d0) f2ccemSH(nsh2)=0.d0
           if (f2ccepSH(nsh2).lt.0.d0) f2ccepSH(nsh2)=0.d0
           if (flccemSH(nsh2).lt.0.d0) flccemSH(nsh2)=0.d0
           if (flccepSH(nsh2).lt.0.d0) flccepSH(nsh2)=0.d0
           if (xf3ccemSH(nsh2).lt.0.d0) xf3ccemSH(nsh2)=0.d0
           if (xf3ccepSH(nsh2).lt.0.d0) xf3ccepSH(nsh2)=0.d0

           write(96,889) x_in, q2_in, f2ccemSH(nsh2), flccemSH(nsh2),
     $          xf3ccemSH(nsh2)
           write(98,889) x_in, q2_in, f2ccepSH(nsh2), flccepSH(nsh2),
     $          xf3ccepSH(nsh2)

 888       format(1x, G14.8, 2x, G14.5,6(2x, G16.8))
 889       format(1x, G14.8, 2x, G14.5,3(2x, G16.8))
           
           enddo
        enddo
      close(74)
      close(84)
      close(94)
      close(96)
      close(98)


      end



