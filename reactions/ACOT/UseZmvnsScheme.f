C----------------------------------------------------------------
C> \brief Calculates F2, FL, xF3, F2gamma and FLgamma using ZMVFNS from QCDNUM
C> \param[out] f2, fl, xf3, f2gamma, flgamma structure functions
C> \param[in] q2, x kinematic bin
C> \param[in] npts total number of points
C> \param[in] polarity of the lepton beam
C> \param[in] charge of the lepton beam
C> \param[in] XSecType DIS process type
C
C  Created by Krzysztof Nowak, 31/01/2012
C---------------------------------------------------------------
      subroutine UseZmvnsScheme(f2, fl, xf3, f2gamma, flgamma,
     $     q2, x, npts, polarity, charge, XSecType, local_hfscheme)

      implicit none
#include "steering.inc"
#include "couplings.inc"
#include "qcdnumhelper.inc"


C Input:
      double precision X(NPMaxDIS),Q2(NPMaxDIS)
      double precision charge, polarity
      integer npts, local_hfscheme
c     ************* PATCH
c      character*(*), intent(in) :: XSecType
      character*(10) XSecType

C Output: 
      double precision F2(NPMaxDIS), FL(NPMaxDIS), xF3(NPMaxDIS)
      double precision F2gamma(NPMaxDIS), FLgamma(NPMaxDIS)
      double precision F2c(NPMaxDIS),FLc(NPMaxDIS)
      double precision F2b(NPMaxDIS),FLb(NPMaxDIS)

C--------------------------------------------------------
C Temporary variables:
      double precision F2m(NPMaxDIS),xF3m(NPMaxDIS),FLm(NPMaxDIS)
      integer i
      double precision ve,ae,au,ad,vu,vd,A_u,A_d,B_u,B_d,pz


C EW param

      double precision sin2th_eff, xkappa, epsilon
      double precision deltar,sweff, sin2thw2
      double precision cau, cad, cvu, cvd
C--------------------------------------------------------
c     Fred tweaks: March 2022
      logical first
      data first /.true./
      save first
C--------------------------------------------------------
      EWFIT=0

C QCDNUM ZMVFNS, caclulate FL, F2 and xF3 for d- and u- type quarks all bins:

      if(XSecType.eq.'CCDIS') then
         if (charge.gt.0) then
            CALL ZMSTFUN(1,CCEP2F,X,Q2,FL,npts,0)
            CALL ZMSTFUN(2,CCEP2F,X,Q2,F2,npts,0)
            CALL ZMSTFUN(3,CCEP3F,X,Q2,XF3,npts,0)    
         else
            CALL ZMSTFUN(1,CCEM2F,X,Q2,FL,npts,0)
            CALL ZMSTFUN(2,CCEM2F,X,Q2,F2,npts,0)      
            CALL ZMSTFUN(3,CCEM3F,X,Q2,XF3,npts,0) 
         endif
      elseif (XSecType.eq.'NCDIS'.or.XSecType.eq.'CHARMDIS'
     $        .or.XSecType.eq.'F2'
     $        .or.XSecType.eq.'FL'.or.XSecType.eq.'BEAUTYDIS') then
C     u-type ( u+c ) contributions 
         CALL ZMSTFUN(1,CNEP2F,X,Q2,FL,npts,0)
         CALL ZMSTFUN(2,CNEP2F,X,Q2,F2,npts,0)
         CALL ZMSTFUN(3,CNEP3F,X,Q2,XF3,npts,0)    
         
C     d-type (d + s + b) contributions
         CALL ZMSTFUN(1,CNEM2F,X,Q2,FLm,npts,0)
         CALL ZMSTFUN(2,CNEM2F,X,Q2,F2m,npts,0)
         CALL ZMSTFUN(3,CNEM3F,X,Q2,XF3m,npts,0) 
      else
         print *, 'UseZmvnsScheme, XSecType',XSecType,
     $        'not supported'
         stop
      endif

c     for NC needs to combine F2p with F2m etc.        

      if(XSecType.eq.'NCDIS'.or.XSecType.eq.'CHARMDIS'.or.
     $     XSecType.eq.'F2'.or.
     $     XSecType.eq.'FL'.or.XSecType.eq.'BEAUTYDIS') then

c ============================================================
c ============================================================
         do i=1, npts

c     ============================================================
            if(EWFIT==0) then
C     
C EW couplings of the electron
C
               ve = -0.5d0 + 2.*sin2thw
               ae = -0.5d0         
            
C
C and quarks
C         
         
               au = 0.5d0
               ad = -0.5d0
            
               vu = au - (4.d0/3.d0)*sin2thw
               vd = ad + (2.d0/3.d0)*sin2thw
               
               PZ = 4.d0 * sin2thw * cos2thw * (1.+Mz**2/Q2(i))
c ============================================================
            else 

               if(first) then
                  write(6,*) " error: wrap_ew()  not implemented "
c                 first=.false.
               endif
c               stop 
c               call wrap_ew(q2(i),sweff,deltar,cau,cad,cvu,cvd,polarity,charge)
               sin2thw2 = 1.d0 - MW**2/MZ**2
               sin2th_eff = sweff
               xkappa = sin2th_eff/sin2thw2
               epsilon = xkappa -1.0
               ve = -0.5d0 + 2.*sin2th_eff
               ae = -0.5d0
               
               vu = cvu - (4.d0/3.d0)*epsilon*sin2thw2
               vd = cvd + (2.d0/3.d0)*epsilon*sin2thw2
               au = cau
               ad = cad
*     
*     Feed the EW parameters to APFEL 
*
               if (mod(local_hfscheme,10).eq.5) then
                  write(6,*) " *** NOT IMPLEMENTED "
                  STOP
               ENDIF               
c$$$               if (mod(local_hfscheme,10).eq.5) then
c$$$                  call SetSin2ThetaW(sin2th_eff)
c$$$                  call SetPropagatorCorrection(deltar)
c$$$                  call SetEWCouplings(vd,vu,ad,au)
c$$$               endif
               
C     Propagator factor PZ
               PZ = 4.d0*sin2thw2*(1.d0 - sin2thw2)*(1.+Mz**2/Q2(i))
               PZ = PZ*(1.d0 - Deltar)
            endif               
c ============================================================
c ============================================================
            PZ = 1./Pz
C     EW couplings of u-type and d-type quarks at the scale Q2
               
            if (charge.gt.0) then
               A_u = e2u        ! gamma
     $              + (-ve-polarity*ae)*PZ*2.*euq*vu !gamma-Z
     $              + (ve**2 + ae**2+2*polarity*ve*ae)*PZ**2*(vu**2+au**2) !Z
               
               A_d = e2d 
     $              + (-ve-polarity*ae)*PZ*2.*edq*vd 
     $              + (ve**2 + ae**2+2*polarity*ve*ae)*PZ**2*(vd**2+ad**2)
               
               B_u = (ae+polarity*ve)*PZ*2.*euq*au !gamma-Z
     $              + (-2.*ve*ae-polarity*(ve**2+ae**2))*(PZ**2)*2.*vu*au !Z
               B_d = (ae+polarity*ve)*PZ*2.*edq*ad 
     $              + (-2.*ve*ae-polarity*(ve**2+ae**2))*(PZ**2)*2.*vd*ad
            else
               A_u = e2u        ! gamma
     $              + (-ve+polarity*ae)*PZ*2.*euq*vu !gamma-Z
     $              + (ve**2 + ae**2-2*polarity*ve*ae)*PZ**2*(vu**2+au**2) !Z
               
               A_d = e2d 
     $              + (-ve+polarity*ae)*PZ*2.*edq*vd 
     $              + (ve**2 + ae**2-2*polarity*ve*ae)*PZ**2*(vd**2+ad**2)
               
               B_u = (-ae+polarity*ve)*PZ*2.*euq*au !gamma-Z
     $              + (2.*ve*ae-polarity*(ve**2+ae**2))*(PZ**2)*2.*vu*au !Z
               B_d = (-ae+polarity*ve)*PZ*2.*edq*ad 
     $              + (2.*ve*ae-polarity*(ve**2+ae**2))*(PZ**2)*2.*vd*ad
               
            endif
            
cv for polarised case should reduce to:
cv         A_u = e2u - ve*PZ*2.*euq*vu +(ve**2 + ae**2)*PZ**2*(vu**2+au**2)
cv         A_d = e2d - ve*PZ*2.*edq*vd +(ve**2 + ae**2)*PZ**2*(vd**2+ad**2)
cv         B_u = -ae*PZ*2.*euq*au + 2.*ve*ae*(PZ**2)*2.*vu*au
cv         B_d = -ae*PZ*2.*edq*ad + 2.*ve*ae*(PZ**2)*2.*vd*ad

            
            F2Gamma(i) = 4.D0/9.D0 * F2(i)  + 1.D0/9.D0 * F2m(i)
            FLGamma(i) = 4.D0/9.D0 * FL(i)  + 1.D0/9.D0 * FLm(i)
            XF3(i)  = B_U*XF3(i)  + B_D*XF3m(i)
            F2(i)   = A_U*F2(i)   + A_D*F2m(i)
            FL(i)   = A_U*FL(i)   + A_D*FLm(i)
            
         enddo
      elseif(XSecType.eq.'CCDIS') then
         do i=1,npts
            F2Gamma(i) = 4.D0/9.D0 * F2(i)  + 1.D0/9.D0 * F2m(i)
            FLGamma(i) = 4.D0/9.D0 * FL(i)  + 1.D0/9.D0 * FLm(i)
         enddo
      else 
         print *, 'UseZmvnsScheme, XSecType',XSecType,
     $        'not supported'
         stop
      endif
      
      end
