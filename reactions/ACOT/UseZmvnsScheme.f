C----------------------------------------------------------------
C> \brief Calculates F2, FL, xF3, F2gamma and FLgamma using ZMVFNS from QCDNUM
C> \param[out] f2, fl, xf3, f2gamma, flgamma structure functions
C> \param[in] q2, x kinematic bin
C> \param[in] npts total number of points
C> \param[in] polarity of the lepton beam
C> \param[in] charge of the lepton beam
C> \param[in] XSecType DIS process type
C
C     Created by Krzysztof Nowak, 31/01/2012
c     modified by Fred Olness: 31 march 2022
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
      double precision F2c(NPMaxDIS),FLc(NPMaxDIS), xF3c(NPMaxDIS)
      double precision F2b(NPMaxDIS),FLb(NPMaxDIS), xF3b(NPMaxDIS)

C--------------------------------------------------------
C Temporary variables:
      double precision F2m(NPMaxDIS),xF3m(NPMaxDIS),FLm(NPMaxDIS)
      double precision F2Mc(NPMaxDIS),FLMc(NPMaxDIS), xF3Mc(NPMaxDIS)
      double precision F2Mb(NPMaxDIS),FLMb(NPMaxDIS), xF3Mb(NPMaxDIS)
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

      
      double precision cnep2fC(-6:6),cnep2fB(-6:6)
      double precision cnem2fC(-6:6),cnem2fB(-6:6)
      double precision cnep3fC(-6:6),cnep3fB(-6:6)
      double precision cnem3fC(-6:6),cnem3fB(-6:6)
c                  tb bb  cb  sb  ub  db  g  d  u  s  c  b  t
      DATA CNEP2Fc/0., 0., 1., 0., 0., 0.,0.,0.,0.,0.,1.,0.,0./
      DATA CNEP2Fb/0., 0., 0., 0., 0., 0.,0.,0.,0.,0.,0.,0.,0./ !*** zero
      DATA CNEM2Fc/0., 0., 0., 0., 0., 0.,0.,0.,0.,0.,0.,0.,0./ !*** zero
      DATA CNEM2Fb/0., 1., 0., 0., 0., 0.,0.,0.,0.,0.,0.,1.,0./
      DATA CNEP3Fc/0., 0.,-1., 0.,-0., 0.,0.,0.,0.,0.,1.,0.,0./
      DATA CNEP3Fb/0., 0.,-0., 0.,-0., 0.,0.,0.,0.,0.,0.,0.,0./ !*** zero
      DATA CNEM3Fc/0.,-0., 0.,-0., 0.,-0.,0.,0.,0.,0.,0.,0.,0./ !*** zero
      DATA CNEM3Fb/0.,-1., 0.,-0., 0.,-0.,0.,0.,0.,0.,0.,1.,0./

      
C--------------------------------------------------------
c      EWFIT=0

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
         
         CALL ZMSTFUN(1,CNEP2Fc,X,Q2,FLc,npts,0)
         CALL ZMSTFUN(2,CNEP2Fc,X,Q2,F2c,npts,0)
         CALL ZMSTFUN(3,CNEP3Fc,X,Q2,XF3c,npts,0)    
         
         CALL ZMSTFUN(1,CNEP2Fb,X,Q2,FLb,npts,0)
         CALL ZMSTFUN(2,CNEP2Fb,X,Q2,F2b,npts,0)
         CALL ZMSTFUN(3,CNEP3Fb,X,Q2,XF3b,npts,0)    
         
C     d-type (d + s + b) contributions
         CALL ZMSTFUN(1,CNEM2F,X,Q2,FLm,npts,0)
         CALL ZMSTFUN(2,CNEM2F,X,Q2,F2m,npts,0)
         CALL ZMSTFUN(3,CNEM3F,X,Q2,XF3m,npts,0) 

         CALL ZMSTFUN(1,CNEM2Fc,X,Q2,FLmC,npts,0)
         CALL ZMSTFUN(2,CNEM2Fc,X,Q2,F2mC,npts,0)
         CALL ZMSTFUN(3,CNEM3Fc,X,Q2,XF3mC,npts,0) 

         CALL ZMSTFUN(1,CNEM2Fb,X,Q2,FLmB,npts,0)
         CALL ZMSTFUN(2,CNEM2Fb,X,Q2,F2mB,npts,0)
         CALL ZMSTFUN(3,CNEM3Fb,X,Q2,XF3mB,npts,0) 
         
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
c            if(EWFIT==0) then
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
c ============================================================
            PZ = 1./Pz
C     EW couplings of u-type and d-type quarks at the scale Q2
               
cv for UNpolarised case should reduce to:
         A_u = e2u - ve*PZ*2.*euq*vu +(ve**2 + ae**2)*PZ**2*(vu**2+au**2)
         A_d = e2d - ve*PZ*2.*edq*vd +(ve**2 + ae**2)*PZ**2*(vd**2+ad**2)
         B_u = -ae*PZ*2.*euq*au + 2.*ve*ae*(PZ**2)*2.*vu*au
         B_d = -ae*PZ*2.*edq*ad + 2.*ve*ae*(PZ**2)*2.*vd*ad

            
            F2Gamma(i) = 4.D0/9.D0 * F2(i)  + 1.D0/9.D0 * F2m(i)
            FLGamma(i) = 4.D0/9.D0 * FL(i)  + 1.D0/9.D0 * FLm(i)
            XF3(i)  = B_U*XF3(i)  + B_D*XF3m(i)
            F2(i)   = A_U*F2(i)   + A_D*F2m(i)
            FL(i)   = A_U*FL(i)   + A_D*FLm(i)
            
            XF3c(i)  = B_U*XF3c(i)  + B_D*XF3mc(i)
            F2c(i)   = A_U*F2c(i)   + A_D*F2mc(i)
            FLc(i)   = A_U*FLc(i)   + A_D*FLmc(i)
            
            XF3b(i)  = B_U*XF3b(i)  + B_D*XF3mb(i)
            F2b(i)   = A_U*F2b(i)   + A_D*F2mb(i)
            FLb(i)   = A_U*FLb(i)   + A_D*FLmb(i)
            
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


c$$$c     =====================================================
c$$$      
c$$$C QCDNUM:
c$$$      double precision cnep2f(-6:6)
c$$$      double precision cnem2f(-6:6)
c$$$      double precision cnep3f(-6:6)
c$$$      double precision cnem3f(-6:6)
c$$$      DATA CNEP2F/0.,0.,1.,0.,1.,0.,0.,0.,1.,0.,1.,0.,0./
c$$$      DATA CNEM2F/0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0./
c$$$      DATA CNEP3F/0.,0.,-1.,0.,-1.,0.,0.,0.,1.,0.,1.,0.,0./
c$$$      DATA CNEM3F/0.,-1.,0.,-1.,0.,-1.,0.,1.,0.,1.,0.,1.,0./
c$$$
c$$$      double precision ccep2f(-6:6)
c$$$      double precision ccem2f(-6:6)
c$$$      double precision ccep3f(-6:6)
c$$$      double precision ccem3f(-6:6)
c$$$      DATA CCEP2F/0.,0.,1.,0.,1.,0.,0.,1.,0.,1.,0.,0.,0./
c$$$      DATA CCEM2F/0.,0.,0.,1.,0.,1.,0.,0.,1.,0.,1.,0.,0./
c$$$      DATA CCEP3F/ 0.,0.,-1.,0.,-1.,0.,0.,1.,0.,1.,0.,0.,0./
c$$$      DATA CCEM3F/0.,0. ,0.,-1.,0.,-1.,0.,0.,1.,0.,1.,0.,0./
c$$$
c$$$      integer NPmaxDIS
c$$$      parameter (NPmaxDIS =  2000)
