cccccccccccccccccccccccccccccccccccccccccccc
      subroutine SF_ACOT_wrap(
     $     x_in,q2_in,
     $     f123l_out,
     $     f123lc_out,
     $     f123lb_out,
     $     hfscheme_in,
     $     icharge_in,
     $     iflag, index, UseKFactors)


cccccccccccccccccccccccccccccccccccccccccccc
      implicit none

!--------      
      include 'steering.inc'
      include 'couplings.inc'
!
      double precision f123l_out(4)
      double precision f123lc_out(4)
      double precision f123lb_out(4)

      double precision f123l_NLO(4)
      double precision f123lc_NLO(4)
      double precision f123lb_NLO(4)

      double precision f123l_LO(4)
      double precision f123lc_LO(4)
      double precision f123lb_LO(4)

      logical UseKFactors
      integer iflag, index

      double precision x_in,q2_in,xmu,x,q
      integer icharge_in, mode_in, hfscheme_in
      
      integer icharge
      integer isch, iset, iflg, ihad
      double precision hmass, xmc,xmb
      double precision sinw2, xmw, xmz
      COMMON /Ischeme/ ISCH, ISET, IFLG, IHAD

      common /fred/ xmc,xmb,HMASS
      common/fredew/ sinw2, xmw, xmz



C
C Local table of k-factors
C
      integer NKfactMax,i
      parameter (NKfactMax=10000)
      double precision akfact(4,NKFACTMAX)


C----------------------------------------------------------------------
C     set "Isch, Iset, Iflg, Ihad" in common block first
      iset =1
      iflg =0
      ihad =1

c      icharge_in=0  !*** photon
c      icharge_in=4  !*** photon + Z boson

      xmc=mch
      xmb=mbt
      sinw2=sin2thw
      xmw=mw
      xmz=mz

      xmu=dsqrt(Q2_in)      !*** mu=Q
      x=x_in
      q=dsqrt(q2_in)
      icharge=icharge_in


      hmass=0.d0
!      hmass=0.938d0             !*** Hadron Mass for target mass corrections

C----------------------------------------------------------------------
C      isch= 0  !*** NLO Massless MS-Bar'                
C      isch= 1  !*** Full ACOT Scheme '                  
C      isch= 2  !*** FFS '                               
C      isch= 3  !*** Simplified ACOT Scheme '            
C      isch= 4  !*** Test Full ACOT Scheme (no NLO Q)'   
C      isch= 5  !*** LO  '                               
C      isch= 6  !*** Massless LO '                       
C      isch= 7  !*** Short-cut2: ACOT w/ Massless NLO-Q '
C      isch= 8  !*** S-ACOT(Chi) [not preferred]'        
C      isch= 9  !*** S-ACOT(Chi) [preferred] '           
C----------------------------------------------------------------------


      if (HFSCHEME_IN.eq.1)   ISCH=6
      if (HFSCHEME_IN.eq.11)  ISCH=5



      if (iflag.eq.1.or..not.UseKFactors) then

! assumes user choice NLO massless or massive
         Call Fgen123L(icharge,   1, X, Q,xmu,F123L_LO) !*** total F
         Call Fgen123L(icharge,   2, X, Q,xmu,F123Lc_LO) !*** F-charm
         Call Fgen123L(icharge,   3, X, Q,xmu,F123Lb_LO) !*** F-bottom
C Store k-factors:
         if (UseKFactors) then

            if (index.gt.NKfactMax) then
               print *,'Error in sf_acot_wrap'
               print *,'Increase NKfactMax from '
     $              ,NKfactMax,' to at least ',Index
               stop
            endif

! make sure the kfactors correspond for massless or massive choice
            if (HFSCHEME.eq.1) then
               isch=0           ! massless NLO
            elseif (HFSCHEME.eq.11) then
               isch=1           ! massive NLO
            endif

            Call Fgen123L(icharge,   1, X, Q,xmu,F123L_NLO) !*** total F
            Call Fgen123L(icharge,   2, X, Q,xmu,F123Lc_NLO) !*** F-charm
            Call Fgen123L(icharge,   3, X, Q,xmu,F123Lb_NLO) !*** F-bottom
            
            do i=1,4 
               akfact(i,index) = F123L_NLO(i)/F123L_LO(i)
            enddo
            print*,'building kfactors for each data point:', index, 
     $           (akfact(i,index),i=1,4)
         endif ! Usekfactors

      else                      ! if iflag.ne.1 and usefactor.eq.false
c     Use k-factor:

         if (HFSCHEME.eq.1) then
            isch=6              ! massless 
         elseif (HFSCHEME.eq.11) then
            isch=5              ! massive 
         endif
         
         Call Fgen123L(icharge,   1, X, Q,xmu,F123L_LO) !*** total F
         Call Fgen123L(icharge,   2, X, Q,xmu,F123Lc_LO) !*** F-charm
         Call Fgen123L(icharge,   3, X, Q,xmu,F123Lb_LO) !*** F-bottom
         do i=1,4
            F123L_out(i)=F123L_LO(i)*akfact(i,index)
         enddo
         
      endif


      RETURN
      END






