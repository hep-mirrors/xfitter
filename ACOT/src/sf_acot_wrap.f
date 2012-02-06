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
c     F1, F2, F3, FL are out via f123l_out 
c     f123l_out(1)=F1
c     f123l_out(2)=F2
c     f123l_out(3)=F3
c     f123l_out(4)=FL
c     same for charm only cotribution: f123lc
c     same for charm only cotribution: f123lb
c     
c     hfscheme_in: for now only NLO massless and massive are possible
c     icharge_in: 0 NC: photon exchange only
c     icharge_in: 4 NC: gamma+gammaZ+Z 
c     icharge_in:-1 CC e-
c     icharge_in:+1 CC e+
c     iflag: flag from FCN (main minimisation routine)
c     index: data index - integer
c     UseKFactors: use of kfactors (NLO/LO) - Logical 
c     (since ACOT takes long time, this is for now set always to TRUE)
cccccccccccccccccccccccccccccccccccccccccccc
      implicit none

!--------      
      include 'steering.inc'
      include 'couplings.inc'
!--------

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

c communication with Fred's code
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
      double precision akfactc(4,NKFACTMAX)
      double precision akfactb(4,NKFACTMAX)



C----------------------------------------------------------------------
C     set "Isch, Iset, Iflg, Ihad" in common block first
      iset =1                   ! dummy
      iflg =0                   ! dummy
      ihad =1                   ! proton


!     taken from couplings.inc
      xmc=mch                   
      xmb=mbt                   
      sinw2=sin2thw
      xmw=mw
      xmz=mz

      xmu=dsqrt(Q2_in)      !*** mu=Q
      x=x_in
      q=dsqrt(q2_in)
      icharge=icharge_in

! Target mass correction!
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

C     get the LO SFs
      ISCH=5    
      
      Call Fgen123L(icharge,   1, X, Q,xmu,F123L_LO) !*** total F
      Call Fgen123L(icharge,   2, X, Q,xmu,F123Lc_LO) !*** F-charm
      Call Fgen123L(icharge,   3, X, Q,xmu,F123Lb_LO) !*** F-bottom

C         

      if (iflag.eq.1.or..not.UseKFactors) then
         
         if (HFSCHEME_IN.eq.1) then 
            ISCH=0
         elseif (HFSCHEME_IN.eq.11) then 
            ISCH=1
         elseif (HFSCHEME_IN.eq.111) then 
            ISCH=9
         endif
         
         Call Fgen123L(icharge,   1, X, Q,xmu,F123L_NLO) !*** total F
         Call Fgen123L(icharge,   2, X, Q,xmu,F123Lc_NLO) !*** total Fc
         Call Fgen123L(icharge,   3, X, Q,xmu,F123Lb_NLO) !*** total Fb
         do i=1,4
            F123L_out(i)=F123L_NLO(i)
            F123Lc_out(i)=F123Lc_NLO(i)
            F123Lb_out(i)=F123Lb_NLO(i)
         enddo
C     Store k-factors:
         if (UseKFactors) then
            
            if (index.gt.NKfactMax) then
               print *,'Error in sf_acot_wrap'
               print *,'Increase NKfactMax from '
     $              ,NKfactMax,' to at least ',Index
               call HF_errlog(104,'F: Error in sf_acot_wrap')
C     stop
            endif
            
            
            do i=1,4 
               
!     protect for the kfactors that return infinity 
               if (F123L_LO(i).eq.0.d0) then
                  akfact(i,index) = 1.d0
               else
                  akfact(i,index) = F123L_NLO(i)/F123L_LO(i)
               endif
               
               if (F123Lc_LO(i).eq.0.d0) then
                  akfactc(i,index) = 1.d0
               else
                  akfactc(i,index) = F123Lc_NLO(i)/F123Lc_LO(i)
               endif
               
               if (F123Lb_LO(i).eq.0.d0) then
                  akfactb(i,index) = 1.d0
               else
                  akfactb(i,index) = F123Lb_NLO(i)/F123Lb_LO(i)
               endif
               
            enddo
            if (index.eq.1) then
               print*,'building kfactors: F1,F2,F3 FL,F2b,F2c,Flb,Flc'
            endif
            print*,
     $           index, 
     $           (akfact(i,index),i=1,4),
     $           akfactb(2,index), akfactc(2,index),
     $           akfactb(4,index), akfactc(4,index)
               
         endif                  ! Usekfactors
         
      else                      ! if iflag.ne.1 and usefactor.eq.false
c     Use k-factor:
         
         
         do i=1,4
            
            F123L_out(i)=F123L_LO(i)*akfact(i,index)
            F123Lc_out(i)=F123Lc_LO(i)*akfactc(i,index)
            F123Lb_out(i)=F123Lb_LO(i)*akfactb(i,index)
            
            
         enddo
      endif
      
      
      
      RETURN
      END






