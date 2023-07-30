cccccccccccccccccccccccccccccccccccccccccccc
      subroutine SF_ACOT_wrap(
     $     x_in,q2_in,
     $     f123l_out,
     $     f123lc_out,
     $     f123lb_out,
     $     hfscheme_in,
     $     icharge_in,
     $     iflag, index, UseKFactors, polar_in)
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
c     icharge_in: 4 NC e+: gamma+gammaZ+Z
c     icharge_in: 5 NC e-: gamma+gammaZ+Z 
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
      

      double precision F123Lxcb(3,4) !*** 3='xcb', 4='123L'

      double precision f123l_out(4)
      double precision f123lc_out(4)
      double precision f123lb_out(4)

      logical UseKFactors
      integer iflag, index,i

      double precision x_in,q2_in,xmu,x,q
      integer icharge_in, mode_in, hfscheme_in

      double precision polar_in, polar
      integer icharge

c communication with Fred's code
      integer isch, iset, iflg, ihad
      double precision hmass, xmc,xmb
      double precision sinw2, xmw, xmz

      COMMON /Ischeme/ ISCH, ISET, IFLG, IHAD
      common /fred/ xmc,xmb,HMASS
      common/fredew/ sinw2, xmw, xmz


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

        if(MASSH.eq.1) then
       xmu=sqrt(hqscale1in*Q2_in+hqscale2in*4*mch*mch)      !*** mu=Q
       elseif (MASSH.eq.2) then
       xmu=dsqrt(hqscale1in*Q2_in+hqscale2in*4*mbt*mbt) 
       endif         
      q=dsqrt(q2_in)
      icharge=icharge_in
      x=x_in
      polar=polar_in

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

         if (HFSCHEME_IN.eq.1) then 
            ISCH=0                  !*** NLO Massless MS-Bar
         elseif (HFSCHEME_IN.eq.11) then 
            ISCH=1                  !*** Full ACOT Scheme
         elseif (HFSCHEME_IN.eq.111) then 
            ISCH=9                  !*** S-ACOT(Chi) [preferred]
         endif



C ---------------------------------------
C     F2, FL, XF3 are already computed by QCDNUM:  (FIO 15 Dec 2012)
c     Pass inside sf_acot_wrap for K-Factor method:
C     Important: This relies on the call to UseZmvnsScheme
         F123Lxcb(1,1)=f123l_out(1)
         F123Lxcb(1,2)=f123l_out(2)
         F123Lxcb(1,3)=f123l_out(3)
         F123Lxcb(1,4)=f123l_out(4)
C ---------------------------------------
         

C----------------------------------------------------------------------
C----------------------------------------------------------------------
c     REWRITE: ENCAPSULATE THE K-FACTORS:  FIO 13 April 2012
C----------------------------------------------------------------------

c     This function returns F-123L, for total-c,b 'xcb', using 'K' factors
c     eventually this will be written in a single pass to speed up calculations 
c     K-factors are stored using data 'index' 
c     BUT no error checking is done on {x,Q}, user must ensure index matches {x,Q}
c     Note: scheme is set inside; this should not change from point to point
c     If index is too big (.gt.NKfactMax) or 0, use full calculation


      if ((.not.UseKFactors).or.(index.eq.0))
     >    then
         write(6,*) ' NOT using K-Factors  (may be slow)'
         Call Fgen123Lxcb(icharge, X, Q,xmu,F123Lxcb, polar) !*** total F no-K-factors, 

      else
         Call Fgen123LxcbK(index,icharge, X, Q,xmu,F123Lxcb, polar) !*** total F
      endif

C----------------------------------------------------------------------
c     UNPACK F123Lxcb ARRAY: 
C         ... IN THE FUTURE WE CAN EDIT THE "subroutine UseAcotScheme" to skip this step
C----------------------------------------------------------------------

      DO I=1,4         
         f123L_out( i)=F123Lxcb(1,i)
         f123Lc_out(i)=F123Lxcb(2,i)
         f123Lb_out(i)=F123Lxcb(3,i)
      enddo
      
      RETURN
      END

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------
c     This function returns F-123L, for total-c,b 'xcb', using 'K' factors
C----------------------------------------------------------------------
      subroutine Fgen123LxcbK(index,icharge, X, Q,xmu,F123Lxcb, polar)

cv      include 'ntot.inc'

      include 'qcdnumhelper.inc'


      double precision F123Lxcb(3,4),F123Lxcb_LO(3,4),F123Lxcb_QCDNUM(3,4), Fsave(4) !*** 3='xcb', 4='123L'
      double precision x, q, xmu, polar
      integer index, icharge
      double precision maxFactor,small
c     data maxFactor,small /1.0d1, 1.e-12/
      data maxFactor,small /999d0, 1.e-12/
      
      COMMON /Ischeme/ ISCH, ISET, IFLG, IHAD
cSG      common /fred/ xmc,xmb,HMASS
cSG      common/fredew/ sinw2, xmw, xmz

C Local table of k-factors
      integer NKfactMax,nTotal,i,j
      parameter (NKfactMax=NPmaxDIS,nTotal=3*4*NKfactMax)
      double precision akFACTxcb(3,4,NKFACTMAX) !*** 3='xcb', 4='123L'
      data  akFACTxcb/nTotal*0.d0/  !*** Zero K-Factor table

c     Create a vector to keep track of which index has k-factors already computed
      logical ifirst(NKfactMax),ihead,isave
      data ifirst, ihead, isave /NKfactMax*.TRUE.,.TRUE.,.TRUE./

c     ------Check index is within bounds
      if(index.gt.NKfactMax) then
         write(6,*) ' Error: index > NKfactMax = ',index,NKfactMax
         stop
      endif


C ---------------------------------------
C     F2, FL, XF3 are already computed by QCDNUM:  (FIO 15 Dec 2012)
c     Save results for K-Factor method:
C     Important: This relies on the call to UseZmvnsScheme
         Fsave(1)=F123Lxcb(1,1)
         Fsave(2)=F123Lxcb(1,2)
         Fsave(3)=F123Lxcb(1,3)
         Fsave(4)=F123Lxcb(1,4)

C-----------------------------------------------------------------------------
c     FIRST TIME THROUGH: FILL K-FACTOR      
C-----------------------------------------------------------------------------
      if(ifirst(index))  then  !***  FIRST TIME THROUGH: FILL K-FACTOR  ===
         ifirst(index)=.FALSE.
         call Fgen123Lxcb(icharge, X, Q,xmu,F123Lxcb, polar)
         IschORIG=Isch
         Isch=5  !*** Massive LO Calculation
         call Fgen123Lxcb(icharge, X, Q,xmu,F123Lxcb_LO, polar)
C---------------------------------!*** OPTION TO USE QCDNUM FOR K-FACTORS
         if(isave) then
            isave=.false.
            iqcdnum=0
            iqcdnum=1
cv            write(6,*) ' enter 1 for QCDNUM K-facs; 0 for LO-Massive ACOT ',iqcdnum
cv            read (5,*) iqcdnum
            write(6,*) ' set: 1 for QCDNUM K-facs; 0 for LO-Massive ACOT  = ',iqcdnum
         endif

         if(iqcdnum.eq.1) then
         if((icharge.eq.0).or.(icharge.eq.4).or.(icharge.eq.5)) THEN !*** NEUTRAL CURRENT ONLY FOR NOW
c         call Fgen123LxcbQCDNUM(index,icharge, X, Q,xmu,F123Lxcb_QCDNUM, polar)
c     !*** 3='xcb', 4='123L'
c         F123Lxcb_LO(1,2)=F123Lxcb_QCDNUM(1,2)  !*** Use QCDNUM F2
c         F123Lxcb_LO(1,4)=F123Lxcb_QCDNUM(1,4)  !*** Use QCDNUM FL
         F123Lxcb_LO(1,2)=Fsave(2)  !*** Use QCDNUM F2
         F123Lxcb_LO(1,4)=Fsave(4)  !*** Use QCDNUM FL
         ENDIF
         endif
C---------------------------------
         Isch=IschORIG  !*** Reset Ischeme
C     Generate K-Factor
         do i=1,3
         do j=1,4
            if(abs(F123Lxcb_LO(i,j)).lt.small) then
               akFACTxcb(i,j,Index)=1.0d0 !**** Default if denom is zero or small
            else
               akFACTxcb(i,j,Index)=F123Lxcb(i,j)/F123Lxcb_LO(i,j)
c
               if(Abs(akFACTxcb(i,j,Index)).gt.maxFactor) then  !*** Sanity Check:
c               write(6,*) ' Warning: K-Fac is large: '
               akFACTxcb(i,j,Index)=sign(maxFactor,akFACTxcb(i,j,Index))
               endif
            endif
         enddo
         enddo

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
C     
C-----------------------------------------------------------------------------
         if(ihead) open(62,file='output/KfactorsACOT.txt')
         if(ihead) write(6,*)   ' K-factors '
         if(ihead) write(6,*) 
     >   ' indx  x    Q    ichg pol   ',
     >   '  F1      F2      F3      FL    ',
     <   '  F1c     F2c     F3c     FLc   ',
     >   '  F1b     F2b     F3b     FLb   '
         ihead=.FALSE.
c     only dump the k-factor info the first time through the loop
            write(62,*)  !**** DUMP COMPLETE INFO 
     $           index,x,q,icharge,polar,
     >        ((F123Lxcb(   i,j),j=1,4),i=1,3),
     >        ((F123Lxcb_LO(i,j),j=1,4),i=1,3),
     >        ((akFACTxcb(i,j,Index),j=1,4),i=1,3)

           write(6,101)   !**** FORMAT FOR SCREEN
     $           index,x,q,icharge,polar,
     >        ((akFACTxcb(i,j,Index),j=1,4),i=1,3)
101        Format(i4,1x,f0.5,1x,f5.1,1x,i3,1x,f5.3,1x,12(f7.2,1x))
c
C-----------------------------------------------------------------------------
         else  !***  NOT FIRST TIME THROUGH: USE K-FACTOR ======================
            IschORIG=Isch
            Isch=5  !*** Massive LO Calculation
            call Fgen123Lxcb(icharge, X, Q,xmu,F123Lxcb_LO, polar)

            Isch=IschORIG  !*** Reset Ischeme
C     Use K-Factor
            do i=1,3
            do j=1,4
               F123Lxcb(i,j)=F123Lxcb_LO(i,j)* akFACTxcb(i,j,Index)
            enddo
            enddo
         endif


c     Adjust F2 and FL if we use QCDNUM K-facs
c     Note: this code is call both first time to re-set F123Lxcb, and 
c     during later calls to overwrite with QCDNUM + K-factor result
         if(iqcdnum.eq.1) then
           if((icharge.eq.0).or.(icharge.eq.4).or.(icharge.eq.5)) THEN !*** NEUTRAL CURRENT ONLY FOR NOW
c               call Fgen123LxcbQCDNUM(index,icharge, X, Q,xmu,F123Lxcb_QCDNUM, polar)
c               F123Lxcb(1,2)=F123Lxcb_QCDNUM(1,2)* akFACTxcb(1,2,Index)
c               F123Lxcb(1,4)=F123Lxcb_QCDNUM(1,4)* akFACTxcb(1,4,Index)

               F123Lxcb(1,2)=Fsave(2)* akFACTxcb(1,2,Index) !*** Use Stored QCDNUM F2
               F123Lxcb(1,4)=Fsave(4)* akFACTxcb(1,4,Index) !*** Use Stored QCDNUM FL


           endif
         ENDIF
  

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
      return
      end

C
