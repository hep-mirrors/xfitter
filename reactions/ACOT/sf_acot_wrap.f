cccccccccccccccccccccccccccccccccccccccccccc
      subroutine SF_ACOT_wrap(
     $     x_in,q2_in,index,
     $     f123l_out,f123lc_out,f123lb_out)
cccccccccccccccccccccccccccccccccccccccccccc
c     F1, F2, F3, FL are out via f123l_out 
c     f123l_out(1..4)=F1,F2,F3,FL
c     same for charm  only cotribution: f123lc
c     same for bottom only cotribution: f123lb
c     
c     hfscheme_in: for now only NLO massless and massive are possible
c     icharge_in: 0 NC: photon exchange only
c     icharge_in: 4 NC e+: gamma+gammaZ+Z
c     icharge_in: 5 NC e-: gamma+gammaZ+Z 
c     icharge_in:-1 CC e-
c     icharge_in:+1 CC e+
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

      integer index,i

      double precision x_in,q2_in,xmu,x,q
      integer icharge_in, mode_in, hfscheme_in

      double precision polar_in, polar
      integer icharge

C----------------------------------------------------------------------
      DOUBLE PRECISION XKMIN,XKMAX,DUM1,DUM2
      INTEGER NORD,ISCHin,KFAC,IQCD

      COMMON /ACOT_SET/ NORD,ISCHin,KFAC,IQCD,XKMIN,XKMAX,DUM1,DUM2      
C----------------------------------------------------------------------
c communication with Fred's code
      integer ISCHout, iset, iflg, ihad
      double precision hmass, xmc,xmb
      double precision sinw2, xmw, xmz

      COMMON /Ischeme/ ISCHout, ISET, IFLG, IHAD
      common /fred/ xmc,xmb,HMASS
      common/fredew/ sinw2, xmw, xmz

      logical ifirst,UseKFactors
      data ifirst /.true./
C----------------------------------------------------------------------
C     set "Isch, Iset, Iflg, Ihad" in common block first
      iset =1                   ! dummy
      iflg =0                   ! dummy
      ihad =1                   ! proton

C     Use ISCH instead of HFSCHEME_IN
      ISCHout=ISCH in

!     taken from couplings.inc
      xmc=mch                   
      xmb=mbt                   
      sinw2=sin2thw
      xmw=mw
      xmz=mz
      UseKFactors=.true.
      if(KFAC.eq.0)  UseKFactors=.false.
      
      q=dsqrt(q2_in)
      xmu=q !*** FIX MU=Q
      icharge=icharge_in
      icharge=4
      x=x_in
      polar=polar_in

! Target mass correction!
C     ************** FIO 31mar2022 ** THIS IS HARD-WIRED FOR NOW:
       hmass=0.938d0             !*** Hadron Mass for target mass corrections

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


      if (.not.UseKFactors)
     >    then
         if(ifirst) write(6,*) ' NOT using K-Factors  (may be slow)'
         ifirst=.false.
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

      include 'qcdnumhelper.inc'
      
      double precision F123Lxcb(3,4),F123Lxcb_LO(3,4),
     >  F123Lxcb_QCDNUM(3,4), Fsave(4) !*** 3='xcb', 4='123L'
      double precision x, q, xmu, polar
      integer index, icharge
      double precision maxFactor,small
c     data maxFactor,small /1.0d1, 1.e-12/
      data maxFactor,small /999d0, 1.e-12/
      
      COMMON /Ischeme/ ISCH, ISET, IFLG, IHAD

C Local table of k-factors
      integer NKfactMax,nTotal,i,j
      parameter (NKfactMax=NPmaxDIS,nTotal=3*4*NKfactMax)  !*** NPmaxDIS defined in qcdnumhelper.inc
      double precision akFACTxcb(3,4,NKFACTMAX) !*** 3='xcb', 4='123L'
      data  akFACTxcb/nTotal*0.d0/  !*** Zero K-Factor table

c     Create a vector to keep track of which index has k-factors already computed
      logical ifirst(NKfactMax),ihead,isave,iprint
      data ifirst, ihead, isave, iprint
     > /NKfactMax*.TRUE.,.TRUE.,.TRUE.,.FALSE./

C----------------------------------------------------------------------
      DOUBLE PRECISION XKMIN,XKMAX,DUM1,DUM2
      INTEGER NORD,ISCHin,KFAC,IQCDNUM

      COMMON /ACOT_SET/ NORD,ISCHin,KFAC,IQCDNUM,XKMIN,XKMAX,DUM1,DUM2      
C----------------------------------------------------------------------

      
c     ------Check index is within bounds
      if(index.gt.NKfactMax) then
         write(6,*) ' Error: index > NKfactMax = ',index,NKfactMax
         stop
      endif


C-----------------------------------------------------------------------------
c     FIRST TIME THROUGH: FILL K-FACTOR      
C-----------------------------------------------------------------------------
      if(ifirst(index))  then  !***  FIRST TIME THROUGH: FILL K-FACTOR  ===
         ifirst(index)=.FALSE.
         call Fgen123Lxcb(icharge, X, Q,xmu,F123Lxcb, polar)
         IschORIG=Isch
         Isch=5  !*** Massive LO Calculation
         call Fgen123Lxcb(icharge, X, Q,xmu,F123Lxcb_LO, polar)
         Isch=IschORIG          !*** Reset Ischeme
         
C==========================================================================
C---------------------------------!*** OPTION TO USE QCDNUM FOR K-FACTORS
         if(iqcdnum.eq.1) then
         if((icharge.eq.0).or.(icharge.eq.4).or.(icharge.eq.5)) THEN  !*** NEUTRAL CURRENT ONLY FOR NOW
         call Fgen123LxcbQCDNUM(index,icharge, X, Q,xmu,F123Lxcb_QCDNUM, polar)
c     !*** 3='xcb', 4='123L'
         F123Lxcb_LO(1,2)=F123Lxcb_QCDNUM(1,2)  !*** Use QCDNUM F2
         F123Lxcb_LO(1,4)=F123Lxcb_QCDNUM(1,4)  !*** Use QCDNUM FL
         ENDIF
         endif

C==========================================================================
C---------------------------------
C     Generate K-Factor
      
      do i=1,3
      do j=1,4
         if(abs(F123Lxcb_LO(i,j)).lt.small) then
            akFACTxcb(i,j,Index)=1.0d0 !**** Default if denom is zero or small
         else
            akFACTxcb(i,j,Index)=F123Lxcb(i,j)/F123Lxcb_LO(i,j)

            if( akFACTxcb(i,j,Index).gt.xkmax ) then
               if(iprint) then
                  write(6,*)'Warning: K-Fac is large:',
     >                 akFACTxcb(i,j,Index),F123Lxcb(i,j),i,j
                endif
                akFACTxcb(i,j,Index)=xkmax
            elseif( akFACTxcb(i,j,Index).lt.xkmin ) then
               if(iprint) then
                  write(6,*)'Warning: K-Fac is small: *** ',
     >                  akFACTxcb(i,j,Index),F123Lxcb(i,j),i,j
               endif
               akFACTxcb(i,j,Index)=xkmin
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
c     END OF  FILL K-FACTOR      
C-----------------------------------------------------------------------------
c===============================================================================           
C-----------------------------------------------------------------------------
c    NOT FIRST TIME THROUGH: USE K-FACTOR 
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
         else  !***  NOT FIRST TIME THROUGH: USE K-FACTOR ======================
           write(76,101)   !**** FORMAT FOR SCREEN
     $           index,x,q,icharge,polar,
     >        ((akFACTxcb(i,j,Index),j=1,4),i=1,3)

            IschORIG=Isch
            Isch=5  !*** Massive LO Calculation
            call Fgen123Lxcb(icharge, X, Q,xmu,F123Lxcb_LO, polar)
            Isch=IschORIG       !*** Reset Ischeme

            if(iqcdnum.eq.1) then
            if((icharge.eq.0).or.(icharge.eq.4).or.(icharge.eq.5)) THEN  !*** NEUTRAL CURRENT ONLY FOR NOW
               call Fgen123LxcbQCDNUM(index,icharge, X, Q,xmu,
     >                  F123Lxcb_QCDNUM, polar)
c     !*** 3='xcb', 4='123L'
               F123Lxcb_LO(1,2)=F123Lxcb_QCDNUM(1,2)  !*** Use QCDNUM F2
               F123Lxcb_LO(1,4)=F123Lxcb_QCDNUM(1,4)  !*** Use QCDNUM FL
            ENDIF
            endif
            
C     Use K-Factor
            do i=1,3
            do j=1,4
               F123Lxcb(i,j)=F123Lxcb_LO(i,j)* akFACTxcb(i,j,Index)
            enddo
            enddo

         endif   !***  END OF: NOT FIRST TIME THROUGH: USE K-FACTOR ======================
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
      return
      end

C
