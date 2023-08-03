C-------------------------------------------------------------------------
C-------------------------------------------------------------------------
      subroutine acotNC_wrapa(x,q2,ipn,f2,f2c,f2b,fl,flc,flb)
C-------------------------------------------------------------------------
C
C Created 29 May 2012. A wraper around sfun allowing for fast k-factors
C
C     Input parameters: x,Q2
C     ipn=1 for proton; nothing else is implemented
C------------------------------------------------------------------------
      implicit none
      double precision x,q2,q,f2,f2c,f2b,fl,flc,flb,
     > f2mstw,f2cmstw,f2bmstw,flmstw,flcmstw,flbmstw
      integer ipn
c----------------------------------------------------------------
c     26 July 2023: Fred: 
c     index keeps track of the sequential data point number
c     this is used to keep track of the k-factors to speed up the ACOT calculation
c     idatapt is only the # for the individual data set; not used for K-facor
      integer index
      data index /0/
      common /acotIndex/  index
c----------------------------------------------------------------
      double precision f123L(4),f123Lc(4),f123Lb(4)
c----------------------------------------------------------------
c     ipn: 1=proton; 2=neutron
         if(ipn.ne.1) stop      !*** FIO: only wired for proton at present  31mar2022
C-------------------------------------------------------------------------
c23456789012345678901234567890123456789012345678901234567890123456789012
C-------------------------------------------------------------------------


c$$$      x=0.1d0
c$$$      q=10.0d0
c$$$      q2=q*q
      
      q=dsqrt(q2)

      
C-------------------------------------------------------------------------
c     THIS IS JUST FOR CROSS-CHECK/DEBUGGING/SANITY CHECK  WITH ACOT
      call mstwnca(x,q,ipn,f2mstw,f2cmstw,f2bmstw,flmstw,
     >         flcmstw,flbmstw)
C-------------------------------------------------------------------------
         index=index+1
         call SF_ACOT_wrap(x,q2,index,f123L,f123Lc,f123Lb)

         f2=f123L(2)
         fL=f123L(4)
         f2c=f123Lc(2)
         fLc=f123Lc(4)
         f2b=f123Lb(2)
         fLb=f123Lb(4)

 811     format(A,1x,i4,1x,f11.6,1x,f6.1,A3,8(1x,f6.2))
         write(6,811) "ACOT/MSTW: x,q,F2,FL,...C..B",
     >   index,x,q," | ",
     >   f2/f2mstw,fL/flmstw,f2c/f2cmstw,fLc/flcmstw,
     >   f2b/f2bmstw,fLb/flbmstw

         write(86,*) index,x,q,
     >   f2,f2mstw,fL,flmstw,f2c,f2cmstw,fLc,flcmstw,
     >   f2b,f2bmstw,fLb,flbmstw



         
cccccccccccccccccccccccccccccccccccccccccccc         
C-------------------------------------------------------------------------
      end
C-------------------------------------------------------------------------
C-------------------------------------------------------------------------


      Subroutine ACOT_SetAlphaS(alphaSzero)
C--------------------------------------------------------------
C  Wraper for INPUT common, to set alphaS
C--------------------------------------------------------------
      implicit none
C Parameter:


      double precision alphaSzero
C Common variables:

      INTEGER alphaSorder,alphaSnfmax
      DOUBLE PRECISION distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ


      COMMON/mstwCommona/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax

C-----------------------------------------
      alphaSQ0 = alphaSZero
C-----------------------------------------
      end






c===================

      Subroutine ACOT_Set_Input(varin,
!     $     distancein,tolerancein,
     $     mCharmin,mBottomin,alphaSQ0in,alphaSMZin,
     $     alphaSorderin,alphaSnfmaxin,iordin,intvarin)
C---------------------------------------------------------------------------
C  Wraper for INPUT common, set parameters
C---------------------------------------------------------------------------
      implicit none
C Input variables:
 
      INTEGER alphaSorderin,alphaSnfmaxin,intvarin(4)
      DOUBLE PRECISION mCharmin,mBottomin,alphaSQ0in,alphaSMZin
      DOUBLE PRECISION var1, var2, var3, var4
      DOUBLE PRECISION varin(4)

      INTEGER iordin
      
C Common variables:
      INTEGER alphaSorder,alphaSnfmax
      DOUBLE PRECISION distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ

      COMMON/mstwCommona/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax

      COMMON/TRprimeCommona/var1,var2,var3,var4 ! G.W. 11/04/2012

      INTEGER iord
      COMMON/iordCommona/iord
!$OMP THREADPRIVATE(/TRprimeCommon/)


C-------------------------------      

      alphaSorder  = alphaSorderin
      alphaSnfmax  = alphaSnfmaxin
      alphaSQ0     = alphaSQ0in
      alphaSMZ     = alphasMZin
      mCharm       = mCharmin
      mBottom      = mBottomin
      iord         = iordin
      tolerance    = 0d0 ! dummy
      distance     = 0d0 ! dummy
      var1         = varin(1)
      var2         = varin(2)
      var3         = varin(3)
      var4         = varin(4)

C--------------------------------------------------------------------------
      end
