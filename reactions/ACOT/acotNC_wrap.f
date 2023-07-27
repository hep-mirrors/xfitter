C-------------------------------------------------------------------------
C-------------------------------------------------------------------------
      subroutine acotNC_wrapa(x,q2,ipn,f2,
     $     f2c,f2b,fl,flc,flb,iflag,idatapt,f2QCDNUM,flQCDNUM
     $     ,usekfactors)
C-------------------------------------------------------------------------
C
C Created 29 May 2012. A wraper around sfun allowing for fast k-factors
C
C  Input parameters: x,Q2
C                   iflag,index        -- iflag=1: get k-factor,
C                                         index: data point global, unique index.
C                   f2QCDNUM,flQCDNUM  -- f2_gamma,fl_gamma from QCDNUM for this x,Q2
C                   UseKFactors        -- use or not kfactors
C------------------------------------------------------------------------
      implicit none
      double precision x,q2,q,f2,f2c,f2b,fl,flc,flb
      integer iflag, ipn, idatapt
      double precision f2QCDNUM,flQCDNUM
      logical usekfactors
c     26 July 2023: Fred: 
c     index keeps track of the sequential data point number
c     this is used to keep track of the k-factors to speed up the ACOT calculation
c     idatapt is only the # for the individual data set; not used for K-facor
      integer index
      common /acotIndex_TMP/  index
      data index /0/


C
C Local table of k-factors
C
      integer NKfactMax
      parameter (NKfactMax=10000)
      double precision akfactF2(NKFACTMAX),akfactFL(NKFACTMAX)
C-------------------------------------------------------------------------
      double precision polar,f123L(4),f123Lc(4),f123Lb(4)
      integer icharge,hfscheme
C-------------------------------------------------------------------------
c23456789012345678901234567890123456789012345678901234567890123456789012
      if(index.ge.NKfactMax) then
         write(6,*) "Error: more data than K-factor array space:",
     >        "index, NKfactMax",index,NKfactMax
      ENDIF
C-------------------------------------------------------------------------
c     ipn: 1=proton; 2=neutron
c     Fred: for UseKFactors=.true.
      UseKFactors=.true.  !*** Fred force

      
      q=dsqrt(q2)


C-------------------------------------------------------------------------
         call mstwnca(x,q,ipn,f2,f2c,f2b,fl,flc,flb)
C-------------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccc
         if(ipn.ne.1) stop      !*** FIO: only wired for proton at present  31mar2022

         hfscheme=111           !*** S-ACOT(Chi) [preferred]  111 -> 9 [need to trace later]
         icharge=4  !*** NC w/ Gamma+Z terms
         iflag=0  !*** not used
         polar=0  !*** only unpolarized 
c         index=0  !*** =0 skip k-factors
c         index=1
         UseKFactors=.false.
         UseKFactors=.true.  !*** Fred force
         
         index=index+1
         call SF_ACOT_wrap(
     $     x,q2,
     $     f123L,
     $     f123Lc,
     $     f123Lb,
     $     hfscheme,
     $     icharge   ,
     $     iflag, index, UseKFactors, polar   )
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


C--   G.W. 12/04/2012 Set these variables via COMMON/TRprimeCommon/.
C      var1=0d0; var2=0d0; var3=0d0; var4=0d0; ! standard TR' scheme
C      var1=0d0; var2=1d0; var3=-2d0/3d0; var4=1d0; ! optimal TR' scheme
C--   Other parameter values can be used to investigate uncertainties
C--   due to the choice of TR' GM-VFNS (see Table 1 of arXiv:1201.6180).
C      var1=0d0; var2=1d0; var3=-1d0; var4=0d0; ! GM-VFNS1
C      var1=0d0; var2=0.5d0; var3=-1d0; var4=0d0; ! GM-VFNS2
C      var1=0d0; var2=0d0; var3=0d0; var4=1d0; ! GM-VFNS3
C      var1=0d0; var2=1d0; var3=0.3d0; var4=0d0; ! GM-VFNS4
C      var1=0.1d0; var2=0d0; var3=0d0; var4=0d0; ! GM-VFNS5
C      var1=-0.2d0; var2=0d0; var3=0d0; var4=0d0; ! GM-VFNS6

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
