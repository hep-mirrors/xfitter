      subroutine mstwnc_wrap(x,q2,ipn,f2,
     $     f2c,f2b,fl,flc,flb,iflag,index,f2QCDNUM,flQCDNUM
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
      integer iflag, index, ipn
      double precision f2QCDNUM,flQCDNUM
      logical usekfactors


C
C Local table of k-factors
C
      integer NKfactMax
      parameter (NKfactMax=10000)
      double precision akfactF2(NKFACTMAX),akfactFL(NKFACTMAX)
C-------------------------------------------------------------------------
     
      q=dsqrt(q2)

      if (iflag.eq.1 .or. .not. UseKFactors) then

         call mstwnc(x,q,ipn,f2,f2c,f2b,fl,flc,flb)

C Store k-factors:
         if (UseKFactors) then
            if (index.gt.NKfactMax) then
               print *,'Error in sfun_wrap'
               print *,'Increase NKfactMax from '
     $              ,NKfactMax,' to at least ',Index
               stop
            endif
            akfactF2(index) = f2/F2QCDNUM
            akfactFL(index) = fl/FlQCDNUM         
         endif
      else
C Use k-factor:
         f2 = F2QCDNUM*akfactF2(index)
         fl = FLQCDNUM*akfactFL(index) 
      endif
      end


      Subroutine RT_SetAlphaS(alphaSzero)
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

      COMMON/mstwCommon/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax

C-----------------------------------------
      alphaSQ0 = alphaSZero
C-----------------------------------------
      end






c===================

      Subroutine RT_Set_Input(
     $     distancein,tolerancein,
     $     mCharmin,mBottomin,alphaSQ0in,alphaSMZin,
     $     alphaSorderin,alphaSnfmaxin,iordin)
C---------------------------------------------------------------------------
C  Wraper for INPUT common, set parameters
C---------------------------------------------------------------------------
      implicit none
C Input variables:
 
      INTEGER alphaSorderin,alphaSnfmaxin
      DOUBLE PRECISION distancein,tolerancein,
     &     mCharmin,mBottomin,alphaSQ0in,alphaSMZin

      INTEGER iordin
      
C Common variables:
      INTEGER alphaSorder,alphaSnfmax
      DOUBLE PRECISION distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ

      COMMON/mstwCommon/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax

      INTEGER iord
      COMMON/iordCommon/iord


C-------------------------------      

      alphaSorder  = alphaSorderin
      alphaSnfmax  = alphaSnfmaxin
      alphaSQ0     = alphaSQ0in
      alphaSMZ     = alphasMZin
      mCharm       = mCharmin
      mBottom      = mBottomin
      iord         = iordin
      tolerance    = tolerancein
      distance     = distancein

C--------------------------------------------------------------------------
      end
