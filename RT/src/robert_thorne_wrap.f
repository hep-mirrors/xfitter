      subroutine sfun_wrap(x,q2,f2p,flp,f1p,rp,f2n,fln,f1n,rn,
     $     f2c,flc,f1c,f2b,flb,f1b,iflag,index,f2QCDNUM,flQCDNUM
     $     ,usekfactors)
C-------------------------------------------------------------------------
C
C Created 19 June 2011. A wraper around sfun allowing for fast k-factors
C
C  Input parameters: x,Q2
C                   iflag,index        -- iflag=1: get k-factor, index: data point global, unique index.
C                   f2QCDNUM,flQCDNUM  -- f2_gamma,fl_gamma from QCDNUM for this x,Q2
C                   UseKFactors        -- use or not kfactors
C------------------------------------------------------------------------
      implicit none
      double precision x,q2,f2p,flp,f1p,rp,f2n,fln,f1n,rn,
     $     f2c,flc,f1c,f2b,flb,f1b
      integer iflag, index
      double precision f2QCDNUM,flQCDNUM
      logical usekfactors

      integer mode 
      data mode/1/
C
C Local table of k-factors
C
      integer NKfactMax
      parameter (NKfactMax=10000)
      double precision akfactF2(NKFACTMAX),akfactFL(NKFACTMAX)
C-------------------------------------------------------------------------

      if (iflag.eq.1 .or. .not. UseKFactors) then
         call sfun(x,q2,mode,f2p,flp,f1p,rp,f2n,fln,f1n,rn,
     $        f2c,flc,f1c,f2b,flb,f1b)
C Store k-factors:
         if (UseKFactors) then
            if (index.gt.NKfactMax) then
               print *,'Error in sfun_wrap'
               print *,'Increase NKfactMax from '
     $              ,NKfactMax,' to at least ',Index
               stop
            endif
            akfactF2(index) = f2p/F2QCDNUM
            akfactFL(index) = flp/FlQCDNUM         
         endif
      else
C Use k-factor:
         f2p = F2QCDNUM*akfactF2(index)
         flp = FLQCDNUM*akfactFL(index) 
      endif
      end


      Subroutine RT_SetAlphaS(alphaSzero)
C--------------------------------------------------------------
C  Wraper for INPUT common, to set alphaS
C--------------------------------------------------------------
      implicit none
C Parameter:
      double precision alphaSzero
c Common:
      double precision alphaS0
      COMMON/INPUT/alphaS0   
C-----------------------------------------
      alphaS0 = alphaSZero
C-----------------------------------------
      end

      Subroutine RT_Set_Input(alphaS0in,alambdain,flavorin,qsctin,
     $     qsdtin,iordin,inullin)
C---------------------------------------------------------------------------
C  Wraper for INPUT common, set parameters
C---------------------------------------------------------------------------
      implicit none
C Input variables:
      double precision alphaS0in,alambdain,flavorin,qsctin,qsdtin
      integer iordin,inullin
      
C Common variables:
      double precision alphaS0,alambda,flavor,qsct,qsdt
      integer iord,inull
      COMMON/INPUT/alphaS0,alambda,flavor,qsct,qsdt,iord,inull
C-------------------------------      

      alphaS0  = alphaS0in
      alambda  = alambdain
      flavor   = flavorin
      qsct     = qsctin
      qsdt     = qsdtin
      iord     = iordin
      inull    = inullin
C--------------------------------------------------------------------------
      end
