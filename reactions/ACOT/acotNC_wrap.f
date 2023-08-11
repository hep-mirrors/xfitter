C-------------------------------------------------------------------------
C-------------------------------------------------------------------------
      subroutine acotNC_wrapa(x,q2,ipn,f2,f2c,f2b,fl,flc,flb)
C-------------------------------------------------------------------------
C     Input parameters: x,Q2
C     ipn=1 for proton; nothing else is implemented
C------------------------------------------------------------------------
      implicit none
      double precision x,q2,f2,f2c,f2b,fl,flc,flb
      integer ipn
c----------------------------------------------------------------
c     26 July 2023: Fred: 
c     index keeps track of the sequential data point number
c     this is used to keep track of the k-factors to speed up the ACOT calculation
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
         index=index+1
         call SF_ACOT_wrap(x,q2,index,f123L,f123Lc,f123Lb)

         f2=f123L(2)
         fL=f123L(4)
         f2c=f123Lc(2)
         fLc=f123Lc(4)
         f2b=f123Lb(2)
         fLb=f123Lb(4)

C-------------------------------------------------------------------------
c     THIS IS JUST FOR CROSS-CHECK/DEBUGGING/SANITY CHECK  WITH ACOT
c      double precision f2mstw,f2cmstw,f2bmstw,flmstw,flcmstw,flbmstw
c      call mstwnca(x,q,ipn,f2mstw,f2cmstw,f2bmstw,flmstw,
c     >         flcmstw,flbmstw)
C-------------------------------------------------------------------------
C-------------------------------------------------------------------------
      end
C-------------------------------------------------------------------------
C-------------------------------------------------------------------------
c===================
      Subroutine ACOT_Set_Input(varinacot,intvarin)
C---------------------------------------------------------------------------
C  Wraper for  common, set parameters
C---------------------------------------------------------------------------
      implicit none
      DOUBLE PRECISION XKMIN,XKMAX,DUM1,DUM2
      INTEGER intvarin(4)
      DOUBLE PRECISION varinacot(4)
      INTEGER NORD,ISCH,KFAC,IQCD
      COMMON /ACOT_SET/ NORD,ISCH,KFAC,IQCD,XKMIN,XKMAX,DUM1,DUM2
C---------------------------------------------------------------------------
C---------------------------------------------------------------------------
      xkmin=varinacot(1)   ! K-factor min
      xkmax=varinacot(2)   ! K-factor max
      dum1 =varinacot(3)   ! unused
      dum2 =varinacot(4)   ! unused

      nord=intvarin(1)     !  NORD for calcululation: = [1,2,3]                                  
      isch=intvarin(2)     !  iScheme: 9= preferred ACOT-CHI                                     
      kfac=intvarin(3)     !  Use K-Factor 1=yes, 0=no                                           
      iqcd=intvarin(4)     !  Use QCDNUM for K-factor denominator: 1=yes, 0=no

C--------------------------------------------------------------------------
      end
