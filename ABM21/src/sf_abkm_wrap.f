      subroutine sf_abkm_wrap(x,q2,f2abkm,flabkm,f3abkm,f2cabkm,
     $   flcabkm,f3cabkm,f2babkm,flbabkm,f3babkm,ncflag,charge,
     $   polar,sin2thw,cos2thw,MZ)
C-------------------------------------------------------------------------
C
C Created by RP 09 Jan 2012. A wraper around abkm functions
C Updated to OPENQCDRAD-2.0b4 version Nov 2016
C
C  Input parameters: x,Q2
C                   iflag,index        -- iflag=1: get k-factor, index: data point global, unique index.
C                   f2QCDNUM,flQCDNUM  -- f2_gamma,fl_gamma from QCDNUM for this x,Q2
C                   UseKFactors        -- use or not kfactors
C------------------------------------------------------------------------
      implicit none
      double precision x,q2,f2abkm,flabkm,f3abkm
      double precision f2cabkm,flcabkm,f3cabkm
      double precision f2babkm,flbabkm,f3babkm
      integer nb,nt,ni,ncflag
      double precision charge,polar

C-------------------------------------------------------------------------
c      b1=f2qcd(3,1,22,xb,q2)
c      f2qcd(nb,nt,ni,xb,q2)
      real*8 f2qcd,flqcd,f3qcd,f2charm_ffn,flcharm_ffn
      real*8 f2nucharm,f3nucharm,ftnucharm
      real*8 facgz,faczz
      real*8 facgzf3,faczzf3
      double precision eleVec,eleAxial,sin2thw,cos2thw,Mz,PZ


      nt = 1
c --------------- Neutral Currents !  ----------------    
      if(ncflag.eq.1) then

c new rewriten version of ABM, now with Z cont available need to calculate full SFs      
c we also take polarisation into account      

c------- Z-ELECTRON COUPLINGS:
       eleVec=   -0.5d0 +  2.0d0 * sin2thw
       eleAxial= -0.5d0

       if (charge.gt.0) then
          facgz = - eleVec-polar*eleAxial
          faczz = eleVec**2 + eleAxial**2 +2.0d0 * polar*eleAxial*eleVec
          facgzf3 = - eleAxial-polar*eleVec
          faczzf3 = 2*eleAxial*eleVec + polar*(eleVec**2+eleAxial**2)
       else
          facgz = - eleVec+polar*eleAxial
          faczz = eleVec**2 + eleAxial**2 -2.0d0 * polar*eleAxial*eleVec
          facgzf3 = - eleAxial+polar*eleVec
          faczzf3 = 2*eleAxial*eleVec - polar*(eleVec**2+eleAxial**2)
       endif
c      print*,'sf_abkm_wrap, charge,polar,sin2thw,facgz,faczz',charge,polar,
c     $   sin2thw,facgz,faczz

C      Propagator factor PZ
       PZ = 4.d0 * sin2thw * cos2thw * (1.+Mz**2/q2)
       PZ = 1./Pz

       f2abkm = f2qcd(3,nt,22,x,q2) + facgz*PZ*f2qcd(3,nt,25,x,q2) 
     $ + faczz*PZ*PZ*f2qcd(3,nt,23,x,q2)
       flabkm = flqcd(3,nt,22,x,q2) + facgz*PZ*flqcd(3,nt,25,x,q2) 
     $ + faczz*PZ*PZ*flqcd(3,nt,23,x,q2)
       f3abkm = facgzf3*PZ*f3qcd(3,nt,25,x,q2) + faczzf3*PZ*PZ
     $ *f3qcd(3,nt,23,x,q2)
c add the negative sign in front of Y_xF3 for neg charge      
       if(charge.gt.0) then
           f3abkm = -1*f3abkm
       endif
c      print*,'ABM f2abkm,gg,gz,zz ',f2abkm,f2qcd(3,1,22,x,q2),
c     $ f2qcd(3,1,23,x,q2),f2qcd(3,1,25,x,q2)

c c quark
       f2cabkm = f2charm_ffn(x,q2,8)
       flcabkm = flcharm_ffn(x,q2,8)
       f3cabkm = 0.0d0
c b quark
       f2babkm = f2charm_ffn(x,q2,10)
       flbabkm = flcharm_ffn(x,q2,10)
       f3babkm = 0.0d0


c --------------- Charged Currents !  ----------------    
      elseif(ncflag.eq.0) then

       ni = 24
       if(charge.gt.0) then
c W+       
        nb = 6
        else
c W-       
        nb = 7
       endif

c divide all SFs by 2 to get e+/-
       f2abkm = f2qcd(nb,nt,ni,x,q2)/2
       flabkm = flqcd(nb,nt,ni,x,q2)/2
       f3abkm = f3qcd(nb,nt,ni,x,q2)/2
 
c c quark
       f2cabkm = f2nucharm(nb,nt,ni,x,q2,8)/2
       flcabkm = f2nucharm(nb,nt,ni,x,q2,8)/2
     &     -ftnucharm(nb,nt,ni,x,q2,8)/2
       f3cabkm = f3nucharm(nb,nt,ni,x,q2,8)/2
c b quark
       f2babkm = 0.0d0
       flbabkm = 0.0d0
       f3babkm = 0.0d0
      endif

      RETURN
      end

      Subroutine ABKM_Set_Input(kschemepdfin,kordpdfin,rmass8in,
     $      rmass10in,msbarmin,hqscale1in,hqscale2in) 
C---------------------------------------------------------------------------
C  Wraper for INPUT common, set parameters
C---------------------------------------------------------------------------
      implicit none
C Input variables:
      double precision rmass8in,rmass10in
      integer kschemepdfin,kordpdfin
      logical msbarmin
      double precision hqscale1in,hqscale2in
      
C Common variables:
      integer npdftot,kordpdf,kschemepdf,kpdfset,kordkernel
c      common /forpdfset/ kschemepdf,kordpdf
      common /forpdfset/ npdftot,kordpdf,kschemepdf,kpdfset,kordkernel

      double precision rmass,rmassp,rcharge
      COMMON /MASSES/ rmass(150),rmassp(50),rcharge(150)

      double precision q20,q2rep,q2s,q20alphas,alphas0,alpsz,alpss
      double precision alpsc,alpsb,alpst,tscale,rscale,fscale,hqscale1
      double precision hqscale2
      integer nfeff,kordalps,kfeff,kordhq,kordf2,kordfl,kordf3
      logical alsmz 

      common /FORALPSRENORM/ q20,q2rep,q2s,q20alphas,alphas0,alpsz,alpss
     , ,alpsc,alpsb,alpst,tscale,rscale,fscale,hqscale1,hqscale2
     , ,nfeff,kordalps,kfeff
     , ,kordhq,kordf2,kordfl,kordf3
     , ,alsmz

      logical msbarm,hqnons,bmsnfopt,bmsnnlo,vloop
      double precision ddnnlohq
      common /forschemedef/ ddnnlohq,msbarm,hqnons
     ,  ,bmsnfopt,bmsnnlo,vloop
C-------------------------------      

      rmass(8)  = rmass8in
      rmass(10) = rmass10in
      kschemepdf= kschemepdfin
c set same order for pdf, light, heavy quarks      
      print*,'kordpdfin, rmass8in,rmass10in  ', kordpdfin, rmass8in,rmass10in
      kordpdf   = kordpdfin
      kordhq    = kordpdfin
      kordf2    = kordpdfin
c follow recommendation of Sergey Alekhin for FL    
      kordfl    = kordpdfin+1
      kordf3    = kordpdfin
      kordalps  = kordpdfin
c run in running m scheme
      msbarm    = msbarmin 
      hqscale1  = hqscale1in
      hqscale2  = hqscale2in
       
C--------------------------------------------------------------------------
      end

