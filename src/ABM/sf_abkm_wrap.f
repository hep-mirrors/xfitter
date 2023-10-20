C
C OZ 10.10.17 updated to openqcdrad-2.1
C New user routines for PDFs and alpha_S have to be provided in this version

      FUNCTION useralphas(q2,kschemepdf,kordpdf,kpdfset)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'APSCOM6.'

      !useralphas = HF_Get_alphas(q2)
      q = dsqrt(q2)
      useralphas = alphas_wrapper(q)
      !print*,'useralphas = ',useralphas

      RETURN
      END


      FUNCTION userpdfs(xb,q2,i,kschemepdf,kordpdf,kpdfset)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer i
      dimension pdfsff(-6:6)

      !CALL HF_Get_PDFs(xb,q2,PDFSFF)
      q = dsqrt(q2)
      CALL pdf_xfxq_wrapper(xb,q,PDFSFF)
      userpdfs = PDFSFF(i)

      RETURN
      END

      DOUBLE PRECISION FUNCTION DiLog(X)
      double precision x
      double precision ddilog
      dilog = ddilog(x)
      return
      end


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

      subroutine sf_abkm_wrap_order(x,q2,f2abkm,flabkm,f3abkm,f2cabkm,
     $   flcabkm,f3cabkm,f2babkm,flbabkm,f3babkm,ncflag,charge,
     $   polar,sin2thw,cos2thw,MZ,kordpdfin)
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

      ! 18.10.2023 SZ change order
      integer kordpdfin
      double precision q20,q2rep,q2s,q20alphas,alphas0,alpsz,alpss
      double precision alpsc,alpsb,alpst,tscale,rscale,fscale,hqscale1
      double precision hqscale2
      integer nfeff,kordalps,kfeff,kordhq,kordf2,kordfl,kordf3
      logical alsmz
      !integer npdftot_saved,kordpdf_saved,kschemepdf_saved,kpdfset_saved,kordkernel_saved
      common /FORALPSRENORM/ q20,q2rep,q2s,q20alphas,alphas0,alpsz,alpss
     , ,alpsc,alpsb,alpst,tscale,rscale,fscale,hqscale1,hqscale2
     , ,nfeff,kordalps,kfeff
     , ,kordhq,kordf2,kordfl,kordf3
     , ,alsmz
      integer kordf2_saved,kordfl_saved,kordf3_saved,kordalps_saved

      ! 18.10.2023 SZ change order
      kordf2_saved = kordf2
      kordfl_saved = kordfl
      kordf3_saved = kordf3
      kordalps_saved = kordalps
      !print*,'kordpdfin  ', kordpdfin
      !kordpdf   = kordpdfin
      !kordhq    = kordpdfin
      kordf2    = kordpdfin
c follow recommendation of Sergey Alekhin for FL
      kordfl    = kordpdfin+1
      kordf3    = kordpdfin
      kordalps  = kordpdfin

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

      ! 18.10.2023 SZ change order
      kordf2 = kordf2_saved
      kordfl = kordfl_saved
      kordf3 = kordf3_saved
      kordalps = kordalps_saved

      RETURN
      end

      Subroutine ABKM_Set_Input(kschemepdfin,kordpdfin,rmass8in,
     $      rmass10in,msbarmin,hqscale1in,hqscale2in,flagthinterface)
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

C OZ 17.10.17 Flags to check that this routine is called only once.
C In the future it should be removed, now it is needed to avoid
C interference with legacy call from init_theory.f
      integer flagthinterface
      integer flaginit
      data flaginit /0/
      save flaginit

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

C OZ 17.10.17 TODO avoid this in the future
C (or stop execution if called not from theory interface)
      if(flagthinterface.eq.0 .and. flaginit.eq.1) then
        return
      endif
      if(flagthinterface.eq.1) then
        print *,'ABKM_init called from theory interface'
        flaginit = 1
      else
        print *,'ABKM_init called not from theory interface'
      endif

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

c 10.10.2017 Discussion with Sergey Alekhin:
c The parameter HQNONS drives the nonsinglet contribution to the charm production.
c It is infrared unsafe in the NNLO therefore there are pro and contra for including it and it is up to user.
c In ABMP16 fit it was set to .false.
c (makes small difference which reaches few % only at highest Q2 of the charm HERA data and is negligible for practical purposes)
      hqnons = .false.

C--------------------------------------------------------------------------
      end

      Subroutine ABKM_Set_Input_OrderFl(flordin)
C  OZ 1.10.2017 set O(alpha_S) for F_L
C---------------------------------------------------------------------------
      implicit none
C Input variables:
      integer flordin

C Common variables:
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

C-------------------------------
      kordfl = kordf2+flordin
C-------------------------------
      end
