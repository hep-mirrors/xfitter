      subroutine sf_abkm_wrap(x,q2,f2abkm,flabkm,f3abkm,f2cabkm,
     $    flcabkm,f3cabkm,f2babkm,flbabkm,f3babkm,nb,nt,ni)
C-------------------------------------------------------------------------
C
C Created by RP 09 Jan 2012. A wraper around abkm functions
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
      integer nb,nt,ni

C-------------------------------------------------------------------------
c      b1=f2qcd(3,1,22,xb,q2)
c      f2qcd(nb,nt,ni,xb,q2)
      real*8 f2qcd,flqcd,f3qcd,f2charm_ffn,flcharm_ffn
      real*8 f2nucharm,f3nucharm,ftnucharm

      f2abkm = f2qcd(nb,nt,ni,x,q2)
      flabkm = flqcd(nb,nt,ni,x,q2)

c in NC no f3 ST 
      if(nb.eq.3.and.ni.eq.22) then
        f3abkm  = 0.0d0
c c quark
        f2cabkm = f2charm_ffn(x,q2,8)
        flcabkm = flcharm_ffn(x,q2,8)
        f3cabkm = 0.0d0
c b quark
        f2babkm = f2charm_ffn(x,q2,10)
        flbabkm = flcharm_ffn(x,q2,10)
        f3babkm = 0.0d0
c in CC there is f3

      elseif(nb.eq.6.or.nb.eq.7.and.ni.eq.24) then
        f3abkm  = f3qcd(nb,nt,ni,x,q2)
c c quark
        f2cabkm =f2nucharm(nb,nt,ni,x,q2,8)
        flcabkm =f2nucharm(nb,nt,ni,x,q2,8)-ftnucharm(nb,nt,ni,x,q2,8)
        f3cabkm =f3nucharm(nb,nt,ni,x,q2,8)
c b quark
        f2babkm =f2nucharm(nb,nt,ni,x,q2,10)
        flbabkm =f2nucharm(nb,nt,ni,x,q2,10)-ftnucharm(nb,nt,ni,x,q2,10)
        f3babkm =f3nucharm(nb,nt,ni,x,q2,10)
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
      integer npdftot,kordpdf,kschemepdf,kpdfset
c      common /forpdfset/ kschemepdf,kordpdf
      common /forpdfset/ npdftot,kordpdf,kschemepdf,kpdfset

      double precision rmass,rmassp,rcharge
      COMMON /MASSES/ rmass(150),rmassp(50),rcharge(150)

      double precision q2rep,q2s,q20alphas,alphas0,alpsz,alpss
      double precision alpsc,alpsb,alpst,tscale,rscale,fscale,hqscale1
      double precision hqscale2
      integer nfeff,kordalps,kfeff,kordhq,kordf2,kordfl,kordf3

      common /FORALPSRENORM/ q2rep,q2s,q20alphas,alphas0,alpsz,alpss
     , ,alpsc,alpsb,alpst,tscale,rscale,fscale,hqscale1,hqscale2
     , ,nfeff,kordalps,kfeff
     , ,kordhq,kordf2,kordfl,kordf3

      logical msbarm,hqnons
      common /forschemedef/ msbarm,hqnons
C-------------------------------      

      rmass(8)  = rmass8in
      rmass(10) = rmass10in
      kschemepdf= kschemepdfin
c set same order for pdf, light, heavy quarks      
      kordpdf   = kordpdfin
      kordhq    = kordpdfin
      kordf2    = kordpdfin
      kordfl    = kordpdfin
      kordf3    = kordpdfin
      kordalps  = kordpdfin
c run in running m scheme
      msbarm    = msbarmin 
      hqscale1  = hqscale1in
      hqscale2  = hqscale2in
       
C--------------------------------------------------------------------------
      end

