************************************************************************
*
*       Initialization routine for QEDEVOL
*
************************************************************************
      subroutine qedevol_ini

      implicit double precision (a-h,o-z)

#include "steering.inc"
#include "alphas.inc"
#include "couplings.inc"
#include "thresholds.inc"

      dimension xf(-6:7)

C     Evolution parameters
      data mxord/4/                                       !maximum order of evolution (QCD+QED)
      data iordqcd/2/                                     !QCD order: 1 - LO, 2 - NLO, 3 - NNLO
      data iordqed/1/                                     !QED order: 0 - no QED, 1 - LO
      data aem0/0.00119306554042/,rem20/3.157729D0/       !alphaem/(2*pi)

      common /aem/ aem0,rem20,q2b,q2t

C     ------------------------------------------------------------------
C     Declarations for the nxn evolution toolbox
C     ------------------------------------------------------------------
      parameter (nstoru = 1000000)                 !size of local store
      dimension storu(nstoru)                               !local store
      dimension iqlim(2)

      dimension idPiju(28,4),idAiju(5),idAlfa(4)      !identifier arrays
      dimension idw1(4,4,4),idf1(4),ida1(4,4,4)       !identifier arrays
      dimension idw2(2,2,4),idf2(2),ida2(2,2,4)
      dimension idw3(1,1,4),idf3(1),ida3(1,1,4)
      dimension idw4(1,1,4),idf4(1)
      dimension idw5(1,1,4),idf5(1)
      dimension idw6(1,1,4),idf6(1)
      dimension idw7(1,1,4),idf7(1)
      dimension idw8(1,1,4),idf8(1)
      dimension idw9(1,1,4),idf9(1)
      dimension idw10(1,1,4),idf10(1)
      dimension itypes(6)                                   !table types
      dimension itypes1(6)                                  !table types
      dimension start1(4,1000)
      dimension start2(2,1000)
      dimension start3(1,1000)
      dimension start4(1,1000)
      dimension start5(1,1000)
      dimension start6(1,1000)
      dimension start7(1,1000)
      dimension start8(1,1000)
      dimension start9(1,1000)
      dimension start10(1,1000)

      data itypes/6*0/                                 !initialise types
      data itypes1/6*0/                                !initialise types

      external AsVal1,AsVal2,AsVal3,AemVal1
      common /qcdqedord/ iordqcd,iordqed

      common /nnevol/ storu,idw1,ida1,idf1,idw2,ida2,idf2,
     $idw3,ida3,idf3,idw4,idf4,idw5,idf5,idw6,idf6,idw7,idf7,idw8,idf8,
     $idw9,idf9,idw10,idf10

      lun = 6 !stdout, -6 stdout w/out banner page

C     Order should not exceed that of the nxn weight calculation
      if(iordqcd+iordqed.gt.mxord) stop 'Evolution order too large'

C     Set evolution parameters
      call setord(iordqcd)
      call setint('nopt',444)
      alphas = hf_get_alphas(mz*mz)
      call setalf(0.1176d0,mz*mz)                     !input alphas
      call grpars(nx, xmi, xma, nq, qmi, qma, iord)

      q0 = starting_scale

C also add threshold values:

      q2c = qc
      q2b = qb
      q2t = qt

      iqc  = iqfrmq(q2c)                          !charm threshold
      iqb  = iqfrmq(q2b)                          !bottom threshold
      iqt  = iqfrmq(q2t)                          !top threshold
      iq0  = iqfrmq(q0)                             !starting scale

C     ------------------------------------------------------------------
C     Do evolution with the nxn evolution toolbox
C     ------------------------------------------------------------------

C     Weight tables
      itypes(1) = 5
      itypes(2) = 84
      itypes1(2) = 28
C     Put 14 pdf and 4 alpha tables in the store
      itypes(5) = 14
      itypes(6) = 4

C     isetw is the table set identifier assigned by QCDNUM
      call MakeTab(storu,nstoru,itypes,0,0,isetw,nwordsu)
      call MakeTab(storu,nstoru,itypes1,0,0,isetw1,nwordsu)

      print *,'+-------------------------------------------------------'
      print *,'Starting scale and flavour thresholds indexes: ',
     $iq0,iqc,iqb,iqt
      print *,'Starting scale and flavour thresholds values: ',
     $qfrmiq(iq0),qfrmiq(iqc),qfrmiq(iqb),qfrmiq(iqt)
      print *,'--------------------------------------------------------'

C     Calculate  evolution weigths
      call FilWTqcd(storu,isetw,idpiju,idaiju,iordqcd)
      call FilWTqed(storu,isetw1,idpiju,iordqed)

C     Fill tables of alphas and alphaem values
      idAlfa(1) = 1000*isetw+601
      idAlfa(2) = 1000*isetw+602
      idAlfa(3) = 1000*isetw+603
      idAlfa(4) = 1000*isetw+604

      call EvFillA(storu,idAlfa(1),AsVal1)  !LO QCD
      call EvFillA(storu,idAlfa(2),AsVal2)  !NLO QCD
      call EvFillA(storu,idAlfa(3),AsVal3)  !NNLO QCD
      call EvFillA(storu,idAlfa(4),AemVal1) !LO QED

C     Setup the identifiers for nxn evolution
      ityp = 0
      do i = 1,4
        do j = 1,4
          ityp = ityp+1
          do k = 1,mxord
            ida1(i,j,k) = idAlfa(k)
            idw1(i,j,k) = idPiju(ityp,k)
          enddo
        enddo
      enddo

      do i = 1,2
        do j = 1,2
          ityp = ityp+1
          do k = 1,mxord
            ida2(i,j,k) = idAlfa(k)
            idw2(i,j,k) = idPiju(ityp,k)
          enddo
        enddo
      enddo

      do k = 1,mxord
        ida3(1,1,k) = idAlfa(k)
        idw3(1,1,k) = idPiju(21,k)
        idw4(1,1,k) = idPiju(22,k)
        idw5(1,1,k) = idPiju(23,k)
        idw6(1,1,k) = idPiju(24,k)
        idw7(1,1,k) = idPiju(25,k)
        idw8(1,1,k) = idPiju(26,k)
        idw9(1,1,k) = idPiju(27,k)
        idw10(1,1,k) = idPiju(28,k)
      enddo

C     PDF table identifiers
      idf1(1) = 1000*isetw+501                                  !Delta_S
      idf1(2) = 1000*isetw+502                                    !Sigma
      idf1(3) = 1000*isetw+503                                    !gluon
      idf1(4) = 1000*isetw+504                                   !photon
      idf2(1) = 1000*isetw+505                                  !Delta_V
      idf2(2) = 1000*isetw+506                                        !V
      idf3(1) = 1000*isetw+507                                 !Delta_ds
      idf4(1) = 1000*isetw+508                                 !Delta_uc
      idf5(1) = 1000*isetw+509                                 !Delta_sb
      idf6(1) = 1000*isetw+510                                 !Delta_ct
      idf7(1) = 1000*isetw+511                                     !V_ds
      idf8(1) = 1000*isetw+512                                     !V_uc
      idf9(1) = 1000*isetw+513                                     !V_sb
      idf10(1) = 1000*isetw+514                                    !V_ct

      do ix = 1,nx
        x = xfrmix(ix)
        call ExternalSetQEDEVOL(x,q0,xf)
        start1(1,ix)=xf(-2)+xf(2)-((xf(-1)+xf(1))+(xf(-3)+xf(3)))
        start1(2,ix)=xf(-2)+xf(2)+((xf(-1)+xf(1))+(xf(-3)+xf(3)))
        start1(3,ix)=xf(0)
        start2(1,ix)=xf(2)-xf(-2)-((xf(1)-xf(-1))+(xf(3)-xf(-3)))
        start2(2,ix)=xf(2)-xf(-2)+((xf(1)-xf(-1))+(xf(3)-xf(-3)))
        start3(1,ix)=xf(1)+xf(-1)-(xf(3)+xf(-3))
        start7(1,ix)=xf(1)-xf(-1)-(xf(3)-xf(-3))
        start1(4,ix)=xf(7)
      enddo



      iqlim(1) = iq0
      iqlim(2) = iq0
      nf = 1
      do while (nf.gt.0)
        iqlim(1) = iqlim(2)
        iqlim(2) = 99999
        call EvDglap(storu,idw1,ida1,idf1,start1,4,4,iqlim,nf,eps)
        call EvDglap(storu,idw2,ida2,idf2,start2,2,2,iqlim,nf,eps)
        call EvDglap(storu,idw3,ida3,idf3,start3,1,1,iqlim,nf,eps)
        call EvDglap(storu,idw7,ida3,idf7,start7,1,1,iqlim,nf,eps)
      enddo

      do ix = 1,nx
        start4(1,ix) = 0.5D0*EvPdfij(storu,idf1(1),ix,iqc,1)
     $                 +0.5D0*EvPdfij(storu,idf1(2),ix,iqc,1)
        start8(1,ix) = 0.5D0*EvPdfij(storu,idf2(1),ix,iqc,1)
     $                 +0.5D0*EvPdfij(storu,idf2(2),ix,iqc,1)
        start5(1,ix) = -0.25D0*EvPdfij(storu,idf1(1),ix,iqb,1)
     $                 +0.25D0*EvPdfij(storu,idf1(2),ix,iqb,1)
     $                 -0.5D0*EvPdfij(storu,idf3(1),ix,iqb,1)
        start9(1,ix) = -0.25D0*EvPdfij(storu,idf2(1),ix,iqb,1)
     $                 +0.25D0*EvPdfij(storu,idf2(2),ix,iqb,1)
     $                 -0.5D0*EvPdfij(storu,idf7(1),ix,iqb,1)
      enddo

      iqlim(1) = iqc
      iqlim(2) = iqc
      nf = 1
      do while (nf.gt.0)
        iqlim(1) = iqlim(2)
        iqlim(2) = 99999
        call EvDglap(storu,idw4,ida3,idf4,start4,1,1,iqlim,nf,eps)
        call EvDglap(storu,idw8,ida3,idf8,start8,1,1,iqlim,nf,eps)
      enddo

      if (iqt.gt.0) then
         do ix = 1,nx
            start6(1,ix) = +0.25D0*EvPdfij(storu,idf1(1),ix,iqt,1)
     $           +0.25D0*EvPdfij(storu,idf1(2),ix,iqt,1)
     $           -0.5D0*EvPdfij(storu,idf4(1),ix,iqt,1)
            start10(1,ix) = +0.25D0*EvPdfij(storu,idf2(1),ix,iqt,1)
     $           +0.25D0*EvPdfij(storu,idf2(2),ix,iqt,1)
     $           -0.5D0*EvPdfij(storu,idf8(1),ix,iqt,1)
         enddo
      endif

      iqlim(1) = iqb
      iqlim(2) = iqb
      nf = 1
      do while (nf.gt.0)
        iqlim(1) = iqlim(2)
        iqlim(2) = 99999
        call EvDglap(storu,idw5,ida3,idf5,start5,1,1,iqlim,nf,eps)
        call EvDglap(storu,idw9,ida3,idf9,start9,1,1,iqlim,nf,eps)
      enddo


      if (iqt.gt.0) then
         iqlim(1) = iqt
         iqlim(2) = iqt
         nf = 1
         do while (nf.gt.0)
            iqlim(1) = iqlim(2)
            iqlim(2) = 99999
            call EvDglap(storu,idw6,ida3,idf6,start6,1,1,iqlim,nf,eps)
            call EvDglap(storu,idw10,ida3,idf10,start10,1,1,iqlim,nf,eps)
         enddo
      endif

c      call dumptab(storu,isetw,11,'qcdweights.wt','')
c      call dumptab(storu,isetw1,12,'qedweights.wt','')

      return
      end

************************************************************************
*
*       Evolution routine for QEDEVOL
*
************************************************************************
      subroutine qedevol_main

      implicit double precision (a-h,o-z)

#include "steering.inc"
#include "alphas.inc"
#include "couplings.inc"
#include "thresholds.inc"

      dimension xf(-6:7)

C     evolution parameters
      data mxord/4/                                          !maximum order of evolution
      data iordqcd/2/                                        !QCD order: 1 - LO, 2 - NLO, 3 - NNLO
      data iordqed/1/                                        !QED order: 0 - no QED, 1 - LO
      data aem0/0.00119306554042/,rem20/3.157729D0/          !alphaem/(2*pi)

      common /aem/ aem0,rem20,q2b,q2t

C     Pdf output
c      data ichk/1/                                  !yes/no check limits

C     ------------------------------------------------------------------
C     Declarations for the nxn evolution toolbox
C     ------------------------------------------------------------------
      parameter (nstoru = 1000000)                 !size of local store
      dimension storu(nstoru)                               !local store
      dimension iqlim(2)

      dimension idPiju(28,4),idAiju(5),idAlfa(4)      !identifier arrays
      dimension idw1(4,4,4),idf1(4),ida1(4,4,4)       !identifier arrays
      dimension idw2(2,2,4),idf2(2),ida2(2,2,4)
      dimension idw3(1,1,4),idf3(1),ida3(1,1,4)
      dimension idw4(1,1,4),idf4(1)
      dimension idw5(1,1,4),idf5(1)
      dimension idw6(1,1,4),idf6(1)
      dimension idw7(1,1,4),idf7(1)
      dimension idw8(1,1,4),idf8(1)
      dimension idw9(1,1,4),idf9(1)
      dimension idw10(1,1,4),idf10(1)
      dimension idf(0:13)
      dimension itypes(6)                                   !table types
      dimension itypes1(6)                                  !table types
      dimension start1(4,1000)
      dimension start2(2,1000)
      dimension start3(1,1000)
      dimension start4(1,1000)
      dimension start5(1,1000)
      dimension start6(1,1000)
      dimension start7(1,1000)
      dimension start8(1,1000)
      dimension start9(1,1000)
      dimension start10(1,1000)
      dimension def(-6:6,12)
      data def /
C--   tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--   -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6..
     + 1.,-1., 1.,-1., 1.,-1., 0.,-1., 1.,-1., 1.,-1., 1.,   !Delta_S
     + 1., 1., 1., 1., 1., 1., 0., 1., 1., 1., 1., 1., 1.,   !Sigma
     +-1., 1.,-1., 1.,-1., 1., 0.,-1., 1.,-1., 1.,-1., 1.,   !Delta_V
     +-1.,-1.,-1.,-1.,-1.,-1., 0., 1., 1., 1., 1., 1., 1.,   !V
     + 0., 0., 0.,-1., 0., 1., 0., 1., 0.,-1., 0., 0., 0.,   !Delta_ds
     + 0., 0., 0., 1., 0.,-1., 0., 1., 0.,-1., 0., 0., 0.,   !V_ds
     + 0., 0.,-1., 0., 1., 0., 0., 0., 1., 0.,-1., 0., 0.,   !Delta_uc
     + 0., 0., 1., 0.,-1., 0., 0., 0., 1., 0.,-1., 0., 0.,   !V_uc
     + 0.,-1., 0., 1., 0., 0., 0., 0., 0., 1., 0.,-1., 0.,   !Delta_sb
     + 0., 1., 0.,-1., 0., 0., 0., 0., 0., 1., 0.,-1., 0.,   !V_sb
     +-1., 0., 1., 0., 0., 0., 0., 0., 0., 0., 1., 0.,-1.,   !Delta_ct
     + 1., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0.,-1. /  !V_ct
      data itypes/6*0/                                 !initialise types
      data itypes1/6*0/                                !initialise types

      external AsVal1,AsVal2,AsVal3,AemVal1
      common /qcdqedord/ iordqcd,iordqed

      common /nnevol/ storu,idw1,ida1,idf1,idw2,ida2,idf2,
     $idw3,ida3,idf3,idw4,idf4,idw5,idf5,idw6,idf6,idw7,idf7,idw8,idf8,
     $idw9,idf9,idw10,idf10

      idf(0) = idf1(3)     !gluon
      idf(1) = idf1(1)     !Delta_S
      idf(2) = idf1(2)     !Sigma
      idf(3) = idf2(1)     !Delta_V
      idf(4) = idf2(2)     !V
      idf(5) = idf3(1)     !Delta_ds
      idf(6) = idf7(1)     !V_ds
      idf(7) = idf4(1)     !Delta_uc
      idf(8) = idf8(1)     !V_uc
      idf(9) = idf5(1)     !Delta_sb
      idf(10) = idf9(1)    !V_sb
      idf(11) = idf6(1)    !Delta_ct
      idf(12) = idf10(1)   !V_ct
      idf(13) = idf1(4)    !photon

      lun = 6 !stdout, -6 stdout w/out banner page

C     Order should not exceed that of the nxn weight calculation
      if(iordqcd+iordqed.gt.mxord) stop 'Evolution order too large'

C     Set evolution parameters
c      call setord(iordqcd)
c      call setint('nopt',444)
c      alphas = hf_get_alphas(mz*mz)
      call grpars(nx, xmi, xma, nq, qmi, qma, iord)
c      call setalf(alphas,mz*mz)
      q0 = starting_scale
      q2c = qc
      q2b = qb
      q2t = qt
      iqc  = iqfrmq(q2c)                           !charm threshold
      iqb  = iqfrmq(q2b)                          !bottom threshold
      iqt  = iqfrmq(q2t)                          !top threshold
      iq0  = iqfrmq(q0)                             !starting scale

C     ------------------------------------------------------------------
C     Do evolution with the nxn evolution toolbox
C     ------------------------------------------------------------------

C     Weight tables
      itypes(1) = 5
      itypes(2) = 84
      itypes1(2) = 28
C     Put 14 pdf and 4 alpha tables in the store
      itypes(5) = 14
      itypes(6) = 4


      do ix = 1,nx
        x = xfrmix(ix)
        if (x.ne.0) then 
           call ExternalSetQEDEVOL(x,q0,xf)
           start1(1,ix)=xf(-2)+xf(2)-((xf(-1)+xf(1))+(xf(-3)+xf(3)))
           start1(2,ix)=xf(-2)+xf(2)+((xf(-1)+xf(1))+(xf(-3)+xf(3)))
           start1(3,ix)=xf(0)
           start2(1,ix)=xf(2)-xf(-2)-((xf(1)-xf(-1))+(xf(3)-xf(-3)))
           start2(2,ix)=xf(2)-xf(-2)+((xf(1)-xf(-1))+(xf(3)-xf(-3)))
           start3(1,ix)=xf(1)+xf(-1)-(xf(3)+xf(-3))
           start7(1,ix)=xf(1)-xf(-1)-(xf(3)-xf(-3))
           start1(4,ix)=xf(7)
        endif
      enddo

      iqlim(1) = iq0
      iqlim(2) = iq0
      nf = 1
      do while (nf.gt.0)
        iqlim(1) = iqlim(2)
        iqlim(2) = 99999
        call EvDglap(storu,idw1,ida1,idf1,start1,4,4,iqlim,nf,eps)
        call EvDglap(storu,idw2,ida2,idf2,start2,2,2,iqlim,nf,eps)
        call EvDglap(storu,idw3,ida3,idf3,start3,1,1,iqlim,nf,eps)
        call EvDglap(storu,idw7,ida3,idf7,start7,1,1,iqlim,nf,eps)
      enddo


      do ix = 1,nx
        start4(1,ix) = 0.5D0*EvPdfij(storu,idf1(1),ix,iqc,1)
     $                 +0.5D0*EvPdfij(storu,idf1(2),ix,iqc,1)
        start8(1,ix) = 0.5D0*EvPdfij(storu,idf2(1),ix,iqc,1)
     $                 +0.5D0*EvPdfij(storu,idf2(2),ix,iqc,1)
        start5(1,ix) = -0.25D0*EvPdfij(storu,idf1(1),ix,iqb,1)
     $                 +0.25D0*EvPdfij(storu,idf1(2),ix,iqb,1)
     $                 -0.5D0*EvPdfij(storu,idf3(1),ix,iqb,1)
        start9(1,ix) = -0.25D0*EvPdfij(storu,idf2(1),ix,iqb,1)
     $                 +0.25D0*EvPdfij(storu,idf2(2),ix,iqb,1)
     $                 -0.5D0*EvPdfij(storu,idf7(1),ix,iqb,1)
      enddo

      iqlim(1) = iqc
      iqlim(2) = iqc
      nf = 1
      do while (nf.gt.0)
        iqlim(1) = iqlim(2)
        iqlim(2) = 99999
        call EvDglap(storu,idw4,ida3,idf4,start4,1,1,iqlim,nf,eps)
        call EvDglap(storu,idw8,ida3,idf8,start8,1,1,iqlim,nf,eps)
      enddo


      if (iqt.gt.0) then
         do ix = 1,nx
            start6(1,ix) = +0.25D0*EvPdfij(storu,idf1(1),ix,iqt,1)
     $           +0.25D0*EvPdfij(storu,idf1(2),ix,iqt,1)
     $           -0.5D0*EvPdfij(storu,idf4(1),ix,iqt,1)
            start10(1,ix) = +0.25D0*EvPdfij(storu,idf2(1),ix,iqt,1)
     $           +0.25D0*EvPdfij(storu,idf2(2),ix,iqt,1)
     $           -0.5D0*EvPdfij(storu,idf8(1),ix,iqt,1)
         enddo
      endif
 
      iqlim(1) = iqb
      iqlim(2) = iqb
      nf = 1
      do while (nf.gt.0)
        iqlim(1) = iqlim(2)
        iqlim(2) = 99999
        call EvDglap(storu,idw5,ida3,idf5,start5,1,1,iqlim,nf,eps)
        call EvDglap(storu,idw9,ida3,idf9,start9,1,1,iqlim,nf,eps)
      enddo


      if (iqt.gt.0) then
         iqlim(1) = iqt
         iqlim(2) = iqt
         nf = 1
         do while (nf.gt.0)
            iqlim(1) = iqlim(2)
            iqlim(2) = 99999
            call EvDglap(storu,idw6,ida3,idf6,start6,1,1,iqlim,nf,eps)
            call EvDglap(storu,idw10,ida3,idf10,start10,1,1,iqlim,nf,eps)
         enddo
      endif

      call EVPCOPY (storu, idf, def, 1, 8)

      return
      end

      Subroutine QEDEVOLsubr(x, qmu2, xf)
C-------------------------------------------------------
C
C External PDF reading for QEDEVOL
C
C--------------------------------------------------------
      implicit double precision (a-h,o-z)
*
#include "steering.inc"
#include "thresholds.inc"
      double precision x,qmu2
      double precision xdelta,xsigma,xgluon,xphoton,xdeltav,xv,xdeltads,
     $ xdeltauc,xdeltasb,xvds,xvuc,xvsb
      dimension xf(-6:7)

      parameter (nstoru = 1000000)                 !size of local store
      dimension storu(nstoru)
      dimension idw1(4,4,4),idf1(4),ida1(4,4,4)
      dimension idw2(2,2,4),idf2(2),ida2(2,2,4)
      dimension idw3(1,1,4),idf3(1),ida3(1,1,4)
      dimension idw4(1,1,4),idf4(1)
      dimension idw5(1,1,4),idf5(1)
      dimension idw6(1,1,4),idf6(1)
      dimension idw7(1,1,4),idf7(1)
      dimension idw8(1,1,4),idf8(1)
      dimension idw9(1,1,4),idf9(1)
      dimension idw10(1,1,4),idf10(1)
     
      common /nnevol/ storu,idw1,ida1,idf1,idw2,ida2,idf2,idw3,
     $ida3,idf3,idw4,idf4,idw5,idf5,idw6,idf6,idw7,idf7,idw8,idf8,
     $idw9,idf9,idw10,idf10

      call evtable(storu,idf1(1),x,1,qmu2,1,xdelta,1)
      call evtable(storu,idf1(2),x,1,qmu2,1,xsigma,1)
      call evtable(storu,idf1(3),x,1,qmu2,1,xgluon,1)
      call evtable(storu,idf1(4),x,1,qmu2,1,xphoton,1)
      call evtable(storu,idf2(1),x,1,qmu2,1,xdeltav,1)
      call evtable(storu,idf2(2),x,1,qmu2,1,xv,1)
      call evtable(storu,idf3(1),x,1,qmu2,1,xdeltads,1)
      call evtable(storu,idf4(1),x,1,qmu2,1,xdeltauc,1)
      call evtable(storu,idf5(1),x,1,qmu2,1,xdeltasb,1)


      iqt  = iqfrmq(qt)                          !top threshold

      if (iqt.gt.0) then
         call evtable(storu,idf6(1),x,1,qmu2,1,xdeltact,1)
      endif
      call evtable(storu,idf7(1),x,1,qmu2,1,xvds,1)
      call evtable(storu,idf8(1),x,1,qmu2,1,xvuc,1)
      call evtable(storu,idf9(1),x,1,qmu2,1,xvsb,1)


      if (iqt.gt.0) then
         call evtable(storu,idf10(1),x,1,qmu2,1,xvct,1)
      endif

      do i = -6,7
        xf(i) = 0d0
      enddo

      xf(0) = xgluon
      xf(7) = xphoton
      xf(1) = 0.125d0*(xsigma-xdelta+xv-xdeltav+2.d0*xdeltads+2.d0*xvds)
      xf(-1)= 0.125d0*(xsigma-xdelta-xv+xdeltav+2.d0*xdeltads-2.d0*xvds)
      xf(2) = 0.25d0*(xsigma+xdelta+xv+xdeltav)
      xf(-2) = 0.25d0*(xsigma+xdelta-xv-xdeltav)
      xf(3) = 0.125d0*(xsigma-xdelta+xv-xdeltav-2.d0*xdeltads-2.d0*xvds)
      xf(-3)= 0.125d0*(xsigma-xdelta-xv+xdeltav-2.d0*xdeltads+2.d0*xvds)

      if (qmu2.gt.hf_mass(1)**2) then
      xf(1) = 0.125d0*(xsigma-xdelta+xv-xdeltav+2.d0*xdeltads+2.d0*xvds)
      xf(-1)= 0.125d0*(xsigma-xdelta-xv+xdeltav+2.d0*xdeltads-2.d0*xvds)
      xf(2) = 0.125d0*(xsigma+xdelta+xv+xdeltav+2.d0*xdeltauc+2.d0*xvuc)
      xf(-2) =0.125d0*(xsigma+xdelta-xv-xdeltav+2.d0*xdeltauc-2.d0*xvuc)
      xf(3) = 0.125d0*(xsigma-xdelta+xv-xdeltav-2.d0*xdeltads-2.d0*xvds)
      xf(-3)= 0.125d0*(xsigma-xdelta-xv+xdeltav-2.d0*xdeltads+2.d0*xvds)
      xf(4) = 0.125d0*(xsigma+xdelta+xv+xdeltav-2.d0*xdeltauc-2.d0*xvuc)
      xf(-4) =0.125d0*(xsigma+xdelta-xv-xdeltav-2.d0*xdeltauc+2.d0*xvuc)
      endif

      if (qmu2.gt.hf_mass(2)**2) then
      xf(1) = 1.d0/12.d0*(xsigma-xdelta+xv-xdeltav
     $+4.d0*xdeltads+4.d0*xvds+2.d0*xdeltasb+2.d0*xvsb)
      xf(-1)= 1.d0/12.d0*(xsigma-xdelta-xv+xdeltav
     $+4.d0*xdeltads-4.d0*xvds+2.d0*xdeltasb-2.d0*xvsb)
      xf(3) = 1.d0/12.d0*(xsigma-xdelta+xv-xdeltav
     $-2.d0*xdeltads-2.d0*xvds+2.d0*xdeltasb+2.d0*xvsb)
      xf(-3)= 1.d0/12.d0*(xsigma-xdelta-xv+xdeltav
     $-2.d0*xdeltads+2.d0*xvds+2.d0*xdeltasb-2.d0*xvsb)
      xf(5) = 1.d0/12.d0*(xsigma-xdelta+xv-xdeltav
     $-2.d0*xdeltads-2.d0*xvds-4.d0*xdeltasb-4.d0*xvsb)
      xf(-5)= 1.d0/12.d0*(xsigma-xdelta-xv+xdeltav
     $-2.d0*xdeltads+2.d0*xvds-4.d0*xdeltasb+4.d0*xvsb)
      endif
     
      if (qmu2.gt.hf_mass(3)**2) then
      xf(2) = 1.d0/12.d0*(xsigma+xdelta+xv+xdeltav
     $+4.d0*xdeltauc+4.d0*xvuc+2.d0*xdeltact+2.d0*xvct)
      xf(-2)= 1.d0/12.d0*(xsigma+xdelta-xv-xdeltav
     $+4.d0*xdeltauc-4.d0*xvuc+2.d0*xdeltact-2.d0*xvct)
      xf(4) = 1.d0/12.d0*(xsigma+xdelta+xv+xdeltav
     $-2.d0*xdeltauc-2.d0*xvuc+2.d0*xdeltact+2.d0*xvct)
      xf(-4)= 1.d0/12.d0*(xsigma+xdelta-xv-xdeltav
     $-2.d0*xdeltauc+2.d0*xvuc+2.d0*xdeltact-2.d0*xvct)
      xf(6) = 1.d0/12.d0*(xsigma+xdelta+xv+xdeltav
     $-2.d0*xdeltauc-2.d0*xvuc-4.d0*xdeltact-4.d0*xvct)
      xf(-6)= 1.d0/12.d0*(xsigma+xdelta-xv-xdeltav
     $-2.d0*xdeltauc+2.d0*xvuc-4.d0*xdeltact+4.d0*xvct)
      endif

      return
      end

      subroutine ExternalSetQEDEVOL(x,q0,xf)
*
      implicit none
#include "steering.inc"
**
*     Input Variables
*
      double precision x
      double precision q0
**
*     Internal Variables
*
      integer ipdf
      double precision gluon
      double precision photon
      double precision pdf_from_text
      double precision qstrange,Ubar,Dbar,H1U,H1D
      double precision sea,dbmub,dval,uval
      double precision dfac,ParDumpFactor
      parameter(ParDumpFactor=1.d-3)
**
*     Output Variables
*
      double precision xf(-6:7)
*
*     Set PDFs to zero
*
      do ipdf=-6,7
         xf(ipdf) = 0d0
      enddo
      if(x.gt.1d0) x = 1d0
c      print *,photon(x)
*
*     Construct PDFs addording to the PDF decomposition
*

      if(PDF_DECOMPOSITION.eq.'LHAPDF')then
c         q0 = sqrt(starting_scale)
         call evolvePDF(x, q0, xf)

      elseif(PDF_DECOMPOSITION.eq.'QCDNUM_GRID')then
         xf(-3) = ( pdf_from_text(x,3) - pdf_from_text(x,6) ) / 2d0
         xf(-2) = pdf_from_text(x,4)
         xf(-1) = pdf_from_text(x,5)
         xf(0)  = pdf_from_text(x,0)
         xf(1)  = pdf_from_text(x,1) - pdf_from_text(x,5)
         xf(2)  = pdf_from_text(x,2) - pdf_from_text(x,4)
         xf(3)  = ( pdf_from_text(x,3) + pdf_from_text(x,6) ) / 2d0

      elseif(Index(PDF_DECOMPOSITION,'D_U_Dbar_Ubar').gt.0)then ! D,U,Dbar,Ubar.
         xf(-3) = qstrange(x)
         xf(-2) = Ubar(x)
         xf(-1) = Dbar(x)
         xf(0)  = gluon(x)
         xf(7)  = photon(x)
         xf(1)  = H1D(x) - xf(-3)
         xf(2)  = H1U(x)
         xf(3)  = xf(-3)

      elseif(Index(PDF_DECOMPOSITION,'Sea').gt.0)then
         xf(-2) = sea(x) / 4d0 - dbmub(x) / 2d0
         xf(-1) = sea(x) / 4d0 + dbmub(x) / 2d0
         xf(0)  = gluon(x)
         xf(7)  = photon(x)
         xf(1)  = dval(x) + xf(-1)
         xf(2)  = uval(x) + xf(-2)

      elseif(PDF_DECOMPOSITION.eq.'Diffractive')then
         dfac = dexp(-ParDumpFactor/(1.00001d0-x))
*
         xf(-3) = dfac * Uval(x)
         xf(-2) = xf(-3)
         xf(-1) = xf(-3)
         xf(0)  = dfac * gluon(x)
         xf(1)  = xf(-3)
         xf(2)  = xf(-3)
         xf(3)  = xf(-3)

      elseif(Index(PDF_DECOMPOSITION,'Dbar_Ubar').gt.0)then
         xf(-3) = qstrange(x)
         xf(-2) = ubar(x)
         xf(-1) = dbar(x) - xf(-3)
         xf(0)  = gluon(x)
         xf(7)  = photon(x)
         xf(1)  = dval(x) + xf(-1)
         xf(2)  = uval(x) + xf(-2)
         xf(3)  = xf(-3)

      else
         print *,'Unknown PDF Decomposition: '//PDF_DECOMPOSITION
         print *,'Stop in evolution'
         call HF_Stop
      endif
*
      return
      end
