C     ------------------------------------------------------------------
      program QEDevol
C     ------------------------------------------------------------------

C     Evolve PDFs using the nxn evolution toolbox of QCDNUM
C
C     We can do upto NNLO QCD and LO QED evolution
C     in FFNS and VFNS

C     ------------------------------------------------------------------
C     Declarations
C     ------------------------------------------------------------------

      implicit double precision (a-h,o-z)

C      include 'partonevolution.inc'
c      include 'mrst.inc'
      include 'apfel.inc'

C     Pdf output
      data ichk/1/                                  !yes/no check limits

C     ------------------------------------------------------------------
C     Declarations for the nxn evolution toolbox
C     ------------------------------------------------------------------
      parameter (nstoru = 20000000)                 !size of local store
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

C     Routines that set alphas and pdf starting values
      external AsVal1,AsVal2,AsVal3,AemVal1
      common /qcdqedord/ iordqcd,iordqed

      open(unit = 101,file='output/xvalues_qcdnum_'//setup//'.dat')
      open(unit = 102,file='output/qvalues_qcdnum_'//setup//'.dat')

      open(unit = 201,file='input/xdeltas_'//setup//'_init.dat')
      open(unit = 202,file='input/xsigma_'//setup//'_init.dat')
      open(unit = 203,file='input/xgluon_'//setup//'_init.dat')
      open(unit = 204,file='input/xphoton_'//setup//'_init.dat')
      open(unit = 205,file='input/xdeltav_'//setup//'_init.dat')
      open(unit = 206,file='input/xv_'//setup//'_init.dat')
      open(unit = 207,file='input/xdeltads_'//setup//'_init.dat')
      open(unit = 208,file='input/xdeltauc_'//setup//'_init.dat')
      open(unit = 209,file='input/xdeltasb_'//setup//'_init.dat')
      open(unit = 210,file='input/xdeltact_'//setup//'_init.dat')
      open(unit = 211,file='input/xvds_'//setup//'_init.dat')
      open(unit = 212,file='input/xvuc_'//setup//'_init.dat')
      open(unit = 213,file='input/xvsb_'//setup//'_init.dat')
      open(unit = 214,file='input/xvct_'//setup//'_init.dat')

      open(unit = 301,file='output/xdeltas_qcdnum_'//setup//'.dat')
      open(unit = 302,file='output/xsigma_qcdnum_'//setup//'.dat')
      open(unit = 303,file='output/xgluon_qcdnum_'//setup//'.dat')
      open(unit = 304,file='output/xphoton_qcdnum_'//setup//'.dat')
      open(unit = 305,file='output/xdeltav_qcdnum_'//setup//'.dat')
      open(unit = 306,file='output/xv_qcdnum_'//setup//'.dat')
      open(unit = 307,file='output/xdeltads_qcdnum_'//setup//'.dat')
      open(unit = 308,file='output/xdeltauc_qcdnum_'//setup//'.dat')
      open(unit = 309,file='output/xdeltasb_qcdnum_'//setup//'.dat')
      open(unit = 310,file='output/xdeltact_qcdnum_'//setup//'.dat')
      open(unit = 311,file='output/xvds_qcdnum_'//setup//'.dat')
      open(unit = 312,file='output/xvuc_qcdnum_'//setup//'.dat')
      open(unit = 313,file='output/xvsb_qcdnum_'//setup//'.dat')
      open(unit = 314,file='output/xvct_qcdnum_'//setup//'.dat')

C     lun = 6 stdout, -6 stdout w/out banner page
      lun = -6
      call qcinit(lun,' ')

C     Make x-mu2 grid
      call gxmake(xmin,iwt,5,nxin,nx,iosp)          !x-grid
      call gqmake(qq,wt,5,nqin,nq)                  !mu2-grid

C     Define renormalisation scale
      call setabr(amu2,bmu2)                        !renor scale

C     Order should not exceed that of the nxn weight calculation
      if(iordqcd+iordqed.gt.mxord) stop 'Evolution order too large'

C     Set evolution parameters
      call setord(mxord-1)                          !LO, NLO, NNLO
      call setalf(as0,r20)                          !input alphas
      iqc  = iqfrmq(q2c)                            !charm threshold
      iqb  = iqfrmq(q2b)                            !bottom threshold
      iqt  = iqfrmq(q2t)                            !top threshold
      call setcbt(nfin,iqc,iqb,iqt)                 !FFNS/VFNS
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

      print *,'--------------------------------------------------------'
      print *,''//scheme//' evolution'
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

      call setord(mxord-1)

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
        read(unit=201,fmt=*) start1(1,ix),start1(1,ix),start1(1,ix)
        read(unit=202,fmt=*) start1(2,ix),start1(2,ix),start1(2,ix)
        read(unit=203,fmt=*) start1(3,ix),start1(3,ix),start1(3,ix)
        read(unit=204,fmt=*) start1(4,ix),start1(4,ix),start1(4,ix)
        read(unit=205,fmt=*) start2(1,ix),start2(1,ix),start2(1,ix)
        read(unit=206,fmt=*) start2(2,ix),start2(2,ix),start2(2,ix)
        read(unit=207,fmt=*) start3(1,ix),start3(1,ix),start3(1,ix)
        read(unit=208,fmt=*) start4(1,ix),start4(1,ix),start4(1,ix)
        read(unit=209,fmt=*) start5(1,ix),start5(1,ix),start5(1,ix)
        read(unit=210,fmt=*) start6(1,ix),start6(1,ix),start6(1,ix)
        read(unit=211,fmt=*) start7(1,ix),start7(1,ix),start7(1,ix)
        read(unit=212,fmt=*) start8(1,ix),start8(1,ix),start8(1,ix)
        read(unit=213,fmt=*) start9(1,ix),start9(1,ix),start9(1,ix)
        read(unit=214,fmt=*) start10(1,ix),start10(1,ix),start10(1,ix)
      enddo

C     Evolve
      print *,''
      print *,'Start evolution...'
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
        if (iqlim(2).lt.nq) then
          do ix = 1,nx
            if (nf.eq.3.or.nf.eq.5) then
              start1(1,ix) = start1(1,ix)+getalfn(iqlim(2)-1,2,ierr)
     $        *(fcrossk(storu,idAiju(1),isetw,idf1(1),ix,iqlim(2)-1)
     $        +fcrossk(storu,idAiju(4),isetw,idf1(2),ix,iqlim(2)-1)
     $        +fcrossk(storu,idAiju(5),isetw,idf1(3),ix,iqlim(2)-1))
            else if (nf.eq.4) then
              start1(1,ix) = start1(1,ix)+getalfn(iqlim(2)-1,2,ierr)
     $        *(fcrossk(storu,idAiju(1),isetw,idf1(1),ix,iqlim(2)-1)
     $       -fcrossk(storu,idAiju(4),isetw,idf1(2),ix,iqlim(2)-1)
     $       -fcrossk(storu,idAiju(5),isetw,idf1(3),ix,iqlim(2)-1))
            endif
            start1(2,ix) = start1(2,ix)+getalfn(iqlim(2)-1,2,ierr)
     $      *(fcrossk(storu,idAiju(1),isetw,idf1(2),ix,iqlim(2)-1)
     $      +fcrossk(storu,idAiju(4),isetw,idf1(2),ix,iqlim(2)-1)
     $      +fcrossk(storu,idAiju(5),isetw,idf1(3),ix,iqlim(2)-1))
            start1(3,ix) = start1(3,ix)+getalfn(iqlim(2)-1,2,ierr)
     $      *(fcrossk(storu,idAiju(2),isetw,idf1(2),ix,iqlim(2)-1)
     $      +fcrossk(storu,idAiju(3),isetw,idf1(3),ix,iqlim(2)-1))
            start2(1,ix) = start2(1,ix)+getalfn(iqlim(2)-1,2,ierr)
     $      *fcrossk(storu,idAiju(1),isetw,idf2(1),ix,iqlim(2)-1)
            start2(2,ix) = start2(2,ix)+getalfn(iqlim(2)-1,2,ierr)
     $      *fcrossk(storu,idAiju(1),isetw,idf2(2),ix,iqlim(2)-1)
            start3(1,ix) = start3(1,ix)+getalfn(iqlim(2)-1,2,ierr)
     $      *fcrossk(storu,idAiju(1),isetw,idf3,ix,iqlim(2)-1)
            start7(1,ix) = start7(1,ix)+getalfn(iqlim(2)-1,2,ierr)
     $      *fcrossk(storu,idAiju(1),isetw,idf7,ix,iqlim(2)-1)
          enddo
        endif
      enddo

      do ix = 1,nx
        start4(1,ix) = 0.5D0*EvPdfij(storu,idf1(1),ix,iqc-1,1)
     $                 +0.5D0*EvPdfij(storu,idf1(2),ix,iqc-1,1)
        start5(1,ix) = -0.25D0*EvPdfij(storu,idf1(1),ix,iqb-1,1)
     $                 +0.25D0*EvPdfij(storu,idf1(2),ix,iqb-1,1)
     $                 -0.5D0*EvPdfij(storu,idf3(1),ix,iqb-1,1)
        start8(1,ix) = 0.5D0*EvPdfij(storu,idf2(1),ix,iqc-1,1)
     $                 +0.5D0*EvPdfij(storu,idf2(2),ix,iqc-1,1)
        start9(1,ix) = -0.25D0*EvPdfij(storu,idf2(1),ix,iqb-1,1)
     $                 +0.25D0*EvPdfij(storu,idf2(2),ix,iqb-1,1)
     $                 -0.5D0*EvPdfij(storu,idf7(1),ix,iqb-1,1)
      enddo

      iqlim(1) = iqc-1
      iqlim(2) = iqc-1
      nf = 1
      do while (nf.gt.0)
        iqlim(1) = iqlim(2)
        iqlim(2) = 99999
        call EvDglap(storu,idw4,ida3,idf4,start4,1,1,iqlim,nf,eps)
        call EvDglap(storu,idw8,ida3,idf8,start8,1,1,iqlim,nf,eps)
        if (iqlim(2).le.nq) then
          do ix = 1,nx
            if (nf.eq.3) then
              start4(1,ix) = start4(1,ix)+getalfn(iqlim(2)-1,2,ierr)
     $        *(fcrossk(storu,idAiju(1),isetw,idf4(1),ix,iqlim(2)-1)
     $        -fcrossk(storu,idAiju(4),isetw,idf1(2),ix,iqlim(2)-1)
     $        -fcrossk(storu,idAiju(5),isetw,idf1(3),ix,iqlim(2)-1))
            else
              start4(1,ix) = start4(1,ix)+getalfn(iqlim(2)-1,2,ierr)
     $        *fcrossk(storu,idAiju(1),isetw,idf4,ix,iqlim(2)-1)
            endif
            start8(1,ix) = start8(1,ix)+getalfn(iqlim(2)-1,2,ierr)
     $      *fcrossk(storu,idAiju(1),isetw,idf8,ix,iqlim(2)-1)
          enddo
        endif
      enddo

      do ix = 1,nx
        start6(1,ix) = +0.25D0*EvPdfij(storu,idf1(1),ix,iqt-1,1)
     $                 +0.25D0*EvPdfij(storu,idf1(2),ix,iqt-1,1)
     $                 -0.5D0*EvPdfij(storu,idf4(1),ix,iqt-1,1)
        start10(1,ix) = +0.25D0*EvPdfij(storu,idf2(1),ix,iqt-1,1)
     $                 +0.25D0*EvPdfij(storu,idf2(2),ix,iqt-1,1)
     $                 -0.5D0*EvPdfij(storu,idf8(1),ix,iqt-1,1)
      enddo

      iqlim(2) = iqb-1
      nf = 1
      do while (nf.gt.0)
        iqlim(1) = iqlim(2)
        iqlim(2) = 99999
        call EvDglap(storu,idw5,ida3,idf5,start5,1,1,iqlim,nf,eps)
        call EvDglap(storu,idw9,ida3,idf9,start9,1,1,iqlim,nf,eps)
        if (iqlim(2).le.nq) then
          do ix = 1,nx
            if (nf.eq.4) then
              start5(1,ix) = start5(1,ix)+getalfn(iqlim(2)-1,2,ierr)
     $        *(fcrossk(storu,idAiju(1),isetw,idf5(1),ix,iqlim(2)-1)
     $        -fcrossk(storu,idAiju(4),isetw,idf1(2),ix,iqlim(2)-1)
     $        -fcrossk(storu,idAiju(5),isetw,idf1(3),ix,iqlim(2)-1))
            else
              start5(1,ix) = start5(1,ix)+getalfn(iqlim(2)-1,2,ierr)
     $        *fcrossk(storu,idAiju(1),isetw,idf5,ix,iqlim(2)-1)
            endif
            start9(1,ix) = start9(1,ix)+getalfn(iqlim(2)-1,2,ierr)
     $      *fcrossk(storu,idAiju(1),isetw,idf9,ix,iqlim(2)-1)
          enddo
        endif
      enddo

      if (iqt.gt.0) then
        iqlim(2) = iqt-1
        nf = 1
        do while (nf.gt.0)
          iqlim(1) = iqlim(2)
          iqlim(2) = 99999
          call EvDglap(storu,idw6,ida3,idf6,start6,1,1,iqlim,nf,eps)
          call EvDglap(storu,idw10,ida3,idf10,start10,1,1,iqlim,nf,eps)
          if (iqlim(2).le.nq) then
            do ix = 1,nx
              if (nf.eq.5) then
                start6(1,ix) = start6(1,ix)+getalfn(iqlim(2)-1,2,ierr)
     $          *(fcrossk(storu,idAiju(1),isetw,idf6(1),ix,iqlim(2)-1)
     $          -fcrossk(storu,idAiju(4),isetw,idf1(2),ix,iqlim(2)-1)
     $          -fcrossk(storu,idAiju(5),isetw,idf1(3),ix,iqlim(2)-1))
              endif
              start10(1,ix) = start10(1,ix)+getalfn(iqlim(2)-1,2,ierr)
     $        *fcrossk(storu,idAiju(1),isetw,idf10,ix,iqlim(2)-1)
            enddo
          endif
        enddo
      endif

      print *,'Evolution finished.'

      do i = 1,nx
        write(unit=101,fmt='(f30.20)') xfrmix(i)
      enddo

      do i = 1,nq
        write(unit=102,fmt='(f30.20)') qfrmiq(i)
      enddo

      j = nq
        do i = 1,nx
          write(unit=301,fmt='(i4.1,f30.20,f30.20)')
     $    i,xfrmix(i),EvPdfij(storu,idf1(1),i,j,1)
          write(unit=302,fmt='(i4.1,f30.20,f30.20)')
     $    i,xfrmix(i),EvPdfij(storu,idf1(2),i,j,1)
          write(unit=303,fmt='(i4.1,f30.20,f30.20)')
     $    i,xfrmix(i),EvPdfij(storu,idf1(3),i,j,1)
          write(unit=304,fmt='(i4.1,f30.20,f30.20)')
     $    i,xfrmix(i),EvPdfij(storu,idf1(4),i,j,1)
          write(unit=305,fmt='(i4.1,f30.20,f30.20)')
     $    i,xfrmix(i),EvPdfij(storu,idf2(1),i,j,1)
          write(unit=306,fmt='(i4.1,f30.20,f30.20)')
     $    i,xfrmix(i),EvPdfij(storu,idf2(2),i,j,1)
          write(unit=307,fmt='(i4.1,f30.20,f30.20)')
     $    i,xfrmix(i),EvPdfij(storu,idf3(1),i,j,1)
          write(unit=308,fmt='(i4.1,f30.20,f30.20)')
     $    i,xfrmix(i),EvPdfij(storu,idf4(1),i,j,1)
          write(unit=309,fmt='(i4.1,f30.20,f30.20)')
     $    i,xfrmix(i),EvPdfij(storu,idf5(1),i,j,1)
          write(unit=310,fmt='(i4.1,f30.20,f30.20)')
     $    i,xfrmix(i),start6(1,i)
          write(unit=311,fmt='(i4.1,f30.20,f30.20)')
     $    i,xfrmix(i),EvPdfij(storu,idf7(1),i,j,1)
          write(unit=312,fmt='(i4.1,f30.20,f30.20)')
     $    i,xfrmix(i),EvPdfij(storu,idf8(1),i,j,1)
          write(unit=313,fmt='(i4.1,f30.20,f30.20)')
     $    i,xfrmix(i),EvPdfij(storu,idf9(1),i,j,1)
          write(unit=314,fmt='(i4.1,f30.20,f30.20)')
     $    i,xfrmix(i),start10(1,i)
        enddo

      close(unit=101)
      close(unit=102)
      close(unit=201)
      close(unit=202)
      close(unit=203)
      close(unit=204)
      close(unit=205)
      close(unit=206)
      close(unit=207)
      close(unit=208)
      close(unit=209)
      close(unit=210)
      close(unit=211)
      close(unit=212)
      close(unit=213)
      close(unit=214)
      close(unit=301)
      close(unit=302)
      close(unit=303)
      close(unit=304)
      close(unit=305)
      close(unit=306)
      close(unit=307)
      close(unit=308)
      close(unit=309)
      close(unit=310)
      close(unit=311)
      close(unit=312)
      close(unit=313)
      close(unit=314)
      end
