C------------------
      SUBROUTINE APEQSOL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'APSCOM6.'
      INCLUDE 'PDFCOM.'
      include 'CONSTCOM.'

      common  /forqgspl/ xb0,q2,xlog,an,an2,an3,kn,iqn,ixn,isn,ihqn

      real*8 dert(-5:2,nflim,-nxmlim:nxplim)
      real*8 fsp(nxtot),bs(nxtot),cs(nxtot),ds(nxtot),xx(nxtot)

      integer kpdfseq(8),kpdfseqi(8),nfcpdfseq (8)
!  A sequence of the PDF evolution in different schemes: 
!  the 3-flavour PDFs with the fixed number of flavours 
!  in the loops (kpdf=0), the 3-flavour PDFs with the 
!  variable number of flavours in the loops (kpdf=-1),  
!  the 4-flavour PDFs (kpdf=1), and the 5-flavour PDFs (kpdf=2). 
      data kpdfseq /-1, 0, 1, 2, -2, -3, -4, -5/
!..........
      data kpdfseqi /0, 0, 0, 1, 0, 0, -2, -3/
!..................
      data nfcpdfseq /3, 3, 4, 5, 4, 4, 5, 5/

      if (2*nf+1.gt.nflim) stop 'NF > NFLIM IN APEQSOL'
      IF (Q2MIN.GT.q20) then 
        print *, ' Q2MIN>Q20 IN APEQSOL',q2min,q20
        stop
      end if

      q2st=q20
      xbmax=0.999d0
      x1=0.2
      xlog1=log(x1)
      xlog2=log(1-x1)

      call pgridini(xbmin)

! clean the main PDF's array Y(...) and the 
! corresponding spline coefficients.
      do jj=-5,2
        do iq=1,nflim
          do i=-nsmgrid,nspgrid
            do ix=-nxmgrid,nxpgrid
              y(jj,iq,ix,i)=0.
              bcoeff(jj,iq,ix,i)=0.
              ccoeff(jj,iq,ix,i)=0.
              dcoeff(jj,iq,ix,i)=0.
            end do
          end do
        end do
      end do

      do ix=-nxmgrid,nxpgrid
        xx(ix+nxmgrid+1)=xgrid(ix)
      end do

!  Loop over the PDF schemes with the 3-flavour schemes (kpdf=0)
!  evolved in any case, and 4-(5-)flavour PDFs if NF=4(5)
      npdfs2=nf-1
      npdfs1=2
      if (.not.bmsnfopt) then 
!  evolve also LO and (N)NLO PDFs for the BMSN prescription of the VFN scheme
        npdfs2=8
      end if 
      if (vloop) then 
!  evolve also  the 3-flavour PDFs with a variable number of loops
!  in the \alpha_s evolution (kpdf=-1)
        npdfs1=1
      end if 

      do kk=npdfs1,npdfs2
!  Set the evolved PDF set 
        kpdf=kpdfseq(kk)
!  Set the input PDF set 
        kpdfi=kpdfseqi(kk)
!  Set number of fermions in the initial state equal 
        nfc=nfcpdfseq(kk)
!  Save status of the evolution order key
        kordkernels=kordkernel
!  Set the order of the matching conditions
!    nominal case
        kome=kordhq   
!    special case of the 4-, 5-flavour PDFs employed in the VFN scheme
        if (kpdf.eq.-2.or.kpdf.eq.-4) then     
          kome=0  
          kordkernel=0            !  LO evolution
        end if
        if (kpdf.eq.-3.or.kpdf.eq.-5) then    
          kome=1  
          if(bmsnnlo) then  
            kordkernel=1         !   NLO evolution
          else
            kordkernel=2         !   NNLO evolution
          end if
        end if          
! fill the grids for 3-flavour initial distributions 
      if (kpdf.eq.0) then 
        CALL APEQSGRID(Q2min,Q2max,xbmin,kpdf,q2st)
        do ix=-nxmgrid,nxpgrid
          xx(ix+nxmgrid+1)=xgrid(ix)
        end do
        DO IQ=1,7
          DO IX=-nxmgrid,nxpgrid-1
            Y(kpdf,IQ,IX,0)=XQG0(0,IQ,XGRID(IX),ix)
            fsp(ix+nxmgrid+1)=y(kpdf,iq,ix,0)
          end do
          fsp(nxpgrid+nxmgrid+1)=y(kpdf,iq,nxpgrid,0)
          call spline (nxmgrid+nxpgrid+1,xx,fsp,bs,cs,ds)
          do ix=-nxmgrid,nxpgrid
            bcoeff(kpdf,iq,ix,0)=bs(ix+nxmgrid+1)
            ccoeff(kpdf,iq,ix,0)=cs(ix+nxmgrid+1)
            dcoeff(kpdf,iq,ix,0)=ds(ix+nxmgrid+1)
          end do
        end do
      end if

      if (kpdf.eq.-1) then 
!  Find a position of m_c^2 in the Q-grid for the 3-flavour PDFs with 
!  the fixed number of flavours in the loops (kpdf=0)
        s=log(log(vfnth(4)**2/0.04)/log(q2ini(0)/0.04))
        if (s.ge.0.) then
          is=int(s/dels1(0))
        else
          is=int(s/dels2(0))-1
        end if
!  Set m_c^2 as a starting scale for the evolution of the 3-flavour PDFs 
!  with the variable number of flavours in the loops (kpdf=-1)
        q2=0.04*exp(exp(sgrid(is,0))*log(q2ini(0)/0.04))
        CALL APEQSGRID(Q2min,Q2max,xbmin,kpdf,q2)
        do ix=-nxmgrid,nxpgrid
          xx(ix+nxmgrid+1)=xgrid(ix)
        end do
        DO IQ=1,7
          DO IX=-nxmgrid,nxpgrid-1
            y(kpdf,iq,ix,0)=y(0,iq,ix,is)
            fsp(ix+nxmgrid+1)=y(kpdf,iq,ix,0)
          end do
          fsp(nxpgrid+nxmgrid+1)=y(kpdf,iq,nxpgrid,0)
          call spline (nxmgrid+nxpgrid+1,xx,fsp,bs,cs,ds)
          do ix=-nxmgrid,nxpgrid
            bcoeff(kpdf,iq,ix,0)=bs(ix+nxmgrid+1)
            ccoeff(kpdf,iq,ix,0)=cs(ix+nxmgrid+1)
            dcoeff(kpdf,iq,ix,0)=ds(ix+nxmgrid+1)
          end do
        end do
      end if

! fill the grids for 4(5)-flavour initial distribution
      if (nfc.ge.4) then 
! find a position of m_c(b)^2 in the 3(4)-flavour Q-grid 
        s=log(log(vfnth(nfc)**2/0.04)/log(q2ini(kpdfi)/0.04))
        if (s.ge.0.) then
          is=int(s/dels1(kpdfi))
        else
          is=int(s/dels2(kpdfi))-1
        end if
! generate the 4(5)-flavour distribution at the scale close to m_c(b)
        q2=0.04*exp(exp(sgrid(is,kpdfi))*log(q2ini(kpdfi)/0.04))
! re-set the boundary scale and re-generate the 4(5)-flavour Q-grid 
        CALL APEQSGRID(Q2min,Q2max,Xbmin,kpdf,q2)
        an=xqg(0,0.1d0,q2,kpdf)/4./pi
        isch0=kpdfi
        isch1=kpdf
! Generate FOPT 4-,5- flavour PDFs,
        call fillvfx(is,2*nfc,kome,isch0,isch1)
! restore alpha_s grid spoiled in fillvfx
        CALL APEQSGRID(Q2min,Q2max,Xbmin,kpdf,q2)
! and set the boundary condition for the 4(5)-flavour PDFs evolution 
        DO IQ=1,2*nfc+1
          DO IX=-nxmgrid,nxpgrid
            fsp(ix+nxmgrid+1)=y(kpdf,iq,ix,is)
            y(kpdf,iq,ix,0)=y(kpdf,iq,ix,is)
          end do
          call spline (nxmgrid+nxpgrid+1,xx,fsp,bs,cs,ds)
          do ix=-nxmgrid,nxpgrid
            bcoeff(kpdf,iq,ix,0)=bs(ix+nxmgrid+1)
            ccoeff(kpdf,iq,ix,0)=cs(ix+nxmgrid+1)
            dcoeff(kpdf,iq,ix,0)=ds(ix+nxmgrid+1)
          end do
        end do
      end if

! high-Q evolution

! first step is performed with a reduced step on Q
      khalf=0
      call dsxqgt(dert,0,dels1(kpdf)*nspgrid,kpdf)
        DO IQ=1,2*nfc+1
          DO IX=NXPGRID-1,-nxmgrid,-1
            yhalf(kpdf,IQ,IX)=Y(kpdf,IQ,IX,0)
     +                  +dert(kpdf,iq,ix)*dels1(kpdf)/2.
          end do
        end do
      do iq=1,2*nfc+1
        do ix=-nxmgrid,nxpgrid
          fsp(ix+nxmgrid+1)=yhalf(kpdf,iq,ix)
        end do
        call spline (nxmgrid+nxpgrid+1,xx,fsp,bs,cs,ds)
        do ix=-nxmgrid,nxpgrid
          bcoeffh(iq,ix)=bs(ix+nxmgrid+1)
          ccoeffh(iq,ix)=cs(ix+nxmgrid+1)
          dcoeffh(iq,ix)=ds(ix+nxmgrid+1)
        end do
      end do

      khalf=1
      call dsxqgt(dert,0,dels1(kpdf)*nspgrid,kpdf)
        DO IQ=1,2*nfc+1
          DO IX=NXPGRID-1,-nxmgrid,-1
            y(kpdf,IQ,IX,1)=Y(kpdf,IQ,IX,0)+dert(kpdf,iq,ix)*dels1(kpdf)
          end do
        end do
      do iq=1,2*nfc+1
        do ix=-nxmgrid,nxpgrid
          fsp(ix+nxmgrid+1)=y(kpdf,iq,ix,1)
        end do
        call spline (nxmgrid+nxpgrid+1,xx,fsp,bs,cs,ds)
        do ix=-nxmgrid,nxpgrid
          bcoeff(kpdf,iq,ix,1)=bs(ix+nxmgrid+1)
          ccoeff(kpdf,iq,ix,1)=cs(ix+nxmgrid+1)
          dcoeff(kpdf,iq,ix,1)=ds(ix+nxmgrid+1)
        end do
      end do

! nominal steps elsewhere 

      DO IS=1,NSPGRID-1
        khalf=0
        call dsxqgt(dert,is,dels1(kpdf)*nspgrid,kpdf)
            DO IQ=1,2*nfc+1
              DO IX=NXPGRID-1,-nxmgrid,-1
                Y(kpdf,IQ,IX,IS+1)=Y(kpdf,IQ,IX,IS-1)
     +                        +2.D0*dert(kpdf,iq,ix)*DELS1(KPDF)
              end do
            END DO
          do iq=1,2*nfc+1
            do ix=-nxmgrid,nxpgrid
              fsp(ix+nxmgrid+1)=y(kpdf,iq,ix,is+1)
            end do
            call spline (nxmgrid+nxpgrid+1,xx,fsp,bs,cs,ds)
            do ix=-nxmgrid,nxpgrid
              bcoeff(kpdf,iq,ix,is+1)=bs(ix+nxmgrid+1)
              ccoeff(kpdf,iq,ix,is+1)=cs(ix+nxmgrid+1)
              dcoeff(kpdf,iq,ix,is+1)=ds(ix+nxmgrid+1)
            end do
          end do
      end do

! low-Q evolution for the 3-flavor PDF only

      if (nfc.eq.3) then 

      IF (Q2MIN.le.q2ini(kpdf)) then 
        khalf=0
        call dsxqgt(dert,0,dels2(kpdf)*nsmgrid,kpdf)
          DO IQ=1,2*nfc+1
            DO IX=NXPGRID-1,-nxmgrid,-1
             yhalf(kpdf,IQ,IX)=Y(kpdf,IQ,IX,0)
     -                     -dert(kpdf,iq,ix)*DELS2(KPDF)/2.
            end do
          END DO
        do iq=1,2*nfc+1
          do ix=-nxmgrid,nxpgrid
            fsp(ix+nxmgrid+1)=yhalf(kpdf,iq,ix)
          end do
          call spline (nxmgrid+nxpgrid+1,xx,fsp,bs,cs,ds)
          do ix=-nxmgrid,nxpgrid
            bcoeffh(iq,ix)=bs(ix+nxmgrid+1)
            ccoeffh(iq,ix)=cs(ix+nxmgrid+1)
            dcoeffh(iq,ix)=ds(ix+nxmgrid+1)
          end do
        end do
        khalf=1
        call dsxqgt(dert,0,dels2(kpdf)*nsmgrid,kpdf)
          DO IQ=1,2*nfc+1
            DO IX=NXPGRID-1,-nxmgrid,-1
              y(kpdf,IQ,IX,-1)=y(kpdf,IQ,IX,0)
     -                         -dert(kpdf,iq,ix)*DELS2(KPDF)
            end do
          END DO
        do iq=1,2*nfc+1
          do ix=-nxmgrid,nxpgrid
            fsp(ix+nxmgrid+1)=y(kpdf,iq,ix,-1)
          end do
          call spline (nxmgrid+nxpgrid+1,xx,fsp,bs,cs,ds)
          do ix=-nxmgrid,nxpgrid
            bcoeff(kpdf,iq,ix,-1)=bs(ix+nxmgrid+1)
            ccoeff(kpdf,iq,ix,-1)=cs(ix+nxmgrid+1)
            dcoeff(kpdf,iq,ix,-1)=ds(ix+nxmgrid+1)
          end do
        end do

        DO IS=-1,-NSMGRID+1,-1
          khalf=0
          call dsxqgt(dert,is,dels2(kpdf)*nsmgrid,kpdf)
              DO IQ=1,2*nfc+1
                DO IX=NXPGRID-1,-nxmgrid,-1
                  Y(kpdf,IQ,IX,IS-1)=Y(kpdf,IQ,IX,IS+1)
     -                          -2.D0*dert(kpdf,iq,ix)*DELS2(KPDF)
                end do
              END DO
            do iq=1,2*nfc+1
              do ix=-nxmgrid,nxpgrid
                fsp(ix+nxmgrid+1)=y(kpdf,iq,ix,is-1)
              end do
              call spline (nxmgrid+nxpgrid+1,xx,fsp,bs,cs,ds)
              do ix=-nxmgrid,nxpgrid
                bcoeff(kpdf,iq,ix,is-1)=bs(ix+nxmgrid+1)
                ccoeff(kpdf,iq,ix,is-1)=cs(ix+nxmgrid+1)
                dcoeff(kpdf,iq,ix,is-1)=ds(ix+nxmgrid+1)
              end do
            end do
        END DO
      end if
      end if

!  Restore status of the evolution order key
      kordkernel=kordkernels
      end do

      khalf=0

      RETURN
      END
c--------------------
      subroutine dsxqgt(dert,is,stot,kpdf)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      INCLUDE 'APSCOM6.'
      INCLUDE 'PDFCOM.'
      include 'CONSTCOM.'
      include 'PRECCOM.'
      
      EXTERNAL QGSPLITlinI,QGSPLITlogI
      COMMON  /FORQGSPL/ XB0,Q2,XLOG,AN,AN2,an3,kn,iqn,ixn,isn,ihqn

      real*8 dert(-5:2,nflim,-nxmlim:nxplim)
     ,  ,p(nflim),pn(nflim),pnn(nflim)
      
      Q2=0.04*exp(exp(sgrid(is,kpdf))*log(q2ini(kpdf)/0.04))

      nfe=nfloops(q2,kpdf)

      BET0=11.-2./3.*nfe
      bet1=102-38/3.*nfe

      an=Y(kpdf,0,0,is)/4./pi

      alr=-log(rscale)
      an2=an**2
      an3=an**3

        DO IQ=1,2*nfc+1
          IQN=IQ
          ISN=IS
          do ix=nxpgrid-1,-nxmgrid,-1
            IXN=IX
            XB0=XGRID(IX)
            DS=0.D0
            kn=kpdf
              if (xb0.lt.x1) then 
                xc=xb0
                CALL GAUSS1(QGSPLITlogI
     ,                      ,log(xb0),log(xbmax),kprecdq,DS,e)
              else 
                CALL GAUSS1(QGSPLITlinI,xb0,xbmax,kprecdq,DS,e)
              end if
              call QGSPLIT0(IQ,p,XB0)
              if (kordkernel.ge.1) call QGSPLITnlo0(IQ,pn,XB0)
              if (kordkernel.ge.2) 
     -                 call QGSPLITnnlo0(IQ,pnn,XB0,kernelmod)
            DO JQ=1,2*nfc+1
              scon=an*2*p(jq)
              if (kordkernel.ge.1) 
     -             scon=scon+2*an2*(2*pn(jq)-p(jq)*bet0*alr)
              if (kordkernel.ge.2) 
     -             scon=scon+an3*(pnn(jq)-8*bet0*alr*pn(jq)
     -                      -2*(bet1*alr-bet0**2*alr**2)*p(jq))
              if (khalf.eq.0) then              
                DS=DS+scon*Y(kpdf,IQ,IX,IS)
              else
                DS=DS+scon*yhalf(kpdf,IQ,IX)
              end if
            end do
            dert(kpdf,iq,ix)=ds*log(q2/0.04)
          end do
        end do

      return
      end
C------------------
      DOUBLE PRECISION FUNCTION QGSPLITlogI(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON  /FORQGSPL/ XB0,Q2,XLOG,AN,AN2,an3,kn,iqn,ixn,isn,ihqn

      xb=exp(x)
      QGSPLITlogI=QGSPLITlinI(xb)*(xb)

      return 
      end
C------------------
      DOUBLE PRECISION FUNCTION QGSPLITlinI(xb)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'APSCOM6.'
      INCLUDE 'PDFCOM.'
      include 'CONSTCOM.'

      COMMON  /FORQGSPL/ XB0,Q2,XLOG,AN,AN2,an3,kn,iqn,ixn,isn,ihqn
      real*8 px(nflim),pnx(nflim),pnnx(nflim),q(nflim)
      real*8 pp(0:2,nflim)

      alr=-log(rscale)
      BET0=11.-2./3.*nfe
      bet1=102.-38./3.*nfe

      QGS=0.D0

      z=XB0/xb
      call pintx(0,iqn,z,pp)

      call QGSPLITX(IQN,px,XB0/XB)
      if (kordkernel.ge.1) then
        call QGSPLITnloX(IQN,pnx,XB0/XB)
      end if
      if (kordkernel.ge.2) then
        call QGSPLITnnloX(IQN,pnnx,XB0/XB,kernelmod)
      end if
      DO JQ=1,2*nfc+1
        pdfc=XQGX1(kn,jq,XB,ISn)
        scon=2*an*pp(0,jq)
        if (kordkernel.ge.1) scon=scon+2*an2*(2*pp(1,jq)
     -                -pp(0,jq)*bet0*alr)
        if (kordkernel.ge.2) scon=scon+an3*(pp(2,jq)
     -                -8*bet0*alr*pp(1,jq)
     -                -2*(bet1*alr-bet0**2*alr**2)*pp(0,jq))
        QGS=QGS+scon*pdfc
        sconx=2*an*px(jq)
        if (kordkernel.ge.1) 
     -         sconx=sconx+2*an2*(2*pnx(jq)-px(jq)*bet0*alr)
        if (kordkernel.ge.2) 
     -         sconx=sconx+an3*(pnnx(jq)-8*bet0*alr*pnx(jq)
     -                    -2*(bet1*alr-bet0**2*alr**2)*px(jq))

        if (khalf.eq.0) then
          QGS=QGS+sconx*(Y(kn,JQ,IXN,ISN)-pdfc)
        else
          QGS=QGS+sconx*(yhalf(kn,JQ,IXN)-pdfc)
        end if
      end do

      QGSPLITlinI=QGS*XB0/XB**2

      RETURN
      END
C-----------------
      subroutine APEQSGRID(q2low,q2high,xblow,kpdf,q2st)
      implicit double precision (a-h,o-z)

      include 'APSCOM6.'
      include 'PDFCOM.'
      include 'CONSTCOM.'

      real*8 p(nflim)
      real*8 fsp(nxtot),bs(nxtot),cs(nxtot),ds(nxtot),xx(nxtot)
!  initial scale for the QCD evolution
      q2ini(kpdf)=q2st
! set up the steps for the x- and Q-grid

! log-log transform for Q 
      smax=log(log(q2high/0.04)/log(q2ini(kpdf)/0.04))
      smin=-log(log(q2low/0.04)/log(q2ini(kpdf)/0.04))
      dels1(kpdf)=smax/nspgrid
      dels2(kpdf)=smin/nsmgrid
! log(x) transform for the low-x region x 
      xlog0=log(xblow)
      delx2=-(xlog0-xlog1)/nxmgrid
! polynomial, (1-x)**2, transform for the high-x region 
      delx1=(1-x1)**2/nxpgrid

! fill the x- and Q-grids 

! small-x region 
      do i=-nxmgrid,0
        xgrid(i)=exp(xlog1+delx2*i)
      end do

! large-x region 
      do i=0,nxpgrid-1
        xgrid(i)=1-((1-x1)**2-delx1*i)**(1./2.)
      end do
      xgrid(nxpgrid)=1.

! small-Q region 
      do i=-nsmgrid,0
        sgrid(i,kpdf)=dels2(kpdf)*i
      end do

! large-Q region 
      do i=0,nspgrid
        sgrid(i,kpdf)=dels1(kpdf)*i
      end do

!  set up the grid for spline calculation
      do ix=-nxmgrid,nxpgrid
        xx(ix+nxmgrid+1)=xgrid(ix)
      end do

!  fill the values of \alpha_s as the PDF with iq=0 
      do is=-nsmgrid,nspgrid
        Q2=0.04*exp(exp(sgrid(is,kpdf))*log(q2ini(kpdf)/0.04))

        if (kpdf.eq.-1) then 
          an=alphas_vfn(rscale*Q2)
        end if
        if (kpdf.eq.0) then 
          an=alphas_ffn(rscale*Q2)
        end if
        if (kpdf.eq.1.or.kpdf.eq.-2.or.kpdf.eq.-3) then 
          an=alphas_ffn4(rscale*Q2)
        end if
        if (kpdf.eq.2.or.kpdf.eq.-4.or.kpdf.eq.-5) then 
          an=alphas_ffn5(rscale*Q2)
        end if

!  calculate the spline coefficients for PDFs and \alpha_s
        do ix=-nxmgrid,nxpgrid
          y(kpdf,0,ix,is)=an
          fsp(ix+nxmgrid+1)=y(kpdf,0,ix,is)
        end do
        call spline (nxmgrid+nxpgrid+1,xx,fsp,bs,cs,ds)
        do ix=-nxmgrid,nxpgrid
          bcoeff(kpdf,0,ix,is)=bs(ix+nxmgrid+1)
          ccoeff(kpdf,0,ix,is)=cs(ix+nxmgrid+1)
          dcoeff(kpdf,0,ix,is)=ds(ix+nxmgrid+1)
        end do
      end do

      RETURN
      END
C------------------
      real*8 function XQGX1(k,iq,XB,IS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'APSCOM6.'

      data xbsave /-1d0/

      save xbsave,ix,aa
      
      if (xb.ne.xbsave) then 
        do ix=-nxmgrid,nxpgrid-1
          if (xb.lt.xgrid(ix+1)) goto 300
        end do
 300    xbsave=xb
      end if

      aa=xb-xgrid(ix)

      if (khalf.eq.0) then 
        xqgx1=y(k,iq,ix,is)+aa*bcoeff(k,iq,ix,is)
     +     +aa**2*ccoeff(k,iq,ix,is)+aa**3*dcoeff(k,iq,ix,is)
      else 
        xqgx1=yhalf(k,iq,ix)+aa*bcoeffh(iq,ix)
     +     +aa**2*ccoeffh(iq,ix)+aa**3*dcoeffh(iq,ix)
      end if

      RETURN
      END
