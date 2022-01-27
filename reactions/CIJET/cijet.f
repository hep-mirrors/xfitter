c CIJET1.1 interface made for xFitter, Jun Gao, 2018.09.10
c Updated for 2021 xFitter release, Toni Makela, 2021.11.19

c performing grid reading and initialization 
      subroutine cijetinit(fname, mufip, murip, oqcd, fut, estat)

! new line
      implicit none

! maximum number of data sets and points per set
      include 'cijetInclude.h'

! input scale choices, CI shape choice and couplings, only color-singlet
      real(8) :: mufip, murip, fut
      real(8) :: cpl(3)=0.d0 ! color singlet couplings, LL, LR, RR
      real(8) :: ash(10) ! individual coupling combinations 

! input QCD order
      integer :: oqcd ! 0 for LO, 1 for NLO

      integer, save :: befirst=1

! store non-zero error status identifying number if need be
      integer :: estat

! beta function needed for scale evaluation
      real(8), parameter :: beta0=(11*3.d0-2*5.d0)/12.d0/3.1415926d0
      real(8), parameter :: scaf(9)=(/0.5,0.5,0.5,1.,1.,1.,2.,2.,2./)
      real(8), parameter :: scar(9)=(/0.5,1.,2.,0.5,1.,2.,0.5,1.,2./)

! individual set workspace
      integer :: nt, cnt, pointor(3*nbk*mpint+1)
      integer :: i, j, k, lhc, ang, scalescheme, Nx, Nq
      real(8) :: chid, chiu, masscutd, masscutu, mu0, Rcone, fct
      integer :: sid, bid1, bid2, nstart, nend, nca
      real(8) :: sqrts, slog, wgt(5), as, sc, tm1, tm2, tm3, tm4
      character*100 :: xx

! shared common block for latter convolution
      common /cigrid/ Nscale, xgd, qgd, smu, ciwgt,
     -       srof, cset, pset, nxgd, nqgd, sig, ffnm

! initialization of commons
      if(befirst==1) then
      sid=0;cset=0;pset=0;srof=1.d0
      ciwgt=0.d0;Nscale=5000.d0
      xgd=0.d0;qgd=0.d0;smu=0.d0
      befirst=0
      endif

      write(*,*) "*****************************************"
      write(*,*) "*  Reading CIJET1.1 table with parameters:"
      write(*,*) "*  mufip, murip, oqcd, fut, fname="
      write(*,*) "*", mufip, murip, oqcd, fut, trim(fname)
      write(*,*) "*  Jun Gao, 2018.09.10"
      write(*,*) "*****************************************"

! end of initialization of commons

! init couplings apart from 1/Lambda^2
      ash(1) =1.d0 !cpl(1)+cpl(3)
      ash(2) =1.d0 !cpl(1)^2+cpl(3)^2
      ash(3) =1.d0 !cpl(2)^2
      ash(4) =1.d0 !cpl(1)+cpl(3)
      ash(5) =1.d0 !cpl(2)
      ash(6) =1.d0 !cpl(1)+cpl(3)
      ash(7) =1.d0 !cpl(2)
      ash(8) =1.d0 !cpl(1)^2+cpl(3)^2
      ash(9) =1.d0 !cpl(2)^2
      ash(10)=1.d0 !cpl(1)^2+cpl(3)^2

! first load to find index of components
      nt=1
      cnt=0
      pointor=0
      open(unit=155, file=fname, status="old", action="read")
      do while(0.lt.1)
      read(155, *, end=200) xx
      cnt=cnt+1
      if(xx.eq."JCC") then
      pointor(nt)=cnt
      nt=nt+1
      endif
      enddo
 200  pointor(nt)=cnt+1

      if(mod(nt-1, 3*nbk).ne.0) then
      print *, "wrong grid file read!"
      stop
      endif
      close(155)

! determine the correct scale choice
      do i=1, 9
      if(mufip==scaf(i)) then
      sid=1+(i-1)/3
      endif
      enddo
      if(sid==0) then
      print *, "scale not available!";stop
      endif
      slog=dlog(murip)     
 
! update record
      cset=cset+1
      if((cset.gt.mset)) then
      estat=21111801
      print *, "F: cset in xfitter.f exceeds mset in xfitterInclude.h."
      print *, "   Increase mset and recompile.",fname;
      return
      endif
      pset(cset)=(nt-1)/(3*nbk)
      ffnm(cset)=fname
      sig(cset)=1
      srof(cset)=murip/mufip
! end first

      open(unit=156, file=fname, status="old", action="read")
! begin of the main loop on the kinematic bins
      do nca=1, pset(cset)

      read(156, *) lhc, xx, sc, as
      read(156, *) sqrts, ang, scalescheme, Rcone
      read(156, *) chid, chiu, masscutd, masscutu, mu0
      read(156, *) Nx, Nq
      
! transfer bin-xsec to double differential
      fct=fut/(chiu-chid)/(masscutu-masscutd)

! storing needed grid information
      if((Nx.gt.xnd).or.(Nq.gt.qnd)) then
      print *, "X or Q grid nodes wrong!"; stop
      endif

      nxgd(cset,nca)=Nx
      nqgd(cset,nca)=Nq

      do i=1, Nx
      read(156, *) xgd(cset,nca,i)
      enddo 
      do i=1, Nq
      read(156, *) qgd(cset,nca,i) 
      enddo
      smu(cset,nca)=mu0
      if(lhc==0) sig(cset)=-1 


! must include all color-singlet case
      if(xx.eq."fitll") then
      print *, "grid file should be for CS at least!"
      stop
      endif

! loop over scale choices
      do j=1, 3

! loop over 42 coefficients
      do k=1, nbk

! skip to the correct scale choice
      bid1=nbk*(3*(nca-1)+j-1)+k
      bid2=nbk*(3*(nca-1)+j-1)+k+1

      nstart=pointor(bid1)
      nend=pointor(bid2)-1

! account for the header lines
      if((j==3).and.(k==nbk).and.(nca.ne.pset(cset))) then
      nend=nend-(pointor(1)-1)
      endif

! drop JCC
      read(156, *) xx

! block by block
      do i=nstart+1, nend
! skipping
      if(j.ne.sid) then
      read(156, *) xx
      endif
! working
      if(j.eq.sid) then
      read(156, *) nx, nq, wgt(1:5)
! tranfer to double differentail
      wgt=wgt*fct
! adding all components with correct couplings
      select case (k)
! O(0)
! inteference
      case (1)
!      ciwgt(1,nx,nq,1:5,cset,nca)=ciwgt(1,nx,nq,1:5,cset,nca)
!     -   +(cpl(1)+cpl(3))*wgt(1:5)
      ciwgt(1,nx,nq,1:5,cset,nca)=ciwgt(1,nx,nq,1:5,cset,nca)
     -   +(ash(1))*wgt(1:5)
! square
      case (5)
!      ciwgt(2,nx,nq,1:5,cset,nca)=ciwgt(2,nx,nq,1:5,cset,nca)
!     -   +(cpl(1)**2+cpl(3)**2)*wgt(1:5)
      ciwgt(2,nx,nq,1:5,cset,nca)=ciwgt(2,nx,nq,1:5,cset,nca)
     -   +(ash(2))*wgt(1:5)
      case (8)
!      ciwgt(2,nx,nq,1:5,cset,nca)=ciwgt(2,nx,nq,1:5,cset,nca)
!     -   +(cpl(2)**2)*wgt(1:5)
      ciwgt(3,nx,nq,1:5,cset,nca)=ciwgt(3,nx,nq,1:5,cset,nca)
     -   +(ash(3))*wgt(1:5)
! O(1)
! inteference, const
      case (11)
      if(oqcd.ne.0) then
!      ciwgt(3,nx,nq,1:5,cset,nca)=ciwgt(3,nx,nq,1:5,cset,nca)
!     -   +(cpl(1)+cpl(3))*wgt(1:5)
      ciwgt(4,nx,nq,1:5,cset,nca)=ciwgt(4,nx,nq,1:5,cset,nca)
     -   +(ash(4))*wgt(1:5)
      endif
      case (15)
      if(oqcd.ne.0) then
!      ciwgt(3,nx,nq,1:5,cset,nca)=ciwgt(3,nx,nq,1:5,cset,nca)
!     -   +(cpl(2))*wgt(1:5)
      ciwgt(5,nx,nq,1:5,cset,nca)=ciwgt(5,nx,nq,1:5,cset,nca)
     -   +(ash(5))*wgt(1:5)
      endif
! inteference, logs
      case (12)
      if(oqcd.ne.0) then
!      ciwgt(4,nx,nq,1:5,cset,nca)=ciwgt(4,nx,nq,1:5,cset,nca)
!     -   +(cpl(1)+cpl(3))*wgt(1:5)
      ciwgt(6,nx,nq,1:5,cset,nca)=ciwgt(6,nx,nq,1:5,cset,nca)
     -   +(ash(6))*wgt(1:5)
      endif
      case (16)
      if(oqcd.ne.0) then
!      ciwgt(4,nx,nq,1:5,cset,nca)=ciwgt(4,nx,nq,1:5,cset,nca)
!     -   +(cpl(2))*wgt(1:5)
      ciwgt(7,nx,nq,1:5,cset,nca)=ciwgt(7,nx,nq,1:5,cset,nca)
     -   +(ash(7))*wgt(1:5)
      endif
! square, const
      case (19)
      if(oqcd.ne.0) then
!      ciwgt(5,nx,nq,1:5,cset,nca)=ciwgt(5,nx,nq,1:5,cset,nca)
!     -   +(cpl(1)**2+cpl(3)**2)*wgt(1:5)
      ciwgt(8,nx,nq,1:5,cset,nca)=ciwgt(8,nx,nq,1:5,cset,nca)
     -   +(ash(8))*wgt(1:5)
      endif
      case (25)
      if(oqcd.ne.0) then
!      ciwgt(5,nx,nq,1:5,cset,nca)=ciwgt(5,nx,nq,1:5,cset,nca)
!     -   +(cpl(2)**2)*wgt(1:5)
      ciwgt(9,nx,nq,1:5,cset,nca)=ciwgt(9,nx,nq,1:5,cset,nca)
     -   +(ash(9))*wgt(1:5)
      endif
! square, logs
      case (20)
      if(oqcd.ne.0) then
!      ciwgt(6,nx,nq,1:5,cset,nca)=ciwgt(6,nx,nq,1:5,cset,nca)
!     -   +(cpl(1)**2+cpl(3)**2)*wgt(1:5)
      ciwgt(10,nx,nq,1:5,cset,nca)=ciwgt(10,nx,nq,1:5,cset,nca)
     -   +(ash(10))*wgt(1:5)
      endif
      case default
      read(156, *) xx
      end select

      endif

      enddo
! end of block

      enddo
! end of coefficients

! adding back mur dependence
      if((j.eq.sid).and.(oqcd.ne.0)) then
      ciwgt(4,1:xnd,1:qnd*xnd,1:5,cset,nca)=
     -    ciwgt(4,1:xnd,1:qnd*xnd,1:5,cset,nca)
     -  +2.d0*slog*beta0*ciwgt(1,1:xnd,1:qnd*xnd,1:5,cset,nca)
      endif

      enddo
! end of scale choices

      enddo
! end of kinematic loop

      close(156)

      return

      end subroutine cijetinit

c perform grid convolution and calculation of xsecs
      subroutine cijetxsec(fname, invlamsq, cpl, npt, pres)

! new line
      implicit none

! maximum number of data sets and points per set
      include 'cijetInclude.h'

! input 1/Lambda^2 [in TeV^-2]
      real(8) :: invlamsq, res(mpint), pres(mpint)
! sign (+ for destructive, -for con.)
      real(8) :: cpl(3) ! color singlet couplings, LL, LR, RR
      real(8) :: ash(10) ! individual coupling combinations 
! number of data points and additional multi-factor
      integer :: npt

! workspace
      real(8) :: XPDF(-6:6, xnd, qnd), xi(xnd), qi(qnd)
      real(8) :: pdfpg(xnd*qnd, xnd, 5)
      integer :: pid, nca, Nx, Nq, sgg
      integer :: i, j, k, ii, jj, kk, iqd
      real(8) :: norm(qnd), tmp, cialphas, tm1, tm2, tm3, tm4
      external :: cialphas, cievolPDF
      

! shared common block for latter convolution
      common /cigrid/ Nscale, xgd, qgd, smu, ciwgt,
     -       srof, cset, pset, nxgd, nqgd, sig, ffnm

! find the correct ID
      pid=0; res=0.d0
      do i=1, cset
      if(fname==ffnm(i)) pid=i
      enddo
      if(pid==0) then
      print *, "CI grid not loaded!",fname; stop
      endif

! find the correct couplings
      ash(1) =cpl(1)+cpl(3)
      ash(2) =cpl(1)**2+cpl(3)**2
      ash(3) =cpl(2)**2
      ash(4) =cpl(1)+cpl(3)
      ash(5) =cpl(2)
      ash(6) =cpl(1)+cpl(3)
      ash(7) =cpl(2)
      ash(8) =cpl(1)**2+cpl(3)**2
      ash(9) =cpl(2)**2
      ash(10)=cpl(1)**2+cpl(3)**2

! loop over the kinematic bins
      do nca=1, pset(pid)

      XPDF=0.d0
! passing grid configuration
      Nx=nxgd(pid,nca); Nq=nqgd(pid,nca); sgg=sig(pid)
      xi=xgd(pid, nca, 1:xnd); qi=qgd(pid, nca, 1:qnd)
 
! calculate the PDFs in each grid points
      do i=1, Nx
      do j=1, Nq
      call cievolPDF(xi(i), qi(j), XPDF(-6, i, j))
      do k=-6, 6
! interpolation kernal is x**2*pdf
      XPDF(k, i, j)=XPDF(k, i, j)*xi(i)
      enddo
      enddo
      enddo
! end of PDF calls

! calculation of the PDF grid quantity
      pdfpg=0.d0

      do i=1, Nx
      do ii=1, Nx
      do j=1, Nq
      do k=1, 4
      do kk=k+1, 5
      pdfpg(j+Nq*(ii-1), i, 1)=
     -   pdfpg(j+Nq*(ii-1), i, 1)
     -   +XPDF(k, i, j)*XPDF(sgg*kk, ii, j)
     -   +XPDF(kk, i, j)*XPDF(sgg*k, ii, j)
     -   +XPDF(-k, i, j)*XPDF(-sgg*kk, ii, j)
     -   +XPDF(-kk, i, j)*XPDF(-sgg*k, ii, j)
      pdfpg(j+Nq*(ii-1), i, 3)=
     -   pdfpg(j+Nq*(ii-1), i, 3)
     -   +XPDF(k, i, j)*XPDF(-sgg*kk, ii, j)
     -   +XPDF(kk, i, j)*XPDF(-sgg*k, ii, j)
     -   +XPDF(sgg*k, i, j)*XPDF(-kk, ii, j)
     -   +XPDF(sgg*kk, i, j)*XPDF(-k, ii, j)
      enddo
      enddo
      do k=1, 5
      pdfpg(j+Nq*(ii-1), i, 2)=
     -   pdfpg(j+Nq*(ii-1), i, 2)
     -   +XPDF(k, i, j)*XPDF(sgg*k, ii, j)
     -   +XPDF(-k, i, j)*XPDF(-sgg*k, ii, j)
      pdfpg(j+Nq*(ii-1), i, 4)=
     -   pdfpg(j+Nq*(ii-1), i, 4)
     -   +XPDF(k, i, j)*XPDF(-sgg*k, ii, j)
     -   +XPDF(sgg*k, i, j)*XPDF(-k, ii, j)
      pdfpg(j+Nq*(ii-1), i, 5)=
     -   pdfpg(j+Nq*(ii-1), i, 5)
     -   +XPDF(k, i, j)*XPDF(0, ii, j)
     -   +XPDF(sgg*k, i, j)*XPDF(0, ii, j)
     -   +XPDF(-k, i, j)*XPDF(0, ii, j)
     -   +XPDF(-sgg*k, i, j)*XPDF(0, ii, j)
      enddo
      enddo
      enddo
      enddo
! end PDF grid

! calculate alphas grid quantity
      do i=1, Nq
      norm(i)=cialphas(srof(pid)*qi(i)) 
      enddo
! end alphas grid

! now ready for cross sections
! inteference  
      tmp=0.d0
      do k=1, 5
      do ii=1, Nx
      do j=1, Nq
      do i=1, Nx
      iqd=j+(ii-1)*Nq
! O(0)
      tmp=tmp+pdfpg(iqd,i,k)*(ash(1)*ciwgt(1,i,iqd,k,pid,nca))
     -    *norm(j)
! O(1) const
      tmp=tmp+pdfpg(iqd,i,k)*(ash(4)*ciwgt(4,i,iqd,k,pid,nca)+
     -    ash(5)*ciwgt(5,i,iqd,k,pid,nca))
     -    *norm(j)**2
! O(1) logs
      tmp=tmp+pdfpg(iqd,i,k)*(ash(6)*ciwgt(6,i,iqd,k,pid,nca)+
     -    ash(7)*ciwgt(7,i,iqd,k,pid,nca))
     -    *(-0.5d0*dlog(abs(invlamsq)*(smu(pid,nca)/1.0d3)**2))
     -    *norm(j)**2
      enddo
      enddo
      enddo
      enddo
      res(nca)=res(nca)+tmp*((Nscale(pid)/1.d3)**2)*invlamsq    
! square 
      tmp=0.d0
      do k=1, 5
      do ii=1, Nx
      do j=1, Nq
      do i=1, Nx
      iqd=j+(ii-1)*Nq
! O(0)
      tmp=tmp+pdfpg(iqd,i,k)*(ash(2)*ciwgt(2,i,iqd,k,pid,nca)+
     -    ash(3)*ciwgt(3,i,iqd,k,pid,nca))
! O(1) const
      tmp=tmp+pdfpg(iqd,i,k)*(ash(8)*ciwgt(8,i,iqd,k,pid,nca)+
     -    ash(9)*ciwgt(9,i,iqd,k,pid,nca))
     -    *norm(j)
! O(1) logs
      tmp=tmp+pdfpg(iqd,i,k)*(ash(10)*ciwgt(10,i,iqd,k,pid,nca))
     -    *(-0.5d0*dlog(abs(invlamsq)*(smu(pid,nca)/1.0d3)**2))
     -    *norm(j)
      enddo
      enddo
      enddo
      enddo
      res(nca)=res(nca)+tmp*(((Nscale(pid)/1.d3)**2)*invlamsq)**2    
! end of cross sections

      enddo
! end of loop on kinematic bins

      npt=pset(pid); pres=0.d0
      pres(1:pset(pid))=res(1:pset(pid))

      return

c end of PDF calculation
      end subroutine cijetxsec


!!!!!!Routines for alphas and PDFs

c link to PDF routines 
      subroutine cievolPDF(x, q, res)

! new line
      implicit none
      real(8) :: x, q, res(-6:6)
      real(8) :: xf(-6:6)

      res=0.d0
! e.g.,
      call pdf_xfxq_wrapper(x, Q, xf)
! test with LHAPDF     call evolvePDF(x, Q, xf)

      res(-6:6)=xf
      return

      end subroutine cievolPDF

c link to PDF routines 
      real(8) function cialphas(q)

! new line
      implicit none
      real(8) :: q
      real(8), external :: alphas_wrapper
! test with LHAPDF      real(8), external :: alphasPDF 

      cialphas=0.118d0
! e.g.,
      cialphas=alphas_wrapper(Q)
! test with LHAPDF      cialphas=alphasPDF(Q)
      return

      end function cialphas
