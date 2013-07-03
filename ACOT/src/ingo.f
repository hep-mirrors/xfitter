c ... v10: combine zm1 and zmchi functions: FIO 
c ... 
c ... 
c ... sfs9.f JUST RELABEL SFS9 TO DISTINGUISH FROM 7: 
C ...    includes updated chi definition
c ... 30.07.2011


c ... sfs7.f TS/IS
c ... 26.07.2011

*************************************************************************
      subroutine zmCHI(nCHI,x,Q2,muf2,mur2,iord,iboson,isf,
     >   ftot,fc,fb,ft,fij)
      implicit none
      double precision x,Q2, nCHI
      double precision muf2,mur2 ! factorization/renormalization scale
      integer iord ! perturbative order: including terms up to O(alpha_s^iord)
      integer iboson ! chose exchange boson: 0: full NC ew, 1:gamma gamma
      integer isf ! choose structure function: 0: FL; 1,2,3: F_1,2,3
      double precision ftot
      double precision fc
      double precision fb
      double precision ft
      double precision fij(0:6,1:6)

      double precision calFij
      integer i,j

c ... init
      do i = 0,6
         do j = 1,6 
            fij(i,j) = 0d0            
         enddo   
      enddo   
      
      ftot = 0d0
      fc = 0d0
      fb = 0d0
c .......................................................................

      do i = 0,6
         do j = 1,5 ! no top in final state here
            fij(i,j) = x*calFij(nCHI,x,Q2,muf2,mur2,i,j,iord,iboson,isf)
            ftot = ftot + fij(i,j)
         enddo   
      enddo   

c .......................................................................

      fc = fij(0,4)+fij(1,4)+fij(2,4)+fij(3,4)+fij(4,4)
      fb = fij(0,5)+fij(1,5)+fij(2,5)+fij(3,5)+fij(4,5)+fij(5,5)
      ft = fij(0,6)+fij(1,6)+fij(2,6)+fij(3,6)+fij(4,6)+fij(5,6) 
     >    +fij(6,6)

      end
      
*************************************************************************

      double precision function calFg(x,Q2,nf,iord,isf)
      implicit none
      double precision x,Q2
      integer nf,iord,isf  

      double precision xmin,xmax,eps
      double precision dinteg

      double precision calFgreg,calFgplus,calFgdelta
      external hcalFgreg ! Integrands

      double precision var(1:10)
      integer ivar(1:10)
      common / tmp_var / var,ivar

      calFg = 0d0

      var(1) = x
      var(2) = Q2
      ivar(1) = nf
      ivar(2) = iord
      ivar(3) = isf

c ... accuracy of integration
      eps = 1d-3     
c ... integration bounds
      xmin = x
c      xmax = 1d0
      xmax = 0.999999999d0
      
c ... regular part      
      calFgreg = dinteg(hcalFgreg,xmin,xmax,eps)

c ... plus distributions
      calFgplus = 0d0 ! no such contribution for gluon-initiated subprocess

c ... local part prop. to delta distribution
      calFgdelta = 0d0 ! no such contribution for gluon-initiated subprocess
      
      calFg = calFgreg + calFgplus + calFgdelta

      end


*************************************************************************
*************************************************************************

      double precision function calFij(nchi,x,Q2,muf2,mur2,
     >   i,j,iord,iboson,isf)
      implicit none
      double precision x,Q2
      double precision muf2,mur2 ! factorization/renormalization scale
      integer i ! which initial parton: i=0: gluon, i=1: u+ubar etc
      integer j ! heaviest final state quark different from i
      double precision mi ! mass of initial parton
      double precision mj ! mass of final quark j
      double precision mips ! phase space mass of i 
      double precision mjps ! phase space mass of j
      integer iord ! perturbative order: including terms up to O(alpha_s^iord)
      integer iboson ! chose exchange boson: 0: full NC ew, 1:gamma gamma
      integer isf ! choose structure function: 0: FL; 1,2,3: F_1,2,3

      double precision chi,nCHI
      double precision aplus,aminus,abarplus,abarminus
      double precision calFg, calFns, calFpsij
      double precision a2j,a2jm1
      integer k
      double precision mass, W, W2

      double precision massps(0:6)  !*** NOTE: include gluon mass i=0
c      data massps / 0d0,0d0,0d0,0d0,1.3d0,4.5d0,175.0d0 /
      data massps / 0d0,0d0,0d0,0d0,1.3d0,4.5d0,0.0d0 /

      double precision xmc,xmb,HMASS
      common /fred/xmc,xmb,HMASS
c 
c     =======================================================================
      massps(4)=xmc
      massps(5)=xmb
c     =======================================================================
      calFij = 0.0d0

      if ((i .gt. 6) .or. (i . lt. 0)) stop 'calFij: i out of range'   !*** Note, includes gluon range
      if ((j .gt. 6) .or. (j . lt. 1)) stop 'calFij: j out of range'

      if (muf2 .ne. Q2  ) stop ' (muf2 .ne. Q2 )  not yet implemented'
      if (muf2 .ne. mur2) stop ' (muf2 .ne. mur2) not yet implemented'

c .......................................................................
      mass = max(massps(i),massps(j))  !*** pick largest mass
      chi = x * (1.0d0 +  (nchi*mass)**2 / Q2)  !*** Note: nchi=0 is massless case
      if (chi .ge. 1d0) return  !*** kinematics out of range: x>1 effective
c .......................................................................
c     CHECK IF W IS ABOVE THRESHOLD:
      W  = SQRT (-Q2 + 1. + Q2 / X)
      W2=Q2*(1./x-1.)+hmass**2
      W=Sqrt(W2)
      If(w.le.mass) then
c         write(6,*) ' W below mass threshold: (W,m)= ',W,mass
c         return
      endif
c .......................................................................
c     =======================================================================
      if (i .eq. 0) then ! gluon initiated contribution

            a2j =0d0
            do k = 1,j
               call get_couplingsINGO(Q2,k,iboson,aplus,aminus,abarplus,
     >                        abarminus)
               a2j = a2j + aplus/j ! average sum of couplings^2 
            enddo

            a2jm1 = 0d0
            do k = 1,j-1
               call get_couplingsINGO(Q2,k,iboson,aplus,aminus,abarplus,
     >                        abarminus)
               a2jm1 = a2jm1 + aplus/(j-1) ! average sum of couplings^2 
            enddo

            calFij = a2j*calFg(chi,Q2,j,iord,isf)
     >             - a2jm1*calFg(chi,Q2,j-1,iord,isf)

c     =======================================================================
      else  ! quark initiated  

            call get_couplingsINGO(Q2,i,iboson,aplus,aminus,abarplus,
     >                        abarminus)         
         if (i .eq. j) then
            calFij = aplus * calFns(chi,Q2,i,0,iord,isf)
         endif
         
         calFij = calFij + aplus * ( calFns(chi,Q2,i,j,iord,isf) -
     >                               calFns(chi,Q2,i,j-1,iord,isf) )

            a2j =0d0
            do k = 1,j
               call get_couplingsINGO(Q2,k,iboson,aplus,aminus,abarplus,
     >                        abarminus)
               a2j = a2j + aplus/j ! average sum of couplings^2 
            enddo

            a2jm1 = 0d0
            do k = 1,j-1
               call get_couplingsINGO(Q2,k,iboson,aplus,aminus,abarplus,
     >                        abarminus)
               a2jm1 = a2jm1 + aplus/(j-1) ! average sum of couplings^2 
            enddo
            
        calFij = calFij + a2j   * calFpsij(chi,Q2,i,j,iord,isf)
     >                  - a2jm1 * calFpsij(chi,Q2,i,j-1,iord,isf)  

      end if  

c     =======================================================================
c     =======================================================================
c     =======================================================================


      end

c .......................................................................

      double precision function calFns(x,Q2,i,nf,iord,isf)
      implicit none
      double precision x,Q2
      integer i,nf,iord,isf  

      double precision xmin,xmax,eps
      double precision dinteg

      double precision calFnsreg,calFnsplus,calFnsdelta
      external hcalFnsreg,hcalFnsplus! Integrands
      double precision calFnsdel

      double precision var(1:10)
      integer ivar(1:10)
      common / tmp_var / var,ivar

      calFns = 0d0

      var(1) = x
      var(2) = Q2
      ivar(1) = nf
      ivar(2) = iord
      ivar(3) = i
      ivar(4) = isf

c ... accuracy of integration
      eps = 1d-3     
c ... integration bounds
      xmin = x
c      xmax = 1d0
      xmax = 0.99999999999d0
      
c ... regular part      
      calFnsreg = dinteg(hcalFnsreg,xmin,xmax,eps)

c ... plus distributions
      calFnsplus = dinteg(hcalFnsplus,xmin,xmax,eps)
c      calFnsplus = 0d0 ! zero for F_L

c ... local part prop. to delta distribution
      calFnsdelta = calFnsdel(x,Q2,i,nf,iord,isf)

      calFns = calFnsreg + calFnsplus + calFnsdelta

      end

c .......................................................................

      double precision function calFnsdel(x,Q2,i,nf,iord,isf)
      implicit none
      double precision x,Q2
      integer i,nf,iord,isf  

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      double precision pdfix      
      double precision Cnsdel,CLnsdelta,C2nsdelta

      double precision pdf(-6:6)

      calFnsdel = 0d0

c ... pdfs at x
      call mypdf(x,Q2,pdf,IsetPDF,iord,A1,Z1)
      pdfix = pdf(i) + pdf(-i)

c ... Wilson coefficient 
      if (isf .eq. 0) then ! C_L
         Cnsdel = CLnsdelta(x,Q2,nf,iord) ! doesn't depend on x
      elseif (isf .eq. 2) then ! C_2
         Cnsdel = C2nsdelta(x,Q2,nf,iord) ! doesn't depend on x
      else 
         stop 'isf value not provided'   
      end if   

      calFnsdel = pdfix * Cnsdel

      end

*************************************************************************  

c ... purely singlet contribution to structure functions
      double precision function calFpsij(x,Q2,i,nf,iord,isf)
      implicit none
      double precision x,Q2
      integer i,nf,iord,isf  

      double precision xmin,xmax,eps
      double precision dinteg

      double precision calFpsijreg,calFpsijplus,calFpsijdelta
      external hcalFpsijreg ! Integrands

      double precision var(1:10)
      integer ivar(1:10)
      common / tmp_var / var,ivar

      double precision calFpsijdel

      calFpsij = 0d0

      var(1) = x
      var(2) = Q2
      ivar(1) = nf
      ivar(2) = iord
      ivar(3) = isf
      ivar(4) = i

c ... accuracy of integration
      eps = 1d-3     
c ... integration bounds
      xmin = x
      xmax = 1d0

c ... regular part      
      calFpsijreg = dinteg(hcalFpsijreg,xmin,xmax,eps)

c ... plus distributions
      calFpsijplus = 0d0 ! no such contribution 

c ... local part prop. to delta distribution
      calFpsijdelta = 0d0 ! initialize

      calFpsijdelta = calFpsijdel(x,Q2,i,nf,iord,isf) 

      calFpsij = calFpsijreg + calFpsijplus + calFpsijdelta

      end


c .......................................................................

      double precision function calFpsijdel(x,Q2,i,nf,iord,isf)
      implicit none
      double precision x,Q2
      integer i,nf,iord,isf  

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      double precision pdfi
      double precision Cpsdel,C2psdelta

      double precision pdf(-6:6)

      calFpsijdel = 0d0

c ... pdfs at x
      call mypdf(x,Q2,pdf,IsetPDF,iord,A1,Z1)

      pdfi = pdf(i) + pdf(-i)

c ... Wilson coefficient 
      if (isf .eq. 0) then ! C_L
         Cpsdel = 0d0
      elseif (isf .eq. 2) then ! C_2
         Cpsdel = C2psdelta(x,Q2,nf,iord) ! doesn't depend on x
      else 
         stop 'isf value not provided'   
      end if   

      calFpsijdel = pdfi * Cpsdel

      end

c .......................................................................

c .......................................................................

      double precision function hcalFgreg(z)
      implicit none
      double precision z

      double precision var(1:10)
      integer ivar(1:10)
      common / tmp_var / var,ivar

      double precision glu,Cgreg
      double precision CLgregular,C2gregular

      double precision x,Q2
      integer nf,iord,isf

      double precision pdf(-6:6)

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      hcalFgreg = 0d0
      glu = 0d0

      x = var(1)
      Q2 = var(2)
      nf = ivar(1)
      iord = ivar(2)
      isf = ivar(3)

c ... gluon distribution at x/z
      call mypdf(x/z,Q2,pdf,IsetPDF,iord,A1,Z1)
      glu = pdf(0)

c ... Wilson coefficient 
      if (isf .eq. 0) then
         Cgreg = CLgregular(z,Q2,nf,iord)
      else if (isf .eq. 2) then   
         Cgreg = C2gregular(z,Q2,nf,iord)
      else 
         stop 'isf value not provided'   
      end if   

      hcalFgreg = glu/z * Cgreg

      end

c .......................................................................

      double precision function hcalFnsplus(z)
      implicit none
      double precision z

      double precision var(1:10)
      integer ivar(1:10)
      common / tmp_var / var,ivar

      double precision Cnspl
      double precision CLnsplus,C2nsplus

      double precision x,Q2
      integer i,nf,iord,isf

      double precision pdf(-6:6)

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      double precision pdfi,pdfix

      hcalFnsplus = 0d0
      pdfi = 0d0

      x = var(1)
      Q2 = var(2)
      nf = ivar(1)
      iord = ivar(2)
      i = ivar(3)
      isf = ivar(4)

c ... pdfs at x/z
      call mypdf(x/z,Q2,pdf,IsetPDF,iord,A1,Z1)
      pdfi = pdf(i) + pdf(-i)

c ... pdfs at x
      call mypdf(x,Q2,pdf,IsetPDF,iord,A1,Z1)
      pdfix = pdf(i) + pdf(-i)

c ... Wilson coefficient 
      if (isf .eq. 0) then ! C_L
         Cnspl = CLnsplus(z,Q2,nf,iord)
      elseif (isf .eq. 2) then ! C_2
         Cnspl = C2nsplus(z,Q2,nf,iord)
      else 
         stop 'isf value not provided'   
      end if   

      hcalFnsplus = Cnspl*(pdfi/z -pdfix)

      end

c .......................................................................

      double precision function hcalFnsreg(z)
      implicit none
      double precision z

      double precision var(1:10)
      integer ivar(1:10)
      common / tmp_var / var,ivar

      double precision Cnsreg
      double precision CLnsregular,C2nsregular

      double precision x,Q2
      integer i,nf,iord,isf

      double precision pdf(-6:6)

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      double precision pdfi

      hcalFnsreg = 0d0
      pdfi = 0d0
      Cnsreg = 0d0

      x = var(1)
      Q2 = var(2)
      nf = ivar(1)
      iord = ivar(2)
      i = ivar(3)
      isf = ivar(4)

c ... gluon distribution at x/z
      call mypdf(x/z,Q2,pdf,IsetPDF,iord,A1,Z1)
      
      pdfi = pdf(i) + pdf(-i)

c ... Wilson coefficient 
      if (isf .eq. 0) then ! C_L
         Cnsreg = CLnsregular(z,Q2,nf,iord)
      elseif (isf .eq. 2) then ! C_2   
         Cnsreg = C2nsregular(z,Q2,nf,iord)
      else 
         stop 'isf value not provided'   
      end if   

      hcalFnsreg = pdfi/z * Cnsreg

      end

c .......................................................................

c ... original integrand
      double precision function hcalFpsijreg(z)
      implicit none
      double precision z

      double precision var(1:10)
      integer ivar(1:10)
      common / tmp_var / var,ivar

      double precision Cpsreg
      double precision CLpsregular,C2psregular

      double precision x,Q2
      integer i,nf,iord,isf

      double precision pdf(-6:6)

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      double precision pdfi
      integer i1

      hcalFpsijreg = 0d0

      x = var(1)
      Q2 = var(2)
      nf = ivar(1)
      iord = ivar(2)
      isf = ivar(3)
      i = ivar(4)

c ... PDFs at x/z
      call mypdf(x/z,Q2,pdf,IsetPDF,iord,A1,Z1)

      pdfi = pdf(i) + pdf(-i)

c ... Wilson coefficient 
      if (isf .eq. 0) then ! C_L
         Cpsreg = CLpsregular(z,Q2,nf,iord)
      elseif (isf .eq. 2) then ! C_2
         Cpsreg = C2psregular(z,Q2,nf,iord)
      else 
         stop 'isf value not provided'   
      endif   

      hcalFpsijreg = pdfi/z * Cpsreg

      end 
