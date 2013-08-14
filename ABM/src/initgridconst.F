      subroutine initgridconst
    
      implicit real*8 (a-h,o-z) 

      real*8 cc(2)

      include 'CONSTCOM.'
      include 'APSCOM6.'
      include 'PDFCOM.'
      include 'PRECCOM.'

! set-up the math constants
      data pi /3.1415927d0/
      data zeta2 /1.6449340668482264365D0/
      data zeta3 /1.2020569031595942854D0/
      data zeta4 /1.0823232337111381916D0/

! set-up the QCD constants
      data nc,nf,nfe,cf,cg,tr /3,6,3,1.33333d0,3.d0,0.5d0/

! set-up the CKM matrix (PDG 2010)
      data ckm / 0.97428d0, 0.2253d0,   3.47d-3
     ,          ,0.2252d0,  0.97345d0,  4.1d-2
     ,          ,8.62d-3,   4.03d-2,    0.999152d0  /

! sum of the quark charges squared for different number of flavours
      data qsum /0.1111d0, 0.5555d0, 0.6666d0 
     ,          ,1.1111d0, 1.2222d0, 1.6666d0/
! sum of the quark charges for different number of flavours
      data qsum0 /-0.3333d0, 0.3333d0, 0d0 
     ,,            0.6666d0, 0.3333d0, 0.3333d0/

! Weak-mixing angle sine squared (PDG 2010)
      data sintw2 /0.23116d0/
! Z-boson mass (PDG 2010)
      rmassp(41)=91.1876
! W-boson mass (PDG 2010)
      rmassp(42)=80.358
! Vector and axial quark couplings
      vqu=0.5-4./3.*sintw2
      aqu=0.5
      vqd=-0.5+2./3.*sintw2
      aqd=-0.5
      vaq2u=vqu**2+aqu**2
      vaq2d=vqd**2+aqd**2
      vaq2sum(3)=vaq2u+2*vaq2d
      vaq2sum(4)=vaq2sum(3)+vaq2u
      vaq2sum(5)=vaq2sum(4)+vaq2d
      vaqsum(3)=4./3.*vqu-4./3.*vqd
      vaqsum(4)=vaqsum(3)+4./3.*vqu
      vaqsum(5)=vaqsum(4)-2./3.*vqd
! Vector and axial lepton couplings
      vlu=0.5
      alu=0.5
      vld=-0.5+2*sintw2
      ald=-0.5
      val2u=vlu**2+alu**2
      val2d=vld**2+ald**2

! A step for calculation of the PDF derivatives
      data delder /0.01/
! A tolerance factor for the numerical solution of the 
! renorm-group equation for the strong coupling constant
      data alphastol /1d-4/

!  Default number of intervals for the integration in the heavy-quark 
!  matching conditions and the DIS structure functions 
      data nmthq, nf2hq, nflhq, nf3hq /1, 1, 1, 1/
!  and the light-parton structure functions
      data nflqcd,nf2qcd1,nf2qcd2,nf3qcd /1, 1, 1, 1/

c set-up of the grids
      q2min=0.8d0
      q2max=0.99d8
c      q2min=1.5d0
c      q2max=1.1d5
      q2ini=9d0
      xbmin=5d-6
c      xbmin=1d-5
      nxpgrid=100
      nxmgrid=100
      nspgrid=60
      nsmgrid=60
      xbmax=0.999d0
      x1=0.3
      xlog1=log(x1)
      xlog2=log(1-x1)

!  Take the 3-flavour NNLO PDFs as a default
      kschemepdf=0
      kordpdf=2

!  Set the number of the PDF members in each set 
!  with no PDF uncertainties by default
      npdftot=0
!  Take the central PDF member as a default
      kpdfset=0

#ifdef ABKM09
      npdftot=25
#endif
#ifdef ABM11
      npdftot=28
#endif
#ifdef MSTW08
      npdftot=40
#endif
#ifdef JR09
      npdftot=26
#endif
#ifdef CT10
      npdftot=0
      if (kschemepdf.eq.2) npdftot=52
#endif
#ifdef NN21
      npdftot=99
#endif
#ifdef HERAPDF1
      npdftot=20
#endif
#ifdef USER
      npdftot=0          !  consider only PDF central values
#endif
#ifdef TOYLH
      npdftot=0          !  consider only PDF central values
#endif

!  Take square of the CKM matrix elements to save computing time
      do i=1,3
      do j=1,3
        ckm2(i,j)=ckm(i,j)**2
      end do
      end do

!  Default values: 

! c-quark mass and charge
      rmass(8)=1.5d0
      rcharge(8)=0.6666666d0
! b-quark mass and charge
      rmass(10)=4.5d0
      rcharge(10)=0.3333333d0
! t-quark mass and charge
      rmass(12)=173d0
      rcharge(12)=0.6666666d0

! NNLO approximation for the strong coupling 
      kordalps=2
! NLO approximation for the heavy-quark contribution to DIS
      kordhq=1
! NNLO approximation for the light-parton contribution to DIS
      kordf2=2
      kordfl=2
      kordf3=2

! The pole mass definition is used as a default
      msbarm=.false.
! The non-singlet contribution is included to the heavy-quark 
! electroproduction by default
      hqnons=.true.
! Variant A of the NNLO term in the heavy-quark electroproduction Wilson
! coefficient of Ref. [1205.5727] is used as a default 
      ddnnlohq=0.

!  Set the factorization scale as sqrt(Q2*hqscale1 + 4m^2*hqscale2) for the 
!  pair heavy-quark DIS production and as sqrt(Q2*hqscale1 + m^2*hqscale2) 
!  for the single heavy-quark DIS production
      hqscale1=1d0
      hqscale2=1d0

! ren.scale=fac.scale as a default
      rscale=1d0

! The interpolation of the OMEs is applied
      omeint=.true.
      call omegridini(xbmin,nxpgrid,nsmgrid)
! The interpolation of the DIS coefficient functions is applied
      lpcint=.true.
      call cgridini(xbmin,nxpgrid,nsmgrid)

      return 
      end
