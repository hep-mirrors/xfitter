      subroutine initgridconst
    
      implicit real*8 (a-h,o-z) 

      real*8 cc(2)

      include 'CONSTCOM.'
      include 'APSCOM6.'
      include 'PDFCOM.'
      include 'PRECCOM.'

! Set-up the math constants
      data pi /3.1415927d0/
      data zeta2 /1.6449340668482264365D0/
      data zeta3 /1.2020569031595942854D0/
      data zeta4 /1.0823232337111381916D0/

! Set-up the QCD constants
      data nc,nf,nfe,cf,cg,tr /3,6,3,1.33333d0,3.d0,0.5d0/

! Z-boson mass (PDG 2012)
      data rmassp(41) /91.1876d0/
! W-boson mass (PDG 2012)
      data rmassp(42) /80.385/
! c-quark mass and charge
      data rmass(8) /1.5d0/
      data rcharge(8) /0.6666666d0/
! b-quark mass and charge
      data rmass(10) /4.5d0/
      data rcharge(10) /0.3333333d0/
! t-quark mass and charge
      data rmass(12) /172d0/
      data rcharge(12) /0.6666666d0/

! The 4-flavour matching is performed by default at the scale of c-quark mass 
      data vfnth(4) /1.5d0/   
! The 5-flavour matching is performed by default at the scale of b-quark mass
      data vfnth(5) /4.5d0/ 

! Set-up the CKM matrix (PDG 2010)
      data ckm / 0.97428d0, 0.2253d0,   3.47d-3
     ,          ,0.2252d0,  0.97345d0,  4.1d-2
     ,          ,8.62d-3,   4.03d-2,    0.999152d0  /

! Sum of the quark charges squared for different number of flavours
      data qsum /0.1111d0, 0.5555d0, 0.6666d0 
     ,          ,1.1111d0, 1.2222d0, 1.6666d0/
! Sum of the quark charges for different number of flavours
      data qsum0 /-0.3333d0, 0.3333d0, 0d0 
     ,,            0.6666d0, 0.3333d0, 0.3333d0/

! Weak-mixing angle sine squared (PDG 2012)
      data sintw2 /0.23116d0/

! A step for calculation of the PDF derivatives
      data delder /0.01/

! A tolerance factor for the numerical solution of the 
! renorm-group equation for the strong coupling constant
      data alphastol /1d-4/

! Default number of intervals for the integration in the heavy-quark 
! matching conditions and the DIS structure functions 
      data nmthq, nf2hq, nflhq, nf3hq /1, 1, 1, 1/
! and the light-parton structure functions
      data nflqcd,nf2qcd1,nf2qcd2,nf3qcd /1, 1, 1, 1/
! Default value for the number intervals for computing convolutions of 
! splitting functions with the PDFs in the evolution equations 
      data kprecdq /1/

! Take the 3-flavour NNLO PDFs as a default
      data kschemepdf /0/
      data kordpdf /2/

! Set the number of the PDF members in each set 
! with no PDF uncertainties by default
      data npdftot /0/
! Take the central PDF member as a default
      data kpdfset /0/

! NNLO approximation for the strong coupling 
      data kordalps /2/
! NLO approximation for the heavy-quark contribution to DIS
      data kordhq /1/
! NNLO approximation for the light-parton contribution to DIS
      data kordf2,kordfl,kordf3 /2,3,2/

! The pole mass definition is used as a default
      data msbarm /.false./
! The non-singlet contribution is included to the heavy-quark 
! electroproduction by default
      data hqnons /.true./
! The strong coupling constant in the 4- and 5-flavour schemes are 
! generated from the 3-flavour input value
      data alsmz /.false./
! The number of loops in running of \alpha_s is fixed to 3. 
      data vloop /.false./

! The BMSN prescriptoin for the VFN scheme is based on the heavy-quark PDFs 
! obtained from the fixed-order-perturbation-theory NLO matching conditions 
      data bmsnfopt /.true./
      data bmsnnlo /.true./

! Variant A of the NNLO term in the heavy-quark electroproduction Wilson
! coefficient of Ref. [1205.5727] is used as a default 
      data ddnnlohq /0d0/

! Set the factorization scale as sqrt(Q2*hqscale1 + 4m^2*hqscale2) for the 
! pair heavy-quark DIS production and as sqrt(Q2*hqscale1 + m^2*hqscale2) 
! for the single heavy-quark DIS production
      data hqscale1, hqscale2 /1d0, 1d0/

! Ren.scale=fac.scale as a default
      data rscale /1d0/

! The interpolation of the OMEs is applied
      data omeint /.true./
! The interpolation of the DIS Wilson coefficients is applied
      data lpcint /.true./

! Set-up of the PDF and interpolating grids
!   the margins of x-grid
      data xbmin, xbmax /5d-6, 0.999d0/
!   midpoint of the x-grid: log(x) interpolation below x1 and log(1-x) above x1
      data x1 /0.2d0/
!   the number of nodes in the x-grid above and below x1, respectively
      data nxpgrid, nxmgrid /100, 100/
!   the margins of Q2-grid
      data q2min, q2max  /0.8d0, 0.99d8/
!   midpoint of the Q2-grid  
      data q2ini /9d0, 9d0, 9d0, 9d0, 9d0, 9d0, 9d0, 9d0/
!   the number of nodes in the Q2-grid above and below q2ini, respectively
      data nspgrid, nsmgrid /60, 60/
!   nominal x-grid is used by default 
      data khalf /0/ 
     
! Create interpolation grids for the OMEs and Wilson coefficients
      xlog1=log(x1)
      xlog2=log(1-x1)
      call omegridini
      call cgridini
c RP this call will be needed for evolution, currently not used
c     call pgridini

!  Take square of the CKM matrix elements to save computing time
      do i=1,3
      do j=1,3
        ckm2(i,j)=ckm(i,j)**2
      end do
      end do

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

      return 
      end
