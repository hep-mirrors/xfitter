
c ... PDF interface
      subroutine mypdf(x,Q2,pdf,iset,iorder,A1,Z1)
C============================================================================= 
C      Front end program for  N3LO pdf's: Fred Olness  27 April 2016 (needs checking)
C
C============================================================================= 
      implicit none

      include 'couplings.inc'  !*** Get quark masses from this:

      double precision x,Q2 ! Input
      integer iset,iorder   ! Input
      double precision A1,Z1           ! Input
      double precision pdf(-6:6),xpdf(-6:6) ! Output

      double precision Ctq6Pdf
      double precision Q      
      integer i,nf
      integer isetold
      save isetold

      integer pdfini, nfix
      save pdfini

      double precision singlet,xm0,xmc,xmb,xmt,HMASS,q2c,q2b,q2t
      integer ierr

c      common /fred/ xm0,xmc,xmb,HMASS
c     =================================================================================================
 
!     taken from couplings.inc
      xmc=mch                   
      xmb=mbt     
      xmt=mtp


      do i=-6,6,1
         pdf(i)=0.0d0
      enddo
      if(x.ge.0.99d0) return

      Q=sqrt(Q2)

      ierr = 0
      if ((x.gt.1d0) .or. (x.lt.0d0)) then 
         ierr = -1
         write(6,*)  ' error: bad x in mypdf = ',x
      end if   

c     =================================================================================================
c     use qcdnum pdfs
c      call fpdfxq(1,x,q2,xpdf,1)                  !interpolate all pdf's
      call hf_get_pdfs(x,q2,xpdf)                  !interpolate all pdf's

c     use qcdnum mass values
c      call getcbt(nfix,q2c,q2b,q2t)  !**** nfix=0 FFS, or nfix =3,4,5,6 VFNS


      do i=-6,6,1
         pdf(i)=xpdf(i)/x
      enddo

c   convert to CTEQ convention: u=1, d=2
      pdf(+2) = xpdf(+1)/x
      pdf(+1) = xpdf(+2)/x
      pdf(-2) = xpdf(-1)/x
      pdf(-1) = xpdf(-2)/x
    



      if(q.lt.xmc) then
         pdf(+4)=0.0d0
         pdf(-4)=0.0d0
      endif

      if(q.lt.xmb) then
         pdf(+5)=0.0d0
         pdf(-5)=0.0d0
      endif

      if(q.lt.xmt) then
         pdf(+6)=0.0d0
         pdf(-6)=0.0d0
      endif

c     UPDATE: IF PDF IS NEGATIVE: SET TO ZERO: FIO 21 Nov. 2011
      do i=-6,6,1
         if(pdf(i).lt.0.0d0) pdf(i)=0.0d0
      enddo


      return
      end

c     =================================================================================================
c     =================================================================================================
c     =================================================================================================

c ... alpha_s interface
      double precision function myalphas(Q2,iorder,iset) ! alpha_s/(4 pi)
      implicit none
      
      double precision Q2
      integer iorder,iset, ierr, nfout,nf

      double precision Q
      double precision ascteq6,as_4pi

      double precision pi
      double precision as,tmp, alfs
      double precision asfunc

      Q = sqrt(Q2)
      pi = dacos(-1d0)
C============================================================================= 
C============================================================================= 

C     OVERRIDE WITH QCDNUM STUFF
      pi = dacos(-1d0)
      alfs= ASFUNC(q2,NF,IERR)
      myalphas=alfs/(4.0d0*pi)
      RETURN
      end
C============================================================================= 
C============================================================================= 

c ... ew_couplings.f TS/IS
c ... 03.07.2011

*************************************************************************

c ... i: choose flavor: i = 1: up quark, 2: down quark, etc
c ... iboson = 0,1: gamma gamma + 2 gamma Z + Z Z, gamma gamma 
c              2,3,4: gamma Z + Z gamma, Z Z, W W
c ... Output: electroweak couplings squared
c ... see Eq. (1.45) in [3]
      subroutine get_couplingsINGO(Q2,i,iboson,aplus,aminus,abarplus,
     #                     abarminus)
      implicit none
      double precision Q2 ! Input
      integer i,iboson ! Input
      double precision aplus,aminus,abarplus,abarminus ! Output

      double precision eq(1:6)

      double precision vq(1:6),vu,vd
      double precision aq(1:6)
      double precision ve,ae

      double precision sw2

      double precision GF,MZ,MW,pi,alphaem
      double precision chiz,chiw
      double precision gf_HF, convfac_HF
      double precision alphaem_HF,sin2thw_HF,cos2thw_HF
      double precision Mz_HF, Mw_HF, Mh_HF


      common/boson_masses/Mz_HF, Mw_HF, Mh_HF
      common/constants/ gf_HF, convfac_HF
      common/ew_couplings/alphaem_HF,sin2thw_HF,cos2thw_HF


      aplus = 0d0
      aminus = 0d0

!      GF = 1.16637d-5 ! Fermi constant
!      MZ = 91.1876d0
!      MW = 80.398d0
!      alphaem = 1d0/137d0
!      sw2 = 0.232d0 ! sin^2 theta_w
      GF = gf_HF ! Fermi constant
      MZ = mz_HF
      MW = mw_HF
      alphaem = alphaem_HF
      sw2 = sin2thw_HF
      pi = dacos(-1d0)
c ... later to be Q2 dependent



      vu = 1d0/2d0 - 4d0/3d0 * sw2
      vd = -1d0/2d0 + 2d0/3d0 * sw2

c ... Z vector and axial-vector coupling to quarks
      vq(1) = vu
      vq(4) = vu
      vq(6) = vu

      vq(2) = vd
      vq(3) = vd
      vq(5) = vd

      aq(1) = 1d0/2d0
      aq(4) = 1d0/2d0
      aq(6) = 1d0/2d0

      aq(2) = -1d0/2d0
      aq(3) = -1d0/2d0
      aq(5) = -1d0/2d0

c ... Z-coupling to charged lepton
      ve = -1d0/2d0 + 2d0 * sw2 
      ae = -1d0/2d0 

      eq(1) = 2d0/3d0
      eq(4) = 2d0/3d0
      eq(6) = 2d0/3d0

      eq(2) = -1d0/3d0      
      eq(3) = -1d0/3d0      
      eq(5) = -1d0/3d0      

      chiz = GF * MZ**2/(2*dsqrt(2d0)*pi*alphaem) * Q2/(Q2+MZ**2)

      chiw = 1d0/2d0*GF* MW**2/(2*dsqrt(2d0)*pi*alphaem) * Q2/(Q2+MW**2)

      if (iboson .eq. 0) then ! gamma^2 + 2 gamma Z + Z^2
         aplus = eq(i)**2 - 2d0*ve*eq(i)*vq(i)*chiz +
     #           (ve**2+ae**2)*(vq(i)**2+aq(i)**2)*chiz**2

         aminus = -2d0*eq(i)*aq(i)*ae*chiz+4d0*ve*ae*vq(i)*aq(i)*chiz**2
      endif   

      if (iboson .eq. 1) then ! gamma gamma
         aplus = eq(i)**2 
         aminus = 0d0
      endif   

      if (iboson .eq. 2) then ! gamma Z + Z gamma
         aplus = - 2d0*ve*eq(i)*vq(i)*chiz 
         aminus = -2d0*eq(i)*aq(i)*ae*chiz
      endif   
      
      if (iboson .eq. 3) then ! Z Z
         aplus = (ve**2+ae**2)*(vq(i)**2+aq(i)**2)*chiz**2
         aminus = 4d0*ve*ae*vq(i)*aq(i)*chiz**2
      endif   

      if (iboson .eq. 4) then ! W W 
         aplus = 1d0
         aminus = 1d0
         stop 'Not yet included'
      endif   

      abarplus = aplus
      abarminus = -aminus

      end

*************************************************************************

*=============================================================================================
* File: li2.f 
*
* ..A routine for the dilogarithm, by Jos Vermaseren
*
       function li2(x)
       implicit real*8  (a-z)
       dimension b(8)
       integer ncall
       data ncall/0/,pi6/1.644934066848226d+00/,een,vier/1.d+00,.25d+00/
       ncall = 0
       if(ncall.eq.0)go to 2
1      if(x.lt.0)go to 3
       if(x.gt.0.5)go to 4
       z=-dlog(1.-x)
7      z2=z*z
       li2=z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een)-z2*vier
       if(x.gt.een)li2=-li2-.5*u*u+2.*pi6
       return
2      b(1)=een
       b(2)=een/36.
       b(3)=-een/3600.
       b(4)=een/211680.
       b(5)=-een/(30.*362880.d+00)
       b(6)=5./(66.*39916800.d+00)
       b(7)=-691./(2730.*39916800.d+00*156.)
       b(8)=een/(39916800.d+00*28080.)
       ncall=1
       go to 1
3      if(x.gt.-een)go to 5
       y=een/(een-x)
       z=-dlog(een-y)
       z2=z*z
       u=dlog(y)
       li2=z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een)-z2*vier-u*(z+.5*u)-pi6
       return
4      if(x.ge.een)go to 10
       y=een-x
       z=-dlog(x)
6      u=dlog(y)
       z2=z*z
       li2=-z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een-u)+z2*vier+pi6
       if(x.gt.een)li2=-li2-.5*z*z+pi6*2.
       return
5      y=een/(een-x)
       z=-dlog(y)
       z2=z*z
       li2=-z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een)-z2*vier
       return
10     if(x.eq.een)go to 20
       xx=1./x
       if(x.gt.2.)go to 11
       z=dlog(x)
       y=1.-xx
       go to 6
11     u=dlog(x)
       z=-dlog(1.-xx)
       go to 7
20     li2=pi6
       return
       end
c ... alphas_interface.f TS/IS
c ... 03.07.2011
