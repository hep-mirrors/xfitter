!------------------
      real*8 function f2cqqpm(nb,nt,ni,xb,q2)
      implicit double precision (a-h,o-z)

      include 'CONSTCOM.'
      include 'PDFCOM.'

      f2cqqpm=0.

! electron-proton NC scattering
      if(nb.eq.3.and.nt.eq.1.and.ni.eq.22) then
! d- and s-quark contributions 
        f2cqqpm=f2cqqpm+1./9.
     *  * (xqg(2,xb,q2,kschemepdf)+xqg(3,xb,q2,kschemepdf)
     +  +  xqg(6,xb,q2,kschemepdf)+xqg(7,xb,q2,kschemepdf))

! u-quark contribution 
        f2cqqpm=f2cqqpm+4./9.
     *  * (xqg(4,xb,q2,kschemepdf)+xqg(5,xb,q2,kschemepdf))
      end if

! electron-neutron NC scattering
      if(nb.eq.3.and.nt.eq.8.and.ni.eq.22) then
! d- and s-quark contributions 
        f2cqqpm=f2cqqpm+1./9.
     *  * (xqgn(2,xb,q2,kschemepdf)+xqgn(3,xb,q2,kschemepdf)
     +  +  xqgn(6,xb,q2,kschemepdf)+xqgn(7,xb,q2,kschemepdf))
! u-quark contribution 
        f2cqqpm=f2cqqpm+4./9.
     *  * (xqgn(4,xb,q2,kschemepdf)+xqgn(5,xb,q2,kschemepdf))
      end if

! neutrino-proton CC scattering
      if(nb.eq.6.and.nt.eq.1.and.ni.eq.24) then
! s-quark contribution 
        f2cqqpm=f2cqqpm+xqg(6,xb,q2,kschemepdf)*ckm(2,2)**2
! d-quark contribution 
        f2cqqpm=f2cqqpm+xqg(2,xb,q2,kschemepdf)*ckm(1,2)**2
      end if

! antineutrino-proton CC scattering
      if(nb.eq.7.and.nt.eq.1.and.ni.eq.24) then
! s-quark contribution 
        f2cqqpm=f2cqqpm+xqg(7,xb,q2,kschemepdf)*ckm(2,2)**2
! d-quark contribution 
        f2cqqpm=f2cqqpm+xqg(3,xb,q2,kschemepdf)*ckm(1,2)**2
      end if

! neutrino-neutron CC scattering
      if(nb.eq.6.and.nt.eq.8.and.ni.eq.24) then
! s-quark contribution 
        f2cqqpm=f2cqqpm+xqgn(6,xb,q2,kschemepdf)*ckm(2,2)**2
! d-quark contribution 
        f2cqqpm=f2cqqpm+xqgn(2,xb,q2,kschemepdf)*ckm(1,2)**2
      end if

! antineutrino-neutron CC scattering
      if(nb.eq.7.and.nt.eq.8.and.ni.eq.24) then
! s-quark contribution 
        f2cqqpm=f2cqqpm+xqgn(7,xb,q2,kschemepdf)*ckm(2,2)**2
! d-quark contribution 
        f2cqqpm=f2cqqpm+xqgn(3,xb,q2,kschemepdf)*ckm(1,2)**2
      end if

      return
      end
!------------------
      real*8 function f3cqqpm(nb,nt,ni,xb,q2)
      implicit double precision (a-h,o-z)

      f3cqqpm=0.

! neutrino-nucleon CC scattering
      if (nb.eq.6) f3cqqpm=f2cqqpm(nb,nt,ni,xb,q2)
! antineutrino-nucleon CC scattering
      if (nb.eq.7) f3cqqpm=-f2cqqpm(nb,nt,ni,xb,q2)

      return
      end
c-----------------------
      real*8 function xpscqqpm(xb,q2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      include 'PDFCOM.'

      xpscqqpm=0.
      do i=2,7
        xpscqqpm=xpscqqpm+xqg(i,xb,q2,kschemepdf)
      end do

      return 
      end
!------------------
      real*8 function f2qpm(nb,nt,ni,xb,q2)
      implicit double precision (a-h,o-z)

      include 'CONSTCOM.'
      include 'PDFCOM.'

      f2qpm=0.

! electron-nucleon NC scattering
      if (ni.eq.22) f2qpm=f2cqqpm(nb,nt,ni,xb,q2)

! neutrino-proton CC scattering
      if(nb.eq.6.and.nt.eq.1.and.ni.eq.24) then
! s-quark contribution 
        f2qpm=f2qpm+2*xqg(6,xb,q2,kschemepdf)*ckm(2,1)**2
! d-quark contribution 
        f2qpm=f2qpm+2*xqg(2,xb,q2,kschemepdf)*ckm(1,1)**2
! u-quark contribution 
        f2qpm=f2qpm+2*xqg(5,xb,q2,kschemepdf)
      end if

! neutrino-neutron CC scattering
      if(nb.eq.6.and.nt.eq.8.and.ni.eq.24) then
! s-quark contribution 
        f2qpm=f2qpm+2*xqgn(6,xb,q2,kschemepdf)*ckm(2,1)**2
! d-quark contribution 
        f2qpm=f2qpm+2*xqgn(2,xb,q2,kschemepdf)*ckm(1,1)**2
! u-quark contribution 
        f2qpm=f2qpm+2*xqgn(5,xb,q2,kschemepdf)
      end if

! antineutrino-proton CC scattering
      if(nb.eq.7.and.nt.eq.1.and.ni.eq.24) then
! s-quark contribution 
        f2qpm=f2qpm+2*xqg(7,xb,q2,kschemepdf)*ckm(2,1)**2
! d-quark contribution 
        f2qpm=f2qpm+2*xqg(3,xb,q2,kschemepdf)*ckm(1,1)**2
! u-quark contribution 
        f2qpm=f2qpm+2*xqg(4,xb,q2,kschemepdf)
      end if

! antineutrino-neutron CC scattering
      if(nb.eq.7.and.nt.eq.8.and.ni.eq.24) then
! s-quark contribution 
        f2qpm=f2qpm+2*xqgn(7,xb,q2,kschemepdf)*ckm(2,1)**2
! d-quark contribution 
        f2qpm=f2qpm+2*xqgn(3,xb,q2,kschemepdf)*ckm(1,1)**2
! u-quark contribution 
        f2qpm=f2qpm+2*xqgn(4,xb,q2,kschemepdf)
      end if

      return 
      end
!------------------
      real*8 function f3qpm(nb,nt,ni,xb,q2)
      implicit double precision (a-h,o-z)

      include 'CONSTCOM.'
      include 'PDFCOM.'

      f3qpm=0.

! neutrino-proton CC scattering
      if(nb.eq.6.and.nt.eq.1.and.ni.eq.24) then
! s-quark contribution 
        f3qpm=f3qpm + 2*xqg(6,xb,q2,kschemepdf)*ckm(2,1)**2
! d-quark contribution 
        f3qpm=f3qpm + 2*xqg(2,xb,q2,kschemepdf)*ckm(1,1)**2
! u-quark contribution 
        f3qpm=f3qpm - 2*xqg(5,xb,q2,kschemepdf)
      end if

! neutrino-neutron CC scattering
      if(nb.eq.6.and.nt.eq.8.and.ni.eq.24) then
! s-quark contribution 
        f3qpm=f3qpm + 2*xqgn(6,xb,q2,kschemepdf)*ckm(2,1)**2
! d-quark contribution 
        f3qpm=f3qpm + 2*xqgn(2,xb,q2,kschemepdf)*ckm(1,1)**2
! u-quark contribution 
        f3qpm=f3qpm - 2*xqgn(5,xb,q2,kschemepdf)
      end if

! antineutrino-proton CC scattering
      if(nb.eq.7.and.nt.eq.1.and.ni.eq.24) then
! s-quark contribution 
        f3qpm=f3qpm - 2*xqg(7,xb,q2,kschemepdf)*ckm(2,1)**2
! d-quark contribution 
        f3qpm=f3qpm - 2*xqg(3,xb,q2,kschemepdf)*ckm(1,1)**2
! u-quark contribution 
        f3qpm=f3qpm + 2*xqg(4,xb,q2,kschemepdf)
      end if

! antineutrino-neutron CC scattering
      if(nb.eq.7.and.nt.eq.8.and.ni.eq.24) then
! s-quark contribution 
        f3qpm=f3qpm - 2*xqgn(7,xb,q2,kschemepdf)*ckm(2,1)**2
! d-quark contribution 
        f3qpm=f3qpm - 2*xqgn(3,xb,q2,kschemepdf)*ckm(1,1)**2
! u-quark contribution 
        f3qpm=f3qpm + 2*xqgn(4,xb,q2,kschemepdf)
      end if

      f3qpm=f3qpm/xb
 
      return 
      end
c-----------------------
      real*8 function xpsqpm(xb,q2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      include 'PDFCOM.'

      xpsqpm=0.
      do i=2,7
        xpsqpm=xpsqpm+xqg(i,xb,q2,kschemepdf)
      end do

      return 
      end
c-----------------------
      real*8 function xpsqpm2(xb,q2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      include 'PDFCOM.'

      xpsqpm2=0.
      do i=2,5
        xpsqpm2=xpsqpm2+xqg(i,xb,q2,kschemepdf)
      end do

      return 
      end
