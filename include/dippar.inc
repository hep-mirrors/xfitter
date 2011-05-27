C ================================================
C === FILE IS CALLED BY SUBROUTINES IN "fcn.f" ===
C ================================================
C === FILE IS REQUIRED BY "DIPOLE" SUBROUTINES ===
C
      real*8 sig0,xlam,x0,xm
      common/pass1/ sig0,xlam,x0,xm
      real*8 cm
      integer icharm
      common/chsteer/ cm,icharm
