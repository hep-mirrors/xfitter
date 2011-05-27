C
C === FILE IS REQUIRED BY "DIPOLE" SUBROUTINES
C
      integer nf2max
      parameter (nf2max=10000)
      integer nn
      real x(nf2max),q2(nf2max),y(nf2max)
      real xsect(nf2max),xsecte(nf2max)
      common/fract/x,q2,y,xsect,xsecte,nn
