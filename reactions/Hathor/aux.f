      real*8 function ttbari(xx)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      real*8 xx(10)

C      ttbari=ff(xx)
c      print *,xx(1),xx(2),ttbari

      return 
      end
c--------------
      subroutine ttbarr(ll,xx,aa,bb)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      real*8 aa(2),bb(2),xx(2)

      del=1d-7

      aa(1)=del
      bb(1)=1-del
      aa(2)=del
      bb(2)=1-del

      return 
      end
