        real function logshift(mmu,ssig,rrunif) 
c        real function logshift(mu,sig,runif) 
C-----------------------------------------------------------------------
C-
C-   Purpose and Methods: 
C-
C-   Inputs  :
C-   Outputs :
C-   Controls:
C-
C-   Created  12-JUN-2008   Voica Radescu
C-
C-----------------------------------------------------------------------
      IMPLICIT NONE

      
      real zeroth, ANS,ex2,runif
      real mu, sig,x2, mu2, sig2,z1,z2
      external zeroth
      real mmu,rrunif, ssig
c      integer iseed,isdrn,i


      COMMON/PARAM/mu,sig,runif

      mu=mmu
      sig=ssig
      runif=rrunif
c      print*,'logshift', mmu, ssig,rrunif

c      call system_clock(iseed)
c      isdrn=iseed
c      call rmarin(isdrn,0,0)
c      call rluxgo(3,isdrn,0,0)
c      print*,'initialize smeering with a seed isdrn',isdrn

c      do i=1,1000
c         call ranlux(runif,1)
c         print*,'uniform random numbers   ', runif

c         z1 = zeroth(0.0001)
c         z2 = zeroth(10000.)

c      print *,'zeroth(0)=',z1,  ' zeroth(10000.)',z2


c   CALL tZERO(A_in,B_in,X0_out,R_out,EPS_in,MAXF_in,F_in)
      call rzero(0.0001,10000.,x2,ex2,0.0001,10000,zeroth)
c         print *,'result:',x2,ex2
c         z1 = zeroth(x2)
      logshift=x2


 



C-----------------------------------------------------------------------
 999  CONTINUE
      END




