
      subroutine Systematics

c 
c cf CTEQ
c

      implicit none 

#include "ntot.inc"
#include "indata.inc"
#include "systematics.inc"
#include "for_debug.inc"

      integer ir(nsys)
      integer i,j,k,ifail
      double precision term


      do i=1,NSYS
       do j=1,NSYS

        sysa(i,j) = 0.d0
        if (j.eq.i) sysa(i,j) = 1.d0

        do k=1,npoints
         term = BETA(i,k)*BETA(j,k)
         term = term /ALPHA(k)**2
         sysa(i,j) = sysa(i,j) + term
        enddo

       enddo
      enddo


      CALL DINV  (NSYS,sysa,NSYSMAX,IR,IFAIL)


      return
      end


