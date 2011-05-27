
      subroutine Systematics

c 
c cf CTEQ
c

      implicit none 

      include 'ntot.inc'
      include 'systematics.inc'
      include 'for_debug.inc'

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

      if (DEBUG) then
        do i=1,nsys
         write(6,'(18f7.2)')
     +            sysa(i,1),sysa(i,2),sysa(i,3),sysa(i,4),
     +            sysa(i,5),sysa(i,6),sysa(i,7),sysa(i,8),
     +            sysa(i,9),sysa(i,10),sysa(i,11),sysa(i,12),
     +            sysa(i,13),sysa(i,14),sysa(i,15),sysa(i,16),
     +            sysa(i,17),sysa(i,18)
        enddo
      endif

      CALL DINV  (NSYS,sysa,NSYS,IR,IFAIL)

      if (DEBUG) then
        write(6,*)
        write(6,*) 'after inversion '
        do i=1,nsys
         write(6,'(18f5.2)') 
     +            sysa(i,1),sysa(i,2),sysa(i,3),sysa(i,4),
     +            sysa(i,5),sysa(i,6),sysa(i,7),sysa(i,8),
     +            sysa(i,9),sysa(i,10),sysa(i,11),sysa(i,12),
     +            sysa(i,13),sysa(i,14),sysa(i,15),sysa(i,16),
     +            sysa(i,17),sysa(i,18)
        enddo
       endif


      return
      end
