
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


      CALL DINV  (NSYS,sysa,NSYSMAX,IR,IFAIL)


      return
      end


      Subroutine CompressSystematics
C-----------------------------------------
C
C-----------------------------------------
      implicit none
C-----------------------------------------
      include 'steering.inc'
      include 'datasets.inc'      
      include 'ntot.inc'
      include 'systematics.inc'
      integer i,j,nsyssav
      double precision anorm
C-----------------------------------------
      nsyssav = nsys
C Diagonal reference
      do i=1,NSYS
         CompressIdx(i) = i
      enddo

      i = 1
      do while ( i.le.NSYS) 
C Check if
         do j=1,Npoints
            if (Beta(i,j).ne.0) then
C Something present, to the next potential candidate for empty line:
               i = i + 1
               goto 17
            endif 
         enddo
C Empty line... Move top line here
         
         CompressIdx(i) = CompressIdx(NSYS)
         do j=1,Npoints
            Beta(i,j) = Beta(NSYS,j)
         enddo
         NSYS = NSYS - 1
         

 17      continue
      enddo

      print '(''Compress the systematcs array'')'
      print '(''Dimension before conversion='',i4,'' after='',i4)',nsyssav,nsys

      end
