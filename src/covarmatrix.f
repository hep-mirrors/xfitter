      subroutine CoVarMatrix

      implicit none

      include 'ntot.inc'
      include 'indata.inc'
      include 'systematics.inc'
      include 'covar_chi2.inc'

      double precision tmp
      dimension tmp(NTOT,NTOT)

      integer ir(ntot)
      integer i,j,k
      integer ifail


      write(6,*)
      write(6,*) '--- Calculate the covariance matrix ---'
      write(6,*)


      do i=1,ntot
       do j=1,ntot
        cov(i,j) = 0.d0
        tmp(i,j) = 0.d0
        if (j.eq.i) tmp(j,j) = 1.d0
       enddo
      enddo

      do i=1,npoints
        cov(i,i) = alpha(i)**2

       do j=1,npoints

        do k=1,nsys
         cov(i,j) = cov(i,j) + beta(k,i)*beta(k,j)
        enddo

        tmp(i,j) = cov(i,j)

       enddo
      enddo


      CALL DINV  (ntot,tmp,ntot,IR,IFAIL)

      do i=1,npoints
       do j=1,npoints
        cov(i,j) = tmp(i,j)
       enddo
      enddo

      write(6,*)
      write(6,*) '--- covariance matrix done ---'
      write(6,*)


      open(unit=22,file='covmat.dat')
       do i=1,npoints
        do j=1,npoints
         write(22,*) cov(i,j)
        enddo
       enddo
      close(22)

      return
      end

