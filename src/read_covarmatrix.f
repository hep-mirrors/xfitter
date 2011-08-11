
      subroutine Read_CoVarMatrix

      implicit none
      include 'ntot.inc'
      include 'indata.inc'
      include 'covar_chi2.inc'
      include 'systematics.inc'
      integer i,j

      write(6,*) 'in Read_CoVarMatrix'
      write(6,*) 'npoints = ',npoints

      open(unit=22,file='covmat.dat')
       do i=1,npoints
        do j=1,npoints
         read(22,*) cov(i,j)
        enddo
       enddo
      close(22)

      return
      end
