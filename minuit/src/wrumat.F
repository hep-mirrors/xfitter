c --------------------------------------------------------------------------------
	subroutine wrumat
c --------------------------------------------------------------------------------
c write the current version of the matrix umat(i,j) to a file
c format of the file is same as that of umat.inp which can be read by reumat: 
c          1 line  with npar = number of active parameters
c       npar lines with the actual parameter numbers of the active parameters
c    npar**2 lines with the values umat(i,j)
c
#include "./d506dp.inc"
#include "./d506cm.inc"

c fluma is the name of the file to be written (other names not used here)
	character*13 flgnu, fluma
	character*10 flname
	common /fileco/ flgnu, fluma, flname

	common /umatco/ hquad(MNI,MNI), umat(MNI,MNI), vmat(MNI,MNI), 
     &	                tvec(MNI), errmat(MNI,MNI)

	data iumat /37/		!unit number for read/write of umat.inp, umat.out


	write(isyswr,50) fluma
50	format(1x,'umat written to file ',a13)

	open(iumat, file=fluma, Status='UNKNOWN')

	write(iumat,100) npar
100	format(1x,i3)

	do i = 1, npar
	   write(iumat,110) nexofi(i)
110	   format(1x,i3)
	enddo

	do i = 1, npar
	   do j = 1, npar
	      write(iumat,120) umat(i,j)
120	      format(1x,e20.13)
	   enddo
	enddo


c **********************************************************************************
c *** print vmat and eigenvalues for testing ***
	write(iumat,200) npar
200	format(1x,i3)

	do i = 1, npar
	   write(iumat,210) 2.d0/tvec(i)**2
210	   format(1x,e22.15)
	enddo

	do i = 1, npar
	   do j = 1, npar
	      write(iumat,220) vmat(i,j)
220	      format(1x,e22.15)
	   enddo
	enddo
c **********************************************************************************

	close(iumat)

	return
	end

