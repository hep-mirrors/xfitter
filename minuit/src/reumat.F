c --------------------------------------------------------------------------------
	subroutine reumat
c --------------------------------------------------------------------------------
c read matrix umat(i,j) from file named umat.inp
c if no such file, just leaves umat alone, which is identity matrix.
c 
c format of the file is same as that written by wrumat: 
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

	common /stepco/ dstep(MNI)
	common /umatco/ hquad(MNI,MNI), umat(MNI,MNI), vmat(MNI,MNI), 
     &	                tvec(MNI), errmat(MNI,MNI)

	data iumat /37/		!unit number for read/write of umat.inp, umat.out

	open(iumat, file='umat.inp', Status='OLD', ERR=200)

	read(iumat,*,ERR=400) np
	if(np .ne. npar) then
	   write(isyswr,100) np, npar
100	   format(1x,'reumat:  fatal error --',
     &	   '  umat.inp has inconsistent npar',i5,1x,i5)
	   goto 500
	endif


	do i = 1, npar
	   read(iumat,*,ERR=400) ntmp
	   if(ntmp .ne. nexofi(i)) then
	      write(isyswr,110) i, ntmp, nexofi(i)
110	      format(1x,'reumat:  fatal error --',
     &	      '  umat.inp has inconsistent nexofi(',i3,')',
     &	             i5,1x,i5)
	      goto 500
	   endif
	enddo

	do i = 1, npar
	   do j = 1, npar
	      read(iumat,*,ERR=400) umat(i,j)
	   enddo
	enddo
	close(iumat)

	write(isyswr,120) 
120	format(1x,'umat(i,j) has been read in from umat.inp')

c set step sizes to current ERRDEF...
	step = sqrt(UP)
	do i = 1, npar
	   dstep(i) = step
	enddo

	write(isyswr,130) UP, step
130	format(1x,'step sizes set to sqrt(UP)= sqrt(',
     &	    g15.8,') =',g15.8)

	goto 999 

200	write(isyswr,202)
202	format(1x,'file umat.inp doesn''t exist:',
     &	 ' leave umat(i,j) = identity matrix')
	goto 999

400	write(isyswr,410)
410	format(1x,'reumat:  fatal error reading umat.inp')

500	continue
	close(iumat)
	stop

999	continue
	return
	end

