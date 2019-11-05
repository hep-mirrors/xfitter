c ---------------------------------------------------------
      SUBROUTINE MYSTUFF(FCN,FUTIL)
c ---------------------------------------------------------
c Switchyard to select user routines that are invoked 
c by the new minuit command 
c           MYStuff nn
c The thousands digit of nn determines which routine is called.
c The other digits can be used as a control parameter for the 
c routine that is called.
c
#include "./d506dp.inc"
#include "./d506cm.inc"
	EXTERNAL FCN,FUTIL

c get the requested number in form of an integer...
	ncode = WORD7(1) + 0.01

	itype = ncode/1000

c --------------------------------------------------
      write(isyswr,100) ncode
100   format(1x,'mystuff called with argument',i6)
c --------------------------------------------------


	if(itype .eq. 1) then
	   call xrange(fcn,futil)	!called if ncode = 1000 to 1999
	   return

	elseif(itype .eq. 2) then
	   call evplot(fcn,futil)	!called if ncode = 2000 to 2999
	   return

	else
	   write(isyswr,110) ncode
110	   format(1x,'MYSTUFF: do nothing --- undefined ncode=',i9)
	endif

	return
	end
