C-------------------------------------------
      subroutine fopen(fnumber, fname)
C---------------------------------
      implicit none
      integer fnumber
      character*(*) fname
      open (fnumber,file=TRIM(fname),status='unknown')
      end

C-------------------------------------------
      subroutine fclose(fnumber)
C---------------------------------
      implicit none
      integer fnumber
      close (fnumber)
      end
