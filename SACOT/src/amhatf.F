      FUNCTION AMHATF(I)
C                                                   -=-=- amhatf
      implicit none
C                Returns the boundary in mass scale between the regions
C                with I & (I-1) effective "light" quarks.
C
      double precision amhatf, alam, amhat, amn
      integer nhq, i, ii
      logical AmHatSet
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet
      integer  NIN, NOUT, NWRT
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      if (.not. AmhatSet) call SetAmhatF
      ii=abs(i)
      IF ((II.LT.0).OR.(II.GT.9)) THEN
         WRITE (NOUT,*) 'I IS OUT OF RANGE IN AMHATF'
         AMHATF = 0
      ELSE
         AMHATF = AMHAT(II)
      ENDIF
      RETURN
      END
C
C********************************************************************
C
