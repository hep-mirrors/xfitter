      SUBROUTINE QWARN (IWRN, NWRT1, MSG, NMVAR, VARIAB,
C                                                   -=-=- qwarn
     >                  VMIN, VMAX, IACT)
 
C     Subroutine to handle warning messages.  Writes the (warning) message
C     and prints out the name and value of an offending variable to SYS$OUT
C     the first time, and to output file unit # NWRT1 in subsequent times.
C
C     The switch IACT decides whether the limits (VMIN, VMAX) are active or
C     not.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
 
      CHARACTER*(*) MSG, NMVAR

      Save Iw
 
      IW = IWRN
      VR = VARIAB
 
      If (Iw .LT. 100) Then
         WRITE (NWRT1,'(I5, 3X,A/ 1X,A,'' = '',1PD16.7)') IW, MSG,
     >                  NMVAR, VR
      Else
         WRITE (NOUT, '(1X, A/1X, A,'' = '', 1PD16.7/A,I4)')
     >      MSG, NMVAR, VR,
     >      ' !! Error # > 100 !! ; better check file unit #', NWRT1
      EndIf   

      IF  (IW .EQ. 0) THEN
         WRITE (NOUT, '(1X, A/1X, A,'' = '', 1PD16.7/A,I4)')
     >      MSG, NMVAR, VR,
     >      ' Complete set of warning messages on file unit #', NWRT1
         IF (IACT .EQ. 1) THEN
         WRITE (NOUT,'(1X,A/2(1PD15.3))')'The limits are: ', VMIN,VMAX
         WRITE (NWRT1,'(1X,A/2(1PD15.3))')'The limits are: ', VMIN,VMAX
         ENDIF
      ENDIF
 
      IWRN = IW + 1
 
      RETURN
C                         *************************
      END

