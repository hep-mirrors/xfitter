      SUBROUTINE QCDIN
C                                                   -=-=- qcdin
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                 Inputs QCD parameters.
C
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      DIMENSION VALMAS(9)
C              ==== MACHINE DEPENDENCE
C              ==== ASSUMES CHARACTER CODES FOR '0' TO '9' ARE CONSECUTIVE
C              NOTE THAT ICHAR('0') IS NON-STANDARD FORTRAN (ACCORDING TO
C              MICROSOFT).
      CHARACTER ONECH*(1)
      ONECH = '0'
      IASC0 = ICHAR(ONECH)
C
      CALL QCDGET ('LAM',ALAM,IRET1)
      CALL QCDGET ('NFL',ANF,IRET2)
      CALL QCDGET ('ORDER',ORDER,IRET3)
      NF = NINT(ANF)
      NORDER = NINT(ORDER)
 1    WRITE (NOUT, *) 'LambdaMSBAR, # Flavors, loop order ?'
      READ (NIN,*, IOSTAT = IRET) ALAM, NF, NORDER
      ORDER = NORDER
      ANF = NF
      IF (IRET .LT. 0) GOTO 22
      IF (IRET .EQ. 0) THEN
        CALL QCDSET ('LAM',ALAM,IRET1)
        CALL QCDSET ('NFL',ANF,IRET2)
        CALL QCDSET ('ORDER',ORDER,IRET3)
      ENDIF
      IF ((IRET.NE.0) .OR. (IRET1.NE.1) .OR. (IRET2.NE.1)
     >  .OR. (IRET3.NE.1)) THEN
        WRITE (NOUT, *) 'Bad value(s), try again.'
        GOTO 1
      ENDIF
      
 3    WRITE (NOUT, *) 'Input on-shell mass (0) or MSbar mass (1)?'
      READ (NIN,*, IOSTAT = IRET) aImass
      IF (IRET .LT. 0) GOTO 22
      IF (IRET .EQ. 0) 
     >  CALL QCDSET ('MASSTYPE',aImass,IRET1)
      IF ((IRET.NE.0) .OR. (IRET1.NE.1)) THEN
        WRITE (NOUT, *) 'Bad value(s), try again.'
        GOTO 3
      ENDIF

      DO I = 1, NF
        if (int(aImass).eq.0) then !set onshell mass
          CALL QCDGET('PLMASS'//CHAR(I+IASC0),VALMAS(I),IRET1)
 10       WRITE (NOUT, '(1X,A,I2,A)') 'On-shell mass of quark', I, '?'
          READ (NIN,*, IOSTAT=IRET) VALMAS(I)
          IF (IRET .LT. 0) GOTO 22
          IF (IRET .EQ. 0)
     >      CALL QCDSET('PLMASS'//CHAR(I+IASC0),VALMAS(I),IRET1)
          IF ((IRET .NE. 0) .OR. (IRET1 .NE. 1)) THEN
            WRITE (NOUT, *) 'Bad value, try again.'
            GOTO 10
          ENDIF
        else                    !set MSbar mass
          CALL QCDGET('MSMASS'//CHAR(I+IASC0),VALMAS(I),IRET1)
 12       WRITE (NOUT, '(1X,A,I2,A)') 'MSbar mass of quark', I, '?'
          READ (NIN,*, IOSTAT=IRET) VALMAS(I)
          IF (IRET .LT. 0) GOTO 22
          IF (IRET .EQ. 0)
     >      CALL QCDSET('MSMASS'//CHAR(I+IASC0),VALMAS(I),IRET1)
          IF ((IRET .NE. 0) .OR. (IRET1 .NE. 1)) THEN
            WRITE (NOUT, *) 'Bad value, try again.'
            GOTO 12
          ENDIF
        endif                   !int(aImass)
      enddo                     !i=1,Nf

 2      WRITE (NOUT, *) 'Switching factor FSWITCH?'
        READ (NIN,*, IOSTAT = IRET) Fswitch
        IF (IRET .LT. 0) GOTO 22
        IF (IRET .EQ. 0) 
     >    CALL QCDSET ('FSWITCH',FSWITCH,IRET1)
        IF ((IRET.NE.0) .OR. (IRET1.NE.1)) THEN
          WRITE (NOUT, *) 'Bad value(s), try again.'
          GOTO 2
        ENDIF
 
      RETURN
C
 22   WRITE (NOUT, *) 'END OF FILE ON INPUT'
      WRITE (NOUT, *)
      RETURN
      END
C
C *****
C
