      SUBROUTINE QCDSET (NAME,VALUE,IRET)
C                                                   -=-=- qcdset
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  Assign the variable whose name is specified by NAME the value VALUE
C             IRET=  0      variable not found.
C                    1      success.
C                    2      variable found, but bad value.
C             NAME is assumed upper-case and VALUE is real.
C             If necessary, VALUE is converted to integer by NINT(VALUE).
C
      CHARACTER*(*) NAME
      Double precision POLEMASS, MSBMASS
      Integer iMass
      Logical MassesSet
      COMMON / COMQMS/ POLEMASS(9), MSBMASS(9), iMass, MassesSet
      LOGICAL AlphaSet, AmHatSet
      COMMON / QCDPAR_CT / FSWITCH, AL, NF, NORDER, AlphaSet
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet

      PARAMETER (PI=3.1415927, EULER=0.57721566)
C
      IVALUE = NINT(VALUE)
      ICODE  = NAMQCD(NAME)

      IF (ICODE .EQ. 0) THEN
        IRET=0; return
      endif                     !icode .eq. 0
      
      IRET = 1
      AlphaSet = .FALSE.

      IF (ICODE .EQ. 1) THEN
        IF (VALUE.LE.0) GOTO 12
        AL=VALUE
      ELSEIF (ICODE .EQ. 2) THEN
        IF ( (IVALUE .LT. 0) .OR. (IVALUE .GT. 9)) GOTO 12
        NF = IVALUE
      ELSEIF ((ICODE .GE. 3) .AND. (ICODE .LE. 11))  THEN
        IF (VALUE .LT. 0) GOTO 12 
        if (iMass .ne. 0) then  !set on-shell mass, only if iMass=0
          print *, 'STOP: iMass=',iMass
          print *,'Can only change MSbar masses, not pole masses'
          stop
        endif                   !iMass .ne. 0
        POLEMASS(ICODE - 2) = VALUE 
c Masses, switching points, alpha_s are not synchronized
        MassesSet=.false.; AmHatSet=.false.
      ELSEIF ((ICODE .GE. 13) .AND. (ICODE .LE. 13+NF))  THEN
        IF (VALUE .LE. 0) GOTO 12
        CALL SETL1 (ICODE-13, VALUE)
      ELSEIF (ICODE .EQ. 24)  THEN
        IF ((IVALUE .LT. 1) .OR. (IVALUE .GT. 3)) GOTO 12
        NORDER = IVALUE
      ELSEIF (ICODE .EQ. 25)  THEN
        FSWITCH=Value           !pn10 set the switching point factor
        AmHatSet=.false.
      ELSEIF (ICODE .EQ. 26)  THEN
        iMass=nint(Value)          !choice between input on-shell or MSbar mass 
c Masses, switching points, alpha_s are not synchronized
        MassesSet=.false.; AmHatSet=.false.
      ELSEIF (ICODE .ge. 30 .and. ICODE .le. 39)  THEN
        stop 'AMHAT must be set by setting the quark mass'
cmarco 
      ELSEIF ((ICODE .GE. 49).AND. (ICODE .LE. 57)) THEN
        IF (VALUE .LT. 0) GOTO 12 
        if (iMass .ne. 1) then  !set MSbar mass, only if iMass=1
          print *, 'STOP: iMass=',iMass
          print *,'Can only change pole masses, not MSbar masses'
          stop
        endif                   !iMass .ne. 1
        MSBMASS(ICODE - 48) = VALUE   
c Masses, switching points, alpha_s are not synchronized
        MassesSet=.false.; AmHatSet=.false.
      ENDIF                     !icode .eq. 1

      RETURN
C
C              Illegal value
 12   IRET=2
      RETURN
      END
C
C************************************************
C
