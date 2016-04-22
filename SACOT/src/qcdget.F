      SUBROUTINE QCDGET(NAME,VALUE,IRET)
C                                                   -=-=- qcdget
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  Sets VALUE to the value of variable named NAME.
C            IRET=  0          variable not found.
C                   1          success.
C
C            NAME is assumed to be an upper-case character variable.
C
      CHARACTER*(*) NAME
      logical set, AmHatSet
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet
      COMMON / QCDPAR_CT / FSWITCH, AL, NF, NORDER, SET

      Double precision POLEMASS, MSBMASS
      Integer iMass
      Logical MassesSet
      COMMON / COMQMS/ POLEMASS(9), MSBMASS(9), iMass, MassesSet

      PARAMETER (PI=3.1415927, EULER=0.57721566)
C
      ICODE = NAMQCD(NAME)
      IRET = 1
      IF (ICODE .EQ. 1) THEN
         VALUE = AL
      ELSEIF (ICODE .EQ. 2) THEN
         VALUE = NF
      ELSEIF ((ICODE .GE. 3) .AND. (ICODE .LE. 11))  THEN
        if (.not.MassesSet) goto 15 !need mass synchronization
         VALUE = POLEMASS(ICODE - 2) !quark masses
      ELSEIF ((ICODE .GE. 13) .AND. (ICODE .LE. 13+NF))  THEN
         VALUE = ALAM(ICODE - 13)
      ELSEIF (ICODE .EQ. 24) THEN
         VALUE = NORDER
       ELSEIF (ICODE .EQ. 25)  THEN
         Value = FSWITCH !pn10 set the switching point factor
       ELSEIF (ICODE .EQ. 26)  THEN
         Value = iMass !Input mass flag
       ELSEIF (ICODE .ge. 30 .and. ICODE .le. 39)  THEN
         if (.not. AmhatSet) goto 16
         Value = AMHAT(icode-30) !pn11 switching points
      ELSEIF ((ICODE .GE. 49).AND. (ICODE .LE. 57)) THEN
        if (.not.MassesSet) goto 15 !need mass synchronization
        VALUE =MSBMASS(ICODE-48)     !MSbar mass
      ELSE
         IRET=0
      ENDIF

      return

 15   stop 'Quark masses are not synchronized; run MassSync to proceed'
 16   stop 'Switching points are not set; run SetAmhatF to proceed'

      END
C
C *****
C
