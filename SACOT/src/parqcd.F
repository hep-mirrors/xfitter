      SUBROUTINE PARQCD(IACT,NAME,VALUE,IRET)
C                                                   -=-=- parqcd
C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: Parqcd
C Description: Input-output routines for setting and reading QCD parameters
C ListOfFiles: parqcd qcdin qcdout qcdset qcdget namqcd
C ====================================================================

C #Header: /Net/d2a/wkt/1hep/2qcd/RCS/Parqcd.f,v 1.1 97/12/21 20:34:45 wkt Exp $ 
C #Log:	Parqcd.f,v $
c Revision 1.1  97/12/21  20:34:45  wkt
c Initial revision
c 

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C             Actions: 0     type list of variables on unit VALUE.
C                      1     set variable with name NAME to VALUE, if
C                             it exists, else set IRET to 0.
C                      2     find value of variable. If it does not exist,
C                             set IRET to 0.
C                      3     request values of all parameters from terminal.
C                      4     type list of all values on unit VALUE
C
C             IRET =   0     variable not found.
C                      1     successful search
C                      2     variable found, but bad value.
C                      3     bad value for IACT.
C                      4     no variable search (i.e., IACT is 0,3,or 4).
C
C             NAME is assumed upper-case.
C             if necessary, VALUE is converted to integer by NINT(VALUE)
C
      INTEGER IACT,IRET
      CHARACTER*(*) NAME
C
      IRET=1
      IF (IACT.EQ.0) THEN
         WRITE (NINT(VALUE), *)  'LAM(BDA), NFL, ORD(ER), Mi, ',
     >               '(i in 1 to 9), LAMi (i in 1 to NFL)'
         IRET=4
      ELSEIF (IACT.EQ.1) THEN
         CALL QCDSET (NAME,VALUE,IRET)
      ELSEIF (IACT.EQ.2) THEN
         CALL QCDGET (NAME,VALUE,IRET)
      ELSEIF (IACT.EQ.3) THEN
         CALL QCDIN
         IRET=4
      ELSEIF (IACT.EQ.4) THEN
         CALL QCDOUT(NINT(VALUE))
         IRET=4
      ELSE
         IRET=3
      ENDIF
 
      RETURN
      END
C
C*******************************************************************
C
