      FUNCTION NAMQCD(NNAME)
C                                                   -=-=- namqcd
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                   Find integer code corresponding to NAME. If no
C                   match,then return zero.
C
      CHARACTER NNAME*(*), NAME*8
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      COMMON /QCDPAR_CT/ FSWITCH, AL, NF, NORDER, ALPHASET
      LOGICAL ALPHASET
C              ==== MACHINE DEPENDENCE
C              ==== ASSUMES CHARACTER CODES FOR '0' TO '9' ARE CONSECUTIVE
      CHARACTER ONECH*(1)
      ONECH = '0'
      IASC0 = ICHAR(ONECH)
C
C                Use local variable to avoid problems with short
C                passed argument:
      NAME = NNAME
      NAMQCD=0
      IF ( (NAME .EQ. 'ALAM') .OR. (NAME .EQ. 'LAMB') .OR.
     >  (NAME .EQ. 'LAM') .OR. (NAME .EQ. 'LAMBDA') )
     >  NAMQCD=1
      IF ( (NAME .EQ. 'NFL') .OR. (NAME(1:3) .EQ. '#FL') .OR.
     >        (NAME .EQ. '# FL') )
     >  NAMQCD=2
      DO 10 I=1, 9                        !M: Pole quark masses
        IF (NAME .EQ. 'PLMASS'//CHAR(I+IASC0))
     >    NAMQCD=I+2
 10   CONTINUE
      DO 20 I= 0, NF
        IF (NAME .EQ. 'LAM'//CHAR(I+IASC0))
     >    NAMQCD=I+13
 20   CONTINUE
      IF (NAME(:3).EQ.'ORD' .OR. NAME(:3).EQ.'NRD') NAMQCD = 24
      IF (NAME.EQ.'FSWITCH') NAMQCD = 25  !Moved Nf switching points
      IF (NAME.EQ.'MASSTYPE') NAMQCD = 26 !Flag: input MS bar or pole mass
      DO 30 I=1, 9                        !Nf switching points
        IF (NAME .EQ. 'MHAT'//CHAR(I+IASC0))
     >    NAMQCD=I+30
 30   CONTINUE
       DO 40 I= 1, 9                      !MSBMASS: MS bar quark masses 
         IF (NAME .EQ. 'MSMASS'//CHAR(IASC0+I)) 
     >     NAMQCD=I+48         
 40    CONTINUE


      RETURN
      END
C
C***************************************************************
C
C                                                          =-=-= Qcdmsc
