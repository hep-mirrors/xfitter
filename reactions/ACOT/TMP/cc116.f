C2345678901234567890123456789012345678901234567890123456789012345678901234567890
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C       PROGRAM CC103
       SUBROUTINE CC103()
C-----------------------------------------------------------------------------
C      Program to COMPUTE K FACTORS
C      
C      
C      03/01/08 FIO  Move S-ACOT(CHI) INTO TOT-F
C      05/03/2005  FIO update
C      02/03/2005  FIO update
C      02/05/01  FIO. 
C      23 Jan 2008: Include Short-Cut-2: 
C      
C      
C----------------------------------------------------
C----------------------------------------------------
C      Ischeme= 0  !*** Massless MS-Bar
C      Ischeme= 1  !*** Full ACOT Scheme 
C      Ischeme= 2  !*** FFS 
C      Ischeme= 3  !*** Simplified ACOT Scheme 
C      Ischeme= 4  !*** Test Full ACOT Scheme (no NLO Q)
C      Ischeme= 5  !*** LO
C      Ischeme= 6  !*** Massless LO
c      Ischeme= 7  !*** Short-cut2: ACOT w/ Massless NLO-Q
C----------------------------------------------------
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       CHARACTER HEADER*78
       Dimension F123L(4)
       PARAMETER (PI=3.14159265359)
       Common /Ischeme/ IschIN, IsetIN, IflgIN, IhadIN
       Common  / ActInt /  AERR, RERR, iActL, iActU
C-----------------------------------------------------------------------------
       DATA SBIG,  XBJ,   Q,    XMU,   ISET,HMASS, Isch, Ihad, Iflg
     >   /  98596.,0.1D0,10.0D0,10.0D0, 1,  0.938D0,  1,   1,    0 /
       DATA SCALE,  icharge, mode
     >   /  -1    ,    +1  , 0  /
C-----------------------------------------------------------------------------
C----------------------------------------------------------------------
C INITIALIZATION
C----------------------------------------------------------------------
      AERR=0.0    
      RERR=1.E-4  !*** Setting this too small will yield unstable results *** 
      iactl=2      
      iactu=2      

C----------------------------------------------------------------------
C SETUP INTEGRATION PARAMETERS
C----------------------------------------------------------------------
      WRITE(6,202)  AERR, RERR, iActL, iActU
202   FORMAT(' AERR= ',1PG14.7,' RERR= ',1PG14.7,
     >       ' IACTL= ',I3,' IACTU= ',I3)
      WRITE(6,*) ' ENTER:  AERR, RERR, iActL, iActU'
      READ (5,*)           AERR, RERR, iActL, iActU
C----------------------------------------------------------------------
C SET THE LOOP
C----------------------------------------------------------------------
      E=(SBIG-HMASS**2)/(2.0*HMASS) 
1     WRITE(6,101) XBJ,Q
      WRITE(6,102) E,ISET,HMASS,SCALE,Ihad, Iflg,Isch
101   FORMAT('   XBJ= ',1PG14.7,'     Q= ',1PG14.7)
102   FORMAT('     E= ',1PG14.7,'  ISET= ',I3, /,
     >       ' HMASS= ',1PG14.7,' SCALE= ',1PG14.7,/,
     >       ' IHAD = ',I3,' Iflg = ',I3,
     >       ' ISCH = ',I3 )
 
      WRITE(6,*) 'ENTER: XBJ,Q,E,ISET,HMASS,SCALE,ISCH,ICHARGE,Ihad'
      READ (5,*)         XBJ,Q,E,ISET,HMASS,SCALE,ISCH,ICHARGE,Ihad
C----------------------------------------------------------------------
C SET THE SCALES
C----------------------------------------------------------------------
      SBIG = 2.0*E*HMASS + HMASS*HMASS  !*** 11/11/97 RESET SBIG FOR SIG CALC
      Q2=Q*Q
      Y=Q2/(2.0*HMASS*E*XBJ)
      Wsq = HMASS**2 + Q**2 * (1./XBJ - 1.)
      W=SQRT(Wsq)

C----------------------------------------------------
C   ***  COMPUTE MU SCALE
C----------------------------------------------------
       IF(SCALE.GT.0.0) THEN
           XMU = SCALE
       ELSEIF(SCALE.EQ.-1.0) THEN
           XMU = Q
       ENDIF

C----------------------------------------------------------------------
C FUNCTION CALL
C----------------------------------------------------------------------

       IschIN=Isch
       IsetIN=Iset
       IflgIN=Iflg          !**** NOT YET USED
       IhadIN=Ihad
       Call Fcc123L(icharge,Mode, XBJ, Q,XMU, F123L)

C----------------------------------------------------
C   ***  PRINT OUT
C----------------------------------------------------

       WRITE(6 ,103) E,XBJ,Y,HMASS,ISET
     >             ,SCALE,Q,W,XMU
     >             ,( F123L(N123),N123=1,4,1)

C----------------------------------------------------
103    FORMAT(/,
     >   '---------------------------------------------------',/,
     >   ' E    =',1PG14.7,' XBJ  =',1PG14.7,' Y   =',1PG14.7,/,
     >   ' PMASS=',1PG14.7,
     >   ' ISET =',i5,/,
     >   ' SCALE=',1PG14.7,' Q    =',1PG14.7,' W    =',1PG14.7,/,
     >   ' FMU  =',1PG14.7,/,
     >   '---------------------------------------------------',/,
     >   ' TOT  (123 ): ',3(1PG14.7,1X),/,
     >   '---------------------------------------------------'
     >   )


      GOTO 1
      END

C----------------------------------------------------------------------
C2345678901234567890123456789012345678901234567890123456789012345678901234567890
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       SUBROUTINE Fcc123L(icharge,Mode, XBJ, Q,XMU, F123L)
C      Mode not used yet
C-----------------------------------------------------------------------------
C      Program to COMPUTE K FACTORS
C      
C      
C      05/07/2007  Include Z-Z, and G-Z terms
C      05/03/2005  FIO update
C      02/03/2005  FIO update
C      02/05/01  FIO. 
C      23 Jan 2008: Include Short-Cut-2: 
C      
C----------------------------------------------------
C----------------------------------------------------
C      Ischeme= 0  !*** Massless MS-Bar
C      Ischeme= 1  !*** Full ACOT Scheme 
C      Ischeme= 2  !*** FFS 
C      Ischeme= 3  !*** Simplified ACOT Scheme 
C      Ischeme= 4  !*** Test Full ACOT Scheme (no NLO Q)
C      Ischeme= 5  !*** LO
C      Ischeme= 6  !*** Massless LO
c      Ischeme= 7  !*** Short-cut2: ACOT w/ Massless NLO-Q
C----------------------------------------------------
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       CHARACTER HEADER*78
       Dimension xout123(3), XXSCL(4), XMARRAY(6), CHARGE(6)
       Dimension F123L(4)
       PARAMETER (PI=3.14159265359)
       Common /Ischeme/ Isch, Iset, Iflg, Ihad
       Common  / ActInt /  AERR, RERR, iActL, iActU
C-----------------------------------------------------------------------------
         INTEGER UPTYPE(3), DNTYPE(3)
         DIMENSION CKM(3,3) 

       DATA  (UPTYPE(I), I=1,3,1) /1,4,6/
       DATA  (DNTYPE(I), I=1,3,1) /2,3,5/
       DATA ((CKM(I,J), J=1,3,1), I=1,3,1)    !*** ABS OF CKM MATRIX
     1   /0.975268, 0.2209986, 0.00350000,
     2    0.2209541, 0.974422, 0.0409997, 
     3    0.00565041, 0.0407591, 0.999153/

C-----------------------------------------------------------------------------
C--- FOR CHARGED CURRENT
       DATA GLQ,GRQ   /  1.0,0.0/
       DATA IPRINT /0/  !*** NO debugging
C       DATA IPRINT /1/ !*** for debugging
C-----------------------------------------------------------------------------
C                          U      D      S      C      B      T
       DATA XMARRAY  /    0.1,   0.1,   0.2,   1.6,   5.0, 175.0/
       DATA CHARGE   /   +2.0,  -1.0,  -1.0,  +2.0,  -1.0,  +2.0/

C      common /fred/ xmc,xmb,Hmass !*** PULL VALUES FROM QCDNUM  fio 14 FEB. 2011
       common /fred/ xmc,xmb,Hmass !*** PULL FROM XFITTER
C----------------------------------------------------------------------
C PULL MC AND MB FROM QCDNUM USING 
C      common /fred/ xmc,xmb
C----------------------------------------------------------------------
       XMARRAY(4)=XMC           !*** PULL VALUES FROM QCDNUM  fio 14 FEB. 2011
       XMARRAY(5)=XMB           !*** PULL VALUES FROM QCDNUM  fio 14 FEB. 2011
C----------------------------------------------------------------------
C ONLY IMODE=1 TOTAL is implemented at present. NOT c and  b yet.
C----------------------------------------------------------------------
       if(mode.ne.1) then
          write(6,*) 
     >   ' error: for CC only imode=1 F123L-tot is implemented '
          stop
       endif

C----------------------------------------------------------------------
C INITIALIZATION
C----------------------------------------------------------------------
       DO I=1,4,1
          F123L( I)  = 0.0d0
       ENDDO
C----------------------------------------------------------------------
C CHOOSE W+ or W- SCATTERING
C----------------------------------------------------------------------
      IUPDOWN  =-1    !*** W- SCATTERING
      IUPDOWN  =+1    !*** W+ SCATTERING
      IUPDOWN = ICHARGE
      UPDOWN=float(ICHARGE)  !*** CHANGE TO FLOAT
C----------------------------------------------------
C   *************************************************
C----------------------------------------------------
C----------------------------------------------------
C   ***  LOOP OVER PARTON FLAVORS: 
C----------------------------------------------------
      DO IQUARKIN  =1,3,1
      DO IQUARKOUT =1,3,1

C----------------------------------------------------
      IF(UPDOWN.GE.0)  THEN  !*** THIS IS FOR DOWN -> UP
           IPARTIN = +DNTYPE(IQUARKIN )
           IPARTOUT= +UPTYPE(IQUARKOUT)

           IUP= IQUARKOUT   !***  USED TO PICK UP CORRECT CKM TERM
           IDN= IQUARKIN

      ELSE                 !*** THIS IS FOR UP -> DOWN
           IPARTIN = +UPTYPE(IQUARKIN )
           IPARTOUT= +DNTYPE(IQUARKOUT)

           IUP= IQUARKIN    !***  USED TO PICK UP CORRECT CKM TERM
           IDN= IQUARKOUT

      ENDIF
      COUPLING = CKM(IUP,IDN)**2
C----------------------------------------------------------------------
C SETUP PARTON MASSES: 
C----------------------------------------------------------------------
      IF(IPRINT.NE.0) WRITE(6,*) ' PARTON ',IPARTIN,' => ',IPARTOUT

C----------------------------------------------------------------------
C SETUP PARTON MASSES: 
C----------------------------------------------------------------------
       F1M=XMARRAY(ABS(IPARTIN))
       F2M=XMARRAY(ABS(IPARTOUT))

C----------------------------------------------------
C NOTE: ASSUME TOP IS ZERO (NOT IMPLEMENTED IN ALL PDF SETS)
C----------------------------------------------------------------------
      if((abs(IPARTIN).NE.6).and.(abs(IPARTOUT).NE.6)) then
         call TOTF(XBJ,Q,F1M,F2M,GLQ,GRQ,GLQ,GRQ,IPARTIN,IPARTOUT,
     >   XMU,ISET,Ihad,HMASS,Isch,Iflg,Icharge,Xout123) 
      else
           DO I=1,3,1
             Xout123( I)  =0.0
           ENDDO
      endif
C----------------------------------------------------

C----------------------------------------------------
C   ***  ADD PARTON CONTRIBUTION TO TOTAL STRUCTURE FUNCTIONS
C----------------------------------------------------
           DO I=1,3,1
cv             F123L( I)  = F123L( I)  + COUPLING * Xout123( I) 
             F123L( I)  = F123L( I)  + COUPLING/2.d0 * Xout123( I) 
           ENDDO

C----------------------------------------------------
C   ***  PICK OUT STRANGE->CHARM AND SAVE THIS
C----------------------------------------------------
       IF((ABS(IPARTIN).EQ.3).AND.(ABS(IPARTOUT).EQ.4)) THEN    
           DO I=1,3,1
             XXSCL(I)  = XXSCL(I) + COUPLING * Xout123(I)    
C                  *** CAREFUL TO PICK UP ONLY STRANGE->CHARM
           ENDDO
        ENDIF
C----------------------------------------------------
C----------------------------------------------------
      ENDDO  !***  IQUARKOUT
      ENDDO  !***  IQUARKIN
C----------------------------------------------------
C   ***  END LOOP OVER PARTON FLAVORS: 
C----------------------------------------------------

C----------------------------------------------------------------------
C COMPUTE  FL 
C----------------------------------------------------------------------
      rho=Sqrt(1.0d0+(2.0d0*hmass*xbj/Q)**2)  !*** Get Hmass from /fred/ common block 
         F123L(4)=rho**2*F123L(2)- 2.0d0*xbj*F123L(1)
         XXSCL(4)=rho**2*XXSCL(2)- 2.0d0*xbj*XXSCL(1)

C----------------------------------------------------------------------

      Return
      END

C----------------------------------------------------------------------
