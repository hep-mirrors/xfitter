C     FUNCTION SD2(P,F,REGION)   COMPUTES INTEGRAL
C
C                                      B        2
C                                SD2 = $ F(X) DX
C                                      A

C  WITH ABSOLUTE PRECISION  P.
C  FUNCTION F  AND SUBROUTINE REGION   SHOULD BE DECLARED AS EXTERNAL.
C  LIMITS A(L),B(L), L<=2, AS FUNCTIONS OF X(L-1),...,X(1) SHOULD BE COMPUTED
C  IN THE SUBPOUTINE REGION(L,X,A,B)

C  OPTIONS --------------------------------------------------------------

C  OUTPUT: LEVEL OF PROTOCOL OUTPUT, ESTIMATE OF ROOT-MEAN-SQUARE ERROR,
C          INSTANT APPER LIMITS AND VALUES OF INTEGRALS, NUMBRS OF QUADRATURES,
C          TOTAL INSTANT VALUES OF INSTNTEGRALS, NUMBER OF LINES OUTPUT, NUMBER
C          OF POINTS USED:
C /SD2OUT/LETOUT,NLINES,NQUAD(2),NPOINT(3)
C       ,RMSERR(2),CLIMIT(2),SUMTOC(2),TOT(2),EXTRA1(2),EXTRA2(2)
C         IF LETOUT>0, PROGRESS REPORT OF LEVEL 'LETOUT' IS OUTPUT ON CHANNEL 1
C  INPUT: LIMITS ON SEGMEMT AND GAP SIZES:/SD2GAP/ STEPMX(2),GAPMAX(2),GAPMIN(2)
C       : 10 HIT POINTS BETWEEN A AND B, WHERE INTEGRAL FROM A TO HIT SHOULD
C         BE COMPUTED: /SD2HIT/NHITS(2),HIT(2,10),YHIT(2,10)

C ONE-STEP MODE: IF P=0, 1 STEP IN EVERY DITECTION IS MADE, USING 8-INTERVAL
C                       QUADRATURE FORMULA

      real*8 FUNCTION SD2(P,F,REGION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /SD2GAP/STEPMX(2),GAPMAX(2),GAPMIN(2)
     *       /SD2HIT/NHITS(2),HIT(2,10),YHIT(2,10)
     *       /SD2OUT/ RMSERR(2),CLIMIT(2),SUMTOC(2),TOT(2),RQA(2),RQB(2)
     *              ,LETOUT,NLINES,NQUAD(2),NPOINT(3)

C --------------------- ARRAYS OF INTERNAL USE -------------------------------
      COMMON /SD214V/
     * X1(2),G(2),PA(2),ESTER,DS, GA(2,2),GB(2,2),Y(2,21)
     *   ,YI(2,16),YT(2,19),errold,L,N(2),INS(2,16),MADE(2,5)
     *                 ,NOLD,MOLD,MODY

      DIMENSION Z(2,200),DW(2),YW(2),
     *          XACCUM(5),YACCUM(5),JQ(2)

      LOGICAL NOMEM(2)

      DIMENSION SUM(3),IN(2),LB(2),LINS(2),M(2),LOST(2),IHIT(2)
     *,SQUERR(2),KIND(2),A(2),B(2),D(2),P03(2),P200(2),X2(2)
      DIMENSION  LEAP(2),LEAP4(2),E1(2),E11(2),PR(2),PR1(2)
     *,LOT(2),LOT0(2),LOT1(2)    ,IODD(2),DS0(2)
      DIMENSION X(2),SGMAX(2),GMIN(2),GMAX(2)
      DIMENSION  W(81),W1(20),W2(29),W3(32)
      EQUIVALENCE (W(1),W1(1)),(W(21),W2(1)),(W(50),W3(1))
      DATA ACCURA/1.E-8/
C     ACCURA: RELATIVE ACCURACY OF FLOATING POINT NUMBERS MULTIPLIED BY 10
      DATA W1/
     *  1.0,
     1  0.5,
     2  0.1666666666666667,  0.3333333333333333,
     3  0.1250000000000000,  0.3750000000000000,
     4  0.0777777777777778,  0.3555555555555556,  0.0666666666666666,
     5  0.0659722222222222,  0.2604166666666667,  0.1736111111111111,
     6  0.0488095238095238,  0.2571428571428571,  0.0321428571428571,
     *  0.1619047619047619, 
     7  0.0434606481481481,  0.2070023148148148,  0.0765625000000000,
     *  0.1729745370370370/
       DATA W2/
     8  0.0348853615520282,  0.2076895943562610, -0.0327336860670194,
     *  0.3702292768959436, -0.0800705467372134,
     9  0.0318861607142857,  0.1756808035714286,  0.0120535714285714,
     *  0.2158928571428571,  0.0644866071428571,
     1  0.0268341483619261,  0.1775359414248303, -0.0810435706269040,
     *  0.4549462882796216, -0.4351551226551227,  0.3568823152156486,
     1  0.0249332309119636,  0.1548553585207231, -0.0371692317937978,
     1  0.2896582547949735, -0.1101780891754850,  0.1779004767416226,
     1  0.0216394874966304,  0.1570361067503925, -0.1203219637505352,
     2  0.5664988979274694, -0.8165056372199229,  1.3877596689025260,
     * -0.6961065601065601/
       DATA W3/
     1  0.0203347191051236,  0.1398760852657854, -0.0777118702884142,
     3  0.3878961542438828, -0.3769238163321176,  0.5136761795561554,
     * -0.1071474515504153,
     1  0.0180344712157984,  0.1420877946927330, -0.1540253470523532,
     4  0.6997489104402685, -1.3239976056465254,  2.5240777544357791,
     * -3.3578644895056932,  1.9519385114199930,
     1  0.0170872997716259,  0.1285073786774052, -0.1127229050595052,
     5  0.5070427082102236, -0.7562931148441402,  1.1913603495067990,
     * -0.9680052114962480,  0.4930234952338398,
     1  0.0153989471166495,  0.1306411914401402, -0.1839764466493380,
     6  0.8518689889045043, -1.9750740358891965,  4.2762649967201477,
     * -6.9673071185989400,  9.5901711055393664, -5.2379876285833335/

C ENTRY
      ND=2
      SD2=0d0

      DO 2 I=1,ND
    2 NPOINT(I)=0
      NACCUM=0
      NLINES=0
      PA(1)=ABS(P)*1.5
      L=0
          GO TO 30
C               30: RECURSIVE CALL OF SDX-INS
C---------------------------------------------------------------------------
C FUNCTION INSERT (L,IN,F): INSERTS POINT INTO THE GAP 'IN'

   10 IF(INS(L,IN(L)).EQ.1) GO TO 14
      INS(L,IN(L))=1
      X(L)=X1(L)+(IN(L)-0.5)*G(L)
          LB(L)=1
             GO TO 30
   12 YI(L,IN(L))=SUM(L+1)
   14 LI=LINS(L)
      GO TO(1001,1002,1003,1004,1005,1006,1007,1008,1009,1010
     *     ,1011,1012),LI
C--------------------------------------------------------------------------

C RECURSIVE CALL OF SDX-INS

   30 L=L+1

      IF(L.LE.ND)GO TO 60
C                      60: START OF INTEGRATION OF L-LEVEL
C-------------- L=ND+1 ------------
      SUM(L)=F(X)
C L-TH RETURN --------------------
   40 NPOINT(L)=NPOINT(L)+1
      IF(L-1.GT.LETOUT)GO TO 49
      IF(L.LE.LETOUT)GO TO 44
      NACCUM=NACCUM+1
      IF(L.GT.1)                XACCUM(NACCUM)=X(L-1)
      YACCUM(NACCUM)=SUM(L)
      IF(NACCUM.LT.5)GO TO 49
      WRITE(1,5040)(XACCUM(I),YACCUM(I),I=1,5)
 5040 FORMAT(5(E16.6,E15.6))
      NLINES=NLINES+1
      NACCUM=0
      GO TO 49
   44 IF(L.EQ.LETOUT.AND.NACCUM.GT.0)WRITE(1,5041)
     *                            (XACCUM(I),YACCUM(I),I=1,NACCUM)
      NACCUM=0
 5041 FORMAT(5(E16.6,E15.6))
      LO=L-1
      XL1=0d0
      IF(L.GT.1)XL1=X(LO)
      WRITE(1,5042)L,LO,XL1,L,SUM(L),RMSERR(L),NQUAD(L),
     *                           (NPOINT(I),I=1,ND+1)
 5042 FORMAT(' L=',I2,'  X(',I1,')=',E17.9,'  INTEGRAL(',I2,')=',E17.9
     *,'  ER=',E10.2,' NQ',I8,' N',9I9)
      NLINES=NLINES+1
   49 L=L-1
      IF(L.EQ.0)GO TO 50
      LBL=LB(L)
      GO TO     (12,90,116),LBL
C EXIT    -----------------------


   50 SD2=SUM(1)
      RETURN
C========================================================================

C START
C      CHECKS  CONSISTENCY OF STEP AND GAP LIMITS
C      MAKES  B-A>=SGMAX>=GMAX>=GMIN
C      NULLIFIES  SDX ETC.
C      SETS   X1=A, Y(3)=F(A)
C ZERO
   60 SUM(L)=0d0
      TOT(L)=0d0
      M(L)=0
      LOST(L)=0
      SQUERR(L)=0d0
      IHIT(L)=0
      NQUAD(L)=0
      KIND(L)=0
         JQ(L)=0
      CALL REGION(L,X,A,B)
      IF(L.GT.LETOUT)GO TO 70
      NLINES=NLINES+4
      IF(L.GT.1)GO TO 61
      WRITE(1,5060)A(1),B(1)
 5060 FORMAT(' 1  ++++++++++++++++      A1 B1',2E20.9/)
      GO TO 70
   61 LO=L-1
      IF(JQ(LO).EQ.0)WRITE(1,5061)LO,NQUAD(LO),X1(LO),SUM(LO),TOT(LO)
     *,RMSERR(LO)
      IF(JQ(LO).EQ.1)WRITE(1,5061)LO,NQUAD(LO),X1(LO),SUM(LO),TOT(LO)
     *,RMSERR(LO),RQA(LO),RQB(LO)
         JQ(LO)=0
 5061 FORMAT(' STATE OF L=',I1,' : NQ',I3
     *,' INTEGRAL (FROM A TO C=',E18.9,') =',E18.9
     *,'  TOTAL',E17.9,' ERR',E10.2,' EXTRA',2E14.6)
      GO TO(62,63,64,65,66,67,68),LO
   62 WRITE(1,5062)
 5062 FORMAT(/20X,' 2  +++++++++++++++++')
      GO TO 69
   63 WRITE(1,5063)
 5063 FORMAT(/40X,' 3  +++++++++++++++++')
      GO TO 69
   64 WRITE(1,5064)
 5064 FORMAT(/60X,' 4  +++++++++++++++++')
      GO TO 69
   65 WRITE(1,5065)
 5065 FORMAT(/80X,' 5  ++++++++++++++++')
      GO TO 69
   66 WRITE(1,5066)
 5066 FORMAT(/100X,' 6  ++++++++++++++++')
      GO TO 69
   67 WRITE(1,5067)
 5067 FORMAT(/120X,' 7  ++++++++++++++++')
      GO TO 69
   68 WRITE(1,5068)
 5068 FORMAT(/140X,' 8  ++++++++++++++++')

   69 WRITE(1,5069)LO,X(LO),L,A(L),L,B(L)
 5069 FORMAT('  X(',I1,')=',E17.9,'     A(',I1,')=',E17.9,
     *                              '  B(',I1,')=',E17.9)

   70       IF(B(L).EQ.A(L))GO TO 40
C                                 40: L-TH RETURN
C LIMITS
      D(L)=B(L)-A(L)
      AD=ABS(D(L))
           R=ABS(STEPMX(L)/D(L))
      SGMAX(L)=AD*R
              IF(R.LT.ACCURA.OR.R.GT.1.)SGMAX(L)=AD
           R=ABS(GAPMAX(L)/AD)
      GMAX(L)=SGMAX(L)*R
              IF(R.LT.ACCURA.OR.R.GT.1.)GMAX(L)=SGMAX(L)
           R=ABS(GAPMIN(L)/GMAX(L))
      GMIN(L)=GMAX(L)*R
              IF(R.LT.ACCURA.OR.R.GT.1.)GMIN(L)=GMAX(L)*ACCURA
          D1=ABS(A(L))+ABS(B(L))
          GMIN4=D(L)*GMIN(L)/(4.*AD)
      IF(A(L)+GMIN4.EQ.A(L).OR.B(L)-GMIN4.EQ.B(L))GMIN(L)=D1*ACCURA
      IF(P.GT.0d0)GO TO 71


C P=0 CASE
      GMIN(L)=SGMAX(L)/8

C PRECISION
   71 IF(L.GT.1) PA(L)= PA(L-1)/ABS( D(L-1 ) )
      P03(L)=PA(L)*0.3d0
      P200(L)=PA(L)*200d0

C ERASE
      DO 72 I=1,2
      GA(L,I)=0d0
   72 GB(L,I)=0d0
C PUT
      X1(L)=A(L)
C----------------------Y(3)=F(A)
      X(L)=A(L)
      LB(L)=2
             GO TO 30
   90 Y(L,3)=SUM(L+1)
C111111111111111111111111111111111111111111111111111111111

C SET                   GO   TO 140:OLD, 120:RECALL, 130:JUMP
  100 CLIMIT(L)=X1(L)
      NOLD=N(L)
      IF(KIND(L).NE.0)GO TO 140
C --------------------------------------------------------------------
C                             STEP WITH NEW BEGINNING X1

      IF(X1(L).EQ.A(L))GO TO 102
C                           PUT LEFT WING BY SUBROUTINE WING(KIND,LOT,...)
      CALL WING2(2,0,DW,YW)
      Y(L,2)=YW(1)
      Y(L,1)=YW(2)
      GA(L,1)=DW(1)
      GB(L,1)=DW(2)-DW(1)
      Y(L,3)=Y(L,3+N(L))
  102 IF(LOST(L).EQ.0.AND.M(L).NE.0)GOTO 120
C                                  120: RECALL
C ---------------------------------------------------------------------


C NEW  .......                        (CASES M=0, LOST=0 AND M>0, LOST>0)

C GAP SIZE AND RIGHT WING SETTING
  110 IF(LOST(L).GT.0)GOTO 118
C          NOTHING IS LOST. GAP SIZE -   IN MOVE TOWARD B (LOST=0,M=0)
          G(L)=SGMAX(L)
          G1=B(L)-X1(L)
            IF(G1/G(L).LT.1.)G(L)=G1

C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    RETURN TO 40

          IF(((X1(L)+0.5*G(L))-X1(L))/D(L).EQ.0d0) GO TO 40

          IF(IHIT(L).GE.NHITS(L))GO TO 114
             HI=HIT(L,IHIT(L)+1)
                IF(HI.EQ.X1(L))GO TO 112
             IF((HI-X1(L))*(X1(L)+G(L)-HI).LT.0d0)GO TO 114
          G(L)=HI-X1(L)
  112 YHIT(L,IHIT(L))=SUM(L)
          IHIT(L)=IHIT(L)+1
C                                             WING-2
  114     X2(L)=X1(L)+G(L)
C---------------------- Y(4)=F(X2)
          X(L)=X2(L)
          LB(L)=3
                 GO TO 30
  116 Y(L,4)=SUM(L+1)
          GA(L,2)=0d0
            GOTO 119
C                      GAP SIZE -  IN MOVE TOWARD X2. CASE LOST>0 (HENCE M>0)
  118     X2(L)=Z(L,M(L)-2)
          G(L) =X2(L)-X1(L)
          IF(((X1(L)+0.5*G(L))-X1(L))/D(L).EQ.0d0)GO TO 120
C                        120: RECALL             WING-2
          N0=Z(L,M(L)-1)+8.1
                       IF(Z(L,M(L)).GT.0d0)N0=7
          Y(L,4)=Z(L,M(L)-N0)
          Y(L,5)=Z(L,M(L)-N0+1)
          Y(L,6)=Z(L,M(L)-N0+2)

          GA(L,2)=Z(L,M(L)-4)
          GB(L,2)=Z(L,M(L)-3)
C                       SET N
  119     N(L)=1
          INS(L,1)=0
           GO TO 190
C                190: PUT YT FOR TRY

C RECALL                FROM MEMORY Z(200)
  120 LOST(L)=0
      IF(Z(L,M(L)).GT.0.5) GO TO 130
C                           130: JUMP
C --------------------- TAKE FROM MEMORY:

      INSN=-Z(L,M(L))+0.1
      N(L)=Z(L,M(L)-1)+0.1
      X1(L)=Z(L,M(L)-2)
      GB(L,2)=Z(L,M(L)-3)
      GA(L,2)=Z(L,M(L)-4)
      G(L)=Z(L,M(L)-5)
C ---------------------- NODES
      NL=N(L)
      DO 121 I=1,NL+2
  121 Y(L,3+I)=Z(L,M(L)-NL-8+I)
C                        ERASE INS
      DO 122 I=1,NL
  122 INS(L,I)=0

      M(L)=M(L)-NL-9-2*INSN
C                        INSERTED
      IF(INSN.EQ.0)GOTO 190
          DO 124  I=1,INSN
             II=Z(L,M(L)+I)+0.1
          YI(L,II)=Z(L,M(L)+INSN+I)
  124  INS(L,II)=1
      GO TO 190
C           190: MAKE YT FOR TRY

C JUMP ============================ JUMP  OVER THE HOLE
  130 X1(L)=Z(L,M(L)-1)
C                         SUM SDX WITH INTEGRAL OVER THE HOLE
      SUM(L)=SUM(L)+Z(L,M(L)-13)
      GB(L,1)=Z(L,M(L)-12)
      GA(L,1)=Z(L,M(L)-11)
       Y(L,1)=Z(L,M(L)-10)
       Y(L,2)=Z(L,M(L)-9)
       Y(L,3)=Z(L,M(L)-8)
        M(L)  =M(L)-14
C                         GO TO 110:NEW, 120:RECALL
      IF(M(L).EQ.0)GO TO 110
                     GO TO 120
C------------------------------------------------------------------------
C          PROCEDURE      OLD(LOT,KIND)
C
C    RECOINS OLD SEGMENT, OR ITS PART, INTO NEW SEGMENT

C   INPUT :
C           N: OF OLD SEGMENT
C         LOT: NUMBER OF GAPS OF OLD SEGMENT THAT WAS USED (SUMMED OR STORED)
C        KIND: FIXED END OF THE PART THAT WAS USED

C       USES  ARRAYS Y,YI,YNAIL  AND INDICATORS AND GAPS RELATED
C             YT(19): WORKING ARRAY OF MERGE
C   VARIABLES:

C               NOUT: NUMB. OF GAPS TAKEN OUT
C               NEW : N OF NEW SEGMENT
C               KO  : OPPOSIT KIND (KIND OF PART TAKEN OUT)

C              FIND KIND,LOT OF RECOINED PART
  140 KO=KIND(L)
      KIND(L)=3-KO
      NOUT=LOT(L)
      LOT(L)=N(L)-NOUT
      IF(NOUT.GT.0)CALL WING2(KIND(L),LOT(L),DW,YW)
      IF(KIND(L).EQ.1)GO TO 144
C SHIFT                      OF  RIGHT PART
      LOTC=LOT(L)
      DO 141 I=1,LOTC+3
  141 Y(L,2+I)=Y(L,2+NOUT+I)
      DO 142 I=1,LOTC
      YI(L,I)=YI(L,I+NOUT)
  142 INS(L,I)=INS(L,I+NOUT)
C                       PLACE WING
      Y(L,2)=YW(1)
      Y(L,1)=YW(2)
      GA(L,1)=DW(1)
      GB(L,1)=DW(2)-DW(1)
      N(L)=LOTC
           GO TO 190
C                190: MAKE YT FOR TRY

C MERGE ----------------  ARRAYS Y AND YI INTO YT
  144 NEW=LOT(L)+LOT(L)
      YT(L,1)=Y(L,3)
      LOTC=LOT(L)
      DO 145 I=1,LOTC
                   YT(L,I+I)=YI(L,I)
  145              YT(L,1+I+I)=Y(L,3+I)
C                                   ERASE  INS
      DO 146 I=1,NEW
  146 INS(L,I)=0

C -------------------------------------------
C             FIXED END
  150 IF (NOUT.EQ.0) GO TO 176

C     ----------------------------------------------
C          PLACE WINGS
C                                PLACE GAPS AND ERASE NAIL
      GA(L,2)=DW(1)
      GB(L,2)=DW(2)-DW(1)
C     ----------------------------------------
C                        PLACE (FREE) RIGHT WING
      YT(L,2+NEW)=YW(1)
      Y(L,5+NEW)=YW(2)
      GO TO 180
C           180:MIDDLE
C     ----------------------------------------
C                        PLACE (FIXED) RIGHT WING
  176 YT(L,2+NEW)=Y(L,4+N(L))
      Y(L,5+NEW)=Y(L,5+N(L))
C     ----------------------------------------

C MIDDLE                 PLACE MIDDLE OF SEGMENT
  180 DO 182 I=1,NEW+2
  182 Y(L,2+I)=YT(L,I)
C             SET NEW N AND G
      G(L)=G(L)*0.5
      N(L)=NEW
            GO TO 200
C           PUT YT FOR TRY
  190 NL=N(L)
      DO 192 I=1,NL+1
  192 YT(L,I)=Y(L,2+I)


C22222222222222222222222222222222222222222222222222222222222222222222

C TRY :
C       SETS (LOT,KIND) AND TRANSFERS CONTROL TO SUM, OR SET, OR HOLE, OR STORE
C       CHANGES ARRAYS  NAIL AND  INS.   USES YT FOR FILTERING

C33333333333333333333333333333333333333333333333333333333333333333333
C22222222222222222222222222222222222222222222222222222222222222222222222

C     PROCEDURE TRY(LOT,KIND,JOB,NOMEM)


C    TRIES SEGMENT AND ITS PARTS, WHETHER THEY ARE OKEY FOR QUADRATURE
C OUTPUTS: (LOT,KIND): OF SECTION

C TRY ========= ERASE E-FUNCTION MEMORY
  200 DS=0d0
      DO 201 I=1,2
      MADE(L,I)=0
  201 MADE(L,I+2)=0

C CONTROL GAP SIZE.1:  IF G IS TOO BIG, MAKE INSERTIONS WITHOUT TRIES
      LOT(L)=0
           IF(ABS(G(L)).GT.GMAX(L))GO TO 250
C                                  N=1 WITH INS-POINT: N=2 WITH G/2
      KIND(L)=2
      IF(N(L).EQ.1.AND.INS(L,1).EQ.1)GO TO 100
C                 2: IF G IS TOO SMALL, SUM WITHOUT CHECKS
      KIND(L)=1
      LOT(L)=N(L)
C ALTERNATION OF YT
      DO 202 I=1,N(L),2
  202 YT(L,1+I)=-YT(L,1+I)
      E1(L)=E2(1,N(L),1)
           IF(ABS(G(L)).LT.2*GMIN(L))GO TO 500
C --------------------------            500: SUM
C                                   INSERT POINT NEAR BARE ENDS
      IF(GA(L,1).NE.0d0)GO TO 3000
        IF(N(L).EQ.2.AND.10d0*E1(L).LE.PA(L))GO TO 500
C====                                CALL INSERT(1,F)
      IN(L)=1
      LINS(L)=1
               GO TO 10

 1001 IF(N(L).GT.1)GO TO 3000
      GO TO 3001
 3000 IF(GA(L,2).NE.0d0)GO TO 3002
C====                             CALL INSERT(LOT,F)
      IN(L)=LOT(L)
      LINS(L)=2
              GO TO 10
 1002 IF(N(L).GT.1)GO TO 3002
 3001 KIND(L)=2
      LOT(L)=0
      GO TO 100

C                      LEAP: LEAP=1 IS A SIGNAL TO SKIP SOME CHECKS
 3002 LEAP(L)=1

C  TOTAL  ------------------- IS TOTAL SEGMENT OKEY?

C FILTERING
            IF(E1(L).GT.P200(L)) GO TO 230
C                                230: PARTS (WITH LEAPS)
C---------------------------------------------------
      EP=E1(L)+P03(L)
      PR(L)=PA(L)
      IF(EP.LT.PA(L))PR(L)=EP
      IF(INS(L,1).EQ.1)PR(L)=PA(L)
      IF (N(L).NE.2*(N(L)/2)) GO TO 210
C EVEN N
C LEFT END CHECK
C                                                 CHECK LEFT END
      IF(E2(3,N(L),1).LE.PR(L) )GO TO 203
C MAYBE 2ND WING POINT IS BAD: MAKE STRINGENT CHECKS WITH 1ST W-POINTS INSTEAD
       IF(E2(2,N(L),1).GT.P03(L))  GO TO 206
        IF(E2(2,N(L),2).GT.P03(L))  GO TO 206
C ================                                 TRY TOTAL
  203 IF(E2(4,N(L),1).LE.PR(L) )GO TO 500
C                                               500: SUM

C==============================

      IF(INS(L,N(L)).EQ.1)GO TO 229
C                          229: PARTS (WITHOUT  LEAPS)
C                                TRY TO MEND RIGHT END BY INSERTION
      IN(L)=N(L)
      LINS(L)=3
               GO TO 10
 1003 GO TO 208
C=====================================================================
  206 IF(INS(L,N(L)).EQ.1)PR(L)=PA(L)
        IF(E2(3,N(L),2).LE.PR(L) )GO TO 207
        GO TO 230
C             230: PARTS WITH LEAPS
C-------------------------------
  207 IF(INS(L,1).EQ.1)GO TO 229
        IN(L)=1
        LINS(L)=4
                 GO TO 10
 1004 CONTINUE
C    ..........................   TRY TOTAL AGAIN
  208 IF(E2(4,N(L),1).LE.PA(L) )GO TO 500
C                                              500: SUM
      GO TO 229

C===========================================================================
C ODD N                                        L-END  CHECK
  210 IF (E2(2,N(L),1).LE.PR(L) )GO TO 211
      GO TO 214
C--------------------------                    R-END CHECK
  211 DS0(L)=DS
      IF(INS(L,N(L)).EQ.1)PR(L)=PA(L)
      IF(E2(2,N(L),2).LE.PR(L) )GO TO 500
      IF (INS(L,N(L)).EQ.1) GO TO 229
C---------------------            TRY TO MAKE R-END GOOD BY INSERTING
  212 IN(L)=N(L)
      LINS(L)=5
               GO TO 10
 1005 IF(E2(2,N(L),2).LE.PA(L)) GO TO 500
      GO TO 229
C=========================================================================

  214 IF(INS(L,N(L)).EQ.1)PR(L)=PA(L)
      IF(E2(2,N(L),2).GT.PR(L) )GO TO 230
C-------------------------                TRY TO SAVE L-END BY INSERTION
  215 DS0(L)=DS
      IN(L)=1
      LINS(L)=6
                GO TO 10
 1006 IF(E2(2,N(L),1).LE.PA(L) ) GO TO 500

C============================================================================
C PARTS OF SEGMENT: SECTION WITH 'LOT' GAPS
  229 LEAP(L)=0
  230 LOT(L)=N(L)
      LOT1(L)=0
C L-CYCLE ------------- L-END FIXED
  240 LOT0(L)=LOT(L)
      LEAP4(L)=LEAP(L)*(int(LOT(L)/4))
      LOT(L)=LOT(L)-1-LEAP4(L)
      IF(LOT(L).EQ.0) GO TO 250
C                        250: R-CYCLE.    LEAP4 : HOW MUCH TO SKIP
C ----------------------  CHECK FIXED END
      IF(LOT(L).EQ.LOT1(L))GO TO 241
C                              MAKE ROUGH CHECK
      E1(L)=E2(1,LOT(L),1)
C                              FILTERING
      IF(E1(L).GT.P200(L))GO TO 240
C ---------------------------   CORRECTED PRECISION
      EP=E1(L)+P03(L)
      PR(L)=PA(L)
      IF(EP.LT.PA(L))PR(L)=EP
      IF(INS(L,1).EQ.1)PR(L)=PA(L)
      IF(LEAP4(L).EQ.0)GO TO 242
C                               STORE E1, PR1
      LEAP(L)=0
      LOT1(L)=LOT(L)
      LOT(L)=LOT0(L)
      E11(L)=E1(L)
      PR1(L)=PR(L)
C             IF E IS NOT BIG, MAKE OMITTED CHECKS
      GO TO 240
  241 E1(L)=E11(L)
      PR(L)=PR1(L)
  242 IODD(L)=LOT(L)-2*(LOT(L)/2)
C-----------------------                          CHECK FIXED L-END
      IF(E2(3-IODD(L) ,LOT(L),1).LE.PR(L))GO TO 244

      IF(INS(L,1).EQ.1)GO TO 240
      IF(E1(L).GT.PR(L))   GO TO 240
C-----------------------------  TRY TO SAVE FIXED L-END
      IN(L)=1
      LINS(L)=7
                GO TO 10
 1007 IF(E2(3-IODD(L),LOT(L),1).GT.PR(L))GO TO 240
C==============================================
  244 DS0(L)=DS
C --------------------    CHECK FREE END, USING OUTSIDE INSERTED POINT
      IN(L)=LOT(L)+1
      LINS(L)=8
               GO TO 10
 1008 IF(E2(6+IODD(L),LOT(L),1).LE.PR(L) )GO TO 500
      IF(LOT(L).GT.1)GO TO 246
C---                     CARE OF DOUBLE-CHECK AT LOT=1
      IF(INS(L,1).EQ.1)GO TO 240
C -----                   TRY TO SAVE FREE END: USE INSIDE INSERTED POINT
  246 IN(L)=LOT(L)
      LINS(L)=9
               GO TO 10
 1009 IF(E2(4+IODD(L),LOT(L),1).LE.PA(L))GO TO 500
      GO TO 240
 1010 CONTINUE

C==========================================================================
C RIGHT PARTS TO TRY                       R-END IS FIXED
  250 KIND(L)=2
      LOT(L)=N(L)-LOT(L)
C----------- CAN POSSIBLE JUMP BE STORED? ------------------

                JUMP0=0
      IF(M(L).NE.0)JUMP0=Z(L,M(L))
             NOMEM(L)=LOST(L).NE.0.OR.(JUMP0.LT.1.AND.M(L).GT.186)
C            NOMEM: TRUE, IF NO MEMORY FOR STORING JUMP IS AVAILABLE
C------------------------------------------------------------------


C RIGHT CYCLE     IF: KEEP 16 GAPS FOR NEW SEGMENT, STORE WHAT RESTS

  260 IF(N(L)-LOT(L).GE.8)GO TO 300
C                         300: STORE
      LOT(L)=LOT(L)-1
C                                          COMPLETE INSERTIONS
      IN(L)=N(L)-LOT(L)
      LINS(L)=11
                GO TO 10
 1011 IF(LOT(L).EQ.0)GO TO 100
C                       100: SET (MAKE NEW SEGMENT OF FULL OLD ONE)
                   IF(NOMEM(L).OR.ABS(G(L)).GT.GMAX(L))GO TO 260
C-----------                   TENTATIVE   CHECK FIXED END
      E1(L)=E2(1,LOT(L),2)
C                             FILTERING
      IF(E1(L).GT.P200(L))GO TO 260
C                CORRECTED PRECISION
      EP=E1(L)+P03(L)
      PR(L)=PA(L)
      IF(EP.LT.PA(L))PR(L)=EP
      IF(INS(L,N(L)).EQ.1)PR(L)=PA(L)
      IODD(L)=LOT(L)-2*(LOT(L)/2)
C================================            CHECK OF FIXED R-END
 1260 IF(E2(3-IODD(L),LOT(L),2).LE.PR(L)  )GO TO 262
      IF(INS(L,N(L)).EQ.1 )GO TO 260
      IF(E1(L).GT.PR(L))GO TO 260
C--------------------------------------- TRY TO SAVE FIXED R-END
  261 IN(L)=N(L)
      LINS(L)=12
                GO TO 10
 1012 IF(E2(3-IODD(L),LOT(L),2).GT.PA(L))GO TO 260
C------------                              CHECK FREE END WITH OUTSIDE POINT
C           IF OKEY, STORE JUMP.IF NOT: NOT TRY TO SAVE, JUST STORE FOR FUTURE
  262 DS0(L)=DS
      IF(E2(6+IODD(L),LOT(L),2).LE.PA(L) )GO TO 400
C                                              400: HOLE
      GO TO 300
C           300: STORE
C33333333333333333333333333333333333333333333333333333333333333333333

C STORE
C              MEMORY LAYOUT
C      N= LOT                        JUMP FROM X1 TO X2
C      M: -INSN                          M: 1
C     -1: N                             -1: X2
C     -2: X1 OF LEFT END                -2: X1       USE AS R-WING
C     -3: GB(2)                         -3: GB(2)
C     -4: GA(2)                         -4: GA(2)
C     -5: G                             -5: Y(N+5)
C     -6: Y(N+5)                        -6: Y(N+4)
C     -7: Y(N+4)                        -7: Y(N+3)
C     -8: Y(N+3) OF R-END               -8: Y(3)     USE AS L-WING
C     ...                               -9: Y(2)
C   -N-8: Y(3)   OF L-END              -10: Y(1)
C                                      -11: GA(1)
C                                      -12: GB(1)
C  -N- 9: YI(IN(INSN))                 -13: S OF JUMP
C     ...
C  -N- 8-INSN: YI(IN(1))
C  -N- 9-INSN: IN(INSN)
C     ...
C  -N- 8-2*INSN: IN(1)
C
C   ----------------                 ----------------
C   M(-1)=M-N-2*INSN-9                M(-1)=M-14
C--------------------------------------------------------

C STORE
  300   IF(LOST(L).NE.0)GO TO 304
        INPUT=LOT(L)+9
        INSN=0
C                 INSN: NUMBER OF INSERTIONS TO STORE
      DO 301 I=N(L)-LOT(L)+1, N(L)
  301 INSN=INSN+INS(L,I)
                      INPUT=INPUT+2*INSN
C                  IS IT POSSIBLE TO STORE?
      LOST(L)=LOT(L)
      IF(M(L)+INPUT.GT.200)GO TO 304
      LOST(L)=0
C                                STORE INSERTIONS
      J=0
      DO 302 I=1,LOT(L)
          IF(INS(L,N(L)-LOT(L)+I).EQ.0)GO TO 302
        J=J+1
        Z(L,M(L)+J)=I
        Z(L,M(L)+J+INSN)=YI(L,N(L)-LOT(L)+I)
  302 CONTINUE
      M(L)=M(L)+INPUT
C          M: NEW M: MEMORY FILLING
C                                STORE NODES
      DO 303  I=1,LOT(L)+3
  303 Z(L,M(L)-LOT(L)-9+I)=Y(L,2+N(L)-LOT(L)+I)
C                               STORE GAPS ETC.
      Z(L,M(L)-5)=G(L)
      Z(L,M(L)-4)=GA(L,2)
      Z(L,M(L)-3)=GB(L,2)
      Z(L,M(L)-2)=X1(L)+G(L)*(N(L)-LOT(L))
      Z(L,M(L)-1)=LOT(L)
      Z(L,M(L))=-INSN
  304    GO TO 100
C              100: SET

C4444444444444444444444444444444444444444444444444444444444444444444444

C HOLE   STORING HOLES IN THE INTERVAL (A,B), WHERE INTEGRATION IS DONE
  400 IF(M(L).EQ.0) GO TO 402
      IF(Z(L,M(L)).GT.0.5)GO TO 410

C FIRST                    410: APPEND
C           SET  INTEGRAL OVER HOLE TO ZERO
  402 Z(L,M(L)+1)=0d0
      M(L)=M(L)+14
C           SET FLAG Z(M) OF THE RECORD
      Z(L,M(L))=1.
      Z(L,M(L)-8)=Y(L,N(L)+3)
C                             MAKE L-WING DATA
      CALL WING2(2,0,DW,YW)
               Z(L,M(L)-11)=DW(1)
               Z(L,M(L)-9)=YW(1)
      Z(L,M(L)-12)=DW(2)-DW(1)
      Z(L,M(L)-10)=YW(2)
      Z(L,M(L)-1)= X1(L)+G(L)*N(L)
C APPEND
  410 LEFT=N(L)-LOT(L)
      Z(L,M(L)-2)=X1(L)+LEFT*G(L)
      Z(L,M(L)-7)=Y(L,LEFT+3)
      CALL WING2(1,LEFT,DW,YW)
      Z(L,M(L)-6)=YW(1)
      Z(L,M(L)-5)=YW(2)
      Z(L,M(L)-4)=DW(1)
      Z(L,M(L)-3)=DW(2)-DW(1)
C555555555555555555555555555555555555555555555555555555555555555555555555

C SUM:
C    SUMS POINTS  Y  WITH WEIGHTS  W  AND ACCUMULATES THE INTEGRALS SDT AND TOT

  500 L2=LOT(L)/2
C              NW: NUMBER OF DIFFERENT WEIGHTS IN THE QUADRATURE FORMULA
      NW=L2+1
C              JJ: POINTER TO SET OF WEIGHTS
      JJ=(LOT(L)-L2)*NW
C                           LEFT: LEFT PART LEFT ASIDE (NUMBER OF GAPS)
      LEFT=(N(L)-LOT(L))*(KIND(L)-1)
C----------------------------------------------------------------------------

C        QUADRATURE           QU: AVERAGE VALUE OF THE FUNCTION IN THE SECTION
      QU=DS
      IF(L2+L2.NE.LOT(L))QU=(DS+DS0(L))*0.5
      L2=LEFT+2
      LL2=L2+LOT(L)+2
      DO 502 I=1,NW
  502 QU=QU+(Y(L,L2+I)+Y(L,LL2-I))*W(JJ+I)
C                             H: LENGTH OF THE INTEGRATED  SECTION
      H=G(L)*LOT(L)
C                  'IF': FIGHTS ACCUMULATION OF ROUNDING-OFF ERRORS
      IF(M(L).GT.0.AND.LOT(L).EQ.N(L).AND.LOST(L).EQ.0)
     *         H=Z(L,M(L)-2)-X1(L)
      S=QU*H
C                S: INTEGRAL OVER THE SECTION WITH  LOT GAPS
C----------------------------------------------------------------------------

C     INTEGRAL ACCUMULATION     SDT: INTEGRAL FROM A TO X1+H*(2-KIND)
C                               TOT: TOTAL INTEGRAL INCLUDING HOLES
      IF(LEFT.EQ.0)GO TO 505
      RQA(L)=X1(L)+LEFT*G(L)
      RQB(L)=X1(L)+N(L)*G(L)
      JQ(L)=1
C                                 PUT S  IN THE HOLE
               Z(L,M(L)-13)=Z(L,M(L)-13)+S
      GOTO 506

C                                 ADD S  TO SDT
  505           SUM(L)=SUM(L)+S
  506 TOT(L)=TOT(L)+S
      SUMTOC(L)=SUM(L)
C------------------------------------------------
C     ERROR ACCUMULATION
C                       NUMBER OF QUADRATURES
      NQUAD(L)=NQUAD(L)+1
C                       ROUNDING-OFF ERROR
      ROUND=S*ACCURA
C                       SQARED ERROR AND ROOT-MEAN-SQUARED ERROR
      SQUERR(L)=SQUERR(L) +ESTER*ESTER +ROUND*ROUND
      RMSERR(L)=SQRT(SQUERR(L))
C------------------------------------------------
             IF(KIND(L).EQ.2)GO TO 100
C                               100: SET
C           STEP FORWARD                           IS DONE, IF KIND=1
      X1(L)=X1(L)+H
              IF(LOT(L).EQ.N(L))KIND(L)=0


      GO TO 100
C              100: SET
      END




      SUBROUTINE WING2(KIND,LOT,DW,YW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /SD214V/
     * X1(2),G(2),PA(2),ESTER,DS, GA(2,2),GB(2,2),Y(2,21)
     *   ,YI(2,16),YT(2,19),errold,L,N(2),INS(2,16),MADE(2,5)
     *                 ,NOLD,MOLD,MODY

      DIMENSION DW(2),YW(2)
C--------------------------------------------------------------
C   WING SETTING
C              VARIABLES:
C                        I: NUMB. OF GAPS PASSED
C                        J: NUMB. OF WING POINTS FOUND
C                       NI: INDEX OF A CONSIDERED GAP
C              YW: VALUES OF POINTS FOUND
C              DW: (-) COORDINATES (FROM FREE END) OF POINTS
C---------------------------------------------------------------
      NOUT=N(L)-LOT
      KO=3-KIND
      I=0
      J=0
    1 IF(KIND-1)  3,3,4

    3   NI=LOT+1+I
      GOTO 5
    4   NI=NOUT-I
C--------------------------------------------------------
C------------------------------------------------------
C                                                INSERTED POINT INTO WING
C     SKIP, IF INSERTED POINT MAY NOT BE (NI<=0) OR
C                                       MAY BE (NI>0) BUT IS NOT (INS=0)
    5 IF (NI.LE.0) GO TO 2
      IF(INS(L,NI).EQ.0) GO TO 2
        J=J+1
        YW(J)=YI(L,NI)
        DW(J)=(I+0.5)*G(L)
    2 IF(J.GT.1) RETURN
C-------------------------------------------------
C--------------------------------------------------
C                                             NODE INTO WING
      I=I+1
      J=J+1
      YW(J)=Y(L,3+NOUT-I)
      IF(KIND.EQ.1) YW(J)=Y(L,3+LOT+I)
      IF(NOUT.GE.I) GOTO 10
                    GOTO 11
   10   DW(J)=I*G(L)
        GOTO 12
   11   DW(J)=G(L)+GA(L,KO)

   12   IF(J.LT.2) GO TO 1
      RETURN
      END



      real*8 FUNCTION E2(MODE,LOT,KIND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C                                E: ERROR ESTIMATE


C   MODE (OF USING CONTROL POINTS):
C                         POSITION OF CONTROL POINT AT FREE END
C
C                         INSIDE  0  OUTSIDE
C   NUMBER OF     $   2           3
C   CONTROL POINTS$   1       4   2   6
C   AT FIXED END  $   0       5   1   7
C
C   LOT: NUMBER OF GAPS TO KEEP IN THE CHECKED PART OF SEGMENT
C
C   KIND (OF FIXED END):    1= LEFT END
C                           2= RIGHT END
C
C   COMMONS:
C           N: NUMBER OF GAPS IN THE (TOTAL) SEGMENT
C       Y(21): FUNCTION VALUES AT NODES
C           G: (GAP) DIFFERENCE OF COORDINATES OF ADJACENT NODS
C       GA(2): CLOSEST GAP ADJACENT TO SEGMENT
C       GB(2): NEXT GAP AFTER THE CLOSEST
C
C              NUMERATION:

C NODS=  1     2     3     4     5   ...  N+3   N+4   N+5    INDEX IN Y
C        .GB(1).GA(1).  G  .  G  .   ...   .GA(2).GB(2).     GAP VALUES
C GAPS=                 1     2      ...                     INDEX IN YI,INS
C
C     INS(16): (INSERTIONS) INDICATORS OF INSERTED POINTS
C                            INS(I)=0: NO INSERTED POINT IN THE GAP I
C                            INS(I)=1: INSERTED POINT IS IN THE GAP I
C      YI(16): (Y INSERTED) FUNCTION VALUES AT INSERTED POINTS

C    MADE(5): INDICATORS OF COMPUTATIONS MADE. MADE(KM)=N   MEANS THAT FOR
C              LOT <=N AND FOR THE CASE KM=(KIND,MODE) INTERMEDIATE VALUES
C              ARE COMPUTED. FOR MODE=2 NEGATIVE SIGN OF MADE(KM) MEANS THAT
C              ONLY NORMALISATION OF POINTS IS MADE.

C OUTPUT PARAMETER DS:

C       DS: CORRECTION TO THE QUADRATURE FORMULA FOR THE CURRENT MODE AND LOT

C INTERNAL VARIABLES:

C     V(80)=   V(5,16): INTERMEDIAT VALUES DEPENDING ON Y
C     Q(80)=   Q(5,16): COEFFICIENTS DEPENDING ON DISTANCES
C     T(42)=   T(2,21): NORMALIZED AND REORDERED FUNCTION VALUES
C     RA     : RELATIVE DISTANCE FROM THE END OF SPAN TO THE
C              WING POINT (IF RA>0), OR NAIL POINT (IF RA<0), OR INSERTED
C              POINT (IF RA=-0.5)
C     CB     : DISTANCE TO ADDITIONAL WING POINT
C     TA,TB  : FUNCTION VALUES IN THE ADDITIONAL (WING OR NAIL) POINTS
      DIMENSION JR(16),RC(64),CF(16)
      DATA JR/0,0,1,2,4,6,9,12,16,20,25,30,36,42,49,56/
      DATA RC/
     *1.d0,3.d0,
     *4.d0,3.d0,
     *5.d0,10.d0,
     *6.d0,15.d0,10.d0,
     *7.d0,21.d0,35.d0,
     *8.d0,28.d0,56.d0,35.d0,
     *9.d0,36.d0,84.d0,126.d0,
     *10.d0,45.d0,120.d0,210.d0,126.d0,
     *11.d0,55.d0,165.d0,330.d0,462.d0,
     *12.d0,66.d0,220.d0,495.d0,792.d0,462.d0,
     *13.d0,78.d0,286.d0,715.d0,1287.d0,1716.d0,
     *14.d0,91.d0,364.d0,1001.d0,2002.d0,3003.d0,1716.d0,
     *15.d0,105.d0,455.d0,1365.d0,3003.d0,5005.d0,6435.d0,
     *16.d0,120.d0,560.d0,1820.d0,4368.d0,8008.d0,11440.d0,6435.d0/
      DATA CF     /0.16666666667d0,-0.06666666667d0,
     3             0.05000000000d0,-0.06349206349d0,
     5             0.02728174603d0,-0.06000000000d0,
     7             0.01804012346d0,-0.05695045695d0,
     9             0.01316456980d0,-0.05436160694d0,
     1             0.01020545107d0,-0.05215584416d0,
     3             0.00824255288d0,-0.05025673744d0,
     5             0.00685759239d0,-0.04860299139d0/
C     JR(I): POINTER TO THE ITH PART OF THE TABLE RC
C     RC : TABLE OF BINOMIAL COEFFICIENTS
C     CF : WEIGHTS USED FOR COMPUTATION  OF CORRECTIONS DS
C                SKETCH OF STRUCTURE
C     'E' REMEMBERS AND REUSES INTERMEDIATE VALUES V,Q,
C     WHILE THEY ARE VALID. IF NAIL OR INSERTED POINT IS
C     AVAILABLE, 'E' PUTS IT ON PLACE OF A WING
C     POINT, AND ERASES ON OBSOLETE PART OF ITS MEMORY.

C     CHAIN OF BLOCKS:
C     FAST MODE 1 - INDEX PACKING - NORMALIZATION AND ORDERING -
C     WING POINTS LISTING -  INSERTED POINT LISTING
C     AND PARTIAL ERASURE-
C     V,Q FOR MODE=2,5,7
C     V,Q FOR MODES 3,4,6 - ATERNATION OF V -
C     DEVIATION AND E - CORRECTIONS DS.
C
      COMMON /SD214V/
     * X1(2),G(2),PA(2),ESTER,DS, GA(2,2),GB(2,2),Y(2,21)
     *   ,YI(2,16),YT(2,19),errold,L,N(2),INS(2,16),MADE(2,5)
     *                 ,NOLD,MOLD,MODY

      DIMENSION V(2,80),Q(2,80),T(2,42),RA(2,2),CB(2,2)
     *   ,TA(2,2),TB(2,2),Y0(2,2)

      save V,Q,T,RA,CB,TA,TB,Y0

C--------------------------------------------------------------
      N2=LOT/2
      JLOT=JR(LOT)
      IF(MODE-2)10,20,20
C-----------------------------------------------------------------
C FAST MODE 1 ++++++++++++  I1,LAST: FIRST AND LAST POINT IN YT TO USE
   10 I1=1+(N(L)-LOT)*(KIND-1)
      LAST=I1+LOT
      DEV=YT(L,I1)+YT(L,LAST)
      IF(N2)14,14,12
   12 DO 13 I=1,N2
   13 DEV=DEV+(YT(L,I1+I)+YT(L,LAST-I))*RC(JLOT+I)
   14 E2=ABS(DEV*G(L))*4.
      GO TO 1000
C------------------------------------------------------------------
C           GENERAL PROCEDURE FOR WING-MODES.  IOUT:INDICATOR OF MODES 6,7
   20 E2=PA(L)*10d0
      NL=N(L)
      IOUT=0
           IF(MODE.GE.6) IOUT=1
C            PACKING INDECES
      IF(KIND.EQ.2) GO TO 22
      KN=1
      KK=LOT+IOUT
           GO TO 24
   22 KN=NL
      KK=NL+1-LOT-IOUT
C            (KIND,MODE) PACKING
   24 KM=KIND+2*(MODE-2)
      KMQ=KM
           IF(MODE.GT.3) KMQ=5
      KAM=16*(KMQ-1)
      KIN=21*(KIND-1)
      KO=3-KIND
        IF(GA(L,KIND).NE.0d0)GO TO 25
        IF(INS(L,KN).NE.0)GO TO 25
        IF(MODE.EQ.2)GO TO 1000
        IF(MODE.EQ.4.OR.MODE.EQ.6)GO TO 1000
   25 CONTINUE
C----------------------------------------------------------------
C         KM=(KIND,MODE)
C         KN:INDEX OF GAP AT FIXED END
C         KK:INDEX OF GAP AT FREE END
C     ONE-DIMENSIONAL CONVERSION OF ARRAYS V,Q,T:
C         T(KIN+I)=T(KIND,I)
C         V(KAM+I)=V(KMQ,I)=V(KIND,MODE,I)
C         KMQ: INDEX KM OF SHORTENED RANGE FOR ARRAYS Q,V
C         KO: KIND OPPOSITE
C--------------------------------------------------------------

C                 NORMALIZE AND ORDER
      IF(ABS(MADE(L,KIND)).GE.LOT) GO TO 42
C         SKIP NORMALIZATION AND WING POINTS HANDLING

      GO TO(31,33),KIND

C                 NORMALIZE TO LEFT END
   31 Y0(L,1)=Y(L,3)
          DO 32 I=1,LOT+5
   32     T(L,I)=Y(L,I)-Y0(L,1)
      GOTO 35
C------
C                 NORMALIZE TO RIGHT END AND ORDER
   33 Y0(L,2)=Y(L,NL+3)
          DO 34 I=1,LOT+5
   34      T(L,21+I)=Y(L,NL+6-I)-Y0(L,2)
C---------------------------------------------------------------
   35 MADE(L,KIND)=-LOT
C                   BOOK-KEEPING: HOW MUCH IS MADE
C---------------------------------------------------------------
C          WING POINTS
      RA(L,KIND)=GA(L,KIND)/G(L)
      TA(L,KIND)=T(L,KIN+2)
      TB(L,KIND)=T(L,KIN+1)
      CB(L,KIND)=GA(L,KIND)+GB(L,KIND)
C---------
C          INS-POINT USE

   42 IF(INS(L,KN).EQ.0)GO TO 47
      IF(RA(L,KIND).EQ.-0.5d0)GO TO 47
C         IF INS-POINT APPEARED
      RA(L,KIND)=-0.5
      TA(L,KIND)=YI(L,KN)-Y0(L,KIND)
      TB(L,KIND)=T(L,KIN+2)
      CB(L,KIND)=GA(L,KIND)
C         ERASE WING - DEPENDENT PART OF STORE
      MADE(L,KIND)=0
      MADE(L,KIND+2)=0
C-----------------------------------------------------------
   47 IF(MODE.GT.3)GO TO 48
      IF(MADE(L,KMQ).GE.LOT) GO TO 600
C                                600: SKIP Q,V CALCUL.

   48 RK=RA(L,KIND)
      MOD=MODE-1
      GO TO(200,300,400,201,412,201),MOD
C         Q AND V BLOCKS
C-----------------------------------------------------------
C           MODES 2,5,7: LINE SUBTRACTION
  200 TAK=TA(L,KIND)
      GO TO 221
C           MODES 5,7
  201 RK=-LOT-IOUT+0.5
      YA=YI(L,KK)
C---------------------------------------------------
C           JOINT PARTS IF MODES 2,5,7: LINE SUBTRACTION
  220 TAK=YA-Y0(L,KIND)
C                      TAR: LINE COEFFICIENT
  221 TAR=TAK/RK
      KIN3=KIN+3
      DO 224 I=1,LOT
      Q(L,KAM+I)=I+RK
  224 V(L,KAM+I)=(T(L,KIN3+I)+TAR*I)/Q(L,KAM+I)
      GO TO 500
C---------------------------------------------------
C           MODE 3: PARABOLA SUBTRACTION. ONE WING USED
  300 R=CB(L,KIND)/G(L)
      IF(GA(L,KIND).EQ.0d0)GO TO 1000
      TBK=TB(L,KIND)
      GO TO 422
C---------------------------------------------------
C            MODE 4: TWO WINGS USED
  400 IF(LOT.LT.NL)GOTO 412
        IF(INS(L,KK).NE.0)GO TO 414

C---------------------------------------------------
  410 IF(GA(L,KO).EQ.0d0)GO TO 412
C                             USE FREE WING POINT, IF AVAILABLE
      YB=Y(L,KK+4)
      R =-NL-GA(L,KO)/G(L)
      GO TO 420
C ---------------------------------------------------
  412 IF(INS(L,KK).EQ.0) STOP ' USE OF NONEXISTENT INS-POINT'
C            MODES 4,6: USE INSERTED POINT
  414 YB=YI(L,KK)
      R=-LOT-IOUT+0.5
C---------------------------------------------------
C            JOINT PART OF MODES 3,4,6: PARABOLA SUBTRACTION
C---------------------------------------------------
C            RK=RA(KIND): FIXED END NAIL
  420 TBK=YB-Y0(L,KIND)
  422 RB=R-RK
      R1=TBK/(RB*R)
      R2=TA(L,KIND)/(RB*RK)
         A1=R1*RK-R2*R
         A2=R1-R2
C            A1,A2: PARABOLA COEFF.
      KIN3=KIN+3
      DO 424 I=1,LOT
      Q(L,KAM+I)=(I+RK)*(I+R)
  424 V(L,KAM+I)=(T(L,KIN3+I)-I*(A1+I*A2))/Q(L,KAM+I)
C---------------------------------------------------
C            ALTERNATING V'S
  500 DO 502 I=1,LOT,2
  502 V(L,KAM+I)=-V(L,KAM+I)
C---------------------------------------------------
C            STORE KEEPING
      MADE(L,KMQ)=LOT
C---------------------------------------------------
C            DEVIATION AND E
  600 KAMLOT=KAM+LOT
      DEV=V(L,KAMLOT)
      IF(N2.EQ.0) GOTO 603
C                 DO SUM
      DO 602 I=1,N2
  602 DEV=DEV+(V(L,KAM+I)+V(L,KAMLOT-I))*RC(JLOT+I)

C----------------  DEV: WEIGHTED DEVIATION
  603  DEVIAT=DEV*Q(L,KAMLOT)
C---------------   DEVIAT: NATURAL DEVIATION OF EXTRAPOLANT FROM THE END POINT
      E2=ABS(DEVIAT*G(L))
C                  CORRECTION TO CLOSENESS OF NAIL
      IF(MODE.GE.6)E2=E2+E2+E2

      IF(MODE.EQ.4.AND.R+LOT.LT.0d0.AND.LOT.EQ.2)E2=E2*3.0

  605 ESTER=E2
C     CORRECTION DS
      DS=0d0
C                NO DS, IF LOT IS EVEN AND ONLY ONE ADDITIONAL POINT USED
      IF(MODE.EQ.2.AND.LOT.EQ.2*N2) GO TO 1000
      DS=DEV*CF(LOT)
C---------------------------------------------------------
 1000 RETURN

      END
