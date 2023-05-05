*---------------------------------------
!>     Controlled program termination.
!>        - print job summary
!>        - print full error summary
!>        - close files
*           ...
*---------------------------------------
      subroutine HF_stop

      implicit none
      integer   nunits
      parameter (nunits=13)
      integer   lun(nunits), i
      data      lun /7, 24,25, 51,52,53, 61, 76,77, 81,85,87, 90/

      write(6,*) 
      write(6,*)' ******** HF_STOP forced program termination ********'
      write(6,*)

*     call HF_jobsum
      call HF_errsum(6)

      do i=1,nunits
        close(lun(i))
      enddo

      stop
      end


*                                                                      
*  S.Levonian  5.02.2012 (adopted from H1util code of S. Egli)         
*                                                                      
*-----------------------------------------------------------------------
*                    PROPOSAL FOR ERROR NUMBERING                          
*                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~                          
*  For xFitter code:
*  ------------------
*  To add a new error number so that it would not clash with already 
*  existing ones it is proposed to put the date of coding into it:        
*                 error number = YYMMDDnn                            
*           (e.g. 120117nn for the date = 17.01.2012)                 
*  This gives for code developers 100 numbers per day which is enough.
*  You will never clash with yourself and quite unlikely with others.
*  For imported code:
*  ------------------
*  Reserve 100 IDs for every imported module 
* 
*-----------------------------------------------------------------------
!> Error logging facility                                      
!> @param[in] ID error identifier (integer)
!> @param[in] TEXT text string connected to ID                            
!>                (ignored after first appearence of error ID)
!>          
!> Output:                         
!>
!> The following arrays are filled for each error message:      
!>       - ERRTXT:    store of error messages (one big string)
!>       - IDERR:     error identifier                                  
!>       - IPSTRT:    pointer to first character of message in ERRTXT   
!>       - IPEND:     pointer to last character                         
!>       - ICNT:      error count                                       
!>
!> Four levels of severity are supported.
!> The level is coded in the first two characters of the TEXT message:
!>       - I: - informational message,
!>       - W: - warning       (continue execution by default)
!>       - E: - error         (terminate program by default, but continue execution if errors are treated as warnings)
!>       - S: - serious error (terminate program by default)
!>       - F: - fatal error   (terminate program always)
C-----------------------------------------------------------------------    
      subroutine HF_errlog(ID,TEXT)

      implicit none

#include "for_debug.inc"

      INTEGER   MAXERR,MAXTXT,NMODS 
      PARAMETER (MAXERR=400, MAXTXT=20000, NMODS=10)

      INTEGER TAB(MAXERR),AUX(2*MAXERR)
      INTEGER IDERR(MAXERR),IPSTRT(MAXERR),IPEND(MAXERR),ICNT(MAXERR)
      LOGICAL FIRST,FIRMES
      CHARACTER*(*) TEXT
      CHARACTER*(MAXTXT) ERRTXT
      CHARACTER*12  MOD_NAME(NMODS)
      CHARACTER*80 ERRLIN
      CHARACTER*2  TSEV(5), ERM_TYPE
      CHARACTER*23 FULLTX(5)
      INTEGER      NLOST,ITOP,ISEVER,MAX_ERR_ALLOWED,IMESS,I,J
      INTEGER      ID,NEN,IND,LTEXT,NCHAR,ITOTAL, LL,LUN, I_MOD_NAM
      INTEGER      LENB, ITLU, JTLU
      INTEGER getparami
      external getparami
      
      DATA TSEV/'I:','W:','E:','S:','F:'/
      DATA FULLTX/
     +   'Informational messages:',
     +   'Warning messages:',
     +   'Error messages:',
     +   'Serious messages:',
     +   'Fatal messages:'/

      DATA MOD_NAME /'xFitter',  ! ID = YYMMDDnn
!     +               'HFitter',   ! ID = YYMMDDnn
     +               'Minuit',    ! ID =   1- 99
     +               'Acot',      ! ID = 100-199
     +               'Dipole',    ! ID = 200-299
     +               'DY',        ! ID = 300-399
     +               'FastNLO',   ! ID = 400-499
     +               'HS',        ! ID = 500-599
     +               'Hathor',    ! ID = 600-699
     +               'RT',        ! ID = 700-799
     +               'Unknown'    ! ID = 800-999999 
     +              /

      DATA FIRST/.TRUE./,NLOST/0/,ITOP/1/
      DATA ICNT/MAXERR*0/

      IF(FIRST)THEN
        FIRST=.FALSE.
        CALL DTLU(MAXERR,0,TAB,AUX)
*       default max severity level for which continuation is still allowed
        MAX_ERR_ALLOWED = 2
      ENDIF

      MAX_ERR_ALLOWED = getparami("MaxErrAllowed")
      IF (MAX_ERR_ALLOWED.LT.1) THEN
         MAX_ERR_ALLOWED = 2
      ENDIF
      IF (MAX_ERR_ALLOWED.GT.5) THEN
         MAX_ERR_ALLOWED = 5
      ENDIF
      
      LTEXT = MAX(MIN(1,LEN(TEXT)), LENB(TEXT))

* print error message if in debug mode

      IF (DEBUG) THEN
        WRITE(6,*)' ==> HF_ERRLOG - ID =',ID,':  ',TEXT(:LTEXT)
      ENDIF

* get index of error message

      NEN=AUX(3)
      IND=ITLU(ID,TAB,AUX)

* is this an old error message ?

      IF(IND.NE.0)THEN
*       old message, just increment counter
        ICNT(IND)=ICNT(IND)+1
        GOTO 999
      ELSE
*       new message,check if sufficient space to store errortext
*       and other information
        NCHAR=MAX(MIN(1,LEN(TEXT)), LENB(TEXT))
        IF(AUX(1).EQ.AUX(3).OR.ITOP-1+NCHAR.GT.MAXTXT)THEN
*         insufficient space, just count lost error messages
          NLOST=NLOST+1
        ELSE
*         enough space, add error message to tables
          IND=JTLU(ID,TAB,AUX)
          IDERR(IND)=ID
          IPSTRT(IND)=ITOP
          ITOP=ITOP+NCHAR
          IPEND(IND)=ITOP-1
          ERRTXT(IPSTRT(IND):IPEND(IND))=TEXT
          ICNT(IND)=1
        ENDIF
      ENDIF

* check severity and terminate program if necessary 

        I=INDEX(TEXT(1:MIN(5,LTEXT)),':')
        IF(I.GT.0) THEN
          ERM_TYPE=TEXT(I-1:I)
          DO J=1,5
            IF(ERM_TYPE.EQ.TSEV(J)) THEN
              IF(J.GT.MAX_ERR_ALLOWED) THEN
                WRITE(6,*)' ==> HF_ERRLOG - terminate due to error =',ID
     +,                  ':  ',TEXT(:LTEXT)
                CALL HF_STOP
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      
      GOTO 999

      ENTRY HF_ERRSUM(LUN)
*--------------------------------------
*     print error summary on unit LUN 
*--------------------------------------

      IF (FIRST) RETURN

* calculate total number of errors

      NEN=AUX(3)
      ITOTAL=0
      DO I=1,NEN
         ITOTAL=ITOTAL+ICNT(I)
      END DO
      IF(ITOTAL.EQ.0)RETURN

* errlog printout to the requested output unit 

      LL=LUN
      IF (LL.LE.0) THEN
        WRITE(6,*) ' ==> HF_ERRSUM: LUN =',LL,' - no summary printed'
        RETURN
      ENDIF 

*  print header lines

      WRITE(LL,900) ITOTAL,NLOST

*  print header line for the list if at least one error is stored

      WRITE(LL,901)

*  loop over all errors, collect messages for same module and print
*  in order of severness: I,W,S,F,unknown

      DO 100 ISEVER=1,5

* for each severity level loop over all messages and
         if (ISEVER.eq.2) then
C Set color on
            write (LL,'(A10)',advance='no') achar(27)//'[34m'
         endif

         if (ISEVER.eq.3) then
C Set color on
            write (LL,'(A10)',advance='no') achar(27)//'[33m'
         endif

         if (ISEVER.ge.4) then
C Set color on
            write (LL,'(A10)',advance='no') achar(27)//'[31m'
         endif


        FIRMES=.TRUE.
        DO 110 IMESS=1,NEN

*  check if message has proper identification among first 5 characters

          J=IPSTRT(IMESS)
*         this is to left justify the text
          ERRLIN=ERRTXT(J:MIN(IPEND(IMESS),J+79))
          IF(INDEX(ERRLIN(1:5),TSEV(ISEVER)).EQ.0.AND.
     +      (ISEVER.NE.1.OR.INDEX(ERRLIN(1:5),':').NE.0))GOTO 110

* for first message of each severity level print header line to
* indicate the severity level:

          IF(FIRMES)THEN
            FIRMES=.FALSE.
            WRITE(LL,904)FULLTX(ISEVER)
          ENDIF

          IF(IDERR(IMESS).GT.9999999) THEN
            I_MOD_NAM=1
          ELSE
            I_MOD_NAM=MIN(NMODS,2+IABS(IDERR(IMESS))/100)
          ENDIF

*         and print:
          WRITE(LL,902)MOD_NAME(I_MOD_NAM),IDERR(IMESS),ICNT(IMESS),
     +                 ERRLIN(1:LENB(ERRLIN))

*         more to print for this message ?
510       J=J+80
          IF(J.LE.IPEND(IMESS))THEN
            ERRLIN=ERRTXT(J:MIN(IPEND(IMESS),J+79))
            WRITE(LL,903)ERRLIN(1:LENB(ERRLIN))
            GOTO 510
          ENDIF

110     CONTINUE

100   CONTINUE

      write (LL,'(A10)',advance='no') achar(27)//'[0m'

      WRITE(LL,905)

*     reset error tables
      do i=1,MAXERR
          ICNT(i) = 0
      enddo
      CALL DTLU(MAXERR,0,TAB,AUX)
      ITOP=1
      NLOST=0

999   RETURN

900   FORMAT(// 12X,36('*')/12X,3('*'),7X,'Messages Summary',7X,3('*')/
     +       12X,36('*')/' '
     +       /' Total number of logged messages:      ',I6,
     +       /' Total number of messages not recorded:',I8)
901   FORMAT(//,' List of messages sorted by severity level:'//
     + ' *',78('-')/
     + ' *','   Module   | Message |Message|'/
     + ' *','     Name   |    Type | Count | Message Description'/
     + ' *',78('-'))
902   FORMAT(2X,A12,I10,I8,1X,A)
903   FORMAT(31X,A)
904   FORMAT(/2X,A/1X,'*',23('-'))
905   FORMAT(/' *',78('-')//'  End of Message Summary'/)
      END

*--------------------------------------------
!> Lenght of non-blank text
!> @param TEXT
*---------------------------------------------
      INTEGER FUNCTION LENB(TEXT)

      implicit none
      CHARACTER*(*) TEXT
      INTEGER I,LE
*
      IF(TEXT.EQ.' ') THEN
         LENB=0
      ELSE
         LE=LEN(TEXT)
   10    IF(LE.GT.4) THEN
            IF(TEXT(LE/2:LE).EQ.' ') THEN
               LE=LE/2-1
               GOTO 10
            END IF
         END IF
         DO 20 I=LE,1,-1
         IF(TEXT(I:I).NE.' ') GOTO 30
   20    CONTINUE
   30    LENB=I
      END IF
*
  100 RETURN
      END


      SUBROUTINE DTLU(NDIM,NEN,TAB,AUX)
* ----------------------------------------------------------------------
*     TLU - package for fast table look-up
*           (look-up time does not depend on the table length)
*
*     The subprograms of this package use the 'hash'-technique.
*     Advantages of this method are:
*     (1) time for one table look-up is constant (no dependence on the
*         length of the table) and about 7 micro seconds for a
*         function including the call/return time (on DESY IBM 1988)
*     (2) the entries of the table can be in any order (no sorting
*         necessary)
*     (3) the code for a table look-up is short (5 statements) and
*         can be included inline in the user program, thus avoiding
*         the call/return time and reducing the time for one table
*         look-up to 3 micro seconds (for any number of entries)
*     (4) the time for one table look-up is the same, whether the entry
*         is found or not.
*
*     There is one disadvantage: an auxiliary array of length 2*NDIM
*     is necessary for a table of length NDIM (maximum of NDIM entries).
*
*     Subprograms: Subroutine DTLU
*                  Integer Functions ITLU and JTLU
*
*     Two arrays are used, TAB(NDIM) and AUX(2*NDIM), both assumed of
*     type integer. The array TAB(NDIM) is the table with a maximum
*     of NDIM entries. The array AUX is an auxiliary array of length
*     2*NDIM. NDIM has to be at least 10 (otherwise the program will
*     stop).
*     The auxilary array AUX is used internally and it has to be
*     initialized by a call of DTLU before any other operation.
*
*     CALL DTLU(NDIM,NEN,TAB,AUX)
*
*               NDIM = dimension parameter of array TAB
*               NEN  = number of entries already stored in table
*               TAB  = table of entries [TAB(NDIM)]
*               AUX  = auxiliary array of at least 2*NDIM words
*
*     The auxiliary array AUX is initialized. IF NEN not equal to zero,
*     the array is prepared for these entries. DTU has to be called,
*     if the table TAB has been modified by the user. The time for
*     initialization is about 3 mikro seconds for each entry in the
*     table [total (10 + 3 * NEN) micro seconds].
*
*     J = JTLU(ITEM,TAB,AUX)
*
*               ITEM = item to be searched for in the table
*               TAB  = table of entries
*               AUX  = auxiliary array
*
*     The returned value is the index J with ITEM = TAB(J). If the
*     entry was not found, it is added to the table and the correspon-
*     ding index is returned. J = 0 is returned, if there is not enough
*     space in the table (if NDIM = NEN).
*
*     The content of AUX may be checked at any time.
*     AUX(1) contains NDIM, and AUX(3) contains NEN. Thus if
*     AUX(1)=AUX(3), then there is no more space to store an
*     additional entry.
*
*     J = ITLU(ITEM,TAB,AUX)
*
*               ITEM = item to be searched for in the table
*               TAB  = table of entries
*               AUX  = auxiliary array
*
*     The returned value is the index J with ITEM = TAB(J), if the
*     entry was found, and J = 0, if no entry ITEM is found in the
*     table.
*
*     Author: V. Blobel Jan 88
*
* ----------------------------------------------------------------------

      implicit none
      INTEGER NDIM,NEN,I,J,K,NPRIM,ITEM
      INTEGER TAB(NDIM),AUX(2*NDIM)

*     NDIM at least 10, otherwise stop
      IF(NDIM.LT.10) THEN
        CALL HF_ERRLOG(12020101,'F: DTLU - too short array')
        CALL HF_STOP
      ENDIF

*     Clear aux array ...
      DO I=1,2*NDIM
         AUX(I)=0
      END DO
*     ... and define first three words
      AUX(1)=NDIM
      AUX(2)=0
*     determine prime number for hash function
      NPRIM=NDIM-3
      IF(MOD(NPRIM,2).EQ.0) NPRIM=NPRIM-1
   20 NPRIM=NPRIM-2
      IF(NPRIM.GT.5) THEN
         DO 30 I=3,INT(SQRT(FLOAT(NPRIM))),2
         J=NPRIM/I
         IF(I*J.EQ.NPRIM) GOTO 20
   30    CONTINUE
      END IF
      AUX(2)=NPRIM
*     Loop to insert index structure in AUX for the NEN existing
*     entries in the table
      DO 50 K=1,MIN(NEN,NDIM)
      ITEM=TAB(K)
*     search and add (like function JTLU)
      J=AUX(1)+MOD(IABS(ITEM),AUX(2))+1
   40 I=J
      J=AUX(I+5)
      IF(J.NE.0) THEN
         IF(ITEM.NE.TAB(J)) GOTO 40
         AUX(3)=AUX(3)+1
      ELSE
         J=AUX(3)+1
         AUX(3)=J
         AUX(I+5)=J
      END IF
   50 CONTINUE
  100 RETURN
      END

      INTEGER FUNCTION ITLU(ITEM,TAB,AUX)
*     -----------------------------------
      implicit none
      INTEGER TAB(*),AUX(*)
      INTEGER ITEM,KTLU

*     search loop in five statements
      KTLU=AUX(1)+MOD(IABS(ITEM),AUX(2))+1
    1 KTLU=AUX(KTLU+5)
      IF(KTLU.NE.0) THEN
         IF(ITEM.NE.TAB(KTLU)) GOTO 1
      END IF
*     end of search loop - result in KTLU
      ITLU=KTLU
      RETURN
      END


      INTEGER FUNCTION JTLU(ITEM,TAB,AUX)
*     -----------------------------------
      implicit none
      INTEGER TAB(*),AUX(*)
      INTEGER ITEM,I,J
  
      J=AUX(1)+MOD(IABS(ITEM),AUX(2))+1
   10 I=J
      J=AUX(I+5)
      IF(J.NE.0) THEN
*        repeat if not equal
         IF(ITEM.NE.TAB(J)) GOTO 10
      ELSE
*        not in list - add new antry in table TAB
         IF(AUX(1).EQ.AUX(3)) GOTO 100
         J=AUX(3)+1
         AUX(3)=J
         TAB(J)=ITEM
*        insert pointer
         AUX(I+5)=J
      END IF
  100 JTLU=J
      RETURN
      END
