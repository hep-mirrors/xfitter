      SUBROUTINE CNVL1 (IRDR, JRDR, NF, VLAM)
C                                                   -=-=- cnvl1
C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: Setaux
C Description: Auxilary functions for lambda and alpha conversions
C ListOfFiles: cnvl1 zcnvlm trnlam zbrlam rtalf
C ====================================================================

C #Header: /Net/d2a/wkt/1hep/2qcd/RCS/Setaux.f,v 1.1 97/12/21 20:34:57 wkt Exp $ 
C #Log:	Setaux.f,v $
c Revision 1.1  97/12/21  20:34:57  wkt
c Initial revision
c 

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL ZCNVLM
C  Auxiliary routine for SETLAM.
C                  Given Lamda (NF) = VLAM at order IRDR, this subroutine
C                  finds the corresponding Lamda at order JRDR so that the
C                  resulting Alpha(Mu, NRDR) remains approximately the same
C                  within the range of Neff = NF.
 
      COMMON / LAMCNV / AMU, ULAM, NFL, IRD, JRD
      COMMON / IOUNIT / NIN, NOUT, NWRT
 
      DATA ALM, BLM, ERR, AMUMIN / 0.001, 2.0, 0.02, 1.5 /
 
      IRD = IRDR
      JRD = JRDR
      ULAM = VLAM
 
      CALL PARQCD(2, 'NFL', ANF, IRT)
      NTL = NFLTOT()
      IF (NF .GT. NTL) THEN
         WRITE (NOUT, *) ' NF .GT. NTOTAL in CNVL1; set NF = NTOTAL'
         WRITE (NOUT, *) ' NF, NTOTAL = ', NF, NTL
         NF = NTL
      ENDIF
C                                                First match at NFth threshold
      NFL = NF
      AMU = AMHATF(NF)
      AMU = MAX (AMU, AMUMIN)
      VLM1 = QZBRNT (ZCNVLM, ALM, BLM, ERR, IR1)
C                                          Match again at the next threshold
      IF (NF .LT. NTL) THEN
        AMU = AMHATF(NF+1)
        AMU = MAX (AMU, AMUMIN)
        VLM2 = QZBRNT(ZCNVLM, ALM, BLM, ERR, IR2)
      ELSE
        VLM2 = VLM1
      ENDIF
C                              Take the average and return new value of VLAM
      VLAM = (VLM1 + VLM2) / 2
 
      RETURN
C                        ****************************
      END
 
