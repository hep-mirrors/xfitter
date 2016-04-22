      SUBROUTINE StEwCpl
C                                                   -=-=- stewcpl

C 05-07-13 Mis-match of W+/W- corrected in the following lines.
C            QCPLN(2, -1, IPRTN) = CW * MOD(IS+1, 2)
C            QCPLN(3, -1, IPRTN) = CW * MOD(IS, 2)
C Mistake throughout these years! No consequences in CCFR experiments
C since always had W+ +/- W- 

C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: EwCpl1
C Description: Function to extract the Ewk couplings for Scattering and Annih.
C ListOfFiles: stewcpl ewcplg 
C ====================================================================
C Entry points (ewcplsc ewcplan)

C #Header: /Net/d2a/wkt/1hep/2ewk/RCS/EwCpl1.f,v 6.1 98/08/16 17:24:17 wkt Exp $
C #Log:	EwCpl1.f,v $
c Revision 6.1  98/08/16  17:24:17  wkt
c Re-organization; rationalization; initialization for DIS & DY
c 

C Set up the basic Electro-weak coupling matrices for the Boson-Fermion
C Yukawa coupling term in the Effective Lagrangian

C     8/16/98  wkt : migrated from DisPac over to EwkPac.

C    Really should be rewritten to use the GEwQx coefficents, 
C    especially if Z-exchanges are to be seriously incorported.
 
C    Leave the way it is for now for expediency!

C Use helicity basis for the boson polarization index in this section.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NSPN = 2, NGN = 3, NBSN = 4, NPOL = 3)
      PARAMETER (MXQRK = 6,  CH1 = D1/D3, CH2 = D2/D3)

      DIMENSION IP(NSPN, NGN), T3(NSPN), CHQ(NSPN), JBSN(NBSN)

      COMMON / WKCPLG / QCPLN(NBSN, -1:1, -MXQRK:MXQRK)

c ip(1,1) = 1, ip(1,2) = 4, ip(1,3) = 6 - charge 2/3 quarks, u,c,t
c ip(2,1) = 2, ip(2,2) = 3, ip(2,3) = 5 - charge -1/3 quarks, d,s,b

      DATA IP / 1, 2, 4, 3, 6, 5 /
      DATA T3, CHQ / 1, -1, 0.66666667, -0.3333333 /
      
      DATA JBSN / 1, 3, 2, 4 /     ! Swap W+/W- when making quark-anti-quark
                                   ! transformation in Do 60.

      CW = SQRT (2.)
ccpy 03/14/2012
      ID_WZ=1
      Call SetEwk(ID_WZ)
c Set up electroweak couplings
      DO 10 IPRTN = 1, MXQRK
C                     Longitudinal coupling vanishes for all bosons
         DO 20 IBSN = 1, 4
            QCPLN(IBSN, 0, IPRTN) = 0
            QCPLN(IBSN, 0,-IPRTN) = 0
 20      CONTINUE
C                     Right-handed coupling vanishes for Charged Bosons
         DO 30 IBSN = 2, 3
 30         QCPLN(IBSN, 1, IPRTN) = 0
 10   CONTINUE

      DO 40 IG = 1, NGN
         DO 50 IS = 1, NSPN
            IPRTN = IP(IS, IG)
C                                  Non-vanishing Photon couplings
            QCPLN(1, -1, IPRTN) = CHQ(IS)
            QCPLN(1,  1, IPRTN) = CHQ(IS)
C                                 Non-vanishing Charged boson couplings
            QCPLN(2, -1, IPRTN) = CW * MOD(IS+1, 2)
            QCPLN(3, -1, IPRTN) = CW * MOD(IS, 2)
C                                 Neutral Z couplings
            QCPLN(4, -1, IPRTN) = T3(IS) - 2.*CHQ(IS)*SWG2F()
            QCPLN(4,  1, IPRTN) = -2.*CHQ(IS)*SWG2F()
c swg2f returns sin^2(theta_weinberg)
C                                 Anti-quark couplings - by CP
            DO 60 ib = 1, 4
c iterate over ivb
               QCPLN(JBSN(IB), -1, -IPRTN) = QCPLN(IB,  1, IPRTN)
               QCPLN(JBSN(IB),  1, -IPRTN) = QCPLN(IB, -1, IPRTN)
 60            CONTINUE

 50      CONTINUE
 40   CONTINUE

      RETURN
C                        ****************************
      END

