C   25/01/91 605221826  MEMBER NAME  DELTAR   (EPRC91.S)    FVS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        DELTA-R AND EFFECTIVE WEAK MIXING ANGLE FROM ZNCVS
C
C   AUTHOR: H.SPIESBERGER,
C           II. INSTITUT FUER THEORETISCHE PHYSIK
C           UNIVERSITAET HAMBURG
C   LAST MODIFIED: 07.02.1991
C
C---->  
C       Date: 13.October 1998
C       Mods: Added Complex*16 SIGFMS in PIMQQ, renamed variables 
C             containing "$" -> "_" 
C---->  
*
*       Date: 27.May 2012
*       Modified for xFitter by HS
*
C---->  
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCC 22. 9. 86   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC H. SPIESBERGER CCC
C

      subroutine EPRC_INIT(doprint)

      IMPLICIT REAL*8(A-H,M,O-Z)
      logical doprint
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
cv
      double precision pi_hf, alphaem_hf, gf_hf, convfac_hf
      double precision mw_hf, mz_hf, mh_hf, mel_hf, mup_hf
      double precision mdn_hf, mst_hf, mch_hf, mbt_hf
      double precision mtp_hf, mta_hf, mmo_hf

      DO 1 I=1,20
    1 LPAR(I)=0
      DO 2 I=1,12
    2 LPARIN(I)=0
C
C...PRINT TITLE
      if (doprint) then
      WRITE(6,2000)
 2000 FORMAT(/,' ******************************************************'
     F        ,' ',/
     F        ,'   PROGRAM: EPRC91.S(DELTAR)               JAN.91',/
     F        ,'                                    H.SPIESBERGER',/
     F        ,'   CALCULATION OF DELTA-R AS FUNCTION OF  ',/
     F        ,'   M-TOP AND M-HIGGS'
     F        ,' ',/
     F        ,' ******************************************************'
     F        ,/)
      endif
      
C...OPTION FOR THE DETERMINATION OF MW
C   LPAR(4)=1 -> W-mass is input, Gmu is calculated 
C                - use this option for Wmass-fit
C   LPAR(4)=2,3 -> fixed Gmu, W-mass is calculated 
C                - use this option for PDF fits
      LPAR(4)=1

C...QCD CORRECTIONS IN DELTA-R
      LPAR(5)=2
C...Jegerlehner's PARAMETRIZATION FOR VACUUM POLARIZATION: hadr5n12.f
      LPAR(7)=3
cv      LPAR(7)=0
C...SELF ENERGY CORRECTIONS
      LPAR(8)=1
      LPAR(9)=1
      LPAR(10)=1
C...NO QED, BUT WEAK CONTRIBUTIONS IN THE ELECTRON AND QUARK VERTEX
      LPAR(11)=0
      LPAR(12)=1
      LPAR(13)=1
      LPAR(14)=0
C...PURELY WEAK CONTRIBUTIONS
      LPAR(15)=1


      call wrap_constants(pi_hf, alphaem_hf, gf_hf, convfac_hf,
     $     mw_hf, mz_hf, mh_hf, mel_hf, mup_hf,
     $     mdn_hf, mst_hf, mch_hf, mbt_hf, mtp_hf, mta_hf, mmo_hf)


      CALL SETPAR(doprint)

      if (doprint) write(6,*) 'end of EPRC_INIT dr = ',deltar

      return
      end



* ---------------------------------------
      subroutine eprc_effective(q2)
* ---------------------------------------
      IMPLICIT REAL*8(A-H,M,O-Z)

      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HDELTR/ DELTAR,AGF0,DRHOT,DALPMZ,XGMT,ALPQCD,BTOP4,DRPIW2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      COMMON /FORMFF/ ALPFFQ,AKAPPA,GMUFFQ,SWEFF2

      COMPLEX*16 SIGMRG,SIGMRM,CG,CM

cv
      COMMON /PARAM/  POLARI,LLEPT,LQUA
      common/eprc_manu/sm
C...EFFECTIVE WEAK MIXING ANGLE
      LPAR(15)=1

      T=-Q2
      CG=SIGMRG(T)
      PIGGG=DREAL(CG)/T
      CM=SIGMRM(T)
      SM=DREAL(CM)/T
      AKAPPA=1D0-CW/SW*SM/(1D0+PIGGG)
      SWEFF2=SW2*AKAPPA

      return
      end



************************************************************************

************************************************************************

      SUBROUTINE SETPAR(doprint)
C---SETTING OF ELECTROWEAK PARAMETERS
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      logical doprint
      COMPLEX*16 SIGMQW,SIGMQZ
      COMPLEX*16 CMW2,CMZ2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      COMMON /PARAM/  POLARI,LLEPT,LQUA
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /CBMASS/ CMW2,CMZ2
      COMMON /HDELTR/ DELTAR,AGF0,DRHOT,DALPMZ,XGMT,ALPQCD,BTOP4,DRPIW2
      COMMON /GSW1/   MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /SMCON/  VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /KNSTCC/ SXNRCC,SX1NCC
      COMMON /RSCALE/ MSC2
      COMMON /DNONSM/ DNSM1,DNSM2,DNSMR
C
      double precision epMz,epsin2thw,cau,cvu,cad,cvd
      common/couplings/epsin2thw,epMz,
     +       cau,cvu,cad,cvd

cv
      double precision pi_hf1, alphaem_hf1, gf_hf1, convfac_hf1
      double precision mw_hf1, mz_hf1, mh_hf1, mel_hf1, mup_hf1
      double precision mdn_hf1, mst_hf1, mch_hf1, mbt_hf1
      double precision mtp_hf1, mta_hf1, mmo_hf1


      call wrap_constants(pi_hf1, alphaem_hf1, gf_hf1, convfac_hf1, 
     $     mw_hf1, mz_hf1, mh_hf1, mel_hf1, mup_hf1,
     $     mdn_hf1, mst_hf1, mch_hf1, mbt_hf1, mtp_hf1, mta_hf1, 
     $     mmo_hf1)


      PI=pi_hf1
      ALPHA=alphaem_hf1

      ALP1PI=ALPHA/PI
      ALP2PI=ALPHA/2D0/PI
      ALP4PI=ALPHA/4D0/PI
      E=DSQRT(4D0*PI*ALPHA)
cv      GF=1.166389D-5
      GF=gf_hf1

      AGF0=PI*ALPHA/GF/DSQRT(2D0)
cv      SXNORM=PI*ALPHA*ALPHA/2D0*3.8938D5
cv      SX1NRM=ALPHA*ALPHA*ALPHA/16D0/PI*3.8938D5
      SXNORM=PI*ALPHA*ALPHA/2D0*convfac_hf1*1.D3
      SX1NRM=ALPHA*ALPHA*ALPHA/16D0/PI*convfac_hf1*1.D3


C---NON-STANDARD PHYSICS
      DNSM1=0D0
      DNSM2=0D0
      DNSMR=0D0
C     DNSM1=-1D-2
C     DNSM2=-1D-2
C     DNSMR=-1D-2

C---DEFINE PARAMETERS OF THE ELECTROWEAK STANDARD MODEL
cv      ME=.51099906D-3
cv      MMY=.105658387D0
cv      MTAU=1.77682D0
cv      MU=.067D0
cv      MD=.089D0
cv      MS=.231D0
cv      MC=1.299D0
cv      MB=4.5D0
cv      mc=1.4d0
cv      mb=4.75d0

      MW=mw_hf1
      MZ=mz_hf1
      MH=mh_hf1
      Me=mel_hf1
      MU=mup_hf1
      MD=mdn_hf1
      MS=mst_hf1
      MC=mch_hf1
      MB=mbt_hf1
      MT=mtp_hf1
      MTAU=mta_hf1
      MMY=mmo_hf1


      MW2=MW*MW
      MZ2=MZ*MZ
      MH2=MH*MH
      ME2=ME*ME
      MMY2=MMY*MMY
      MTAU2=MTAU*MTAU
      MU2=MU*MU
      MD2=MD*MD
      MS2=MS*MS
      MC2=MC*MC
      MB2=MB*MB
      MT2=MT*MT

C...RENORMALIZATION SCALE
      MSC2=MZ2
C
      IF (LPAR(4).GT.1) THEN
C...Calculate MW from alpha_em, G_mu and MZ
        MW=PARGFX()
      ELSE
        Mw=mw_hf1
C...Call this to get Delta-r, but keep MW from input
        MW=PARGFX()
        Mw=mw_hf1
      ENDIF
cv      Mw = 80.3980d0

      MW2=MW*MW
      CW=MW/MZ
      CW2=CW*CW
cv this should produce 0.2224 for SW2
      SW2=1D0-CW2
      SW=DSQRT(SW2)
      WWIDTH=DIMAG(SIGMQW(MW2))/MW
      ZWIDTH=DIMAG(SIGMQZ(MZ2))/MZ
      CMW2=MW*DCMPLX(MW,-WWIDTH)
      CMZ2=MZ*DCMPLX(MZ,-ZWIDTH)

C---NORMALIZATION OF NC AND CC CROSS SECTIONS
C---ON-MASS SHELL SCHEME
      IF (LPAR(4).EQ.1) THEN
        B=1D0/4D0/CW/SW
        SXNRCC=SXNORM/SW2/SW2
        SX1NCC=SX1NRM/SW2/SW2*8D0
        ELSE
C---MODIFIED ON-MASS SHELL SCHEME (NORMALIZATION TO G-MU)
        B=MZ/SQRT(AGF0)/4D0
        SXNRCC=SXNORM*MW2*MW2/AGF0/AGF0
        SX1NCC=SX1NRM*8D0*MW2*MW2/AGF0/AGF0
        IF (LPAR(4).EQ.3) B=B*SQRT(1D0-DELTAR)
      ENDIF
C---DEFINE FERMION GAUGE BOSON COUPLING CONSTANTS
      LBOSON=LPAR(17)
      VAFI(2,1,1)=0D0
      VAFI(2,2,1)=0D0
      VAFI(2,3,1)=0D0
      IF (LBOSON.LE.2) THEN
        VAFI(1,1,1)=1D0
        VAFI(1,2,1)=-2D0/3D0
        VAFI(1,3,1)=1D0/3D0
        ELSE
        VAFI(1,1,1)=0D0
        VAFI(1,2,1)=0D0
        VAFI(1,3,1)=0D0
      ENDIF
      IF (LBOSON.GE.2.OR.LBOSON.EQ.0) THEN
cep        VAFI(2,1,2)=-B
cep        VAFI(2,2,2)=B
cep        VAFI(2,3,2)=-B
cep        VAFI(1,1,2)=B*(4D0*SW2-1D0)
cep        VAFI(1,2,2)=B*(1D0-8D0*SW2/3D0)
cep        VAFI(1,3,2)=B*(4D0*SW2/3D0-1D0)
cep
        VAFI(2,1,2)=-B
        VAFI(2,2,2)=B*2*cau
        VAFI(2,3,2)=B*2*cad
        VAFI(1,1,2)=B*(4D0*SW2-1D0)
        VAFI(1,2,2)=B*2*cvu
        VAFI(1,3,2)=B*2*cvd
        ELSE
        VAFI(2,1,2)=0D0
        VAFI(2,2,2)=0D0
        VAFI(2,3,2)=0D0
        VAFI(1,1,2)=0D0
        VAFI(1,2,2)=0D0
        VAFI(1,3,2)=0D0
      ENDIF

C----------------------------------------------------------------------
C---MASSES USED FOR HARD BREMSSTRAHLUNG KINEMATIC
      MEI = ME
      MEF = ME
C---QUARK MASSES USED AS REGULATORS IN THE LEPTONIC BREMSSTRAHLUNG
C---DO NOT USE THIS VERSION FOR VERY SMALL X
      MPRO=938.28D-3
      MPRO2=MPRO*MPRO
C     MQI = MU
C     MQF = MU
      MQI = 2D0
      MQF = 2D0
      MEF2 = MEF*MEF
      MQF2 = MQF*MQF
      MEI2 = MEI*MEI
      MQI2 = MQI*MQI
C---COUPLING CONSTANTS
      IGAMMA = 1
      IZ = 2
      IEL = 1
      IFU = 2
      IFD = 3
      INDV = 1
      INDA = 2
C
      DO 1 IF=IEL,IFD
        DO 1 IB1=IGAMMA,IZ
          DO 1 IB2=IGAMMA,IZ
          IF (LBOSON.EQ.2.AND.IB1.EQ.IB2) THEN
             FLIND(INDV,IF,IB1,IB2)=0D0
             GOTO 1
          ENDIF
          FLIND(INDV,IF,IB1,IB2)=
     *     2D0*(VAFI(INDV,IF,IB1)*VAFI(INDV,IF,IB2)
     *         +VAFI(INDA,IF,IB1)*VAFI(INDA,IF,IB2))
    1 CONTINUE
      DO 2 IF=IEL,IFD
        DO 2 IB1=IGAMMA,IZ
          DO 2 IB2=IGAMMA,IZ
          IF (LBOSON.EQ.2.AND.IB1.EQ.IB2) THEN
             FLIND(INDA,IF,IB1,IB2)=0D0
             GOTO 2
          ENDIF
          FLIND(INDA,IF,IB1,IB2)=
     *     2D0*(VAFI(INDV,IF,IB1)*VAFI(INDA,IF,IB2)
     *         +VAFI(INDA,IF,IB1)*VAFI(INDV,IF,IB2))
    2 CONTINUE
C
      DO 3 IVB1 = IGAMMA, IZ
       DO 3 IVB2 = IGAMMA, IZ
        DO 3 IFERM = IFU, IFD
        IF (LBOSON.EQ.1.AND.(IVB1.EQ.2.OR.IVB2.EQ.2)) THEN
          AFIJ(IFERM,IVB1,IVB2)=0D0
          BFIJ(IFERM,IVB1,IVB2)=0D0
          GOTO 3
        ELSEIF (LBOSON.EQ.3.AND.(IVB1.EQ.1.OR.IVB2.EQ.1)) THEN
          AFIJ(IFERM,IVB1,IVB2)=0D0
          BFIJ(IFERM,IVB1,IVB2)=0D0
          GOTO 3
        ELSEIF (LBOSON.EQ.2.AND.IVB1.EQ.IVB2) THEN
          AFIJ(IFERM,IVB1,IVB2)=0D0
          BFIJ(IFERM,IVB1,IVB2)=0D0
          GOTO 3
        ENDIF
        AFIJ(IFERM,IVB1,IVB2)=FLIND(INDV,IFERM,IVB1,IVB2)
     *  *(FLIND(INDV,IEL,IVB1,IVB2) - POLARI*FLIND(INDA,IEL,IVB1,IVB2))
        BFIJ(IFERM,IVB1,IVB2)=FLIND(INDA,IFERM,IVB1,IVB2)
     *  *(FLIND(INDA,IEL,IVB1,IVB2) - POLARI*FLIND(INDV,IEL,IVB1,IVB2))
    3 CONTINUE
C
      LUNOUT=6
C
C---PRINT TITLE
      if (doprint) then
      WRITE(LUNOUT,'(///15(A/))')
     * '##############################################################',
     * '#                                                            #',
     * '#                          EPRC93                            #',
     * '#                                                            #',
     * '#              ELECTROWEAK RADIATIVE CORRECTIONS             #',
     * '#     FOR DEEP INELASTIC LEPTON PROTON SCATTERING AT HERA    #',
     * '#                                                            #',
     * '#                      OCTOBER 1991                          #',
     * '#            (REV.: MARCH 1992 AND OCT. 1995)                #',
     * '#                                                            #',
     * '#                            BY                              #',
     * '#                      H. SPIESBERGER                        #',
     * '#                                                            #',
     * '##############################################################'
C
C---Print theory parameters
      WRITE(LUNOUT,'(///A/)')
     *    ' *****  PARAMETERS FOR EL-WEAK THEORY  *****'
      WRITE(LUNOUT,'(10X,A,I2,10X,A/29X,A/29X,A/30X,A/)')
     *    ' IFIX =', LPAR(4), ' IFIX=1 : MW FIXED',
     *                        '     =2 : GF FIXED',
     *                        '     =3 : GF FIXED,',
     *                        '         AND DELTA-R IN BORN INCLUDED '
      WRITE(LUNOUT,212) MW,MZ,SW2
      WRITE(LUNOUT,213) WWIDTH,ZWIDTH
      WRITE(LUNOUT,214) MU,MC,MD,MB,MS,MT
      WRITE(LUNOUT,2104) ME,MMY,MTAU
      WRITE(LUNOUT,215) MH
      WRITE(LUNOUT,2016) DELTAR
      endif

212   FORMAT(10X,' BOSON MASSES:   MW = ',F10.4,' GEV,',/
     F       10X,'                 MZ = ',F10.4,' GEV,   SW2 = ',F10.4)
213   FORMAT(10X,' BOSON WIDTHS:   GW = ',F10.4,' GEV,    GZ = ',F10.4)
214   FORMAT(10X,
     F       ' FERMION MASSES: MU = ',F10.4,' GEV,  MC = ',F10.4,' GEV',
     F /,10X,'                 MD = ',F10.4,' GEV,  MB = ',F10.4,' GEV',
     F /,10X,'                 MS = ',F10.4,' GEV,  MT = ',F10.4,' GEV')
2104  FORMAT(10X,
     F       ' LEPTON MASSES:  ME = ',F17.12,' GEV',
     F /,10X,'                MMY = ',F17.12,' GEV',
     F /,10X,'               MTAU = ',F17.12,' GEV')
215   FORMAT(10X,' HIGGS MASS:     MH = ',F10.4,' GEV')
2016  FORMAT(10X,' DELTA-R            = ',F10.6)
C
C---OPTIONS FOR VIRTUAL&SOFT CONTRIBUTION
        if (doprint) then
        WRITE(LUNOUT,'(//A/)')
     *              ' *****  OPTIONS FOR VIRT&SOFT CONTRIBUTION  *****'
        WRITE(LUNOUT,216) LPAR
        endif

216     FORMAT(/,' PARAMETER LIST',/
     F,' *****************************************************',/
     F,' **      BORN CROSS SECTION:    LPAR( 1) = ',I4,'     **',/
     F,' **      1-LOOP CORRECTIONS:    LPAR( 2) = ',I4,'     **',/
     F,' **      EXPONENTIATION:        LPAR( 3) = ',I4,'     **',/
     F,' **      MW OR GMU FIXED:       LPAR( 4) = ',I4,'     **',/
     F,' **      QCD CORR IN DELTA-R:   LPAR( 5) = ',I4,'     **',/
     F,' **      PARTON DISTRIB.:       LPAR( 6) = ',I4,'     **',/
     F,' **      SIGMA-GAMMA:           LPAR( 7) = ',I4,'     **',/
     F,' **      SIGMA-GAMMA-Z:         LPAR( 8) = ',I4,'     **',/
     F,' **      SIGMA-Z:               LPAR( 9) = ',I4,'     **',/
     F,' **      SIGMA-W:               LPAR(10) = ',I4,'     **',/
     F,' **      QED CORRECTIONS:       LPAR(11) = ',I4,'     **',/
     F,' **      LEPTONIC QED:          LPAR(12) = ',I4,'     **',/
     F,' **      HADRONIC QED:          LPAR(13) = ',I4,'     **',/
     F,' **      LEPTON-QUARK-INTRF:    LPAR(14) = ',I4,'     **',/
     F,' **      WEAK CORRECTIONS:      LPAR(15) = ',I4,'     **',/
     F,' **      WEAK BOXES:            LPAR(16) = ',I4,'     **',/
     F,' **      GAMMA OR/AND Z         LPAR(17) = ',I4,'     **',/
     F,' **      LEPTONIC LLA INCL.     LPAR(18) = ',I4,'     **',/
     F,' **      X/Y DEFINITION         LPAR(19) = ',I4,'     **',/
     F,' **      FIXED Q2 IN HARD BS    LPAR(20) = ',I4,'     **',/
     F,' *****************************************************'/)
C
C---NON-STANDARD PHYSICS
      if (doprint) then
      WRITE(LUNOUT,'(///A/)')
     *    ' *****  NON-STANDARD CONTRIBUTIONS TO DELTAS *****'
      WRITE(LUNOUT,'(10X,A,F8.4)')
     *      ' DELTA-RHO(0) =',DNSMR
      WRITE(LUNOUT,'(10X,A,F8.4)')
     *      ' DELTA-1      =',DNSM1
      WRITE(LUNOUT,'(10X,A,F8.4)')
     *      ' DELTA-2      =',DNSM2
C
      WRITE(LUNOUT,'(///A///)')
     *' ******************* END OF INITIALIZATION *********************'
      endif

      RETURN
      END
 
************************************************************************
 
C...Running fine structure constant using Burkhard's parametrization
      FUNCTION AEMRUN(Q2)
      IMPLICIT NONE
      REAL*8 AEMRUN, Q2, HADRQQ
      COMPLEX*16 HSSRGG,F
      REAL*8 PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      REAL*8 SW,CW,SW2,CW2,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *      ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
 
      IF (Q2.EQ.0D0) THEN
        AEMRUN=ALPHA
        RETURN
       ELSE
        HSSRGG=ALPHA/PI*(
     &   +(+(Q2+2D0*ME2  )*F(Q2,ME  ,ME  )
     &     +(Q2+2D0*MMY2 )*F(Q2,MMY ,MMY )
     &     +(Q2+2D0*MTAU2)*F(Q2,MTAU,MTAU) - Q2)/3D0
     &   +(+(Q2+2D0*MT2  )*F(Q2,MT  ,MT  ) - Q2/3D0)/2.25D0     )
     &   - DCMPLX(HADRQQ(Q2)*Q2,0D0)
        AEMRUN=ALPHA/(1D0+DREAL(HSSRGG)/Q2)
      ENDIF
 
      RETURN
      END

************************************************************************

      FUNCTION PARGFX()
C
C   DETERMINATION OF MW FOR GIVEN MZ FROM G_MU
C   (LAST MOD HS 20.11.92): HIGHER ORDERS IMPROVED
C
C   INPUT:   Z0 MASS = MZ, HIGGS MASS = MH, TOP MASS= MT
C   OUTPUT:  W MASS = MW,
C   OUTPUT ON COMMONS: DELTA-R = DELTAR, WEAK MIXING ANGLE =SW2
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      DIMENSION DR2(40)
      COMPLEX*16 SIGMQG,SIGFWS,SIGFZS,SIGFMS
      COMPLEX*16 PIZMZ,PIMMZ,PIW0
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      COMMON /HDELTR/ DELTAR,AGF0,DRHOT,DALPMZ,XGMT,ALPQCD,BTOP4,DRPIW2
      COMMON /DLTASM/ DALMZ,DRHOTT,DRHO0,DLTA1,DLTA2,DLTAR
      COMMON /CTMASS/ RDMZ2,RDMW2,PIZMZ,PIMMZ,PIW0

      SQ=DSQRT(1D0-4D0*AGF0/MZ2)
      S02=(1D0-SQ)/2D0
      C02=1D0-S02
      DALPMZ=0.0602D0
C     DALPMZ=-REAL(SIGMQG(MZ2))/MZ2
      XGMT=GF*MT2/SQRT(2D0)/8D0/PI/PI
      DRHOTT=3D0*XGMT
      ALAMB2=0.14D0**2D0
      NF=5
C     IF (LPAR(5).GE.1) DRHOTT=DRHOTT*(1D0+XGMT*(19D0-2D0*PI*PI))
C...NEW: DELTA-RHO OF BARBIERI ET AL. FOR ANY MHIGGS
C   THE RHO PARAMETER
      XHT=MH/MT
      R=1/XHT**2
      ALR=DLOG(R)
      H=49D0/4D0+PI**2+27D0/2D0*ALR+1.5D0*ALR**2
     .  +R/3D0*(2D0-12D0*PI**2+12D0*ALR-27D0*ALR**2)
     .  +R*R/48D0*(1613D0-240D0*PI**2-1500D0*ALR-720D0*ALR**2)
      A0=19D0-2D0*PI**2
      X4=4D0
      T=1D0/X4**2
      ALT=DLOG(T)
      H4=49D0/4D0+PI**2+27D0/2D0*ALT+1.5D0*ALT**2
     .  +T/3D0*(2D0-12D0*PI**2+12D0*ALT-27D0*ALT**2)
     .  +T*T/48D0*(1613D0-240D0*PI**2-1500D0*ALT-720D0*ALT**2)
      DH=(27D0/2D0+3D0*ALT)/T
     .   +4D0-18D0*ALT
     .  +1D0/3D0*(2D0-12D0*PI**2+12D0*ALT-27D0*ALT**2)
     .  +T/24D0*(1613D0-240D0*PI**2-1500D0*ALT-720D0*ALT**2)
     .   -T/48D0*(1500D0+2D0*720D0*ALT)
      DH4=-DH*2D0/X4**3
      A1=-(2.5D0+A0)*5D0
      A3=(H4-2D0*DH4-A0-2D0*A1)/(-32D0)
      A2=(DH4-A1-48D0*A3)/8D0
      IF(XHT.GT.2D0) THEN
      R2=H
             ELSE
      R2=A0+A1*XHT+A2*XHT**2+A3*XHT**3
         ENDIF
C  R2 ERSETZT  19-2*PI**2
      IF (LPAR(5).GE.1) DRHOTT=DRHOTT*(1D0+XGMT*R2)
      DRHOT2=DRHOTT
      IF (LPAR(5).GE.2) THEN
        ALPQCD=4D0*PI/(11D0-2D0/3D0*NF)/DLOG(MT2/ALAMB2)
        DRHOTT=DRHOTT-ALPQCD/PI*(PI*PI/3D0+1D0)
     *       *ALP4PI/2D0/S02/C02*MT2/MZ2
      ENDIF
      IF (LPAR(5).GE.1) THEN
        DELTAR=1D0-(1D0-DALPMZ)*(1D0+C02/S02*DRHOTT)
        ELSE
        DELTAR=DALPMZ-C02/S02*DRHOTT
      ENDIF
C
      I=1
C---SW2 = SIN**2 THETA-W, START WITH APPROXIMATE VALUE
      IF (LPAR(4).EQ.1) GOTO 50
      SQ=DSQRT(1D0-4D0*AGF0/MZ2/(1D0-DELTAR))
      SW2=(1D0-SQ)/2D0
      DR2(1)=0D0
      CW2=1D0-SW2
      SW=DSQRT(SW2)
      CW=DSQRT(CW2)
      MW=MZ*CW
      MW2=MW*MW
      GOTO 51
50    CW=MW/MZ
      CW2=CW*CW
      SW2=1D0-CW2
      SW=DSQRT(SW2)
      MW2=MW*MW

C---DELTAR MEANS THE QUANTITY 'DELTA R' IN MU LIFETIME FORMULA
51    CONTINUE
      LPR15K=LPAR(15)
      LPAR(15)=1
      CALL SIGMR0
      LPAR(15)=LPR15K
      IF (LPAR(5).GE.2) THEN
        ALPQCD=4D0*PI/(11D0-2D0/3D0*NF)/DLOG(MT2/ALAMB2)
        DRHOTT=DRHOT2-ALPQCD/PI*(PI*PI/3D0+1D0)
     *       *ALP4PI/2D0/SW2/CW2*MT2/MZ2

      ENDIF
      IF (LPAR(5).GE.1) THEN
        DELTAR=1D0-(1D0-DALPMZ)*(1D0+CW2/SW2*DRHOTT)
C       DELTAR=DALPMZ-CW2/SW2*DRHOTT
     *        +dREAL(SIGFWS(0D0)-SIGFWS(MW2))/MW2
     *        +2D0*CW/SW*SIGFMS(0D0)/MZ2+FSGS(0D0)-DALPMZ
     *        -CW2/SW2*(dREAL(SIGFZS(MZ2)/MZ2-SIGFWS(MW2)/MW2)
     *                  -3D0/SW2*ALP4PI/4D0*MT2/MW2)
     *            *DSQRT(2D0)*GF*MW2*(1D0-DALPMZ)/PI/ALPHA*SW2
     *        +ALP4PI/SW2*(6D0+(3.5D0/SW2-2D0)*DLOG(CW2))
        IF (LPAR(7).GE.2) DELTAR=DELTAR-DSGMRG(MZ2)
        ELSE
        DELTAR=-dREAL(PIW0)+ALP4PI/SW2*(6D0+(3.5D0/SW2-2D0)*DLOG(CW2))
      ENDIF
      SQ=DSQRT(1D0-4D0*AGF0/MZ2/(1D0-DELTAR))
      IF (LPAR(4).EQ.1) GOTO 60

C---THE CORRECTED VALUE FOR SIN**2 THETA-W
      SW2=(1D0-SQ)/2D0
      DR2(I+1)=DELTAR
      DS2=DABS(DR2(I+1)-DR2(I))
      CW2=1D0-SW2
      CW=DSQRT(CW2)
      SW=DSQRT(SW2)
C---THE CORRECTED VALUE FOR THE W MASS
      MW=MZ*CW
      MW2=MW*MW
C---DELTA-RHO
      BTOP4=1D0/(1D0-DRHOTT)/(1D0-DELTAR)
      BTOP4=SQRT(BTOP4)

C     WRITE(6,200) I,SW2,DELTAR
  200 FORMAT(I3,2F20.10)
      IF(DS2.LT.1D-8) THEN
         PARGFX=MW
         GOTO 60
      END IF
      I=I+1
      IF(I.LE.39) GOTO  51
      PARGFX=MW

60    CONTINUE
      DRHOT=3D0*ALP4PI/4D0*MT2/SW2/MW2
      DRPIW2=ALP4PI/SW2*(6D0+(3.5D0/SW2-2D0)*DLOG(CW2))
C     WRITE(6,99) DRHOT,DRPIW2
   99 FORMAT(/,' DRHOT = ',F15.5,/,' DRPIW = ',F15.5,/)

C...CALCULATION FROM UNRENORMALIZED SELF ENERGIES INCLUDING
C   SINGULAR CONTRIBUTIONS (INDEPENDENCE ON SCALE MSC2)
C     LPRK15=LPAR(15)
C     LPAR(15)=0
      DALMZ=-dREAL(SIGMQG(MZ2))/MZ2
      ALPMZ=ALPHA*(1D0+DALMZ)
      DRHO0=dREAL(SIGFZS(0D0))/MZ2-dREAL(SIGFWS(0D0))/MW2
     *        -SW/CW*2D0*dREAL(SIGFMS(0D0))/MZ2
      DLTA1=CW2*(FSGS(MZ2)-FSZS(MZ2))-CW/SW*(CW2-SW2)*FSMS(MZ2)
      DLTA2=CW2*FSZS(MZ2)-FSWS(MW2)+SW2*FSGS(MZ2)-2D0*SW*CW*FSMS(MZ2)
      DLTAR=DALMZ-CW2/SW2*DRHO0-(CW2-SW2)/SW2*DLTA2+2D0*DLTA1
     *       +ALP4PI/SW2*(6D0+(7D0-4D0*SW2)/2D0/SW2*DLOG(CW2))
C     WRITE(6,200) DRHO0,DLTA1,DLTA2
C 200 FORMAT(8X,F10.5,10X,2F10.5)
C     LPAR(15)=LPRK15
C
C...CHECK 2: CALCULATION FROM RENORMALIZED SELF ENERGIES
C     DPIZ=dREAL(PIZF(MZ2))-MZ2*FSZD(MZ2)
C     DR0REN=-dREAL(SIGMRW(0D0))/MW2-DPIZ
C     DL2REN=SW2*dREAL(SIGMRG(MZ2))/MZ2+CW2*DPIZ
C    *       -2D0*CW*SW*dREAL(SIGMRM(MZ2))/MZ2+REAL(SIGMRW(0D0))/MW2
C     DL1REN=CW2*dREAL(SIGMRG(MZ2)/MZ2-DPIZ)
C    *       -CW/SW*(CW2-SW2)*dREAL(SIGMRM(MZ2))/MZ2
C     WRITE(6,200) DR0REN,DL1REN,DL2REN
C
C...COMPARISON WITH ALTARELLI, BARBIERI, JADACH:
C     S02=(1D0-DSQRT(1D0-4D0*PI*ALPMZ/DSQRT(2D0)/GF/MZ2))/2D0
C     S02=0.23146D0
C     C02=1D0-S02
C     S0=DSQRT(S02)
C     C0=DSQRT(C02)
C     DELTRW=DLTAR-DALMZ
C     DELTRW=1D0-(1D0-DLTAR)/(1D0-DALMZ)
C     DKAPSZ=-C0/S0*dREAL(SIGMQM(MZ2))/MZ2+C02/(C02-S02)*DELTRW
C     EPS2=C02*DRHO0+S02/(C02-S02)*DELTRW-2D0*S02*DKAPSZ
C     EPS3=C02*DRHO0+(C02-S02)*DKAPSZ
C     WRITE(6,201)        EPS3,-EPS2
C 201 FORMAT(8X,20X,2F10.5)

      END

************************************************************************

      FUNCTION DLTRMM(MWA,MZA,MTA,MHA)
C
C   DETERMINATION OF DELTA_R FOR FIXED MASSES
C   (HS 07.08.98)
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 SIGFWS,SIGFZS,SIGFMS
      COMPLEX*16 PIZMZ,PIMMZ,PIW0
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      COMMON /HDELTR/ DELTAR,AGF0,DRHOT,DALPMZ,XGMT,ALPQCD,BTOP4,DRPIW2
      COMMON /DLTASM/ DALMZ,DRHOTT,DRHO0,DLTA1,DLTA2,DLTAR
      COMMON /CTMASS/ RDMZ2,RDMW2,PIZMZ,PIMMZ,PIW0

      MW=MWA
      MZ=MZA
      MT=MTA
      MH=MHA
      MW2=MW*MW
      MZ2=MZ*MZ
      MT2=MT*MT
      MH2=MH*MH
      C0=MW/MZ
      C02=C0*C0
      S02=1D0-C02
      S0=DSQRT(S02)
      CW=C0
      SW=S0
      SW2=S02
      CW2=C02
      DALPMZ=0.0602D0

C...iterate on GF for the calculation of 2nd order terms
      IGFC=0
 1    CONTINUE
      IGFC=IGFC+1

      XGMT=GF*MT2/SQRT(2D0)/8D0/PI/PI
      DRHOTT=3D0*XGMT
      ALAMB2=0.14D0**2D0
      NF=5
C     IF (LPAR(5).GE.1) DRHOTT=DRHOTT*(1D0+XGMT*(19D0-2D0*PI*PI))
C...NEW: DELTA-RHO OF BARBIERI ET AL. FOR ANY MHIGGS
C   THE RHO PARAMETER
      XHT=MH/MT
      R=1/XHT**2
      ALR=DLOG(R)
      H=49D0/4D0+PI**2+27D0/2D0*ALR+1.5D0*ALR**2
     .  +R/3D0*(2D0-12D0*PI**2+12D0*ALR-27D0*ALR**2)
     .  +R*R/48D0*(1613D0-240D0*PI**2-1500D0*ALR-720D0*ALR**2)
      A0=19D0-2D0*PI**2
      X4=4D0
      T=1D0/X4**2
      ALT=DLOG(T)
      H4=49D0/4D0+PI**2+27D0/2D0*ALT+1.5D0*ALT**2
     .  +T/3D0*(2D0-12D0*PI**2+12D0*ALT-27D0*ALT**2)
     .  +T*T/48D0*(1613D0-240D0*PI**2-1500D0*ALT-720D0*ALT**2)
      DH=(27D0/2D0+3D0*ALT)/T
     .   +4D0-18D0*ALT
     .  +1D0/3D0*(2D0-12D0*PI**2+12D0*ALT-27D0*ALT**2)
     .  +T/24D0*(1613D0-240D0*PI**2-1500D0*ALT-720D0*ALT**2)
     .   -T/48D0*(1500D0+2D0*720D0*ALT)
      DH4=-DH*2D0/X4**3
      A1=-(2.5D0+A0)*5D0
      A3=(H4-2D0*DH4-A0-2D0*A1)/(-32D0)
      A2=(DH4-A1-48D0*A3)/8D0
      IF(XHT.GT.2D0) THEN
      R2=H
             ELSE
      R2=A0+A1*XHT+A2*XHT**2+A3*XHT**3
         ENDIF
C  R2 ERSETZT  19-2*PI**2
      IF (LPAR(5).GE.1) DRHOTT=DRHOTT*(1D0+XGMT*R2)
      DRHOT2=DRHOTT

C---DELTAR MEANS THE QUANTITY 'DELTA R' IN MU LIFETIME FORMULA
      LPR15K=LPAR(15)
      LPAR(15)=1
      CALL SIGMR0
      LPAR(15)=LPR15K
      IF (LPAR(5).GE.2) THEN
        ALPQCD=4D0*PI/(11D0-2D0/3D0*NF)/DLOG(MT2/ALAMB2)
        DRHOTT=DRHOT2-ALPQCD/PI*(PI*PI/3D0+1D0)
     *       *ALP4PI/2D0/SW2/CW2*MT2/MZ2
      ENDIF
      IF (LPAR(5).GE.1) THEN
        DELTAR=1D0-(1D0-DALPMZ)*(1D0+CW2/SW2*DRHOTT)
     *        +dREAL(SIGFWS(0D0)-SIGFWS(MW2))/MW2
     *        +2D0*CW/SW*SIGFMS(0D0)/MZ2+FSGS(0D0)-DALPMZ
     *        -CW2/SW2*(dREAL(SIGFZS(MZ2)/MZ2-SIGFWS(MW2)/MW2)
     *                  -3D0/SW2*ALP4PI/4D0*MT2/MW2)
     *            *DSQRT(2D0)*GF*MW2*(1D0-DALPMZ)/PI/ALPHA*SW2
     *        +ALP4PI/SW2*(6D0+(3.5D0/SW2-2D0)*DLOG(CW2))
        IF (LPAR(7).GE.2) DELTAR=DELTAR-DSGMRG(MZ2)
        ELSE
        DELTAR=-dREAL(PIW0)+ALP4PI/SW2*(6D0+(3.5D0/SW2-2D0)*DLOG(CW2))
      ENDIF

c      GF=PI*ALPHA/DSQRT(2D0)/SW2/MW2/(1D0-DELTAR)
c      WRITE(*,*) ' DELTAR, GF = ',DELTAR,GF
      IF (IGFC.LT.10) GOTO 1

      DLTRMM=DELTAR

      END

************************************************************************
C..UNRENORMALIZED SELF ENERGIES, INCLUDING SINGULAR CONTRIBUTIONS

      FUNCTION SIGFZS(Q2)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGFZS,F,NYANT,ANTZH
     *          ,FQWW,FQZH,FQEE,FQMYMY,FQTATA
     *          ,FQUU,FQCC,FQTT,FQDD,FQSS,FQBB
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /DNONSM/ DNSM1,DNSM2,DNSMR
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      COMMON /RSCALE/ MSC2
        DMHZ2=MH2-MZ2
        SMHZ2=MH2+MZ2
        FQWW=F(Q2,MW,MW)
        FQZH=F(Q2,MZ,MH)
        FQEE=F(Q2,ME,ME)
        FQMYMY=F(Q2,MMY,MMY)
        FQTATA=F(Q2,MTAU,MTAU)
        FQUU=F(Q2,MU,MU)
        FQCC=F(Q2,MC,MC)
        FQTT=F(Q2,MT,MT)
        FQDD=F(Q2,MD,MD)
        FQSS=F(Q2,MS,MS)
        FQBB=F(Q2,MB,MB)
        DLHW=DLOG(MH/MW)
        DLZW=-DLOG(CW)
        DLZH=DLOG(MZ/MH)
        DLEW2=DLOG(ME2/MW2)
        DLYW2=DLOG(MMY2/MW2)
        DLAW2=DLOG(MTAU2/MW2)
        NYANT=0D0
        IF(Q2 .NE. 0D0) GOTO 10
         ANTZH=0.5D0*SMHZ2+MH2*MZ2/DMHZ2*2D0*DLZH
        GOTO 20
10      ANTZH=FQZH/Q2*DMHZ2*DMHZ2
        IF(Q2 .LT. 0D0) NYANT=5D0/3D0-DLOG(-Q2/MW2)
        IF(Q2 .GT. 0D0) NYANT=DCMPLX(5D0/3D0-DLOG(Q2/MW2),PI)
20      S2C224=1D0/(24D0*SW2*CW2)
        S2C208=3D0*S2C224
        S2_L=8D0*SW2*(SW2-0.5D0)+1D0
        S2_2=SW2/0.28125D0*(SW2-0.75D0)+1D0
        S2_1=SW2/1.125D0*(SW2-1.5D0)+1D0
C...SINGULAR PART
        AZ2=1D0/16D0/SW2/CW2
        VZN=( 0.5D0        )/2D0/SW/CW
        VZE=(-0.5D0+2D0*SW2)/2D0/SW/CW
        VZU=( 0.5D0-4D0*SW2/3D0)/2D0/SW/CW
        VZD=(-0.5D0+2D0*SW2/3D0)/2D0/SW/CW
        AVN2=AZ2+VZN*VZN
        AVE2=AZ2+VZE*VZE
        AVU2=AZ2+VZU*VZU
        AVD2=AZ2+VZD*VZD
C...SINGULAR PART
        SIGFZS=ALP4PI*Q2*(4D0/3D0
     L *(-AVN2*(DLOG(ME2/MSC2)+DLOG(MMY2/MSC2)+DLOG(MTAU2/MSC2))
     L   -AVE2*(DLOG(ME2/MSC2)+DLOG(MMY2/MSC2)+DLOG(MTAU2/MSC2))
     Q   -AVU2*(DLOG(MU2/MSC2)+DLOG(MC2/MSC2)+DLOG(MT2/MSC2))*3D0
     Q   -AVD2*(DLOG(MD2/MSC2)+DLOG(MS2/MSC2)+DLOG(MB2/MSC2))*3D0)
     G   -(3D0-19D0/6D0/SW2+1D0/6D0/CW2)*DLOG(MW2/MSC2))
     G   -ALP4PI*MZ2*(4D0+1D0/CW2-1D0/SW2)*DLOG(MW2/MSC2)
        SIGFZS=SIGFZS
     F   +ALP4PI/2D0/CW2/SW2*
     F  (ME2*DLOG(ME2/MSC2)+MMY2*DLOG(MMY2/MSC2)+MTAU2*DLOG(MTAU2/MSC2)
     F  +(MU2*DLOG(MU2/MSC2)+MC2*DLOG(MC2/MSC2)+MT2*DLOG(MT2/MSC2)
     F  +MD2*DLOG(MD2/MSC2)+MS2*DLOG(MS2/MSC2)+MB2*DLOG(MB2/MSC2))*
     F  3D0                                                      )
        SIGFZS=SIGFZS+ALP1PI*(
     G  +DFLOAT(LPAR(15))/(12D0*SW2)*
     G   (-(Q2/1.5D0+(2D1*MW2+1D1*Q2)*FQWW                     )*CW2
     G    +(3D0*MW2*FQWW+Q2/6D0-MH2*DLHW-MZ2*DLZW
     G      +0.25D0*(1D1*MZ2-2D0*MH2+Q2)*
     G         (1D0+SMHZ2/DMHZ2*DLZH-DLZW-DLHW+FQZH)
     G      +0.25D0*ANTZH                                      )/CW2
     G    +(Q2/6D0+(2D0*MW2+0.25D0*Q2)*FQWW)*(CW2-SW2)*(CW2-SW2)/CW2))
        SIGFZS=SIGFZS+ALP1PI*(
     L  +S2C224*(S2_L*((2D0*ME2+Q2)*FQEE-Q2/3D0)    -3D0*ME2*FQEE    )
     L  +S2C224*(S2_L*((2D0*MMY2+Q2)*FQMYMY-Q2/3D0) -3D0*MMY2*FQMYMY )
     L  +S2C224*(S2_L*((2D0*MTAU2+Q2)*FQTATA-Q2/3D0)-3D0*MTAU2*FQTATA)
     N  +S2C224*Q2*(NYANT+DLEW2)
     N  +S2C224*Q2*(NYANT+DLYW2)
     N  +S2C224*Q2*(NYANT+DLAW2) )
        SIGFZS=SIGFZS+ALP1PI*(
     2  +S2C208*(S2_2*((2D0*MU2+Q2)*FQUU-Q2/3D0)-3D0*MU2*FQUU    )
     2  +S2C208*(S2_2*((2D0*MC2+Q2)*FQCC-Q2/3D0)-3D0*MC2*FQCC    )
     2  +S2C208*(S2_2*((2D0*MT2+Q2)*FQTT-Q2/3D0)-3D0*MT2*FQTT    )
     1  +S2C208*(S2_1*((2D0*MD2+Q2)*FQDD-Q2/3D0)-3D0*MD2*FQDD    )
     1  +S2C208*(S2_1*((2D0*MS2+Q2)*FQSS-Q2/3D0)-3D0*MS2*FQSS    )
     1  +S2C208*(S2_1*((2D0*MB2+Q2)*FQBB-Q2/3D0)-3D0*MB2*FQBB    ) )

C...NON-STANDARD CONTRIBUTION
      SIGFZS=SIGFZS
     *   +DCMPLX(Q2,0D0)*(-2D0*DNSM1+(CW2-SW2)/SW2*(DNSM2+DNSMR))
     *   +DCMPLX(MZ2,0D0)*DNSMR

        RETURN
        END

************************************************************************

      FUNCTION SIGFWS(Q2)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGFWS,F,FQ0W,FQZW,FQHW,FQUD,FQCS,FQTB
     *          ,FQ0E,FQ0MY,FQ0TA
     *          ,ANT1,ANT2,ANT3,ANT4,ANT5,ANT6,ANT7,ANT8,ANT9
     *          ,ANT71,ANT81
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /DNONSM/ DNSM1,DNSM2,DNSMR
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      COMMON /RSCALE/ MSC2
        FQ0W=F(Q2,0D0,MW)
        FQZW=F(Q2,MZ,MW)
        FQHW=F(Q2,MH,MW)
        ANT71=1D0-DLOG(MZ2/MW2)/SW2+FQZW
        IF (MH2.EQ.MW2) THEN
          ANT81=1D0+MH2/MW2+FQHW
          ELSE
          ANT81=1D0-DLOG(MH2/MW2)*MH2/(MH2-MW2)+FQHW
        ENDIF
        IF (MD2.EQ.MU2) THEN
          DLMUD=0D0
          ELSE
          DLMUD=1D0-(MU2+MD2)/(MU2-MD2)*DLOG(MU/MD)
        ENDIF
        FQ0E=F(Q2,0D0,ME)
        FQ0MY=F(Q2,0D0,MMY)
        FQ0TA=F(Q2,0D0,MTAU)
        FQUD=F(Q2,MU,MD)
        FQCS=F(Q2,MC,MS)
        FQTB=F(Q2,MT,MB)
        IF(Q2 .NE. 0D0) GOTO 10
30      ANT9=MW2
        ANT1=ME2
        ANT2=MMY2
        ANT3=MTAU2
        ANT7=0.5D0*(MW2+MZ2)+MW2*MZ2/(MZ2-MW2)*DLOG(MW2/MZ2)
        IF (MH2.EQ.MW2) THEN
          ANT8=0.5D0*(MW2-MH2)
          ELSE
          ANT8=0.5D0*(MW2+MH2)+MW2*MH2/(MH2-MW2)*DLOG(MW2/MH2)
        ENDIF
        IF (MD2.EQ.MU2) THEN
          ANT4=0.5D0*(MU2-MD2)
          ELSE
          ANT4=0.5D0*(MU2+MD2)+MU2*MD2/(MD2-MU2)*DLOG(MU2/MD2)
        ENDIF
        ANT5=0.5D0*(MC2+MS2)+MC2*MS2/(MS2-MC2)*DLOG(MC2/MS2)
        ANT6=0.5D0*(MT2+MB2)+MT2*MB2/(MB2-MT2)*DLOG(MT2/MB2)
        GOTO 100
10      ANT9=2D0*MW2*MW2*FQ0W/Q2
        ANT1=2D0*ME2*ME2*FQ0E/Q2
        ANT2=2D0*MMY2*MMY2*FQ0MY/Q2
        ANT3=2D0*MTAU2*MTAU2*FQ0TA/Q2
        ANT7=(MZ2-MW2)*(MZ2-MW2)*FQZW/Q2
        ANT8=(MH2-MW2)*(MH2-MW2)*FQHW/Q2
        ANT4=(MD2-MU2)*(MD2-MU2)*FQUD/Q2
        ANT5=(MS2-MC2)*(MS2-MC2)*FQCS/Q2
        ANT6=(MB2-MT2)*(MB2-MT2)*FQTB/Q2
100     CONTINUE
C...SINGULAR PART
        SIGFWS=ALP4PI/6D0/SW2*
     L  (-DLOG(ME2/MSC2)*  (Q2            -0.5D0*ME2  )
     L   -DLOG(ME2/MSC2)*  (Q2-2.5D0*ME2              )
     L   -DLOG(MMY2/MSC2)* (Q2            -0.5D0*MMY2 )
     L   -DLOG(MMY2/MSC2)* (Q2-2.5D0*MMY2             )
     L   -DLOG(MTAU2/MSC2)*(Q2            -0.5D0*MTAU2)
     L   -DLOG(MTAU2/MSC2)*(Q2-2.5D0*MTAU2            )  )
        SIGFWS=SIGFWS+ALP4PI/6D0/SW2*(
     Q   -DLOG(MU2/MSC2)*(Q2-2.5D0*MU2-0.5D0*MD2)*3D0
     Q   -DLOG(MD2/MSC2)*(Q2-2.5D0*MD2-0.5D0*MU2)*3D0
     Q   -DLOG(MC2/MSC2)*(Q2-2.5D0*MC2-0.5D0*MS2)*3D0
     Q   -DLOG(MS2/MSC2)*(Q2-2.5D0*MS2-0.5D0*MC2)*3D0
     Q   -DLOG(MT2/MSC2)*(Q2-2.5D0*MT2-0.5D0*MB2)*3D0
     Q   -DLOG(MB2/MSC2)*(Q2-2.5D0*MB2-0.5D0*MT2)*3D0   )
     G   +ALP4PI/SW2*DLOG(MW2/MSC2)*(MW2*(1D0-SW2/CW2)+19D0/6D0*Q2)
        SIGFWS=SIGFWS+ALP1PI*(
     G +DFLOAT(LPAR(15))/(12D0*SW2)*
     G  (+(-(7D0*(MZ2+MW2)+1D1*Q2)*ANT71
     G     -4D0*MZ2*DLOG(MZ2/MW2)-Q2/1.5D0+2D0*ANT7         )*CW2
     G   -((4D0*MW2+1D1*Q2)*FQ0W+4D0*MW2+Q2/0.09375D0-ANT9  )*SW2
     G   +(1D0-DLOG(MZ2/MW2)/SW2+FQZW       )*3D0*MW2*SW2*SW2/CW2
     G   +( ANT81                                       )*3D0*MW2
     G   -(MZ2*DLOG(MZ2/MW2)+MH2*DLOG(MH2/MW2))*0.5D0+Q2/3D0
     G   -(MW2+MZ2-0.5D0*Q2)*ANT71*0.5D0+ANT7*0.25D0
     G   -(MH2+MW2-0.5D0*Q2)*ANT81*0.5D0+ANT8*0.25D0             ) )
        SIGFWS=SIGFWS+ALP1PI*(
     L +1D0/(12D0*SW2)*
     L  ((Q2-ME2*0.5D0)*(F(Q2,0D0,ME)+1D0)-Q2/3D0-0.25D0*ANT1)
     L +1D0/(12D0*SW2)*
     L  ((Q2-MMY2*0.5D0)*(F(Q2,0D0,MMY)+1D0)-Q2/3D0-0.25D0*ANT2)
     L +1D0/(12D0*SW2)*
     L  ((Q2-MTAU2*0.5D0)*(F(Q2,0D0,MTAU)+1D0)-Q2/3D0-0.25D0*ANT3) )
       SIGFWS=SIGFWS+ALP1PI*(
     Q +1D0/(4D0*SW2)*
     Q  (+(Q2-(MU2+MD2)*0.5D0)*(FQUD+DLMUD)-Q2/3D0-0.5D0*ANT4)
     Q +1D0/(4D0*SW2)*
     Q  (+(Q2-(MC2+MS2)*0.5D0)*(FQCS+1D0
     Q      -(MC2+MS2)/(MC2-MS2)*DLOG(MC/MS))-Q2/3D0-0.5D0*ANT5)
     Q +1D0/(4D0*SW2)*
     Q  (+(Q2-(MT2+MB2)*0.5D0)*(FQTB+1D0
     Q      -(MT2+MB2)/(MT2-MB2)*DLOG(MT/MB))-Q2/3D0-0.5D0*ANT6))

C...NON-STANDARD CONTRIBUTION
      SIGFWS=SIGFWS+
     *    DCMPLX(Q2,0D0)*(-2D0*DNSM1+(CW2-SW2)/SW2*DNSM2+CW2/SW2*DNSMR)

      RETURN
      END

************************************************************************

      FUNCTION SIGFMS(Q2)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGFMS,F,FQEE,FQMYMY,FQTATA
     *          ,FQUU,FQCC,FQTT,FQDD,FQSS,FQBB
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /DNONSM/ DNSM1,DNSM2,DNSMR
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      COMMON /RSCALE/ MSC2
        FQEE=F(Q2,ME,ME)
        FQMYMY=F(Q2,MMY,MMY)
        FQTATA=F(Q2,MTAU,MTAU)
        FQUU=F(Q2,MU,MU)
        FQCC=F(Q2,MC,MC)
        FQTT=F(Q2,MT,MT)
        FQDD=F(Q2,MD,MD)
        FQSS=F(Q2,MS,MS)
        FQBB=F(Q2,MB,MB)
C...SINGULAR PART
        VZE=(-0.5D0+2D0*SW2)/2D0/SW/CW
        VZU=( 0.5D0-4D0*SW2/3D0)/2D0/SW/CW
        VZD=(-0.5D0+2D0*SW2/3D0)/2D0/SW/CW
        SIGFMS=ALP4PI*Q2*(-4D0/3D0
     L *(         VZE*(DLOG(ME2/MSC2)+DLOG(MMY2/MSC2)+DLOG(MTAU2/MSC2))
     Q   -2D0    *VZU*(DLOG(MU2/MSC2)+DLOG(MC2/MSC2)+DLOG(MT2/MSC2))
     Q   +1D0    *VZD*(DLOG(MD2/MSC2)+DLOG(MS2/MSC2)+DLOG(MB2/MSC2)))
     G   -(3D0*CW2+1D0/6D0)/CW/SW*DLOG(MW2/MSC2))
     G   -ALP4PI*2D0*MW2/CW/SW*DLOG(MW2/MSC2)
C...FINITE PART
        SIGFMS=SIGFMS-ALP1PI*( +DFLOAT(LPAR(15))/4D0/SW2*
     G  (-(+(3D0*SW*CW+SW/(CW*6D0))*Q2
     G     +(4D0*SW*CW+SW/(CW*0.75D0))*MW2)*F(Q2,MW,MW)-Q2*SW/CW/9D0)
     L  +1D0/(3D0*SW*CW)*
     L   (0.25D0-SW2)*  ((2D0*ME2+Q2)*FQEE-Q2/3D0)
     L  +1D0/(3D0*SW*CW)*
     L   (0.25D0-SW2)*  ((2D0*MMY2+Q2)*FQMYMY-Q2/3D0)
     L  +1D0/(3D0*SW*CW)*
     L   (0.25D0-SW2)*  ((2D0*MTAU2+Q2)*FQTATA-Q2/3D0) )
        SIGFMS=SIGFMS-ALP1PI*(
     2  +1D0/(2.25D0*SW*CW)*(0.375D0-SW2)*((2D0*MU2+Q2)*FQUU-Q2/3D0)
     2  +1D0/(2.25D0*SW*CW)*(0.375D0-SW2)*((2D0*MC2+Q2)*FQCC-Q2/3D0)
     2  +1D0/(2.25D0*SW*CW)*(0.375D0-SW2)*((2D0*MT2+Q2)*FQTT-Q2/3D0)
     1  +1D0/(9D0*SW*CW)*(0.75D0-SW2)*((2D0*MD2+Q2)*FQDD-Q2/3D0)
     1  +1D0/(9D0*SW*CW)*(0.75D0-SW2)*((2D0*MS2+Q2)*FQSS-Q2/3D0)
     1  +1D0/9D0/SW/CW*(0.75D0-SW2)*((2D0*MB2+Q2)*FQBB-Q2/3D0))

C...NON-STANDARD CONTRIBUTION
      SIGFMS=SIGFMS+
     *      DCMPLX(Q2,0D0)/CW/SW*(SW2*DNSM1-CW2*(DNSM2+DNSMR))

      RETURN
      END

************************************************************************

      FUNCTION SIGFGS(Q2)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGFGS,F
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      COMMON /RSCALE/ MSC2

C...SINGULAR PART
        SIGFGR=ALP4PI*Q2*(4D0/3D0
     L        *(-(DLOG(ME2/MSC2)+DLOG(MMY2/MSC2)+DLOG(MTAU2/MSC2))
     Q          -4D0/3D0*(DLOG(MU2/MSC2)+DLOG(MC2/MSC2)+DLOG(MT2/MSC2))
     Q          -1D0/3D0*(DLOG(MD2/MSC2)+DLOG(MS2/MSC2)+DLOG(MB2/MSC2)))
     G       +3D0*DLOG(MW2/MSC2))
        SIGFGS=DCMPLX(SIGFGR,0D0)
      IF (LPAR(7).EQ.1) THEN
C...HADRONIC CONTRIBUTION FROM EFFECTIVE QUARK LOOPS
        SIGFGS=SIGFGS+ALP1PI*(
     G   -DFLOAT(LPAR(15))*((3D0*Q2+4D0*MW2)*F(Q2,MW,MW)           )/4D0
     L   +(+(Q2+2D0*ME2  )*F(Q2,ME  ,ME  )
     L     +(Q2+2D0*MMY2 )*F(Q2,MMY ,MMY )
     L     +(Q2+2D0*MTAU2)*F(Q2,MTAU,MTAU) - Q2)/3D0
     2   +(+(Q2+2D0*MU2  )*F(Q2,MU  ,MU  )
     2     +(Q2+2D0*MC2  )*F(Q2,MC  ,MC  )
     2     +(Q2+2D0*MT2  )*F(Q2,MT  ,MT  ) - Q2)/2.25D0
     1   +(+(Q2+2D0*MD2  )*F(Q2,MD  ,MD  )
     1     +(Q2+2D0*MS2  )*F(Q2,MS  ,MS  )
     1     +(Q2+2D0*MB2  )*F(Q2,MB  ,MB  ) - Q2)/9D0        )
      ELSEIF(LPAR(7).GE.2) THEN
C...HADRONIC CONTRIBUTION FROM PARAMETRIZATION
        SIGFGS=SIGFGS+ALP1PI*(
     G   -DFLOAT(LPAR(15))*((3D0*Q2+4D0*MW2)*F(Q2,MW,MW)           )/4D0
     L   +(+(Q2+2D0*ME2  )*F(Q2,ME  ,ME  )
     L     +(Q2+2D0*MMY2 )*F(Q2,MMY ,MMY )
     L     +(Q2+2D0*MTAU2)*F(Q2,MTAU,MTAU) - Q2)/3D0
     2   +(+(Q2+2D0*MT2  )*F(Q2,MT  ,MT  ) - Q2/3D0)/2.25D0     )
     H   - DCMPLX(HADRQQ(Q2)*Q2,0D0)
      ELSEIF(LPAR(7).EQ.0) THEN
       SIGFGS=DCMPLX(0D0,0D0)
      ELSE
      ENDIF
      RETURN
      END

************************************************************************

      FUNCTION FSGS(Q2)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGFGS
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      COMMON /RSCALE/ MSC2

      IF (Q2.EQ.0D0) THEN
        FSGS=ALP4PI*(4D0/3D0
     L        *(-(DLOG(ME2/MSC2)+DLOG(MMY2/MSC2)+DLOG(MTAU2/MSC2))
     Q          -4D0/3D0*(DLOG(MU2/MSC2)+DLOG(MC2/MSC2)+DLOG(MT2/MSC2))
     Q          -1D0/3D0*(DLOG(MD2/MSC2)+DLOG(MS2/MSC2)+DLOG(MB2/MSC2)))
     G       +3D0*DLOG(MW2/MSC2))
        FSGS=FSGS-ALP1PI*DFLOAT(LPAR(15))/6D0
      ELSE
        FSGS=dREAL(SIGFGS(Q2))/Q2
      ENDIF

      RETURN
      END

C***********************************************************************

      FUNCTION FSZS(Q2)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGFZS
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      FSZS=dREAL(SIGFZS(Q2)-SIGFZS(0D0))/Q2+DSGMRG(MZ2)
      RETURN
      END

      FUNCTION FSWS(Q2)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGFWS
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      FSWS=dREAL(SIGFWS(Q2)-SIGFWS(0D0))/Q2+DSGMRG(MZ2)
      RETURN
      END

      FUNCTION FSMS(Q2)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGFMS
      FSMS=dREAL(SIGFMS(Q2)-SIGFMS(0D0))/Q2
      RETURN
      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FUNCTION FSZD(Q2)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      Q21=Q2*(1D0+1D-4)
      Q22=Q2*(1D0-1D-4)
      FSZD=(FSZS(Q22)-FSZS(Q21))/(Q22-Q21)
      RETURN
      END

************************************************************************

      SUBROUTINE SIGMR0
C
C   RENORMALIZED NC SELF ENERGIES, DETERMINATION OF COUNTERTERMS
C
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGFZS,SIGFWS,SIGFMS,SIGMQM,SIGMQG
     *          ,PIZMZ,PIMMZ,PIW0
     *          ,CAFINZ,CAFINW,CAFING,SMM1,SMM2
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /CTMASS/ RDMZ2,RDMW2,PIZMZ,PIMMZ,PIW0
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
C
      CAFINW=SIGFWS(MW2)
      CAFINZ=SIGFZS(MZ2)
      RDMW2=dREAL(CAFINW)/MW2
      RDMZ2=dREAL(CAFINZ)/MZ2
      RDMZ20=RDMZ2
      CAFING=SIGMQG(MZ2)
      I=0
    1 I=I+1
      PIMMZ=(SIGFMS(MZ2)+SIGFMS(0D0))/MZ2-CW/SW*(RDMZ20-RDMW2)
      IF (LPAR(5).LT.1) GOTO 2
      RDMZ21=RDMZ2-dREAL(PIMMZ*PIMMZ/(1D0+CAFING/MZ2))
      DRDMZ2=DABS(RDMZ21-RDMZ20)
      IF (DRDMZ2.GT.1D-6.AND.I.LT.30) THEN
        RDMZ20=RDMZ21
        GOTO 1
      ENDIF
      RDMZ2=RDMZ21
    2 CONTINUE

      S1=MZ2*(1D0+1D-4)
      S2=MZ2*(1D0-1D-4)
      SMM1=SIGMQM(S1)
      SMM2=SIGMQM(S2)
      PIZMZ=(SIGFZS(S1)-SIGFZS(S2)
     *       -SMM1*SMM1/(S1+SIGMQG(S1))
     *       +SMM2*SMM2/(S2+SIGMQG(S2)))/(S1-S2)-FSGS(0D0)
     *     +(CW2/SW2-1D0)*(RDMZ2-RDMW2-2D0*SW/CW*SIGFMS(0D0)/MZ2)
      IF (LPAR(7).GE.2) THEN
        DDALPP=DSGMRG(MZ2)
        PIZMZ=PIZMZ+DCMPLX(DDALPP,0D0)
      ENDIF

      PIW0=SIGFWS(MW2)/MW2-RDMW2+CW2/SW2*(RDMZ2-RDMW2)
     *    -2D0*CW/SW*SIGFMS(0D0)/MZ2-FSGS(0D0)
      IF (LPAR(7).GE.2) THEN
        DDALPP=DSGMRG(MZ2)
        PIW0=PIW0+DCMPLX(DDALPP,0D0)
      ENDIF

      RETURN
      END

************************************************************************

      FUNCTION SIGMQG(Q2)
C
C   RENORMALIZED PHOTON SELF ENERGY
C
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGMQG,SIGFGS

      SIGMQG=SIGFGS(Q2)-Q2*FSGS(0D0)

      RETURN
      END

************************************************************************

      FUNCTION SIGMQM(Q2)
C
C   RENORMALIZED PHOTON Z MIXING
C
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGMQM,PIZMZ,PIMMZ,PIW0,SIGFMS
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /CTMASS/ RDMZ2,RDMW2,PIZMZ,PIMMZ,PIW0

      SIGMQM=SIGFMS(Q2)-SIGFMS(0D0)
     *      +Q2*(2D0*SIGFMS(0D0)/MZ2-CW/SW*(RDMZ2-RDMW2))

      RETURN
      END

************************************************************************

      FUNCTION SIGMQZ(Q2)
C
C   RENORMALIZED Z SELF ENERGY
C
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGMQZ,PIZMZ,PIMMZ,PIW0,SIGFZS,SIGFMS
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      COMMON /CTMASS/ RDMZ2,RDMW2,PIZMZ,PIMMZ,PIW0

      SIGMQZ=SIGFZS(Q2)-RDMZ2*MZ2
     *      +(Q2-MZ2)*(-FSGS(0D0)-2D0*(CW2-SW2)/SW/CW*SIGFMS(0D0)/MZ2
     *                 +(CW2-SW2)/SW2*(RDMZ2-RDMW2))
      IF (LPAR(7).GE.2) THEN
        DDALPP=DSGMRG(MZ2)
        SIGMQZ=SIGMQZ+DCMPLX(DDALPP,0D0)*(Q2-MZ2)
      ENDIF

      RETURN
      END

************************************************************************

      FUNCTION SIGMQW(Q2)
C
C   RENORMALIZED W SELF ENERGY
C
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGMQW,PIZMZ,PIMMZ,PIW0,SIGFWS,SIGFMS
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      COMMON /CTMASS/ RDMZ2,RDMW2,PIZMZ,PIMMZ,PIW0

      SIGMQW=SIGFWS(Q2)-RDMW2*MW2
     *      +(Q2-MW2)*(-FSGS(0D0)-2D0*CW/SW*SIGFMS(0D0)/MZ2
     *                 +CW2/SW2*(RDMZ2-RDMW2))
      IF (LPAR(7).GE.2) THEN
        DDALPP=DSGMRG(MZ2)
        SIGMQW=SIGMQW+DCMPLX(DDALPP,0D0)*(Q2-MW2)
      ENDIF

      RETURN
      END

************************************************************************

      FUNCTION PIGQQ(Q2)
C
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 PIGQQ,SIGMQG

      IF (Q2.EQ.0D0) THEN
        PIGQQ=DCMPLX(0D0,0D0)
      ELSE
        PIGQQ=SIGMQG(Q2)/Q2
      ENDIF
      RETURN
      END

************************************************************************

      FUNCTION PIMQQ(Q2)
C
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 PIMQQ,PIGQQ,PIZQQ,SIGMQM,SIGFMS
      COMPLEX*16 PIZMZ,PIMMZ,PIW0
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      COMMON /CTMASS/ RDMZ2,RDMW2,PIZMZ,PIMMZ,PIW0
      COMMON /RSCALE/ MSC2

      IF (Q2.EQ.0D0) THEN
        VZE=(-0.5D0+2D0*SW2)/2D0/SW/CW
        VZU=( 0.5D0-4D0*SW2/3D0)/2D0/SW/CW
        VZD=(-0.5D0+2D0*SW2/3D0)/2D0/SW/CW
        PIMQQ=-ALP4PI*4D0/3D0
     L *(     VZE*(DLOG(ME2/MSC2)+DLOG(MMY2/MSC2)+DLOG(MTAU2/MSC2))
     Q   -2D0*VZU*(DLOG(MU2/MSC2)+DLOG(MC2/MSC2)+DLOG(MT2/MSC2))
     Q   +1D0*VZD*(DLOG(MD2/MSC2)+DLOG(MS2/MSC2)+DLOG(MB2/MSC2)))
     G +DFLOAT(LPAR(15))*ALP4PI
     G   *(-(3D0*CW2+1D0/6D0)/CW/SW*DLOG(MW2/MSC2)
     G     +1D0/SW2*((4D0*SW*CW+SW/(CW*0.75D0))/6D0+SW/CW/9D0) )
        PIMQQ=PIMQQ+2D0*SIGFMS(0D0)/MZ2-CW/SW*(RDMZ2-RDMW2)
      ELSE
        PIMQQ=SIGMQM(Q2)/Q2/(1D0+PIGQQ(Q2))/(1D0+PIZQQ(Q2))
      ENDIF

      RETURN
      END

************************************************************************

      FUNCTION PIZQQ(Q2)
C
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 PIZQQ,SIGMQZ,SIGMQM,SIGMQG
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2

      SM=SIGMQM(Q2)
      PIZQQ=(SIGMQZ(Q2)-SM*SM/(Q2+SIGMQG(Q2)))/(Q2-MZ2)

      RETURN
      END

************************************************************************

      FUNCTION PIZQQR(Q2)
C
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 PIZQQR,SIGFZS,SIGFMS,SIGFWS
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HDELTR/ DELTAR,AGF0,DRHOT,DALPMZ,XGMT,ALPQCD,BTOP4,DRPIW2
      COMMON /DNONSM/ DNSM1,DNSM2,DNSMR

      PIZQQR=(SIGFZS(Q2)-SIGFZS(MZ2))/(Q2-MZ2)-FSGS(0D0)
     *   -2D0*(CW2-SW2)/CW/SW*SIGFMS(0D0)/MZ2
     *   +(CW2-SW2)/SW2*(dREAL(SIGFZS(MZ2))/MZ2-dREAL(SIGFWS(MW2))/MW2)
     *   +DALPMZ-(CW2-SW2)/SW2*DRHOT
      IF (LPAR(7).GE.2) THEN
        DDALPP=DSGMRG(MZ2)
        PIZQQR=PIZQQR+DCMPLX(DDALPP,0D0)
      ENDIF
      DPIZNS=-2D0*DNSM1+(CW2-SW2)/SW2*(DNSM2+DNSMR)
      PIZQQR=PIZQQR-DCMPLX(DPIZNS,0D0)-DCMPLX(DNSMR,0D0)

      RETURN
      END

************************************************************************

      FUNCTION PIWQQ(Q2)
C
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 PIWQQ,SIGMQW
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2

      PIWQQ=SIGMQW(Q2)/(Q2-MW2)

      RETURN
      END

************************************************************************

      FUNCTION PIWQQR(Q2)
C
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 PIWQQR,SIGFWS,SIGFZS,SIGFMS
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HDELTR/ DELTAR,AGF0,DRHOT,DALPMZ,XGMT,ALPQCD,BTOP4,DRPIW2
      COMMON /DNONSM/ DNSM1,DNSM2,DNSMR

      PIWQQR=(SIGFWS(Q2)-SIGFWS(MW2))/(Q2-MW2)
     *      -FSGS(0D0)-2D0*CW/SW*SIGFMS(0D0)/MZ2
     *      +CW2/SW2*(dREAL(SIGFZS(MZ2))/MZ2-dREAL(SIGFWS(MW2))/MW2)
     *      +DALPMZ-CW2/SW2*DRHOT+DRPIW2
      IF (LPAR(7).GE.2) THEN
        DDALPP=DSGMRG(MZ2)
        PIWQQR=PIWQQR+DCMPLX(DDALPP,0D0)
      ENDIF
      DPIWNS=-2D0*DNSM1+(CW2-SW2)/SW2*DNSM2+CW2/SW2*DNSMR
      PIWQQR=PIWQQR-DCMPLX(DPIWNS,0D0)

      RETURN
      END

************************************************************************

      SUBROUTINE SETFFQ(T)
C---FORM FACTORS
C   RECALCULATE EFFECTIVE COUPLING CONSTANTS INCLUDING SELF ENERGIES
C   RUNNING ALPHA, EFFECTIVE SW2, FORM FACTOR KAPPA (FROM PI_Z)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 PIGQQ,PIMQQ,PIZQQR,CFHFB
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      COMMON /PARAM/  POLARI,LLEPT,LQUA
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HDELTR/ DELTAR,AGF0,DRHOT,DALPMZ,XGMT,ALPQCD,BTOP4,DRPIW2
      COMMON /GSW1/   MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /SMCON/  VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /SMCON1/ VAFI1(2,3,2)
     &               ,AFIJ1(3,2,2),BFIJ1(3,2,2),FLIND1(2,3,2,2)
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /KNSTCC/ SXNRCC,SX1NCC
      COMMON /FORMFF/ ALPFFQ,AKAPPA,GMUFFQ,SWEFF2
      COMMON /DLTASM/ DALMZ,DRHOTT,DRHO0,DLTA1,DLTA2,DLTAR
C
C  SELF ENERGIES
C
      IF (LPAR(7).GE.1) THEN
        PIGGG=dREAL(PIGQQ(T))
        ELSE
        PIGGG=0D0
      ENDIF
      ALPFFQ=1D0/(1D0+PIGGG)
      EFFQ=SQRT(ALPFFQ)
      IF (LPAR(8).EQ.1) THEN
        SM=dREAL(PIMQQ(T))
        ELSE
        SM=0D0
      ENDIF
      AKAPPA=1D0-CW/SW*SM
      SWEFF2=SW2*AKAPPA
      IF (LPAR(9).EQ.1) THEN
        PIGZRM=dREAL(PIZQQR(T))
        ELSE
        PIGZRM=0D0
      ENDIF
      GMUFFQ=1D0-PIGZRM
      BFFQ=SQRT(GMUFFQ)
      IF (LPAR(5).GE.1) BFFQ=BFFQ*BTOP4
C
C---DEFINE FERMION GAUGE BOSON COUPLING CONSTANTS
      IF (LPAR(4).EQ.1) THEN
        B0=1D0/4D0/CW/SW
        B=B0
        ELSE
C---NORMALIZED TO G-MU
        B0=MZ/SQRT(AGF0)/4D0
        B=B0
        IF (LPAR(9).GE.1) B=B0*SQRT(1D0-DELTAR)
      ENDIF
      IF (LPAR(2).EQ.1.AND.LPAR(9).EQ.1) B=B*BFFQ
      RHO1=B/B0

      VAFI1(2,1,1)=0D0
      VAFI1(2,2,1)=0D0
      VAFI1(2,3,1)=0D0
      VAFI1(1,1,1)=EFFQ
      VAFI1(1,2,1)=-2D0/3D0*EFFQ
      VAFI1(1,3,1)=1D0/3D0*EFFQ
      VAFI1(2,1,2)=-B
      VAFI1(2,2,2)=B
      VAFI1(2,3,2)=-B
      VAFI1(1,1,2)=B*(4D0*SWEFF2-1D0)
      VAFI1(1,2,2)=B*(1D0-8D0*SWEFF2/3D0)
      VAFI1(1,3,2)=B*(4D0*SWEFF2/3D0-1D0)
C
C..VERTEX CORRECTIONS
C..ELECTRON VERTEX
      IF ((LPAR(12).EQ.1).AND.(LPAR(11).EQ.0)) THEN
        LF=1
        DO 1 IVA=1,2
        VAFI1(IVA,LF,1)=VAFI1(IVA,LF,1)
     *                   +EFFQ*dREAL(CFHFB(T,IVA,LF,1,ME2))
        VAFI1(IVA,LF,2)=VAFI1(IVA,LF,2)
     *                   +RHO1*dREAL(CFHFB(T,IVA,LF,2,ME2))
    1   CONTINUE
      ENDIF
      IF ((LPAR(12).EQ.1).AND.(LPAR(11).EQ.1)) THEN
        LF=1
        DO 2 IVA=1,2
        VAFI1(IVA,LF,1)=VAFI1(IVA,LF,1)+dREAL(CFHFB(T,IVA,LF,1,ME2))
        VAFI1(IVA,LF,2)=VAFI1(IVA,LF,2)+dREAL(CFHFB(T,IVA,LF,2,ME2))
    2   CONTINUE
      ENDIF
C..QUARK VERTEX
      IF ((LPAR(13).EQ.1).AND.(LPAR(11).EQ.0)) THEN
        DO 3 LF=2,3
        DO 3 IVA=1,2
        VAFI1(IVA,LF,1)=VAFI1(IVA,LF,1)
     *                   +EFFQ*dREAL(CFHFB(T,IVA,LF,1,MQI2))
        VAFI1(IVA,LF,2)=VAFI1(IVA,LF,2)
     *                   +RHO1*dREAL(CFHFB(T,IVA,LF,2,MQI2))
    3   CONTINUE
      ENDIF
      IF ((LPAR(13).EQ.1).AND.(LPAR(11).EQ.1)) THEN
        DO 4 LF=2,3
        DO 4 IVA=1,2
        VAFI1(IVA,LF,1)=VAFI1(IVA,LF,1)+dREAL(CFHFB(T,IVA,LF,1,MQI2))
        VAFI1(IVA,LF,2)=VAFI1(IVA,LF,2)+dREAL(CFHFB(T,IVA,LF,2,MQI2))
    4   CONTINUE
      ENDIF
C
      LBOSON=LPAR(17)
      IF (LBOSON.EQ.1) THEN
        DO 11 LF=1,3
        DO 11 IVA=1,2
   11   VAFI1(IVA,LF,2)=0
      ENDIF
      IF (LBOSON.EQ.3) THEN
        DO 12 LF=1,3
        DO 12 IVA=1,2
   12   VAFI1(IVA,LF,1)=0
      ENDIF
C
      IGAMMA = 1
      IZ = 2
      IEL = 1
      IFU = 2
      IFD = 3
      INDV = 1
      INDA = 2
C
      DO 5 IF=IEL,IFD
        DO 5 IB1=IGAMMA,IZ
          DO 5 IB2=IGAMMA,IZ
          IF (LBOSON.EQ.2.AND.IB1.EQ.IB2) THEN
             FLIND1(INDV,IF,IB1,IB2)=0D0
             GOTO 5
          ENDIF
          FLIND1(INDV,IF,IB1,IB2)=
     *     2D0*(VAFI1(INDV,IF,IB1)*VAFI1(INDV,IF,IB2)
     *         +VAFI1(INDA,IF,IB1)*VAFI1(INDA,IF,IB2))
    5 CONTINUE
      DO 6 IF=IEL,IFD
        DO 6 IB1=IGAMMA,IZ
          DO 6 IB2=IGAMMA,IZ
          IF (LBOSON.EQ.2.AND.IB1.EQ.IB2) THEN
             FLIND1(INDA,IF,IB1,IB2)=0D0
             GOTO 6
          ENDIF
          FLIND1(INDA,IF,IB1,IB2)=
     *     2D0*(VAFI1(INDV,IF,IB1)*VAFI1(INDA,IF,IB2)
     *         +VAFI1(INDA,IF,IB1)*VAFI1(INDV,IF,IB2))
    6 CONTINUE
C
      DO 7 IVB1 = IGAMMA, IZ
       DO 7 IVB2 = IGAMMA, IZ
        DO 7 IFERM = IFU, IFD
          IF (LBOSON.EQ.2.AND.IVB1.EQ.IVB2) THEN
            AFIJ1(IFERM,IVB1,IVB2)=0D0
            BFIJ1(IFERM,IVB1,IVB2)=0D0
            GOTO 7
          ENDIF
        AFIJ1(IFERM,IVB1,IVB2)=FLIND1(INDV,IFERM,IVB1,IVB2)
     *  *(FLIND1(INDV,IEL,IVB1,IVB2)-POLARI*FLIND1(INDA,IEL,IVB1,IVB2))
        BFIJ1(IFERM,IVB1,IVB2)=FLIND1(INDA,IFERM,IVB1,IVB2)
     *  *(FLIND1(INDA,IEL,IVB1,IVB2)-POLARI*FLIND1(INDV,IEL,IVB1,IVB2))
    7 CONTINUE
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      FUNCTION HADRQQ(S)
C  HADRONIC IRREDUCIBLE QQ SELF-ENERGY: TRANSVERSE
* THIS IS A SLIGHTLY MODIFIED VERSION OF BURKHARDTS CODE. DB. TR.,01/91
C     parametrize the real part of the photon self energy function
C     by  a + b ln(1+C*:S:) , as in my 1981 TASSO note but using
C     updated values, extended using RQCD up to 100 TeV
C     for details see:
C     H.Burkhardt, F.Jegerlehner, G.Penso and C.Verzegnassi
C     in CERN Yellow Report on "Polarization at LEP" 1988
C     H.BURKHARDT, CERN/ALEPH, AUGUST 1988
C     negative values mean t - channel (spacelike)
C     positive values mean s - channel (timelike )
C     in the space like values around 1 GeV are typical for luminosity
C     the values at 92 GeV ( Z mass ) give the light quark contribution
C     to delta r
C     take care of the sign of REPI when using this in different
C     programs
C     Here REPI was chosen to
C     be positive (so that it corresponds directly to delta alpha)
C     often its assumed to be negative.
C
C     the imaginary part is proportional to R (had / mu cross section)
C     and is therefore 0 below threshold ( in all the spacelike region)
C     note also that alpha_s usually has been derived from the measured
C     values of R.
C     Changing alpha_s to values incompatible with current data
C     would imply to be also inconsistent with RE,IM PI
C     defined here
C
C     H.BURKHARDT
C
      IMPLICIT REAL*8(A-H,M,O-Z)
C
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      real q,st2,der,errder,deg,errdeg
C
      DATA A1,B1,C1/   0.0   ,  0.00835,  1.0  /
      DATA A2,B2,C2/   0.0   ,  0.00238,  3.927 /
      DATA A3,B3,C3/ 0.00165 ,  0.00300,  1.0  /
      DATA A4,B4,C4/ 0.00221 ,  0.00293,  1.0  /
      DATA PI/3.141592653589793D0/,ALFAIN/137.0359895D0/,INIT/0/
C
      hadrqq=0d0
      IF (LPAR(7).EQ.2) THEN
        IF(INIT.EQ.0) THEN
          INIT=1
          ALFA=1./ALFAIN
          ALFAPI=1./PI/ALFAIN
        ENDIF
        T=ABS(S)
        IF(T.LT.0.3**2) THEN
          REPIAA=A1+B1*LOG(1.+C1*T)
        ELSEIF(T.LT.3.**2) THEN
          REPIAA=A2+B2*LOG(1.+C2*T)
        ELSEIF(T.LT.100.**2) THEN
          REPIAA=A3+B3*LOG(1.+C3*T)
        ELSE
          REPIAA=A4+B4*LOG(1.+C4*T)
        ENDIF
C     as imaginary part take -i alfa/3 Rexp
        HADRQQ=REPIAA

C...HADRONIC CONTRIBUTION FROM JEGERLEHNER'S hadr5n12 (2012)
      ELSEIF(LPAR(7).EQ.3) THEN
        st2=sw2
c...Q2 is space-like for DIS
        q=-sqrt(abs(s))
        call hadr5n12(q,st2,der,errder,deg,errdeg)
        hadrqq=der
      ENDIF

      END

      FUNCTION DSGMRG(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       DIFFERENCE OF THE RENORMALIZED PHOTON SELF ENERGY DIVIDED BY Q2
C       OF THE PERTURBATIVE RESULT WITH EFFECTIVE QUARK MASSES AND
C       THE PARAMETRIZATION OF BURKHARD
C       14.05.91 HS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SGEFFQ,F
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART

      IF (LPAR(7).LT.2) THEN
C...NO CORRECTION
        DSGMRG=0D0
      ELSE
        SGEFFQ=+ALP1PI*(
     2   +(+(Q2+2D0*MU2  )*F(Q2,MU  ,MU  )
     2     +(Q2+2D0*MC2  )*F(Q2,MC  ,MC  ) - 2D0*Q2/3D0)/2.25D0
     1   +(+(Q2+2D0*MD2  )*F(Q2,MD  ,MD  )
     1     +(Q2+2D0*MS2  )*F(Q2,MS  ,MS  )
     1     +(Q2+2D0*MB2  )*F(Q2,MB  ,MB  ) - Q2        )/9D0        )
        DSGMRG=-HADRQQ(Q2)-DREAL(SGEFFQ)/Q2
      ENDIF
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        FUNCTION F(Q2,RM,RN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SKALARES EINSCHLEIFENINTEGRAL, DOPPELTLANG, 'REGULAERER ANTEIL'
C       F(Q2,RM,RN)=B0(Q2,RM,RN)-B0(0D0,RM,RN)  'SUBTRAHIERTES F'
C       Q2=QUADRAT DES DIE SCHLEIFE DURCHLAUFENDEN IMPULSES
C       RM,RN: MASSEN DER TEILCHEN AUF BEIDEN ARMEN
C       D O U B L E P R E C I S I O N
C-----------------------------------------------------------------------
C       19.10.83
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMPLEX*16 F
        DOUBLEPRECISION M,N,PI,A,S,T,B,C,D,Q2,RM,RN,U,V,W
        DATA PI/3.1415926535897932384626438D0/
        M=RM
        N=RN
505     IF (M .EQ. N ) GOTO 30
506     IF (N .EQ. 0.D0) GOTO 310
507     IF (M .EQ. 0.D0) GOTO 300
C---->  ALLGEMEINER FALL
510     IF (Q2 .NE. 0D0) GOTO 520
        B=0D0
        A=0D0
        GOTO 560
520     U=M*M+N*N
        V=M*M-N*N
        W=M*N
        S=M+N
        T=M-N
        C=DSQRT(DABS(S*S-Q2))
        D=DSQRT(DABS(T*T-Q2))
        B=1D0+(V/Q2-U/V)*DLOG(N/M)
540     IF (2D0*W .LE. DABS(Q2-U)) GOTO 550
        B=B-2D0*C*D/Q2*DATAN(D/C)
        A=0D0
        GOTO 560
550     A=C*D/Q2
        B=B-DSIGN(1D0,Q2-U)*A*DLOG((C+D)*(C+D)/(4D0*W))
        A=PI*A
        IF (Q2 .GE. U) GOTO 560
        A=0D0
560     CONTINUE
570     F=DCMPLX(B,A)
        RETURN
C---->  GLEICHE MASSEN
30      IF (Q2 .NE. 0D0) GOTO 40
        B=0D0
        A=0D0
        GOTO 560
40      U=4D0*M*M
        IF (DABS(Q2).LT.U/1D4) THEN
          B=Q2/6D0/M/M
          A=0D0
          GOTO 560
        ENDIF
        V=DSQRT(DABS(1D0-U/Q2))
        IF ((Q2 .GE. 0D0) .AND. (Q2 .LT. U)) GOTO 50
        B=2D0-V*DLOG((V+1D0)*(V+1D0)/U*DABS(Q2))
        A=PI*V
        IF (Q2 .GE. U) GOTO 560
        A=0D0
        GOTO 560
50      B=2D0-2D0*V*DATAN(1D0/V)
        A=0D0
        GOTO 560
C---->  EINE MASSE NULL
300     M=N
310     IF (Q2 .NE. M*M) GOTO 320
        A=0D0
        B=1D0
        GOTO 560
320     B=1D0
        IF(Q2 .EQ. 0D0) B=0D0
        A=B*(1D0-M*M/(Q2+(1D0-B)))
        B=B-A*DLOG(DABS(1D0-Q2/M/M))
        A=PI*A
        IF (Q2 .GT. M*M) GOTO 560
        A=0D0
        GOTO 560
        END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     VERTEX CORRECTIONS
C     IR-FINITE PARTS, DOUBLE-LOG CONSTRIBUTIONS SUBTRACTED
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION CFHFB(T,IVA,LF,LB,MF2)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 CFHFB,CLAMB1,CLAMB2,CLAMB3,CMW2,CMZ2
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /SMCON/  VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /CBMASS/ CMW2,CMZ2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
      IF ((LPAR(11).EQ.1).AND.
     I    (LPAR(12).EQ.1.OR.LPAR(13).EQ.1)) THEN
C---PHOTON-EXCHANGE
        CFHFB=VAFI(IVA,LF,LB)*VAFI(1,LF,1)**2 * ALP4PI*CLAMB1(T,MF2)
        ELSE
        CFHFB=DCMPLX(0D0,0D0)
      ENDIF
      IF (LPAR(15).EQ.0) RETURN
C---Z-EXCHANGE
      IF (IVA.EQ.1) THEN
        IVAS=2
        ELSE
        IVAS=1
      ENDIF
      FCZ=VAFI(IVA,LF,LB)*(VAFI(1,LF,2)**2+VAFI(2,LF,2)**2)
     *   +VAFI(IVAS,LF,LB)*2D0*VAFI(1,LF,2)*VAFI(2,LF,2)
      CFHFB=CFHFB + ALP4PI*FCZ*CLAMB2(T,CMZ2)
C---W-EXCHANGE
      IF (LF.EQ.1.AND.LB.EQ.1) GOTO 10
      LFS=3
      IF (LF.EQ.3) LFS=2
      FCW=(VAFI(1,LFS,LB)+VAFI(2,LFS,LB))/4D0/SW2
      IF (LF.EQ.1.AND.LB.EQ.2) FCW=1D0/8D0/CW/SW/SW2
      CFHFB=CFHFB + ALP4PI*FCW*CLAMB2(T,CMW2)
C---TRIPLE GAUGE BOSON VERTEX
   10 CONTINUE
      IF (LB.EQ.1) THEN
        FC3B=3D0/4D0/SW2
        ELSE
        FC3B=-3D0/4D0/SW2*CW/SW
      ENDIF
      IF (MOD(LF,2).EQ.0) FC3B=-FC3B
      CFHFB=CFHFB + ALP4PI*FC3B*CLAMB3(T,CMW2)
      RETURN
      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
        FUNCTION CLN(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       KOMPLEXER NATUERLICHER LOGARITHMUS MIT SCHNITT= NEG. REELLE ACHS
C-----------------------------------------------------------------------
C       20.07.83 UE 20.08.85
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMPLEX*16 CLN,Z
        REAL*8 X,Y,R,PHI,SIGN
        EXTERNAL SIGN
C---->  SIGN(A,B)=|A|*SGN(B), NICHT DIE INTRINSIC FUNCTION VERWENDEN!
C       PI=3.1415926535897932384
        X=DREAL(Z)
        Y=DIMAG(Z)
        R=ABS(Z)


cv        IF (X)  10,20,30
cv20      PHI=SIGN(1.5707963267948966192D0,Y)
cv        GOTO 100
cv10      IF (Y .EQ. 0D0) GOTO 40
cv        PHI=DATAN(Y/X)+SIGN(3.1415926535897932384D0,Y)
cv        GOTO 100
cv40      PHI=3.1415926535897932384D0
cv        GOTO 100
cv30      PHI=DATAN(Y/X)

        IF (X.eq.0) then
           PHI=SIGN(1.5707963267948966192D0,Y)
           GOTO 100
        elseif (x.lt.0) then
           IF (Y .EQ. 0D0) GOTO 40
           PHI=DATAN(Y/X)+SIGN(3.1415926535897932384D0,Y)
           GOTO 100
 40        PHI=3.1415926535897932384D0
           GOTO 100
        elseif (x.gt.0) then
           PHI=DATAN(Y/X)
        endif


100     CLN=DCMPLX(DLOG(R),PHI)
        RETURN
        END
        FUNCTION SIGN(X,Y)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SIGN(X,Y)=|X|*SGN(Y) IM GEGENSATZ ZUR INTRINSIC FUNCTION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        REAL*8 SIGN,X,Y
cv        IF (Y)  10,20,30
cv10      SIGN=-DABS(X)
cv        RETURN
cv20      SIGN=0D0
cv        RETURN
cv30      SIGN=DABS(X)
        if (Y.lt.0D0) then
           SIGN=-DABS(X)
        elseif (Y.eq.0D0) then
           SIGN=0D0
        elseif (Y.ge.0D0) then
           SIGN=DABS(X)
        endif

        RETURN
        END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        FUNCTION CSPEN(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK
C-----------------------------------------------------------------------
C       20.07.83    UE 20.08.85
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMPLEX*16 CSPEN,CLN,W,SUM,Z,U
        REAL*8 RZ,AZ,A1
        REAL*8 B(9)/
     1   0.1666666666666666666666666667D0,
     2  -0.0333333333333333333333333333D0,
     3   0.0238095238095238095238095238D0,
     4  -0.0333333333333333333333333333D0,
     5   0.0757575757575757575757575758D0,
     6  -0.2531135531135531135531135531D0,
     7   1.1666666666666666666666666667D0,
     8  -7.09215686274509804D0         ,
     9  54.97117794486215539D0         /
C     BEACHTE:                 B(N)=B2N
C     B(1)=1./6.
C     B(2)=-1./30.
C     B(3)=1./42.
C     B(4)=-1./30.
C     B(5)=5./66.
C     B(6)=-691./2730.
C     B(7)=7./6.
C     B(8)=-3617./510.
C     B(9)=43867./798.
C     B(10)=-174611./330.
C     B(11)=854513./138.
C     PI=3.1415926535897932384
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...
C
      RZ=DREAL(Z)
      AZ=ABS(Z)
      A1=ABS(1D0-Z)
CH
CH FOR VS FORTRAN:
      IF (AZ.LT.1D-15) THEN
         CSPEN=DCMPLX(0D0,0D0)
         RETURN
      ENDIF
CH
      IF((RZ .EQ. 1D0) .AND. (DIMAG(Z) .EQ. 0D0)) GOTO 40
      IF(RZ.GT.5D-1) GOTO 20
      IF(AZ.GT.1D0) GOTO 10
      W=-CLN(1D0-Z)
      SUM=W-0.25D0*W*W
      U=W
      DO 1 K=1,9
      U=U*W*W/FLOAT(2*K*(2*K+1))
      SUM=SUM+U*B(K)
 1    CONTINUE
      CSPEN=SUM
      RETURN
10    W=-CLN(1D0-1D0/Z)
      SUM=W-0.25D0*W*W
      U=W
      DO 11 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      SUM=SUM+U*B(K)
11    CONTINUE
      CSPEN=-SUM-1.64493406684822643D0-.5D0*CLN(-Z)**2
      RETURN
20    IF(A1.GT.1D0) GOTO 30
      W=-CLN(Z)
      SUM=W-0.25D0*W*W
      U=W
      DO 21 K=1,9
      U=U*W*W/FLOAT(2*K*(2*K+1))
      SUM=SUM+U*B(K)
21    CONTINUE
      CSPEN=-SUM+1.64493406684822643D0-CLN(Z)*CLN(1D0-Z)
      RETURN
30    W=CLN(1D0-1D0/Z)
      SUM=W-0.25D0*W*W
      U=W
      DO 31 K=1,9
      U=U*W*W/FLOAT(2*K*(2*K+1))
      SUM=SUM+U*B(K)
31    CONTINUE
      CSPEN=SUM+3.28986813369645287D0+.5D0*CLN(Z-1D0)**2-
     *         CLN(Z)*CLN(1D0-Z)
50    CONTINUE
      RETURN
40    CSPEN=1.64493406684822643D0
      RETURN
      END
C=======================================================================
C=======================================================================
      FUNCTION CLAMB1(Q2,MF2)
C       PHOTONIC VERTEX CORRECTION
C       IR-FINITE PART
C                                 OHNE LOG**2
C
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 CLAMB1
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
cv      IF (Q2) 10,20,30
cv10    AMB1 = DLOG(-Q2/MF2)
cv     1         +4D0*(PI*PI/12D0 - 1D0)
cv      CLAMB1=DCMPLX(AMB1)
cv      GOTO 40
cv20    CLAMB1 = (0D0,0D0)
cv      GOTO 40
cv30    AMB1 = DLOG(Q2/MF2)
cv     1         +4D0*(PI*PI/3D0 - 1D0)
cv      CLAMB1=DCMPLX(AMB1)
      if (Q2.lt.0.d0) then
         AMB1 = DLOG(-Q2/MF2)
     $        +4D0*(PI*PI/12D0 - 1D0)
         CLAMB1=DCMPLX(AMB1)
         GOTO 40
      elseif (Q2.eq.0.d0) then
         CLAMB1 = (0D0,0D0)
         GOTO 40
      elseif (Q2.gt.0.d0) then
         AMB1 = DLOG(Q2/MF2)
     $        +4D0*(PI*PI/3D0 - 1D0)
         CLAMB1=DCMPLX(AMB1)
      endif

40    RETURN
      END
C=======================================================================
      FUNCTION CLAMB2(Q2,CM2)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 CLAMB2,W,CLN,CSPEN,CM2
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      W=CM2/DCMPLX(Q2,0D0)
      IF(dREAL(W).GT.0) GOTO 100
      CLAMB2=.5D0-2D0*(2D0+W)-(2D0*W+3D0)*CLN(-W)+
     1        2D0*(1D0+W)**2*(CSPEN(1D0+1D0/W)-PI*PI/6D0)
      GOTO 200
100   CLAMB2=.5D0-2D0*(2D0+W)-(2D0*W+3D0)*CLN(W)+
     1  2D0*(1D0+W)**2*(CLN(W)*CLN((W+1D0)/W)-CSPEN(-1D0/W))
     2  -DCMPLX(0D0,1D0)*PI*(3D0+2D0*W-2D0*(1D0+W)**2*CLN((1D0+W)/W))
200   CONTINUE
      RETURN
      END
C=======================================================================
      FUNCTION CLAMB3(Q2,CM2)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 CLAMB3,CM2,W,CLN,CHI,CLH
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      W=CM2/DCMPLX(Q2,0D0)
      IF(dREAL(W).GT.0) GOTO 100
      CHI=SQRT(1D0-4D0*W)
      CLH=CLN((1D0+CHI)/(CHI-1D0))
      CLAMB3=5D0/6D0-2D0/3D0*W+(2D0*W+1D0)/3D0*CHI*CLH
     1        +2D0/3D0*W*(W+2D0)*CLH*CLH
      GOTO 400
cv100   IF (Q2 - 4D0*dREAL(CM2)) 200,300,300
cv200   CHI=SQRT(4D0*W-1D0)
cv      CLH=ATAN(1D0/DREAL(CHI))
cv      CLAMB3=5D0/6D0-2D0/3D0*W+2D0/3D0*(2D0*W+1D0)*CHI*CLH
cv     1        -8D0/3D0*W*(W+2D0)*CLH*CLH
cv      GOTO 400
cv300   CHI=SQRT(1D0-4D0*W)
cv      CLH=CLN( (1D0+CHI)/(1D0-CHI) )
cv      CLAMB3=5D0/6D0-2D0*W/3D0+(2D0*W+1D0)/3D0*CHI*CLH
cv     1       +2D0/3D0*W*(W+2D0)*(CLH*CLH - PI*PI)
cv     2       -DCMPLX(0D0,1D0)*PI*((2D0*W+1D0)/3D0*CHI+
cv     3                            2D0/3D0*W*(W+2D0)*CLH)


 100  IF ((Q2 - 4D0*dREAL(CM2)).lt.0D0) then
         CHI=SQRT(4D0*W-1D0)
         CLH=ATAN(1D0/DREAL(CHI))
         CLAMB3=5D0/6D0-2D0/3D0*W+2D0/3D0*(2D0*W+1D0)*CHI*CLH
     1        -8D0/3D0*W*(W+2D0)*CLH*CLH
         GOTO 400
      elseif ((Q2 - 4D0*dREAL(CM2)).ge.0D0) then
         CHI=SQRT(1D0-4D0*W)
         CLH=CLN( (1D0+CHI)/(1D0-CHI) )
         CLAMB3=5D0/6D0-2D0*W/3D0+(2D0*W+1D0)/3D0*CHI*CLH
     1        +2D0/3D0*W*(W+2D0)*(CLH*CLH - PI*PI)
     2           -DCMPLX(0D0,1D0)*PI*((2D0*W+1D0)/3D0*CHI+
     3        2D0/3D0*W*(W+2D0)*CLH)
         
      endif
      

400   CONTINUE
      RETURN
      END
C=======================================================================
      FUNCTION CLAMB4(Q2,CM1,CM2,MF2)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 CLAMB4,CM1,CM2
     1          ,W1,W2,X1,X2,C12,CLN,CSPEN
      RM1=SQRT(DREAL(CM1))
      RM2=SQRT(DREAL(CM2))
      SM2=(RM1+RM2)*(RM1+RM2)
      DM2=(RM1-RM2)*(RM1-RM2)
      W1=CM1/DCMPLX(Q2,0D0)
      W2=CM2/DCMPLX(Q2,0D0)
      IF (dREAL(CM1).EQ.0) GOTO 401
      IF (dREAL(CM2).EQ.0) GOTO 402
      IF ( (Q2.GT.DM2).AND.(Q2.LT.SM2) ) GOTO 100
      X1=(1D0-W1+W2)/2D0+SQRT((1D0-W1+W2)*(1D0-W1+W2)-4D0*W2)/2D0
      X2=(1D0-W1+W2)/2D0-SQRT((1D0-W1+W2)*(1D0-W1+W2)-4D0*W2)/2D0
      GOTO 200
100   X1=(1D0-W1+W2)/2D0+DCMPLX(0D0,1D0)/2D0
     1                 *SQRT(4D0*W2-(1D0-W1+W2)*(1D0-W1+W2) )
      X2=(1D0-W1+W2)/2D0-DCMPLX(0D0,1D0)/2D0
     1                 *SQRT(4D0*W2-(1D0-W1+W2)*(1D0-W1+W2) )
200   C12=CLN(CM1/CM2)/2D0
      CLAMB4=1D0/6D0 + (W1+W2)/(W1-W2)*C12 - (W1-W2)/3D0*C12
     1  +(W1+W2+1D0)/3D0*(C12-1D0)
     2  +(W1+W2+1D0)/3D0*(X1*CLN(X1/(X1-1D0)) + X2*CLN(-X2/(1D0-X2)) )
     3  -2D0/3D0*(W1+W2+W1*W2)*CLN(X1/(X1-1D0))*CLN(-X2/(1D0-X2))
      GOTO 500
401   W1=W2
      X1=CLN(CM2/MF2)
      X2=CLN((CM2-Q2)/MF2)
      GOTO 404
402   X1=CLN(CM1/MF2)
      X2=CLN((CM1-Q2)/MF2)
404   CLAMB4=.5D0 - W1/3D0 + (1D0-W1*W1)/3D0*CLN((W1-1D0)/W1)
     1       +2D0/3D0*X1 + W1/3D0*(X2*X2-X1*X1+2D0*CSPEN(1D0/(1D0-W1)))
500   RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        FUNCTION SIGFZ(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SIGFZ(Q2): BEITRAEGE ZUR Z0-BOSONSELBSTENERGIE
C                   'REGULAERER' ANTEIL
C       D O U B L E P R E C I S I O N
C-----------------------------------------------------------------------
C       09.11.83
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGFZ
      COMPLEX*16 F,NYANT,ANTZH
     *          ,FQWW,FQZH,FQEE,FQMYMY,FQTATA
     *          ,FQUU,FQCC,FQTT,FQDD,FQSS,FQBB
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
        DMHZ2=MH2-MZ2
        SMHZ2=MH2+MZ2
        FQWW=F(Q2,MW,MW)
        FQZH=F(Q2,MZ,MH)
        FQEE=F(Q2,ME,ME)
        FQMYMY=F(Q2,MMY,MMY)
        FQTATA=F(Q2,MTAU,MTAU)
        FQUU=F(Q2,MU,MU)
        FQCC=F(Q2,MC,MC)
        FQTT=F(Q2,MT,MT)
        FQDD=F(Q2,MD,MD)
        FQSS=F(Q2,MS,MS)
        FQBB=F(Q2,MB,MB)
        DLHW=DLOG(MH/MW)
        DLZW=-DLOG(CW)
        DLZH=DLOG(MZ/MH)
        DLEW2=DLOG(ME2/MW2)
        DLYW2=DLOG(MMY2/MW2)
        DLAW2=DLOG(MTAU2/MW2)
        DLUW2=DLOG(MU2/MW2)
        DLCW2=DLOG(MC2/MW2)
        DLTW2=DLOG(MT2/MW2)
        DLDW2=DLOG(MD2/MW2)
        DLSW2=DLOG(MS2/MW2)
        DLBW2=DLOG(MB2/MW2)
        NYANT=0D0
        IF(Q2 .NE. 0D0) GOTO 10
         ANTZH=0.5D0*SMHZ2+MH2*MZ2/DMHZ2*2D0*DLZH
        GOTO 20
10      ANTZH=FQZH/Q2*DMHZ2*DMHZ2
        IF(Q2 .LT. 0D0) NYANT=5D0/3D0-DLOG(-Q2/MW2)
        IF(Q2 .GT. 0D0) NYANT=DCMPLX(5D0/3D0-DLOG(Q2/MW2),PI)
20      S2C224=1D0/(24D0*SW2*CW2)
        S2C208=3D0*S2C224
        S2_L=8D0*SW2*(SW2-0.5D0)+1D0
        S2_2=SW2/0.28125D0*(SW2-0.75D0)+1D0
        S2_1=SW2/1.125D0*(SW2-1.5D0)+1D0
        SIGFZ=+ALP1PI*(
     G  +DFLOAT(LPAR(15))/(12D0*SW2)*
     G   (-(Q2/1.5D0+(2D1*MW2+1D1*Q2)*FQWW                     )*CW2
     G    +(3D0*MW2*FQWW+Q2/6D0
     G      -MH2*DLHW-MZ2*DLZW
     G      +0.25D0*(1D1*MZ2-2D0*MH2+Q2)*
     G         (1D0+SMHZ2/DMHZ2*DLZH
     G             -DLZW-DLHW      +FQZH   )
     G      +0.25D0*ANTZH                                      )/CW2
     G    +(Q2/6D0+(2D0*MW2+0.25D0*Q2)*FQWW)*(CW2-SW2)*(CW2-SW2)/CW2)
     L  +S2C224*
     L   (S2_L*((2D0*ME2+Q2)*FQEE-Q2/3D0)    -3D0*ME2*FQEE    )
     L  +S2C224*
     L   (S2_L*((2D0*MMY2+Q2)*FQMYMY-Q2/3D0) -3D0*MMY2*FQMYMY )
     L  +S2C224*
     L   (S2_L*((2D0*MTAU2+Q2)*FQTATA-Q2/3D0)-3D0*MTAU2*FQTATA)
     N  +S2C224*Q2*(NYANT+DLEW2)
     N  +S2C224*Q2*(NYANT+DLYW2)
     N  +S2C224*Q2*(NYANT+DLAW2) )
        SIGFZ=SIGFZ+ALP1PI*(
     2  +S2C208*
     2   (S2_2*((2D0*MU2+Q2)*FQUU-Q2/3D0)    -3D0*MU2*FQUU    )
     2  +S2C208*
     2   (S2_2*((2D0*MC2+Q2)*FQCC-Q2/3D0)    -3D0*MC2*FQCC    )
     2  +S2C208*
     2   (S2_2*((2D0*MT2+Q2)*FQTT-Q2/3D0)    -3D0*MT2*FQTT    )
     1  +S2C208*
     1   (S2_1*((2D0*MD2+Q2)*FQDD-Q2/3D0)    -3D0*MD2*FQDD    )
     1  +S2C208*
     1   (S2_1*((2D0*MS2+Q2)*FQSS-Q2/3D0)    -3D0*MS2*FQSS    )
     1  +S2C208*
     1   (S2_1*((2D0*MB2+Q2)*FQBB-Q2/3D0)    -3D0*MB2*FQBB    ) )
        RETURN
        END
        FUNCTION SIGFW(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SIGFW(Q2): BEITRAEGE ZUR W+-BOSONSELBSTENERGIE
C                   NICHT RENORMIERT
C       D O U B L E P R E C I S I O N
C-----------------------------------------------------------------------
C       09.11.83  BER 16.08.85HS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGFW
      COMPLEX*16 F,FQ0W,FQZW,FQHW,FQUD,FQCS,FQTB
     *          ,FQ0E,FQ0MY,FQ0TA
     *          ,ANT1,ANT2,ANT3,ANT4,ANT5,ANT6,ANT7,ANT8,ANT9
     *          ,ANT71,ANT81
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
        FQ0W=F(Q2,0D0,MW)
        FQZW=F(Q2,MZ,MW)
        FQHW=F(Q2,MH,MW)
        ANT71=1D0-DLOG(MZ2/MW2)/SW2+FQZW
        IF (MH2.EQ.MW2) THEN
          ANT81=1D0+MH2/MW2+FQHW
          ELSE
          ANT81=1D0-DLOG(MH2/MW2)*MH2/(MH2-MW2)+FQHW
        ENDIF
        IF (MD2.EQ.MU2) THEN
          DLMUD=0D0
          ELSE
          DLMUD=1D0-(MU2+MD2)/(MU2-MD2)*DLOG(MU/MD)
        ENDIF
        FQ0E=F(Q2,0D0,ME)
        FQ0MY=F(Q2,0D0,MMY)
        FQ0TA=F(Q2,0D0,MTAU)
        FQUD=F(Q2,MU,MD)
        FQCS=F(Q2,MC,MS)
        FQTB=F(Q2,MT,MB)
        IF(Q2 .NE. 0D0) GOTO 10
30      ANT9=MW2
        ANT1=ME2
        ANT2=MMY2
        ANT3=MTAU2
        ANT7=0.5D0*(MW2+MZ2)+MW2*MZ2/(MZ2-MW2)*DLOG(MW2/MZ2)
        IF (MH2.EQ.MW2) THEN
          ANT8=0.5D0*(MW2-MH2)
          ELSE
          ANT8=0.5D0*(MW2+MH2)+MW2*MH2/(MH2-MW2)*DLOG(MW2/MH2)
        ENDIF
        IF (MD2.EQ.MU2) THEN
          ANT4=0.5D0*(MU2-MD2)
          ELSE
          ANT4=0.5D0*(MU2+MD2)+MU2*MD2/(MD2-MU2)*DLOG(MU2/MD2)
        ENDIF
        ANT5=0.5D0*(MC2+MS2)+MC2*MS2/(MS2-MC2)*DLOG(MC2/MS2)
        ANT6=0.5D0*(MT2+MB2)+MT2*MB2/(MB2-MT2)*DLOG(MT2/MB2)
        GOTO 100
10      ANT9=2D0*MW2*MW2*FQ0W/Q2
        ANT1=2D0*ME2*ME2*FQ0E/Q2
        ANT2=2D0*MMY2*MMY2*FQ0MY/Q2
        ANT3=2D0*MTAU2*MTAU2*FQ0TA/Q2
        ANT7=(MZ2-MW2)*(MZ2-MW2)*FQZW/Q2
        ANT8=(MH2-MW2)*(MH2-MW2)*FQHW/Q2
        ANT4=(MD2-MU2)*(MD2-MU2)*FQUD/Q2
        ANT5=(MS2-MC2)*(MS2-MC2)*FQCS/Q2
        ANT6=(MB2-MT2)*(MB2-MT2)*FQTB/Q2
100     SIGFW=+ALP1PI*(
     G +DFLOAT(LPAR(15))/(12D0*SW2)*
     G  (+(-(7D0*(MZ2+MW2)+1D1*Q2)*ANT71
     G     -4D0*MZ2*DLOG(MZ2/MW2)-Q2/1.5D0+2D0*ANT7         )*CW2
     G   -((4D0*MW2+1D1*Q2)*FQ0W+4D0*MW2+Q2/0.09375D0-ANT9  )*SW2
     G   +(1D0-DLOG(MZ2/MW2)/SW2+FQZW       )*3D0*MW2*SW2*SW2/CW2
     G   +( ANT81                                       )*3D0*MW2
     G   -(MZ2*DLOG(MZ2/MW2)+MH2*DLOG(MH2/MW2))*0.5D0+Q2/3D0
     G   -(MW2+MZ2-0.5D0*Q2)*ANT71*0.5D0+ANT7*0.25D0
     G   -(MH2+MW2-0.5D0*Q2)*ANT81*0.5D0+ANT8*0.25D0             )
     L +1D0/(12D0*SW2)*
     L  ((Q2-ME2*0.5D0)*(F(Q2,0D0,ME)+1D0)
     L        -Q2/3D0-0.25D0*ANT1)
     L +1D0/(12D0*SW2)*
     L  ((Q2-MMY2*0.5D0)*(F(Q2,0D0,MMY)+1D0)
     L      -Q2/3D0-0.25D0*ANT2)
     L +1D0/(12D0*SW2)*
     L  ((Q2-MTAU2*0.5D0)*(F(Q2,0D0,MTAU)+1D0)
     L    -Q2/3D0-0.25D0*ANT3) )
       SIGFW=SIGFW+ALP1PI*(
     Q +1D0/(4D0*SW2)*
     Q  (+(Q2-(MU2+MD2)*0.5D0)*
     Q     (+FQUD+DLMUD)
     Q   -Q2/3D0         -0.5D0*ANT4)
     Q +1D0/(4D0*SW2)*
     Q  (+(Q2-(MC2+MS2)*0.5D0)*
     Q     (+FQCS+1D0
     Q      -(MC2+MS2)/(MC2-MS2)*DLOG(MC/MS))
     Q   -Q2/3D0          -0.5D0*ANT5)
     Q +1D0/(4D0*SW2)*
     Q  (+(Q2-(MT2+MB2)*0.5D0)*
     Q     (+FQTB+1D0
     Q      -(MT2+MB2)/(MT2-MB2)*DLOG(MT/MB))
     Q   -Q2/3D0          -0.5D0*ANT6))
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION SIGFM(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SIGFM(Q2): BEITRAEGE ZUR Z-GAMMA-MISCHUNGSENERGIE
C                   'REGULAERER ANTEIL'
C       D O U B L E P R E C I S I O N
C-----------------------------------------------------------------------
C       09.11.83  BER. 16.08.85 HS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGFM
      COMPLEX*16 F,FQEE,FQMYMY,FQTATA
     *          ,FQUU,FQCC,FQTT,FQDD,FQSS,FQBB
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
        FQEE=F(Q2,ME,ME)
        FQMYMY=F(Q2,MMY,MMY)
        FQTATA=F(Q2,MTAU,MTAU)
        FQUU=F(Q2,MU,MU)
        FQCC=F(Q2,MC,MC)
        FQTT=F(Q2,MT,MT)
        FQDD=F(Q2,MD,MD)
        FQSS=F(Q2,MS,MS)
        FQBB=F(Q2,MB,MB)
        SIGFM=+ALP1PI*( +DFLOAT(LPAR(15))/4D0/SW2*
     G  (-(+(3D0*SW*CW+SW/(CW*6D0))*Q2
     G     +(4D0*SW*CW+SW/(CW*0.75D0))*MW2)*F(Q2,MW,MW)-Q2*SW/CW/9D0)
     L  +1D0/(3D0*SW*CW)*
     L   (0.25D0-SW2)*  ((2D0*ME2+Q2)*FQEE-Q2/3D0)
     L  +1D0/(3D0*SW*CW)*
     L   (0.25D0-SW2)*  ((2D0*MMY2+Q2)*FQMYMY-Q2/3D0)
     L  +1D0/(3D0*SW*CW)*
     L   (0.25D0-SW2)*  ((2D0*MTAU2+Q2)*FQTATA-Q2/3D0)
     2  +1D0/(2.25D0*SW*CW)*
     2   (0.375D0-SW2)*  ((2D0*MU2+Q2)*FQUU-Q2/3D0)
     2  +1D0/(2.25D0*SW*CW)*
     2   (0.375D0-SW2)*  ((2D0*MC2+Q2)*FQCC-Q2/3D0)
     2  +1D0/(2.25D0*SW*CW)*
     2   (0.375D0-SW2)*  ((2D0*MT2+Q2)*FQTT-Q2/3D0)
     1  +1D0/(9D0*SW*CW)*
     1   (0.75D0-SW2)*  ((2D0*MD2+Q2)*FQDD-Q2/3D0)
     1  +1D0/(9D0*SW*CW)*
     1   (0.75D0-SW2)*  ((2D0*MS2+Q2)*FQSS-Q2/3D0)
     1  +1D0/9D0/SW/CW*(0.75D0-SW2)*((2D0*MB2+Q2)*FQBB-Q2/3D0))
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  RENORMIERTE SELBSTENERGIEN
C  IN DEN UNTERPROGRAMMEN
C     SIGMRG, SIGMRM, SIGMRZ, SIGMRW
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION SIGMRG(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SIGMRG(Q2): BEITRAEGE ZUR PHOTONSELBSTENERGIE  (RENORMIERT)
C              Q2 = Q**2 VIERERIMPULSQUADRAT
C-----------------------------------------------------------------------
C       09.11.83
C       14.05.91 HS: BURKHARD'S PARAMETRIZATION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGMRG,F
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART

      IF (LPAR(7).EQ.1) THEN
C...HADRONIC CONTRIBUTION FROM EFFECTIVE QUARK LOOPS
        SIGMRG=+ALP1PI*(
     G   -DFLOAT(LPAR(15))*((3D0*Q2+4D0*MW2)*F(Q2,MW,MW) - Q2/1.5D0)/4D0
     L   +(+(Q2+2D0*ME2  )*F(Q2,ME  ,ME  )
     L     +(Q2+2D0*MMY2 )*F(Q2,MMY ,MMY )
     L     +(Q2+2D0*MTAU2)*F(Q2,MTAU,MTAU) - Q2)/3D0
     2   +(+(Q2+2D0*MU2  )*F(Q2,MU  ,MU  )
     2     +(Q2+2D0*MC2  )*F(Q2,MC  ,MC  )
     2     +(Q2+2D0*MT2  )*F(Q2,MT  ,MT  ) - Q2)/2.25D0
     1   +(+(Q2+2D0*MD2  )*F(Q2,MD  ,MD  )
     1     +(Q2+2D0*MS2  )*F(Q2,MS  ,MS  )
     1     +(Q2+2D0*MB2  )*F(Q2,MB  ,MB  ) - Q2)/9D0        )
      ELSEIF(LPAR(7).GE.2) THEN
C...HADRONIC CONTRIBUTION FROM PARAMETRIZATION
        SIGMRG=+ALP1PI*(
     G   -DFLOAT(LPAR(15))*((3D0*Q2+4D0*MW2)*F(Q2,MW,MW) - Q2/1.5D0)/4D0
     L   +(+(Q2+2D0*ME2  )*F(Q2,ME  ,ME  )
     L     +(Q2+2D0*MMY2 )*F(Q2,MMY ,MMY )
     L     +(Q2+2D0*MTAU2)*F(Q2,MTAU,MTAU) - Q2)/3D0
     2   +(+(Q2+2D0*MT2  )*F(Q2,MT  ,MT  ) - Q2/3D0)/2.25D0     )
     H   - DCMPLX(HADRQQ(Q2)*Q2,0D0)
      ELSEIF(LPAR(7).EQ.0) THEN
       SIGMRG=DCMPLX(0D0,0D0)
      ELSE
      ENDIF
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      FUNCTION SIGMRM(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SIGMRM(Q2): BEITRAEGE ZUR GAMMA-Z-MISCHUNG (RENORMIERT)
C              Q2 = Q**2 VIEREIMPULSQUADRAT
C-----------------------------------------------------------------------
C       21.10.83
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGMRM,SIGFM,SIGFZ,SIGFW,CAFINW,CAFINZ
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
        SQ=ALP4PI*(DLOG(MU2/MD2)+DLOG(MC2/MS2)+DLOG(MT2/MB2))
        SP=ALP4PI*DLOG(MT2/MB2)*(MT2-MB2)/4D0/SW2/MW2
        CAFINZ=SIGFZ(MZ2)
        CAFINW=SIGFW(MW2)
        RDMZ2=DREAL(CAFINZ)
        RDMW2=DREAL(CAFINW)
        MRENK1=CW/SW*(RDMZ2/MZ2 - RDMW2/MW2) +SQ/6D0/SW/CW + SP*CW/SW
        SIGMRM=
     *   -(SIGFM(Q2)+ Q2 *MRENK1)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      FUNCTION SIGMRZ(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SIGMRZ(Q2): BEITRAEGE ZUR Z-BOSONSELBSTENERGIE (RENORMIERT)
C              Q2 = Q**2 VIEREIMPULSQUADRAT
C-----------------------------------------------------------------------
C       21.10.83
C       14.05.91 HS (HADRONIC VACUUM POLARIZATION)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 CJFINZ,CAFINZ,CAFINW
     *          ,SIGMRZ,SIGFZ,SIGFW
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
        SQ=ALP4PI*(DLOG(MU2/MD2)+DLOG(MC2/MS2)+DLOG(MT2/MB2))
        SP=ALP4PI*DLOG(MT2/MB2)*(MT2-MB2)/4D0/SW2/MW2
        CAFINZ=SIGFZ(MZ2)
        CAFINW=SIGFW(MW2)
        RDMZ2=DREAL(CAFINZ)
        RDMW2=DREAL(CAFINW)
        MRENK0=(CW2/SW2-1D0)*(RDMZ2/MZ2 - RDMW2/MW2) + ALP2PI/3D0
     *         + SQ/3D0/SW2 + SP*(CW2/SW2-1D0)
        CJFINZ=DCMPLX(RDMZ2,0D0)
        SIGMRZ=
     *    (SIGFZ(Q2)-CJFINZ+(Q2-MZ2)*MRENK0)
        IF (LPAR(7).GE.2) THEN
          DDALPP=DSGMRG(MZ2)*(Q2-MZ2)
          SIGMRZ=SIGMRZ+DCMPLX(DDALPP,0D0)
        ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      FUNCTION SIGMRW(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SIGMRW(Q2): BEITRAEGE ZUR Z-BOSONSELBSTENERGIE (RENORMIERT)
C              Q2 = Q**2 VIEREIMPULSQUADRAT
C-----------------------------------------------------------------------
C       21.10.83
C       14.05.91 HS (HADRONIC VACUUM POLARIZATION)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 SIGMRW,CJFINW,CAFINZ,CAFINW
     *          ,SIGFZ,SIGFW
      COMMON /KONST/  PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /GSW/    SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PARLIS/ LPAR(20),LPARIN(12),IPART
        SQ=ALP4PI*(DLOG(MU2/MD2)+DLOG(MC2/MS2)+DLOG(MT2/MB2))
        SP=ALP4PI*DLOG(MT2/MB2)*(MT2-MB2)/4D0/SW2/MW2
        CAFINZ=SIGFZ(MZ2)
        CAFINW=SIGFW(MW2)
        RDMW2=DREAL(CAFINW)
        RDMZ2=DREAL(CAFINZ)
        MRENKW=CW2/SW2*(RDMZ2/MZ2 - RDMW2/MW2) + ALP2PI/3D0
     *         + SQ/3D0/SW2 + SP*CW2/SW2
        CJFINW=DCMPLX(RDMW2,0D0)
        SIGMRW=
     *    (SIGFW(Q2)-CJFINW+(Q2-MW2)*MRENKW)
        IF (LPAR(7).GE.2) THEN
          DDALPP=DSGMRG(MZ2)*(Q2-MW2)
          SIGMRW=SIGMRW+DCMPLX(DDALPP,0D0)
        ENDIF
      RETURN
      END
