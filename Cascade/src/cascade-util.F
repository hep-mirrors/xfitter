      SUBROUTINE CASCADE
      IMPLICIT None
c common for bases/spring 5.1
      Integer Mxdim
      PARAMETER (MXDIM = 50 )
      Double precision xlb,xub
      Integer ndimb,nwildb,ig,ncallb
      COMMON /BPARM1/ XLB(MXDIM),XUB(MXDIM),NDIMB,NWILDB,
     +   IG(MXDIM),NCALLB
      Double precision xlmc,xumc
      Integer ndimmc,ncallmc
      Common /MCINT1/XLMC(MXDIM),XUMC(MXDIM),NDIMMC,NCALLMC
      Double Precision acc1b,acc2b
      Integer itmx1b,itmx2b
      COMMON /BPARM2/ ACC1B,ACC2B,ITMX1B,ITMX2B
#include "caslujets.inc"
#include "casluco.inc"
#include "caspart.inc"
#include "caspara.inc"
#include "caspara1.inc"
#include "caseffic.inc"
#include "casgki.inc"
#include "casludat1.inc"
#include "casludat3.inc"
#include "casshwr.inc"
#include "casprkt.inc"
#include "casvalues.inc"
#include "casscal.inc"
#include "caspara5.inc"
#include "casmult.inc"
#include "casherup.inc"
#include "caspsipol.inc"
	Integer kf
	Common /cbran1/ kf
      REAL ACC1,ACC2
      Integer iint,ncb
      COMMON /INTEGR/ ACC1,ACC2,IINT,NCB
      Integer ndim
      COMMON/INTDIM/ NDIM
      REAL ANS,ERROR
      Double Precision PYPAR,PYVAR
      Integer IPY
      COMMON /MYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      Integer Ndimc
      COMMON /DIMEN/ NDIMC
      double Precision XD(20)
      Double Precision XGR(20)
      CHARACTER CHIN*100
      CHARACTER * 16 CNAM
      CHARACTER * 7 CINT,CPSI,CCHI
      CHARACTER*20 VERSQQ
      Integer Lerr
      COMMON/ERR/LERR(100)
      Integer IGENFL,IMEPHI
      COMMON/GENWEI/IGENFL,IMEPHI
      Integer IGLU
	Common/CAGLUON/Iglu
      REAL SNGL
      double precision DFUN,FXN1,FXNB
      EXTERNAL DFUN,FXN1,FXNB
      Integer I,J,KPCH,KPAO,IT1,IT2
      Double Precision SQRTS,DDUM,CTIME
      Double Precision DOT
	Integer MSTP,MSTI
	Double Precision PARP,PARI
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      Integer ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON /CASCOL/ ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      Integer pycomp,kc
      External pycomp
      Integer Lun6
      Data Lun6/16/
#include "casvers.inc"
      IGENFL = 0
      DO 10  I=1,100
   10 LERR(I) = 0
      ERROR = 0.0
ccheck      CHIN='PARU(11)=0.010;'
ccheck      CALL PYGIVE(CHIN)
      ALPH=7.299D-3
      ALPHS=0.3D0
      PI=4.D0*DATAN(1.D0)
      IWEI = 0
C...  GIVE ELECTRON FOUR VECTOR
      K(1,1)=21
      K(1,2)=KE
      P(1,1) = 0.0D0
      P(1,2) = 0.0D0
      P(1,3) = DBLE(PLEPIN)
      P(1,5) = PYMASS(KE)
c      P(1,5) = 0
      P(1,4) = DSQRT(P(1,1)**2+P(1,2)**2+P(1,3)**2+P(1,5)**2)
C...  GIVE PROTON FOUR VECTOR
      K(2,1) = 21
      K(2,2) = KP
      P(2,1) = 0.0D0
      P(2,2) = 0.0D0
      P(2,3) = DBLE(PPIN)
      P(2,5) = PYMASS(KP)
c      P(2,5) = 0
      P(2,4) = DSQRT(P(2,1)**2+P(2,2)**2+P(2,3)**2+P(2,5)**2)
      N=2
ccheck      CALL PYLIST(1)
C... CALCULATE CMS ENERGY
      DO 20 I=1,4
         CM(I) =P(1,I)+P(2,I)
   20 CONTINUE

C BOOST TO EP CMS
      CALL PYROBO(0,N,0.D0,0.D0,-CM(1)/CM(4),-CM(2)/CM(4),-CM(3)/CM(4))
      DO 30 I=1,2
         DO 30 J=1,5
            PBEAM(I,J)=P(I,J)
   30 KBEAM(I,J)=K(I,J)
      KINT(1,2) = 22
      KINT(2,2)= KP
      KEB = KPH
      SSS = DOT(CM,CM)
      SQRTS=DSQRT(SSS)
      IF(IPRO.EQ.10.OR.IPRO.EQ.11.OR.IPRO.EQ.10.OR.IPRO.EQ.15.OR.
     & IPRO.EQ.2.or.IPRO.eq.3.OR.ipro.eq.4.or.
     & IPRO.eq.102.or.
     & IPRO.ge.500.and.IPRO.le.515) THEN
      ELSE
         write(lun6,*) ' wrong subprocess selected ', IPRO
         write(lun6,*) ' PROGRAM STOPPED '
         STOP
      ENDIF
      CALL PYNAME(KINT(2,2),CNAM)
      CINT = CNAM(1:7)
      KPCH = 213
      CALL PYNAME(KPCH,CNAM)
      write(lun6,10000)

      write(lun6,10100)

      write(lun6,10200)

      write(lun6,10300) VERSQQ
	
	kf = 2212
	
      if(Iabs(K(1,2)).eq.2212) then
c proton proton
         write(lun6,*) '*       proton proton version selected            '
     +   //'  *'
         KEB=21
	   IRES(1)=1
	   IRES(2)=1
         QMI = 0.
         YMI = 0.
         YMA = 1. 
         
         IF(IPRO.EQ.2) THEN
            if(ihfla.eq.3) then
               kpsi = 333
            elseif(ihfla.eq.4) then
               kpsi = 443
            elseif(ihfla.eq.5) then
               kpsi = 553
            else
               kpsi = 443
            endif

            CALL PYNAME(Kpsi,CNAM)
            Cpsi = CNAM(1:7)
         
            IF(INTER.EQ.0) write(lun6,11001) Cpsi
            if(Ipsipol.eq.1) then
               write(lun6,*) '*      including lepton polarisation ',
     &         Ipsipol,'  *' 
               write(lun6,*) '*      only leptonic decay calculated  ',
     &         '             *' 
     	         write(lun6,*) '*    No final state PS ',
     &         ' with lepton polarisation    *'
     	         write(lun6,*) '*      Final state PS ',
     &         'from gluon switched off       *'
               if(IFPS.eq.2) IFPS=0
               if(IFPS.eq.3) IFPS=1
c               write(lun6,*) ' original decay table '
c               call pylist(12)
               KC = pycomp(kpsi)
c           switch off all decay channels               
               do i=mdcy(kc,2),mdcy(kc,2)+mdcy(kc,3)-1 
                  mdme(i,1) = 0
               end do
               mdme(mdcy(kc,2),1) = 1
               mdme(mdcy(kc,2)+1,1) = 1
c               write(lun6,*) ' modified decay table ',mdcy(kc,2),mdcy(kc,3)
c               call pylist(12)
            endif
         ELSEIF(IPRO.EQ.3) THEN
            if(ihfla.eq.3) then
               kchi = 10331
            elseif(ihfla.eq.4) then
               kchi = 10441
               kchi = 445
            elseif(ihfla.eq.5) then
               kchi = 10551
            else
               kchi = 10441
            endif

            CALL PYNAME(Kchi,CNAM)
            Cchi = CNAM(1:7)
         
            IF(INTER.EQ.0) write(lun6,11002) Cchi
         ENDIF

         IF(IPRO.EQ.10) Then
            IHFLA = 3
c            IF(IRPA.ne.0) write(lun6,11300)
            IF(IRPA.NE.0) THEN
               IF(IHFLA.LT.4) write(lun6,11205)
               IF(IHFLA.GE.4) write(lun6,11305)
            ENDIF
            IF(IRPB.NE.0) write(lun6,11405)
            IF(IRPC.NE.0) THEN
               IF(IHFLA.LT.4) write(lun6,11505)
               IF(IHFLA.GE.4) write(lun6,11605)
            ENDIF
            IF(IRPD.NE.0) write(lun6,11705)
            IF(IRPE.NE.0) THEN
               IF(IHFLA.LT.4) write(lun6,11805)
               IF(IHFLA.GE.4) write(lun6,11905)
            ENDIF
            IF(IRPF.NE.0) THEN
               IF(IHFLA.LT.4) write(lun6,12005)
               IF(IHFLA.GE.4) write(lun6,12105)
            ENDIF
         Endif
         
         IF(IPRO.EQ.11) THEN
            IRPA = 1
            write(lun6,11700)
            write(lun6,11800) MAX(4,IHFLA)
            KPA=4
         ENDIF
c         IF(IPRO.EQ.15) write(lun6,11304)
         IF(IPRO.EQ.102) write(lun6,11302)
         IF(IPRO.EQ.500) write(lun6,14500)
         IF(IPRO.EQ.501) write(lun6,14501)
         IF(IPRO.EQ.503) write(lun6,14503)
         IF(IPRO.EQ.504) write(lun6,14504)
         IF(IPRO.EQ.510) write(lun6,14510)
         IF(IPRO.EQ.511) write(lun6,14511)
         IF(IPRO.EQ.513) write(lun6,14513)
         IF(IPRO.EQ.514) write(lun6,14514)
         IF(IPRO.EQ.504.OR.IPRO.EQ.514) THEN
            write(lun6,11800) MAX(4,IHFLA)
            KPA=MAX(4,IHFLA)
         ENDIF
c	   write(lun6,11801) Mult
      elseif(Iabs(K(2,2)).eq.22.or.
     & (Iabs(K(1,2)).eq.11.and.Iabs(K(2,2)).eq.11)) then
	if(Iabs(K(2,2)).eq.22) then 
c gamma gamma 
         write(lun6,*) '*     gamma gamma version selected                '
     +   //'  *'
         else
         write(lun6,*) '*     e+ e-   version selected                    '
     +   //'  *'
         endif
         KEB=22
	   kf=keb
         KINT(1,2) = 22
         KINT(2,2) = 22 
         IRPA = 1 
         IF(IPRO.EQ.10) write(lun6,11301)
         IF(IPRO.EQ.11) THEN
            write(lun6,11701)
            write(lun6,11800) MAX(4,IHFLA)
            KPA=4
         ENDIF
         IF(IPRO.EQ.102) write(lun6,11302)
        write(lun6,*) '*     beam 1:  direct/resolved (0/1): ',Ires(1),
     +   '            *'
        write(lun6,*) '*     beam 2:  direct/resolved (0/1): ',Ires(2),
     +   '            *'
	   
      else
	   IRES(2)=1
         IRPA = 1
         IF(INTER.EQ.0) write(lun6,10400)
         IF(INTER.EQ.2) write(lun6,10500)
         IF(IPRO.EQ.10) THEN
            IF(NFLAV.LE.3) THEN
               IF(INTER.EQ.0) write(lun6,10600) CINT
               IF(INTER.EQ.2) write(lun6,10700) CINT
            ELSE
               IF(INTER.EQ.0) write(lun6,10800) CINT,NFLAV
               IF(INTER.EQ.2) write(lun6,10900) CINT,NFLAV
            ENDIF
         ENDIF
         if(ihfla.eq.3) then
         kpsi = 333
         elseif(ihfla.eq.4) then
         kpsi = 443
         elseif(ihfla.eq.5) then
         kpsi = 553
         else
         kpsi = 443
         endif

         CALL PYNAME(Kpsi,CNAM)
         Cpsi = CNAM(1:7)
         
         IF(IPRO.EQ.2) THEN
            IF(INTER.EQ.0) write(lun6,11000) CINT,Cpsi
            IF(INTER.EQ.2) write(lun6,11100) CINT
            if(Ipsipol.eq.1) then
               write(lun6,*) '*      including lepton polarisation ',
     &         Ipsipol,'  *' 
               write(lun6,*) '*      only leptonic decay calculated  ',
     &         '             *' 
     	         write(lun6,*) '*    No final state PS ',
     &         ' with lepton polarisation    *'
     	         write(lun6,*) '*      Final state PS ',
     &         'from gluon switched off       *'
               if(IFPS.eq.2) IFPS=0
               if(IFPS.eq.3) IFPS=1
c               write(lun6,*) ' original decay table '
c               call pylist(12)
               KC = pycomp(kpsi)
c           switch off all decay channels               
               do i=mdcy(kc,2),mdcy(kc,2)+mdcy(kc,3)-1 
                  mdme(i,1) = 0
               end do
               mdme(mdcy(kc,2),1) = 1
               mdme(mdcy(kc,2)+1,1) = 1
c               write(lun6,*) ' modified decay table ',mdcy(kc,2),mdcy(kc,3)
c               call pylist(12)
            endif
         ENDIF
        IF(IPRO.EQ.4) then
           write(lun6,11516)
        endif
         IF(IPRO.EQ.10) write(lun6,11200)

         IF(ISEMIH.EQ.1) THEN
            IF(IPRO.EQ.2) then
               write(lun6,11500)
            else
               write(lun6,11400)
            endif
         ENDIF
         IF(IPRO.EQ.11) THEN
            IF(INTER.EQ.0) THEN
               write(lun6,11600) CINT
               write(lun6,11800) MAX(4,IHFLA)
            ENDIF
            IF(INTER.EQ.2) write(lun6,11900) CINT
            KPA=4
         ENDIF
         IF(IPRO.EQ.11) write(lun6,12000)

         IHF = 0
         IF(THEMA.NE.180.0D0) THEN
            write(lun6,12100) THEMA

         ELSE
            write(lun6,*) '*  no cut on max angle of scattered electron   '
     +      //'     *'
         ENDIF
         IF(THEMI.NE.0.0D0) THEN
            write(lun6,12200) THEMI

         ELSE
            write(lun6,*) '*  no cut on min angle of scattered electron   '
     +      //'     *'
         ENDIF
         IF(QMI.GT.0.0) THEN
            write(lun6,12300) QMI

         ELSE
            write(lun6,*) '*  Q^2 _min according to kinematics            '
     +      //'     *'
         ENDIF
         IF(QMA.LT.10D8) THEN
            write(lun6,12400) QMA

         ELSE
            write(lun6,*) '*  Q^2 _max according to kinematics            '
     +      //'     *'
         ENDIF
         IF(YMI.GT.0.0) THEN
            write(lun6,12500) YMI

         ELSE
            write(lun6,*) '*  y_min according to kinematics               '
     +      //'     *'
         ENDIF
         IF(YMA.LT.1.0) THEN
            write(lun6,12600) YMA

         ELSE
            write(lun6,*) '*  y_max according to kinematics               '
     +      //'     *'
         ENDIF
      endif
      IF(IPRO.eq.3.or.IPRO.eq.102) then
      else
      IF(PT2CUT(IPRO).GT.0D0) write(lun6,17100) DSQRT(PT2CUT(IPRO)),IPRO
      endif 
      IF(Ipro.eq.10.or.IPRO.eq.11) then
         IF(PT2CUT(IPRO).GT.0D0) then
            PT2CUT(10) = PT2CUT(IPRO)
            PT2CUT(11) = PT2CUT(10)
c            write(lun6,11888)
         Endif
      endif
      IPY(8) = NFLAV

      write(lun6,*) '*****************************************************'
      IF(IFPS.EQ.0) THEN
         IPY(13) = 0
         IPY(14) = 0
      ENDIF
      IF(IFPS.NE.0) THEN
         write(lun6,*) '##################################################'
     +   //'###'
         write(lun6,*) '#     parton shower selection:                    '
     +   //'  #'
         IF(IFPS.EQ.0) THEN
            IPY(13) = 0
            IPY(14) = 0
         ELSEIF(IFPS.EQ.1) THEN
            IPY(13) = 0
            IPY(14) = 1		
            if(iglu.eq.1.or.iglu.ge.6) then
		   else
               IPY(14) = 0
               write(lun6,*) '############################################'
     +         //'#########'
               write(lun6,*) '#     no initial state PS possible for IGLU>'
     +         //'1       #'
               write(lun6,*) '#        initial state PS switched off      '
     +         //'        #'
               write(lun6,*) '#     see manual.                           '
     +         //'        #'
            endif
         ELSEIF(IFPS.EQ.2) THEN
            IPY(13) = 1
            IPY(14) = 0
         ELSEIF(IFPS.EQ.3) THEN
            IPY(13) = 1
            IPY(14) = 1
            if(ires(1).eq.1.and.iabs(K(1,2)).eq.11) then
               IPY(14) = 0
               write(lun6,*) '############################################'
     +         //'#########'
               write(lun6,*) '#     no initial state PS possible for resol'
     +         //'ved     #'
               write(lun6,*) '#     initial state PS switched off         '
     +         //'        #'
               write(lun6,*) '#     see manual.                           '
     +         //'        #'
            endif
            if(iglu.eq.1.or.iglu.ge.6) then
		   else
               IPY(14) = 0
               write(lun6,*) '############################################'
     +         //'#########'
               write(lun6,*) '#     no initial state PS possible for IGLU>'
     +         //'1 or 6    #'
               write(lun6,*) '#        initial state PS switched off      '
     +         //'        #'
               write(lun6,*) '#     see manual.                           '
     +         //'        #'
            endif
         ENDIF
         if(IPY(14).EQ.1.and.Iccfm.eq.1) Then
            write(lun6,*) '#     CCFM inital state parton showers '
     +      //'including    #'
            write(lun6,*) '#     angular ordering                         '
     +      //'     #'
         endif
         if(IPY(14).EQ.1.and.Iccfm.eq.0) Then
            write(lun6,*) '#     DGLAP inital state parton showers '
     +      //'including   #'
            write(lun6,*) '#     angular ordering                         '
     +      //'     #'
         endif
         if(itimshr.eq.1)
     +   write(lun6,*) '#     timelike inital state partons may shower    '
     +   //'  #'
         if(itimshr.eq.0)
     +   write(lun6,*) '#     no timelike inital state partons            '
     +   //'  #'
         IF(IPY(13).EQ.1) write(lun6,*) '#     final state parton shower  '
     +   //'                   #'
         write(lun6,*) '##################################################'
     +   //'###'
         write(lun6,*) ' '
      ENDIF
	if(ILHA.le.0) then
	write(lun6,*)  '##################################################'
	write(lun6,*)  '#           proton remnant parameters            #'
	write(lun6,*)  '#           treated inside CASCADE               #'
	write(lun6,12701) IRspl
12701	format(' #    energy sharing IREM = ',I2,'                    #')
	write(lun6,*)  '##################################################'
	elseif(ILHA.eq.1) then
	write(lun6,*)  '##################################################'
	write(lun6,*)  '#             output of events in                #' 
	write(lun6,*)  '#          Les Houges Accord Format              #'
	write(lun6,*)  '##################################################'
	elseif(ILHA.gt.1) then
         MSTP(61)=0 ! initial state PS
         MSTP(81)=0 ! multiple interactions
         MSTP(71)=0 
         IPY(13) = 0		
         IPY(14) = 0		
         if(ifps.eq.1) then 
c            MSTP(71)=0 
            ifps = 0
            IPY(14) = 1		
            endif
         if(ifps.eq.2) then 
c            MSTP(71)=1 
            ifps = 0
            IPY(14) = 0		
            IPY(13) = 1		
            endif
         if(ifps.eq.3) then
c            MSTP(71)=1 
            IFPS = 1
            IPY(14) = 1
            IPY(13) = 1		
         Endif
	write(lun6,*)  '##################################################'
	write(lun6,*)  '#              proton remnant and                #' 
	write(lun6,*)  '#                final state PS                  #'
	write(lun6,*)  '#               treated in PYTHIA                #'
	write(lun6,*)  '##################################################'
	
	endif
      write(lun6,12700) SQRTS
c      IF(IRUNAEM.EQ.0) THEN
c         write(lun6,*) ' alpha_em fixed; alpha_em = ',ALPH
c      ELSE
c         write(lun6,*) ' running alpha_em selected '
c      ENDIF
      IF(IRUNA.EQ.0) THEN
         MSTU(111) = 0
         PARU(111) = ALPHS
         write(lun6,*) ' alpha_s fixed; alpha_s = ',ALPHS
      ENDIF
      IF(IQ2.EQ.1) write(lun6,*) ' scale for alpha_s: 4*m_q**2 '
      IF(IQ2.EQ.2) write(lun6,*) ' scale for alpha_s: shat '
      IF(IQ2.EQ.3) write(lun6,*) ' scale for alpha_s: m_q**2 + p_t **'
     +//'2 '
      IF(IQ2.EQ.4) write(lun6,*) ' scale for alpha_s: Q2 '
      IF(IQ2.EQ.5) write(lun6,*) ' scale for alpha_s: Q2 + p_t **2'
      IF(IQ2.EQ.6) write(lun6,*) ' scale for alpha_s: kt**2 '
c      IF(IQ2.EQ.7) write(lun6,*) ' scale for alpha_s:'
c     +//'two*(sL*rt*ru)/((sL*rt*ru)**2+rt**2+ru**2)'
      IF(IQ2.EQ.7) write(lun6,*) ' scale for alpha_s:'
     +//' max(stmax1,stamx2)'

      write(lun6,*) ' scale or alpha_s is multiplied by: ',SCALFA
c      write(lun6,13300) PARU(112),MSTU(112)
c      PYPAR(21) = PARU(112)



      KPAO = KPA


      IF(IPRO.EQ.10.OR.IPRO.EQ.11) NDIM=8
      NDIMC = 9999
   40 CALL dcasrnv(XGR,20)
      DO 50 I=1,20
   50 XD(I)=XGR(I)
      DDUM=DFUN(20,XD)
      NDIMC = 4
      if(IPRO.eq.10) NDIMC=6  ! this is only for using in fit 
      if(IPRO.eq.11) NDIMC=6  ! this is only for using in fit 
c         write(lun6,*) ' NDIM dynamically calculated ',NDIM,NDIMC,IPRO
      IF(NDIMC.GT.20) goto 40
      NDIM = NDIMC
      NIN=0
      NOUT=0
      IF(iint.eq.0) then
********************************************************************
*     Initialization of BASES/SPRING V5.1
********************************************************************
*===> Initialization of BASES by calling BSINIT
c
         CALL BSINIT
         NDIMB=NDIM
         NWILDB = NDIM
         NCALLB = NCB
         ACC1B = DBLE(ACC1)
         ITMX1B = 30
         ITMX2B = 100
c         ITMX1B = 10
c         ITMX2B = 10
         ACC2B = DBLE(ACC2)
         DO 60 I=1,NDIMB
           IG(I) = 1
           XLB(I)=1.0D-12
   60    XUB(I)=1.D0
********************************************************************
*     Nimerical Integration by BASES V5.1
********************************************************************

         CALL BASES( FXNB, AVGI, SD, CTIME, IT1, IT2 )
         ANS = SNGL(AVGI)
         ERROR = SNGL(SD)
         SD=DBLE(ERROR)
         AVGI=DBLE(ANS)
      elseif(iint.eq.1) then
         NDIMMC=NDIM
         DO I=1,NDIMMC
           IG(I) = 1
           XLMC(I)=1.0D-12
           XUMC(I)=1.D0
         End do
         NCALLMC = NCB
         Call MCINT(FXNB,AVGI,SD)
         ANS = SNGL(AVGI)
         ERROR = SNGL(SD)
         SD=DBLE(ERROR)
         AVGI=DBLE(ANS)
      else
      write(lun6,*)' CASCADE: IINT = ',IINT,' not implemented: STOP '
      endif


      KPA = KPAO
      RETURN

10000 FORMAT(' *****************************************************')
10100 FORMAT(' *                                                   *')
10200 FORMAT(' *     You are using the CASCADE MC generator        *')
10300 FORMAT(' *           version  ',A20,'           *')
10400 FORMAT(' *    neutral current interaction selected           *')
10500 FORMAT(' *    charged current interaction selected           *')
10600 FORMAT(' *    gamma + gluon_',A6,' --> q q_bar(light quarks) *')
10700 FORMAT(' *      W + gluon_',A6,' --> q q_bar(light quarks)   *')
10800 FORMAT(' *    gamma + gluon_',A6,' --> q q_bar max flav =',I2,
     +'  *')
10900 FORMAT(' *      W + gluon_',A6,' --> q q_bar max flav =',I2,
     +'    *')
11000 FORMAT(' *    gamma + gluon_',A6,' --> ',A5,' + gluon',
     +'         *')
11001 FORMAT(' *        gluon + gluon   --> ',A5,' + gluon',
     +'          *')
11002 FORMAT(' *        gluon + gluon   --> ',A5,'        ',
     +'          *')
11205 FORMAT(' *               g + g     --> q + q_bar             *')
11305 FORMAT(' *               g + g     --> Q + Q_bar             *')
11405 FORMAT(' *               g + g     --> g + g                 *')
11505 FORMAT(' *               g + q     --> g + q                 *')
11515 FORMAT(' *              g* + q*    --> q                     *')
11516 FORMAT(' *           gamma* + q*   --> q                     *')
11605 FORMAT(' *               g + Q     --> g + Q                 *')
11705 FORMAT(' *               q + q_bar --> g + g                 *')
11805 FORMAT(' *               q + q_bar --> q + q_bar             *')
11905 FORMAT(' *               q + q_bar --> Q + Q_bar             *')
12005 FORMAT(' *               q + q     --> q + q                 *')
12105 FORMAT(' *               q + q     --> Q + Q                 *')
11100 FORMAT(' *      W + gluon_',A6,' --> q q_bar max flav =',I2,
     +'    *')
11200 FORMAT(' *    EPA + gamma* gluon --> q q_bar  used           *')
11300 FORMAT(' *    gluon* gluon* --> q q_bar  used                  *')
11301 FORMAT(' *     gamma gamma --> q q_bar  used                   *')
11302 FORMAT(' *       g g --> h0                                    *')
14503 FORMAT(' *       g g --> Z q qbar                            *')
14504 FORMAT(' *       g g --> Z Q Qbar                            *')
14513 FORMAT(' *       g g --> W q qbar                            *')
14514 FORMAT(' *       g g --> W Q Qbar                            *')
14500 FORMAT(' *       q_i qbar_i --> Z                            *')
14510 FORMAT(' *       q_i q_j --> W                               *')
14501 FORMAT(' *       q g --> Z q                                 *')
14511 FORMAT(' *       q g --> W q                                 *')
11404 FORMAT(' *       q g --> q g                                 *')
11400 FORMAT(' *  semihard approach for BGF Catani et al           *')
11500 FORMAT(' *  semihard approach for BGF Baranov/Zotov et al    *')
11600 FORMAT(' *     gamma + gluon_',A4,
     +   ' --> Q Q_bar                *')
11700 FORMAT(' *       gluon* + gluon*',
     +   ' --> Q Q_bar                 *')
11701 FORMAT(' *       gamma  + gamma',
     +   ' --> Q Q_bar                 *')
11800 FORMAT(' *       heavy flavor produced is :',
     +    I3,'               *')
11801	   format(' *       multiple scattering selected: ',I2,
     +'            *')
11900 FORMAT(' *      W + gluon_',A4,
     +   ' --> c c_bar                   *')
12000 FORMAT(' *    EPA + gamma gluon --> Q Q_bar (massive)   used   *')
12100 FORMAT(' *    scattered electron theta_max = ',F7.3,'          *')
12200 FORMAT(' *    scattered electron theta_min = ',F7.3,'          *')
12300 FORMAT(' *    Q^2 _min = ',F8.3,'                            *')
12400 FORMAT(' *    Q^2 _max = ',F8.3,'                            *')
12500 FORMAT(' *    y _min = ',F7.3,'                               *')
12600 FORMAT(' *    y _max = ',F7.3,'                               *')
12700 FORMAT('  cm energy   ',F10.3,' GeV')
17100 FORMAT(' *  cut on p_t = ',f10.2,' for IPRO = ',I5,
     +       '         *')
11888 FORMAT(' *  cut on pt applied to process 10 and 11           *')
      END
      SUBROUTINE PARTI(IP1,YX,YVAL,WEIGHT,IFL,IST)
      IMPLICIT NONE
      Integer IP1,IFL,IST
      Double Precision YX,YVAL,WEIGHT
#include "caslujets.inc"
#include "casluco.inc"
#include "caspara1.inc"
#include "caspara.inc"
#include "caspart.inc"
#include "casgki.inc"
#include "casvalues.inc"
#include "caspara5.inc"
	Double Precision xg1,xg1_bar,yx1,yx1_bar 
	Double Precision xg2,xg2_bar,yx2,yx2_bar 
	Common /resolv/xg1,xg1_bar,yx1,yx1_bar,xg2,xg2_bar,yx2,yx2_bar
      Double Precision PHI,DBGAM
      COMMON/DIFFA/ PHI
      DOUBLE PRECISION ME
      DIMENSION DBGAM(4)
      Double Precision ALPH_EM,FGAM,X11,COST2,SPHE,CPHE
      Double Precision PEP,PEZ,PDLE,XPP,BOCHCK,PHIE,XGP
      Double Precision SPHI,STHETA,PT,DOT1
      Double Precision yg,yg_bar
      Integer NIPH
C IST = 0, ONLY PARTICLE MOMENTA, BUT SCALE Q2Q IN STRUCTURE FUNCTION
C NOT YET DEFINED.
C IST = 1, CALCULATE ONLY STRUCTURE FUNCTION SINCE Q2Q NOW DEFINED
C
C IF PARTICLE IP1 = 11 (ELECTRON) USE EQUIVALENT PHOTON APPROXIMATION
C
c      write(6,*) 'parti  Q2 = ',Q2,' yx = ',yx
c      write(6,*) 'parti IRES(1) = ',IRES(1)
      YVAL =  -999999.D0
      WEIGHT = -999999.D0
      ME = PYMASS(IP1)
      ALPH_EM = ALPH
	if(IFL.EQ.1) Then
	  XGP=XG1
	  ELSE
	  XGP=XG2
	  ENDIF
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)

C IF 1ST INCOMING PARTICLE IS A ELECTRON AND LOOKING FOR DIRECT PHOTON
      IF((IABS(IP1).EQ.11.OR.IABS(IP1).EQ.13)) THEN
C
C...          WEIZAECKER WILLIAMS APPROXIMATION FOR BREMSTRAHLUNG (EPA)
C...          WITH Q2 DEPENDENCE
         FGAM = (1.D0 - YX +YX*YX/2.D0)/Q2/YX - ME*ME*YX/Q2/Q2

         YVAL = ALPH_EM * FGAM / PI
c         write(6,*) 'parti ',q2min,q2max
         WEIGHT = YX*Q2*DLOG(YMAX/YMIN)*DLOG(Q2MAX/Q2MIN)
C...

C...      CALCULATE COS OF ELECTRON VIA FORMULA GIVEN IN CAMPBELL ET AL
C...      CAN.J.PHYS. 39 1981 (1742)
C...      VIA X Q2 Y
C...      X11= Q2/(2P.Q)
         X11= Q2/YX/SSS
         PHIE = PHI
         IF(PHIE.GT.7.D0) write(6,*) 'fatal PHIE = ',PHIE
         SPHE=DSIN(PHIE)
         CPHE=DCOS(PHIE)
C...      PEL1 = 4 VECTOR OF SCATTERED ELECTRON IN EP CMS SYSTEM
         N=N+1
         K(N,1)=1
         IF(INTER.LT.2) THEN
            K(N,2)=K(1,2)
         ELSEIF(INTER.EQ.2) THEN
            K(N,2)=ISIGN(1,K(1,2))*12
         ENDIF
         K(N,3)=IFL
c. 2*P.l_e

         P(N,5)=PYMASS(K(N,2))
         YG = YX
         YG_BAR = -Q2/(1.d0 - YG)/SSS
         YG_BAR = 0 ! hannes hack

         P(N,4) =  - (yg*P(1,4) + yg_bar*P(2,4)) + P(1,4)
         P(N,3) =  - (yg*P(1,3) + yg_bar*P(2,3)) + P(1,3)
         P(N,1) = - DSQRT(Q2)*SIN(PHIE)          + P(1,1)
         P(N,2) = - DSQRT(Q2)*COS(PHIE)          + P(1,2)

C...      PPH = 4 VECTOR OF (VIRTUAL) PHOTON IN EP CMS SYSTEM
         N=N+1
         K(N,1)=21
         K(N,2)=KEB
         K(N,3)=IFL
         P(N,1)= P(1,1) - P(3,1)
         P(N,2)= P(1,2) - P(3,2)
         P(N,3)= P(1,3) - P(3,3)
         P(N,4)= P(1,4) - P(3,4)
         P(N,5)= -SQRT(ABS(DOT1(N,N)))
c            write(6,*) 'parti: ',P(N,5)**2,q2
         NIPH = N
         IF(IFL.EQ.1) THEN
            NIA1 = N
         ELSE
	      NIA2 = N
	   ENDIF
         IF(IRES(IFL).EQ.0) THEN
	      
	      IF(IFL.EQ.1) THEN
               NIR1 = -99999
		   XGP = yx
		ELSE
		   NIR2 = -99999
		   XGP = yx
		ENDIF
         ELSEIF(IRES(IFL).EQ.1) THEN
	      IF(IFL.EQ.1) THEN
		   XGP = xg1
		ELSE
		   XGP = xg2
		ENDIF
C         XGP = X_GLUON OF THE PHOTON = E_GLUON/E_ELECTRON
C         YX  = E_PHOTON/E_ELECTRON
C         XPP = E_GLUON/E_PHOTON
            XPP=XGP/YX
C...      PPH = NOW 4 VECTOR OF HADR COMP OF PHOTON IN EP CMS SYSTEM
            N = N + 1
            K(N+1,1)=21
            K(N+1,2)=21
            K(N+1,3)=N-1
            P(N+1,1)= P(N-1,1)*XPP
            P(N+1,2)= P(N-1,2)*XPP
            P(N+1,3)= P(N-1,3)*XPP
            P(N+1,4)= DABS(P(N-1,4))*XPP
            P(N+1,5)= PYMASS(21)
		IF(IFL.EQ.1) THEN
               NIA1 = N+1
		ELSE
		   NIA2 = N+1
		ENDIF
C...      PHA = NOW 4 VECTOR OF PHOTON REMNANT  IN EP CMS SYSTEM
            K(N,1)=1
            K(N,2)=21
            K(N,3)=N-1
            P(N,1)= P(1,1) - P(N+1,1)- P(3,1)
            P(N,2)= P(1,2) - P(N+1,2)- P(3,2)
            P(N,3)= P(1,3) - P(N+1,3)- P(3,3)
            P(N,5)= PYMASS(21)
            P(N,4)= DSQRT(P(N,1)**2+P(N,2)**2+P(N,3)**2+P(N,5)**2)
            P(N+1,4)= P(1,4) - P(3,4) - P(N,4)
		IF(IFL.EQ.1) THEN
               NIR1 = N 
		ELSE
		   NIR2 = N 
		ENDIF

c            write(6,*) ' parti yx,xgp,xpp',yx,xgp,xpp
c            call PYLIST(1)
C NOW BOOST TO GAMMA PROTON FRAME
            DBGAM(1) = P(NIPH,1) + P(2,1)
            DBGAM(2) = P(NIPH,2) + P(2,2)
            DBGAM(3) = P(NIPH,3) + P(2,3)
            DBGAM(4) = P(NIPH,4) + P(2,4)
            BOCHCK = (DBGAM(1)/DBGAM(4))**2 + (DBGAM(2)/DBGAM(4))**2 +
     +      (DBGAM(3)/DBGAM(4))**2
            BOCHCK = DSQRT(BOCHCK)
            IF(BOCHCK.GT.0.99999999D0) goto 10
            CALL PYROBO(0,N,0.D0,0.D0,-DBGAM(1)/DBGAM(4),-DBGAM(2)/
     +      DBGAM(4), -DBGAM(3)/DBGAM(4))
            SPHI = PYANGL(P(NIPH,1),P(NIPH,2))
            call PYROBO(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
            STHETA = PYANGL(P(NIPH,3),P(NIPH,1))
            call PYROBO(0,0,-STHETA,0.D0,0.d0,0.d0,0.d0)
            P(N+1,1) = 0.D0
            P(N+1,2) = 0.D0
            P(N+1,3) = XPP*P(NIPH,3)
            P(N+1,4) = XPP*(P(NIPH,4)*P(2,4)-P(NIPH,3)*P(2,3))/
     +                 (P(2,4)-P(2,3))
            P(N+1,3) = P(N+1,4)
            P(N,1)= P(1,1) - P(N+1,1)- P(3,1)
            P(N,2)= P(1,2) - P(N+1,2)- P(3,2)
            P(N,3)= P(1,3) - P(N+1,3)- P(3,3)
            P(N,4)= P(1,4) - P(N+1,4)- P(3,4)
            P(N,5)= DSQRT(DABS(P(N,1)**2+P(N,2)**2+P(N,3)**2-P(N,4)**2))
            N = N + 1
c            call PYLIST(1)
            call PYROBO(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
            call PYROBO(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
            CALL PYROBO(0,N,0.D0,0.D0,DBGAM(1)/DBGAM(4),DBGAM(2)/
     +      DBGAM(4), DBGAM(3)/DBGAM(4))
c            call PYLIST(1)
         ENDIF
C IF 1ST INCOMING PARTICLE IS A PHOTON
      ELSEIF(IABS(IP1).EQ.22.AND.IRES(IFL).EQ.0) THEN
C...      4 VECTOR OF (REAL) PHOTON IN CMS SYSTEM
         N=N+1
         K(N,1)=21
         K(N,2)=KPH
         K(N,3)=IFL
         P(N,1)= P(1,1)
         P(N,2)= P(1,2)
         P(N,3)= P(1,3)
         P(N,4)= P(1,4)
         P(N,5)=PYMASS(KPH)
	   IF(IFL.EQ.1) THEN 
            NIA1 = N
            NIR1 = -99999
	   ELSE
            NIA2 = N
            NIR2 = -99999
	   ENDIF
         YVAL = 1.D0
         WEIGHT = 1.D0
c	   write(6,*) ' in parti ',ires(1),(p(Nia1,i),i=1,5)
      ELSEIF(IABS(IP1).EQ.22.AND.IRES(1).EQ.1) THEN
         WEIGHT = 1.0D0
         YVAL = 1.0D0

C         XGp = X_GLUON OF THE PHOTON = E_GLUON/E_ELECTRON
C         YX  = E_PHOTON/E_ELECTRON
C         XPp = E_GLUON/E_PHOTON
         XPP=XGP/YX
C...      PPH = NOW 4 VECTOR OF HADR COMP OF PHOTON IN EP CMS SYSTEM
         N = N + 1
         K(N+1,1)=21
         K(N+1,2)=21
         K(N+1,3)=IFL
         P(N+1,1)= P(N-1,1)*XPP
         P(N+1,2)= P(N-1,2)*XPP
         P(N+1,3)= P(N-1,3)*XPP
         P(N+1,4)= DABS(P(N-1,4))*XPP
         P(N+1,5)= PYMASS(21)
	   IF(IFL.EQ.1) THEN 
            NIA1 = N+1
	   ELSE
            NIA2 = N+1
	   ENDIF
C...      PHA = NOW 4 VECTOR OF PHOTON REMNANT  IN EP CMS SYSTEM
         K(N,1)=1
         K(N,2)=21
         K(N,3)=IFL
         P(N,1)= P(1,1) - P(N+1,1)- P(3,1)
         P(N,2)= P(1,2) - P(N+1,2)- P(3,2)
         P(N,3)= P(1,3) - P(N+1,3)- P(3,3)
         P(N,5)= PYMASS(21)
         P(N,4)= DSQRT(P(N,1)**2+P(N,2)**2+P(N,3)**2+P(N,5)**2)
         P(N+1,4)= P(1,4) - P(3,4) - P(N,4)
	   IF(IFL.EQ.1) THEN 
            NIR1 = N
	   ELSE
            NIR2 = N 
	   ENDIF
         N = N + 1
c            write(6,*) ' parti yx,xgp,xpp',yx,xgp,xpp
c            call PYLIST(1)
      ELSEIF(IABS(IP1).EQ.2212.AND.IRES(IFL).EQ.1) THEN
c no factor of pi here, this is accounted for in the structure function routine
         WEIGHT = YX*Q2*DLOG(YMAX/YMIN)*DLOG(Q2MAX/Q2MIN)
         YVAL = 1.D0/YX
         PHIE = PHI
c	   write(6,*) yx,q2,phi
         IF(PHIE.GT.7.D0) write(6,*) 'fatal PHIE = ',PHIE
         N=N+1
         K(N,1)=1
         K(N,2)=K(1,2)
         K(N,3)=IFL
C...      PPH = 4 VECTOR OF gluon IN PP CMS SYSTEM
         N=N+1
         K(N,1)=21
         K(N,2)=KEB
         K(N,3)=IFL
         YG = YX
         YG_BAR = -Q2/(1.d0 - YG)/SSS

         P(N,4) =  yg*P(1,4) + yg_bar*P(2,4)
         P(N,3) =  yg*P(1,3) + yg_bar*P(2,3)
         PHIE = PHI
         P(N,1) = DSQRT(Q2)*SIN(PHIE)
         P(N,2) = DSQRT(Q2)*COS(PHIE)
         P(N,5)= -SQRT(ABS(DOT1(N,N)))
c            write(6,*) 'parti: ',P(N,5)**2,q2
c         write(6,*) ' parti ',KEB,K(1,2),IFL
	   IF(IFL.EQ.1) THEN 
            NIA1 = N
	   ELSE
            NIA2 = N
	   ENDIF
         P(N-1,1)= P(1,1) - P(N,1)
         P(N-1,2)= P(1,2) - P(N,2)
         P(N-1,3)= P(1,3) - P(N,3)
         P(N-1,4)= P(1,4) - P(N,4)
         P(N-1,5)= -SQRT(ABS(DOT1(N-1,N-1)))
         NIPH = N
	   IF(IFL.EQ.1) THEN 
         NIR1 = -99999
	   ELSE
         NIR2 = -99999
	   ENDIF  

      ENDIF
      RETURN
   10 write(6,*) ' PARTI boost error '
      RETURN
      END

      SUBROUTINE p_semih(X,WPART)
      IMPLICIT None
#include "casgki.inc"
#include "caslujets.inc"
#include "casluco.inc"
#include "caspart.inc"
#include "caspara1.inc"
#include "caspara.inc"
#include "casvalues.inc"
#include "caspara5.inc"
#include "casshwr.inc"
#include "casmult.inc"
      COMMON/DIFFA/ PHI
      DOUBLE PRECISION ME,MP
      Double Precision XY
      COMMON/XVAR/ XY(10)
      DOUBLE PRECISION STHETA,SPHI
      DIMENSION X(20)
      Integer ndimc
      COMMON /DIMEN/ NDIMC
      Integer IMIX
      COMMON /OALPHAS/ IMIX
      Integer IDEBUG
      COMMON/ INTERN/IDEBUG
      Integer NDIM
      COMMON/INTDIM/ NDIM
      Integer IGENFL,IMEPHI
      COMMON/GENWEI/IGENFL,IMEPHI
      DOUBLE PRECISION QG21,QG22
      COMMON/SEMIH/ QG21,QG22
      REAL SNGL
      Double Precision dcasrn
      EXTERNAL dcasrn
	Integer  ikincut,Ipgg,ns_sel
	Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Double Precision QF,W12
      DATA QF/0.0D0/,W12/0.0D0/

      Integer NDIMEN
      Integer I,J,IN,NP,NPP,KPA2,NRN,IST,NPFIN,KI,KPZ
      Integer nb2,kpf,kstart,kpao,kpao2,nacc
      Double Precision X,WPART,W02,QFT,YX,QT2
      Double Precision QG2MAX,QG2MIN,WQG22,SMIN
      Double Precision SHAT1,PHI,FGAM,FWEI,THETE,PR,CPHP,SPHP

      Double Precision STHP,CTHP
      Double precision xsi_hard,phi_k2,phi_k1,BOCHCK,ECM,WT,PT2,xp
      Double Precision wtglu,xpint,flux,wmat,XZOHRN
      Double Precision DOT,DOT1
	Double Precision XP1MIN,xp1,fwei1,wqg21
	Double Precision xg1,xg1_bar,yx1,yx1_bar 
	Double Precision xg2,xg2_bar,yx2,yx2_bar 
	Common /resolv/xg1,xg1_bar,yx1,yx1_bar,xg2,xg2_bar,yx2,yx2_bar

	double precision xsim,ykt
	common/hzbval/xsim,ykt
	Double Precision sumkt,k1x,k1y,k2x,k2y,xmin1,ratio 

      double precision cut
	logical firstt
      double precision au 
      common/f2fit/au(50),firstt

      real qg22t,shat1t
      common/mytest/qg22t,shat1t
c      write(6,*) ' p_semih ',ipro

c      Idebug = 1
      cut=pt2cut(ipro)
      cut = 0.25
cc      cut = 0.001
cc      cut = 1.
      NDIMEN = NDIM
      Do  IN=1,200
	 Do J=1,5
         K(IN,J) = 0
         P(IN,J) = 0
	 Enddo
      Enddo 
C...  GIVE BEAM  FOUR VECTORS
      DO 20 I=1,2
         DO 20 J=1,5
            K(I,J)=KBEAM(I,J)
   20 P(I,J) = PBEAM(I,J)
      If(Mult.ne.0) then
C... CALCULATE new CMS ENERGY
         DO 21 I=1,4
            CM(I) =P(1,I)+P(2,I)
   21    CONTINUE
      SSS = DOT(CM,CM)
      Endif
   
      ME =P(1,5)
      MP =P(2,5)
      NP=2
      NPP=3*NP
      N=2
      KT2_1=-99999d0
      KT2_2=-99999d0
      QF = 0.D0
      XG1SF=-9999d0
	XG2SF=-9999d0
      if(ihfla.eq.3) then
         kpsi = 333
         elseif(ihfla.eq.4) then
         kpsi = 443
         elseif(ihfla.eq.5) then
         kpsi = 553
         else
         kpsi = 443
         endif
      IF(IPRO.EQ.10) THEN
c         AM(1) = 0.140d0
         AM(1) = 0.250d0
         AM(1) = 0.35d0
         AM(1) = max(AM(1),0.0001d0)
         AM(1) = au(8)
         AM(2) = AM(1)
      ELSEIF(IPRO.EQ.11) THEN
         IF(IWEI.EQ.0) THEN
c            KPA = 4
	      KPA=IHFLA
            IF(IHFLA.GE.4) KPA=IHFLA
         ELSEIF(IWEI.EQ.1) THEN
         ENDIF
         AM(1) = PYMASS(KPA)
         if(IHFLA.EQ.4) then
           AM(1) = au(9)
         endif
ccc         AM(2) = PYMASS(KPA)
         AM(2) = AM(1)
      ELSEIF(IPRO.EQ.15) THEN
         AM(1) = 0.250d0
         AM(1) = 0.0d0
c         AM(1) = max(AM(1),0.0001d0)
c         AM(2) = 2d0
      ELSEIF(IPRO.EQ.2) THEN
         AM(1) = PYMASS(Kpsi)
         AM(2) = 0.d0
      ELSEIF(IPRO.EQ.3) THEN
         KPA=Kchi
         AM(1) = PYMASS(KPA)
         AM(1) = au(8)
         AM(2) = 0.d0
      ELSEIF(IPRO.EQ.4) THEN
         KPA=2
         AM(1) = PYMASS(KPA)
         AM(2) = 0.d0
      ELSEIF(IPRO.EQ.102) THEN
	   KPA = 25
         AM(1) = PYMASS(KPA)
         AM(2) = 0.d0
      ELSEIF(IPRO.EQ.500) THEN
         KPZ = 23
         AM(1) = PYMASS(KPA)
      ELSEIF(IPRO.EQ.501) THEN
         KPZ = 23
         KPA = KPZ
         KPA2= 1
         AM(1) = PYMASS(KPA)
         AM(2) = 0.d0
      ELSEIF(IPRO.EQ.503) THEN
         NP=3
         NPP=3*NP
	   KPA = 1
         KPZ = 23
         AM(1) = PYMASS(KPA)
         AM(2) = PYMASS(KPA)
         AM(3) = PYMASS(KPZ)
      ELSEIF(IPRO.EQ.504) THEN
         NP=3
         NPP=3*NP
	   KPA = IHFLA
         KPZ = 23
         AM(1) = PYMASS(KPA)
         AM(2) = PYMASS(KPA)
         AM(3) = PYMASS(KPZ)
      ELSEIF(IPRO.EQ.510) THEN
      ELSEIF(IPRO.EQ.513) THEN
         NP=3
         NPP=3*NP
c        select charge of W
         KPZ = 24
	   KPA = -2
	   KPA2 = 1
         if(dcasrn().ge.0.5) then 
            KPZ=-24         
	      KPA = 2
            KPA2 = -1
         Endif   
         AM(1) = PYMASS(KPA)
         AM(2) = PYMASS(KPA2)
         AM(3) = PYMASS(KPZ)
      ELSEIF(IPRO.EQ.514) THEN
         NP=3
         NPP=3*NP
c        select charge of W
         KPZ = 24
         if(IHFLA.eq.5) then
            write(6,*) ' p_semih: cabbibo supressed combination selected'
            write(6,*) ' p_semih: heaviest flavor = ',IHFLA
            write(6,*) ' p_semih: not implemented - Program stopped '
            stop
            endif
	   KPA = -IHFLA
         KPA2 = IHFLA-1
         if(dcasrn().ge.0.5) then 
            KPZ=-24         
	      KPA = IHFLA
            KPA2 = -(IHFLA-1)
         Endif   
c         write(6,*) ' KPZ, KPA KPA2 ',KPZ,KPA,KPA2
         AM(1) = PYMASS(KPA)
         AM(2) = PYMASS(KPA2)
         AM(3) = PYMASS(KPZ)
      ELSE
         WRITE(6,*) 'wrong subprocess selected: IPRO = ',IPRO
         WRITE(6,*) '**** PROGRAM STOP ****'
         STOP
      ENDIF
      
      
c      write(6,*) ' p_smeih mass, ipro ',AM(1),ipro
      
	Fwei1 = 1d0
      
C.. HERE THE LIMITS ON Y( PHOTON ENERGY) AND Q**2 ARE CALCULATED
C... ACCORDING TO PROCEEDINGS OF HERA WORKSHOP 1987
C... ALI ET AL
      W02=(AM(1)+AM(2)+MP)**2
c      IF(AM(1).LT.1.0D0) W02=(1.D0 + MP)**2
      IF(AM(1).LT.1.0D0.and.pt2cut(ipro).ne.0.) then
           W02=(DSQRT(PT2CUT(IPRO)) + MP)**2
           else
           W02=(2*DSQRT(cut) + MP)**2
      Endif
      W12=W02-MP*MP
      YMAX=SSS+W12+DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMAX=YMAX/(2.D0*(SSS+ME*ME))
      YMIN=SSS+W12-DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMIN=YMIN/(2.D0*(SSS+ME*ME))
      IF(YMI.GT.YMIN) YMIN=YMI
      IF(YMA.LT.YMAX) YMAX=YMA
c         WRITE(6,10500) YMIN,YMAX
c10500 FORMAT(' limits on y ',/,
c     +' YMIN = ',E10.5,' YMAX = ',E10.5)

C ... select particle code for light flavour production according
C ... to charges
      IF(IPRO.EQ.10) THEN
         IF(INTER.LT.2) THEN
            QF=DFLOAT(PYCHGE(1))**2 + DFLOAT(PYCHGE(2))**2 +
     +      DFLOAT(PYCHGE(3))**2
            IF(IABS(K(1,2)).eq.2212.and.IABS(K(2,2)).eq.2212) QF=3d0
            QF = 2.D0*QF
            KPA = -4
            QFT = - DBLE(dcasrn())*QF
   30       KPA=KPA+1
            IF(IABS(K(1,2)).eq.2212.and.IABS(K(2,2)).eq.2212) Then
               QFT = QFT + 1d0
            else
               QFT = QFT + DFLOAT(PYCHGE(KPA))**2
            endif
            IF(QFT.LT.0.0D0.OR.KPA.EQ.0) GOTO 30
            IF(KPA.GT.3) write(6,*) 'fatal light quark = charm!!!!!! ',
     +      KPA
            IF(KPA.EQ.0) write(6,*) 'fatal no quark !!!!!! ',
     +      KPA,QFT
            KPA2 = - KPA
         ENDIF
      ELSEIF(IPRO.EQ.11) THEN
         IF(IWEI.EQ.0) THEN

         ELSEIF(IWEI.EQ.1) THEN
            IF(dcasrn().GT.0.5) KPA = -KPA
         ENDIF
         KPA2=-KPA
c         AM(1)=PYMASS(KPA)
c         AM(2)=AM(1)
      ELSEIF(IPRO.eq.501) THEN
      ELSEIF(IPRO.ge.503.and.IPRO.le.505) THEN
         IF(IWEI.EQ.0) THEN

         ELSEIF(IWEI.EQ.1) THEN
            IF(dcasrn().GT.0.5) KPA = -KPA
         ENDIF
         KPA2=-KPA
         AM(2)=PYMASS(KPA2)
         AM(1)=PYMASS(KPA)
      ELSEIF(IPRO.ge.513.and.IPRO.le.515) THEN
         AM(2)=PYMASS(KPA2)
         AM(1)=PYMASS(KPA)
      ENDIF

C... YX IS THE PHOTON ENERGY
C... Q2 FROM PHOTON
C... xg2 IS XGLUON ( MOMENTUM FRACTION OF THE GLUON)
C... XMIN = MIN XGLUON TO PRODUCE THE INV. MASS OF GAMMA GLUON SYSTEM
C... XMAX=1.
CCC
C... GENERATE YX,Q2,xg2 ACCORDING TO 1/X SPECTRUM
C... FGAM IS THE WEIGHT OF EPA
      NRN = 0
      XMAX=0.999d0

      IF(Iabs(KE).eq.11.or.KE.eq.22) THEN
c.......
         IF(KE.NE.22) THEN
            IF(YMA.EQ.YMI) then
              YX = YMIN
              Q2 = QMI
              Q2min=q2
              q2max=q2
c              write(6,*) ' here yx,q2',yx,q2
            else
              NRN = NRN + 1

              YX = YMIN*((YMAX/YMIN)**X(NRN))
              IF(YX.GT.YMAX.OR.YX.LT.YMIN) GOTO 90
              Q2MIN=ME*ME*YX*YX/(1.D0-YX)
              IF(QMI.GT.Q2MIN) Q2MIN = QMI
              Q2MAX=YX*SSS - W12
              IF(QMA.LT.Q2MAX) Q2MAX = QMA
              IF(Q2MAX.LT.Q2MIN) GOTO 100
              NRN = NRN + 1
              Q2 = Q2MIN*((Q2MAX/Q2MIN)**X(NRN))
            endif
            XMIN = 0.0D0
            QG22 = 0.D0
            WQG22 = 1.D0
            IF(ISEMIH.EQ.1) THEN
               QT2 = Q2*(1.d0-yx)
               QG2MAX = SSS*YX - Q2 - (AM(1) + AM(2))**2
c               write(6,*) ' p_semih 1 Qg2max '
               QG2MIN = 0.0001d0
ccc               QG2MIN = 0.000001d0
               QG2MAX=DMAX1(QG2MIN,QG2MAX)
               NRN = NRN + 1     
               XZOHRN = X(NRN)  
c 1/qt**2
               QG22 = QG2MIN*((QG2MAX/QG2MIN)**XZOHRN)
               WQG22 = QG22*DLOG(QG2MAX/QG2MIN)
c               write(6,*) 'p_semih ',WQG22
c... end semihard 
            ENDIF

            SMIN = DMAX1(4*cut,2.D0*(AM(1)**2+AM(2)**2))
            XMIN=(SMIN)/(YX*SSS)
            IF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
               XMIN = DMAX1(XMIN,Q2/YX/SSS)
            ENDIF
         ELSEIF(KE.EQ.22) THEN
            Q2 = 0.D0
            YX = 1.D0
            QG22 = 0.D0
            WQG22 = 1.D0
            IF(ISEMIH.EQ.1) THEN
C... here is new for semihard
               QG2MAX = SSS - (AM(1)+AM(2))**2

               QG2MIN = 0.0001d0
               NRN = NRN + 1
               XZOHRN = X(NRN)
c 1/qt**2
               QG22 = QG2MIN*((QG2MAX/QG2MIN)**XZOHRN)
               WQG22 = QG22*DLOG(QG2MAX/QG2MIN)

c... end new semihard
            ENDIF
c               write(6,*) ' p_semih yx,q2 ',yx,q2

            SMIN = DMAX1(4*cut,2.D0*(AM(1)**2+AM(2)**2))
            XMIN= (SMIN)/(YX*SSS)
            

         ELSE
            write(6,*) ' wrong KF selected; program stop '
            YX = 0.D0
         ENDIF
         XG1=YX
	   XG1SF = XG1
	   yx1 = yx
	   IF(XMIN.GE.XMAX) GOTO 110
C ... xg2 = E_parton /E_proton
         IF(IMIX.EQ.0) THEN
           NRN = NRN + 1
           xg2= XMIN*((XMAX/XMIN)**X(NRN))
         ELSEIF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
           xg2= XMIN*((XMAX/XMIN)**X(NDIMEN+1))
         ELSE
           xg2 = 0.D0
         ENDIF
         XG2SF = XG2
         IF(IRES(1).EQ.1) THEN
	      If(Iabs(KE).EQ.11) then
               NRN = NRN + 1

               XP1MIN =DMAX1(4*cut,2.D0*(AM(1)**2+AM(2)**2))
     +           /yx/sss
               xp1 = XP1MIN*(1.D0/XP1MIN)**X(NRN)
               xg1 = xp1*yx
		   XG1SF = xp1  
               QG2MAX = sss*xg1 - Q2 -Qg22 - 4.D0 * AM(1)**2

               QG2MIN = 0.0001d0
               QG2MAX=DMAX1(QG2MIN,QG2MAX)
               NRN = NRN + 1
               QG21 = QG2MIN*((QG2MAX/QG2MIN)**X(NRN))
               WQG21 = QG21*DLOG(QG2MAX/QG2MIN)
	         Fwei1 = WQG21 * DLOG(1.D0/XP1MIN)
               IF(yx1.ne.1d0) THEN
                  xg1_bar = (-qg21)/(1.d0 - xg1)/sss
                  xg1_bar=0d0
                ELSE
                  xg1_bar=0d0
               ENDIF
		ELSEIF(KE.EQ.22) THEN
		
		ENDIF
         ENDIF
	ELSEIF(Iabs(KE).eq.2212) then
c only p-p included
         NRN = NRN + 1
         YX = YMIN*((YMAX/YMIN)**X(NRN))
         XG1=YX
	   XG1SF = XG1
         QG2MIN = 0.0001d0
         Q2MIN = QG2MIN
         Q2MAX=YX*(SSS-MP**2) - W12
                  
         IF(Q2MAX.LT.Q2MIN) GOTO 100
         NRN = NRN + 1
         Q2 = Q2MIN*((Q2MAX/Q2MIN)**X(NRN))
         XMIN = 0.0D0
         XMIN = 0.0D0
         QG22 = 0.D0
         WQG22 = 1.D0
         QT2 = Q2 
         QG2MAX = SSS*YX - Q2 - (AM(1)+AM(2))**2
                  
         QG2MAX=DMAX1(QG2MIN,QG2MAX)
         NRN = NRN + 1
         XZOHRN = X(NRN)
c 1/qt**2
         QG22 = QG2MIN*((QG2MAX/QG2MIN)**XZOHRN)
         WQG22 = QG22*DLOG(QG2MAX/QG2MIN)
c        write(6,*) ' p_semih Q2,QG2: ',Q2,QG22
         IF(IPRO.EQ.3.or.ipro.eq.4.or.ipro.eq.102.or.ipro.eq.500) then
	     Smin = AM(1)**2
           NRN = NRN + 1
           PHI_K1 = 2*PI*X(NRN)
           k1x = DSQRT(Q2)*SIN(PHI_K1)
           k1y = DSQRT(Q2)*COS(PHI_K1)
           NRN = NRN + 1
           PHI_K2 = 2*PI*X(NRN)
           k2x = DSQRT(QG22)*SIN(PHI_K2)
           k2y = DSQRT(QG22)*COS(PHI_K2)
	     
	     sumkt = (k1x+k2x)**2 + (k1y+k2y)**2
	     xmin = (SMIN+sumkt)/(YX*SSS)
	     xg2 = xmin
         Elseif(IPRO.ge.501.and.IPRO.lt.503) then

           SMIN = DMAX1(4*cut,(AM(1)+AM(2))**2)
           
           XMIN=(SMIN+QT2+QG22)/(YX*SSS)
C ... xg2 = E_parton /E_proton
           IF(IMIX.EQ.0) THEN
              NRN = NRN + 1
              xg2= XMIN*((XMAX/XMIN)**X(NRN))
           ELSEIF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
              xg2= XMIN*((XMAX/XMIN)**X(NDIMEN+1))
           ELSE
             xg2 = 0.D0
           ENDIF
         Elseif(IPRO.ge.503.and.ipro.le.515) then

           SMIN = DMAX1(4*cut,(AM(1)+AM(2)+AM(3))**2)
           
           XMIN=(SMIN+QT2+QG22)/(YX*SSS)
C ... xg2 = E_parton /E_proton
           IF(IMIX.EQ.0) THEN
              NRN = NRN + 1
              xg2= XMIN*((XMAX/XMIN)**X(NRN))
           ELSEIF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
              xg2= XMIN*((XMAX/XMIN)**X(NDIMEN+1))
           ELSE
             xg2 = 0.D0
           ENDIF
	   Else

           SMIN = DMAX1(4*cut,(AM(1)+AM(2))**2)
           XMIN=(SMIN+QT2+QG22)/(YX*SSS)
C ... xg2 = E_parton /E_proton
           IF(IMIX.EQ.0) THEN
              NRN = NRN + 1
              xg2= XMIN*((XMAX/XMIN)**X(NRN))
           ELSEIF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
              xg2= XMIN*((XMAX/XMIN)**X(NDIMEN+1))
           ELSE
             xg2 = 0.D0
           ENDIF
	   Endif
	   XG2SF = XG2
      Else
	   Write(6,*) ' IRES(1)=1 for KE = ',KE,' not implemented '
         Stop
	Endif


      IF(SNGL(xg2).GE.0.9999.OR.xg2.LT.XMIN) GOTO 120
      SHAT1 = SSS*XG1*xg2
      IF(SHAT1.LT.(AM(1)+AM(2))**2) GOTO 130

      IST = 0
      YY=YX
      XG=xg2
      PHI = 99999.D0
c	If(ipro.ne.102) then 
c         NRN=NRN+1
c         PHI = 2.D0*PI*X(NRN)
c	   else
c	   phi=phi_k1
c	   endif
	If(ipro.eq.3..or.ipro.eq.4.or.ipro.eq.102.or.ipro.eq.500) then 
	   phi=phi_k1
	   else
         NRN=NRN+1
         PHI = 2.D0*PI*X(NRN)
	   endif
      CALL PARTI(KE,YX,FGAM,FWEI,1,IST)
      IF(YMA.EQ.YMI) then
         FGAM = 1
         FWEI =1 
         Endif
      XP=0.D0
      IRES(2)=1
C here do some gymnastics to make cd .. on theta angle of electron
      IF(THEMA.LT.180.D0.OR.THEMI.GE.0.D0) THEN
C go first to ep LAB system
         CALL PYROBO(3,3,0.D0,0.D0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
         THETE = PYP(3,14)
         CALL PYROBO(3,3,0.D0,0.D0,-CM(1)/CM(4),
     +                -CM(2)/CM(4),-CM(3)/CM(4))
c         CALL PYLIST(1)
         IF(THETE .GT. 180.01D0) THEN
            WRITE(6,*) ' FATAL: theta_electron > 180 deg ',THETE
         ENDIF
         IF(THETE .GT. THEMA) GOTO 140
         IF(THETE .LT. THEMI) GOTO 140
      ENDIF
C end of these gymnastics

C FINAL STATE PROTON

      NPFIN=NIA1+2+NP
      N=NPFIN
      K(NPFIN,1)=1
      K(NPFIN,2)=KP
      K(NPFIN,3)=2

      PR = P(2,3)*(1.D0 - xg2)
      P(NPFIN,4) = PR
      CPHP=1.D0
      SPHP=DSQRT(1.D0 - CPHP**2)
      IF(ISEMIH.EQ.1) THEN
        STHP = DSQRT(QG22)/PR  
      ENDIF
      IF(DABS(STHP).GT.1.D0) goto 190 ! hannes hack
      CTHP=DSQRT(1.D0 - STHP**2)
      P(NPFIN,1)= PR*STHP*CPHP
      P(NPFIN,2)= PR*STHP*SPHP
      P(NPFIN,3)= PR*CTHP
c     write(6,*) ' lin 372 '
c      CALL PYLIST(1)
C MOMENTA OF PARTON
      K(NIA1+1,1)=21
      K(NIA1+1,2)=KGL
      IF(IPRO.EQ.15) K(NIA1+1,2)=2
      K(NIA1+1,3)=2
      DO 40  KI=1,4
         P(NIA1+1,KI)=P(2,KI)-P(NPFIN,KI)
   40 CONTINUE
      IF(ISEMIH.EQ.0) THEN
         P(NIA1+1,4)=ABS(P(NIA1+1,3))
         P(NIA1+1,5)=0.0D0
      ELSE
         P(NIA1+1,4)=DSQRT(P(NIA1+1,1)**2+P(NIA1+1,2)**2+P(NIA1+1,3)**2
     +               - QG22)
         P(NIA1+1,5)=-sqrt(ABS(DOT1(NIA1+1,NIA1+1)))
c here define new momenta: introduce maximum virtuality mu**2 = Q2 + 4M**2
c for the matrix element and kinematics we set Virt_max=kt**2
c     define lightcone variable xg=xg2
c         XG = xg2
c     define lightcone variable xg_bar=(-kt**2)/(1-xg)/SSS
         XG2_bar = (-QG22)/(1.d0-xg2)/SSS
c here check on virt. gamma
         QT2 = PYP(NIA1,9)
         IF(XG1.NE.1D0) THEN
            XG1_BAR = (-QT2)/(1.d0 - XG1)/SSS
         ELSE
            XG1_BAR=0D0
         ENDIF
         
         xg2_bar = 0. ! hannes hack
         xg1_bar = 0. ! hannes hack
         
         
         KT2_1 = ABS(QT2)
                  
         xsi_hard=(xg1+xg2_bar)/(xg2+xg1_bar)

	   xsim = xsi_hard
         ykt = sqrt(abs(QG22))
         YY=YX
         if(xsi_hard.le.0) goto 200 ! hannes hack
         P(NIA1+1,4) =  xg2_bar*P(1,4) + xg2*P(2,4)
         P(NIA1+1,3) =  xg2_bar*P(1,3) + xg2*P(2,3)
	   if(ipro.eq.3.or.ipro.eq.4.or.ipro.eq.102.or.ipro.eq.500) then
            IF(Iabs(KE).eq.11.or.KE.eq.22) THEN
	       NRN = NRN + 1
             PHI_K2 = 2*PI*X(NRN)
	      endif
            else 
            NRN = NRN + 1
            PHI_K2 = 2*PI*X(NRN)
	   Endif
         P(NIA1+1,1) = DSQRT(QG22)*SIN(PHI_K2)
         P(NIA1+1,2) = DSQRT(QG22)*COS(PHI_K2)
         P(NIA1+1,5)=-sqrt(ABS(DOT1(NIA1+1,NIA1+1)))
         KT2_2=ABS(QG22)
c	   write(6,*) KT2_2,QG22
c now we also change the photon: set plus component to zero
         P(NIA1,4) =  xg1*P(1,4) + xg1_bar*P(2,4)
         P(NIA1,3) =  xg1*P(1,3) + xg1_bar*P(2,3)
         P(NIA1,5)=-sqrt(ABS(DOT1(NIA1,NIA1)))
         IF(IRES(1).EQ.1.AND.IABS(KE).NE.2212) Then
            NRN = NRN + 1
            PHI_K1 = 2*PI*X(NRN)
            P(NIA1,1) = DSQRT(QG21)*SIN(PHI_K1)
            P(NIA1,2) = DSQRT(QG21)*COS(PHI_K1)
		DO I=1,4
		P(NIA1-1,I) = P(4,I) - P(NIA1,I)
		ENDDO
		YY = XG1/YX
	   ENDIF
	   
         do i=1,5
            p(n+1,i) = P(NIA1,i)
            p(n+2,i) = P(NIA1+1,i)
         enddo
         P(N+2,3) =  0.5d0*(P(N+2,3)+P(N+2,4))
         P(N+2,4) =  DABS(P(N+2,3))
         P(N+1,3) =  0.5d0*(P(N+1,3)-P(N+1,4))
         P(N+1,4) =  DABS(P(N+1,3))
         K(N+1,1)=13
         K(N+2,1)=13
         N=N+2
c        write(6,*) ' p_semih xg: ',xg

         DO 50   KI=1,4
            P(NPFIN,KI)=P(2,KI)-P(NIA1+1,KI)
   50    CONTINUE
         P(NPFIN,5) = -sqrt(abs(dot1(npfin,npfin)))
      ENDIF
      NIA2 = NIA1+1
c      call PYLIST(1)
      NF1=NIA1+2
      NF2=NIA1+3
      if(IPRO.ge.503.and.IPRO.le.515) NF2=NIA1+4
      KPF=KPA
      K(NF1,1)=2
      K(NF1,2)=KPF
      K(NF1,3)=NIA2
      P(NF1,5)=AM(1)
      K(NF2,1)=1
      K(NF2,2)=KPA2
      K(NF2,3)=NIA2
      P(NF2,5)=AM(2)
      if(IPRO.ge.503.and.IPRO.le.515) then
        K(NF1+1,1)=1
        K(NF1+1,2)=KPZ
        K(NF1+1,3)=NIA2
        P(NF1+1,5)=AM(3)
        elseif(ipro.eq.15) then
        K(NF1,2)=21
        K(NF2,2)=2 
        elseif(ipro.eq.501) then
        K(NF1,2)=KPA
        K(NF2,2)=KPA2 
        
      endif   
      if(kpa.eq.0) then
         write(6,*) 'p-semih kpa ',kpa
      endif
      NB2 = NIA2
C...   VECTOR OF GAMMA GLUON CMS SYSTEM

      DBCMS(1)=  P(NIA1,1) + P(NB2,1)
      DBCMS(2)=  P(NIA1,2) + P(NB2,2)
      DBCMS(3)=  P(NIA1,3) + P(NB2,3)
      DBCMS(4)=  P(NIA1,4) + P(NB2,4)

      DO 60  I=1,4
         P(NF1,I)=0.0D0
         P(NF2,I)=0.0D0
         if(IPRO.ge.503.and.IPRO.le.515) P(NF1+1,I)=0.0D0
   60 CONTINUE
      SHAT=DOT(DBCMS,DBCMS)
      
      
      IF(SHAT.LE.0.0) THEN
c         call pylist(1)
         GOTO 160
      ENDIF
      if(xg2.lt.Q2/(YX*SSS)) then
cc        write(6,*) ' xg2, xmin ',xg2,q2/(yx*sss),sqrt(shat)
c        write(6,*) ' q2 ',q2,dbcms(1)**2+DBCMS(2)**2,dbcms(3)**2-DBCMS(4)**2
cc        write(6,*) PYP(NIA1,9),PYP(NIA2,9),DOT1(NIA1,NIA1),DOT1(NIA2,NIA2)
      endif

      IF(IPRO.EQ.10) THEN
         IF(INTER.LT.2) THEN
            KSTART = -4
            QF=DFLOAT(PYCHGE(1))**2 + DFLOAT(PYCHGE(2))**2 +
     +      DFLOAT(PYCHGE(3))**2
            IF(IABS(K(1,2)).eq.2212.and.IABS(K(2,2)).eq.2212) QF=3d0
            QF = 2.D0*QF
            KPA = KSTART
            QFT = - DBLE(dcasrn())*QF
   70       KPA=KPA+1
            IF(IABS(K(1,2)).eq.2212.and.IABS(K(2,2)).eq.2212) then
               QFT = QFT + 1.d0
            else
               QFT = QFT + DFLOAT(PYCHGE(KPA))**2
            endif
            IF(QFT.LT.0.0D0.OR.KPA.EQ.0) GOTO 70
            IF(KPA.GT.IABS(KSTART)-1) THEN
               write(6,*) 'fatal light quark = charm!!!!!! ',KPA
            ENDIF
            IF(KPA.EQ.0) write(6,*) 'fatal no quark !!!!!! ',
     +      KPA,QFT
            KPA2 = - KPA
            IF(IGENFL.EQ.0) THEN
               KPAO = KPA
               KPAO2 = KPA2
            ELSE
               KPA=KPAO
               KPA2 = KPAO2
            ENDIF
            K(NF1,2)=KPA
            K(NF2,2)=KPA2
         ENDIF
      ENDIF

      IF(IPRO.EQ.2) THEN
         K(NF1,2)=KPSI
         K(NF2,2)=21
      ENDIF
c      call PYLIST(1)
C NOW BOOST TO GAMMA GLUON
      BOCHCK = (DBCMS(1)/DBCMS(4))**2 + (DBCMS(2)/DBCMS(4))**2 +
     +(DBCMS(3)/DBCMS(4))**2
      BOCHCK = DSQRT(BOCHCK)
      IF(BOCHCK.GT.0.99999999D0) goto 150
      CALL PYROBO(0,N,0.D0,0.D0,-DBCMS(1)/DBCMS(4),
     +-DBCMS(2)/DBCMS(4),-DBCMS(3)/DBCMS(4))
      SPHI = PYANGL(P(NIA1,1),P(NIA1,2))
      call PYROBO(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
      STHETA = PYANGL(P(NIA1,3),P(NIA1,1))
      call PYROBO(0,0,-STHETA,0.D0,0.d0,0.d0,0.d0)
c      write(6,*) ' after boost '
c      call PYLIST(1)
C NOW  LOOK THAT WE REALLY HAVE ENOUGH ENERGY IN GAMMA GLUON CMS SYSTEM
C...  ECM = CMS ENERGY OF GAMMA GLUON SYSTEM
c         write(6,*) ECM
      ECM =DSQRT(SHAT)

      If(ipro.eq.3..or.ipro.eq.4.or.ipro.eq.102.or.ipro.eq.500)  then 
      else
c cut only for a 2-> 2 process, if resonance, EMC = AM(1) by construction
c if cutted, only precision problems are encountered
         if(IPRO.ge.503.and.IPRO.le.515) then
            IF(ECM.LE.(AM(1)+AM(2)+AM(3))) GOTO 170
         else
            IF(ECM.LE.(AM(1)+AM(2))) GOTO 170
         endif

	endif 
      IF(IPRO.EQ.3.or.IPRO.eq.4.or.IPRO.EQ.102.or.ipro.eq.500) Then
         K(NF2,2)=0
         K(NF2,1)=-1
         K(NF2,3)=0
         P(NF2,5)=0
         K(NF1,1)=1
         K(NF1,3)=4
        DO I=1,4
            P(NF1,I)=0
            P(NF2,I)=0
         Enddo
         P(NF1,4)=ECM
         P(NF1,5)=AM(1)
c        WT to corrrect for it afterwards	   
	   WT = (2.D0*PI)**NPP
c	   call pylist(1)
	ELSE
c         write(6,*) ' IPRO = 13 ,NRN= ',NRN
         IF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
            XY(1)=X(NDIMEN+3)
            XY(2)=X(NDIMEN+4)
            if(IPRO.ge.503.and.IPRO.le.515) then
              XY(3)=X(NDIMEN+5)
              XY(4)=X(NDIMEN+6)
              XY(5)=X(NDIMEN+7)
           endif
         ELSE
            NRN = NRN + 1
            XY(1) = X(NRN)
            NRN = NRN + 1
            XY(2) = X(NRN)
            if(IPRO.ge.503.and.IPRO.le.515) then
              NRN = NRN + 1
              XY(3) = X(NRN)
              NRN=NRN+1
              XY(4)=X(NRN)
              NRN=NRN+1
              XY(5)=X(NRN)
           endif
         ENDIF
c         write(6,*) am(1),am(2),ipro
c         write(6,*) 'p_semih ',XY(1),XY(2),NDIM,NDIMEN,IMIX,IWEI
         CALL PHASE(NP,ECM,AM,PCM,WT)
         IF(WT.LE.0.D0) GOTO 180
         DO 80  I=1,4
            P(NF1,I)=PCM(I,1)
            P(NF2,I)=PCM(I,2)
            if(IPRO.ge.503.and.IPRO.le.515) then
              P(NF1+1,I)=PCM(I,3)
            endif
   80    CONTINUE 
c         write(6,*) ' p_semih ', KPA,KPA2,K(NF1,2),K(NF2,2)
	Endif

C BOOST BACK TO OVERALL CMS 
      call PYROBO(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
      call PYROBO(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
      CALL PYROBO(0,N,0.D0,0.D0, DBCMS(1)/DBCMS(4),DBCMS(2)/DBCMS(4),
     +DBCMS(3)/DBCMS(4))
cc remove cuts to meoffsh routines


      XG = xg2
      if(xg2.GT.1.d0) write(6,*) 'p_semih xg2>1',xg2,ipro
c      call pylist(1)
C... WTGLU IS THE WEIGHT FOR XGLUON GENERATION
      XPINT=DLOG(XMAX/XMIN)
      WTGLU = XPINT
	if(ipro.eq.3.or.ipro.eq.4.or.ipro.eq.102.or.ipro.eq.500) wtglu=1
      WTGLU = WTGLU*WQG22
      WMAT=WT *WTGLU*FGAM*FWEI * Fwei1
c	write(6,*) 'p_semih: FGAM,FWEI,FWEI1 ',FGAM,FWEI,FWEI1
      WMAT=WMAT/((2.D0*PI)**NPP)
C...  FLUX = FLUX USED FOR X-SECTION CALCULATION
      FLUX=(2.D0*PI)**4
      FLUX=FLUX/(2.D0*xg2*xg1*sss)
      WPART = WMAT *FLUX
	IF(IPRO.EQ.3.or.ipro.eq.4.or.ipro.eq.102.or.ipro.eq.500) then
c        write(6,*) 'p_semih WQG22 ',WQG22,Fgam,Fwei,Fwei1
	  Wpart = WTGLU*FGAM*FWEI * Fwei1
c        IF(Iabs(KE).eq.11.and.IPRO.eq.4) then
c	     Wpart = WTGLU*FWEI
c        endif
        if(Ipro.eq.3) then
          Wpart=Wpart*Flux/((2.D0*PI)**3)
c 
c xsection for resonance production gg-> chi_c
c sig = \int dx1 dx2 g(x1) g(x2) sigma_hat \delta(x1x2 * sss - M2)
c     = \int dx1 dx2 g(x1) g(x2) sigma_hat 1/(x1*sss) \delta(x2 - M2/(x1*sss))
c     = \int dx1 g(x1) g(x2=M2/x1*sss) sigma_hat x2/M2 
c     = \int dx1/x1  x1g(x1) x2g(x2) sigma_hat/M2
          Wpart=Wpart/AM(1)**2
        endif
	Endif
      qg22t = qg22   ! hannes hack
      shat1t =shat1  ! hannes hack

      PT2H = PT2
      SHH = SHAT
      NDIMC = NRN
      CALL PYEDIT(14)
c      write(6,*)  'p_semih: WPART ',WPART
      IF(IGENFL.EQ.0) THEN
         KPAO = KPA
         KPAO2 = KPA2
      ELSE
         KPA = KPAO
         KPA2 = KPAO2
      ENDIF
	
c      write(6,*) ' end of p_semih ',NDIMC,NRN
c      CALL PYLIST(1)
	
	
      if(wpart.ne.wpart) write(6,*) ' p_semih wpart = ',wpart,xsi_hard
c      write(6,*) ' p_semih end ,kpa,kpa1,kpa2,kpao2',kpa,kpa1,kpa2,kpao2
      RETURN
   90 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: ylimit ; RETURN ',yx,ymin,ymax,x(1)
      ENDIF
      NDIMC = 9999
      RETURN
  100 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: Q2limit ; RETURN ',Q2,Q2MIN,Q2MAX
         NDIMC = 9999
      ENDIF
      RETURN
  110 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: xlimit ; RETURN xmin,xmax',xmin,xmax
         write(6,*) ' p_semih: xlimit ; RETURN smin,y,sss',SMIN,YX,SSS
      ENDIF
      NDIMC = 9999
      RETURN
  120 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: xg2imit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  130 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: slimit ; RETURN ',shat1,am(1),am(2),IPRO
         write(6,*) ' sss,xg1,xg2,q2,qg22 ',sss,xg1,xg2,q2,qg22
         write(6,*) ' p_semih q2,y,sss',q2,yx,sss
         write(6,*) ' p_semih 2*pq ',2*dot1(nia1,2),sss*yx
      ENDIF
      NDIMC = 9999
      RETURN
  140 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: theta limit ; RETURN ',themi,thema,thete
      ENDIF
      NDIMC = 9999
      RETURN
  150 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: bochck limit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  160 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: shat ; RETURN ',shat,shat1,kpa
        write(6,*) ' p_semih: QG2MAX,QG22 ',QG2MAX,QG22
        write(6,*) ' p_semih: yx,xg2,q2,s ',yx,xg2,q2,sss
c        write(6,*) ' p_semih 2*pq ',2*dot1(nia1,2),sss*yx
        write(6,*) ' p_semih nia1,nia2 ',nia1,nia2
c        write(6,*) ' p_semih 2*ktgam*ktglu ',2*(P(Nia1,1)*P(nia2,1)+
c     &  P(Nia1,2)*P(NIA2,2))
      ENDIF
      NDIMC = 9999
      RETURN
  170 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: ECM limit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  180 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: PHASE WT=0 ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  191 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: PTCUT limit ; RETURN ',PT2,KPA
      ENDIF
      NDIMC = 9999
      RETURN

  190 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: STHP limit ; RETURN ',DSQRT(QG22),PR,STHP
      ENDIF
      NDIMC = 9999
      RETURN
  200 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: xsi_hard limit ; RETURN ',xsi_hard
      ENDIF
      NDIMC = 9999
      RETURN


      END

      SUBROUTINE meoffsh(WT1)
      IMPLICIT None
C
c
c use the momenta ala CCH
c
c   P1 (Kin1)
c   ---------------------------
c      \
c    K1 \ KIA1     / P4 (KF1)
c        \        /
c         --------
c        /        \
c    K2 / KIA2     \ P3 (KF2)
c      / 
c   --------------------------
c   P2 (Kin2)
c
c

C

#include "caslujets.inc"
#include "caspara1.inc"
#include "caspara.inc"
#include "casluco.inc"
#include "casludat1.inc"
#include "caspara5.inc"
#include "caspart.inc"
#include "casgki.inc"
#include "casshwr.inc"
#include "casmult.inc"
#include "casscal.inc"
      Integer IGENFL,IMEPHI
      COMMON/GENWEI/IGENFL,IMEPHI
      DOUBLE PRECISION K1K2,K1P2,K1P3,K1P4,K2P1,K2P3,K2P4,K12,K22
      DOUBLE PRECISION P1P2,P1P3,P1P4,P2P3,P2P4,P3P4,M2,PT2Q
      Double Precision k1p3t,k12t,qq1
      Double Precision XPQ1(-6:6),XPQ2(-6:6),SFSCAL1,SFSCAL2 
      Double Precision STMAX1,STMAX2
      common/myscale/STMAX1,STMAX2
      Double Precision Xdum,SFSdum,STMdum,XPQdum(-6:6)
      DOUBLE PRECISION DSHN(4),DSH(4),dshnf(4),pfin(4)
      Double Precision DScale
      Double Precision scale
      Double Precision CUT
c      DATA CUT /1.d-10/
c      DATA CUT /1.d-3/
      DATA CUT /1.d-6/ ! as smallx
c      DATA CUT /1.d-8/ ! as smallx
      Double Precision Sfcut
      DATA Sfcut /1.d-3/ 
c      DATA Sfcut /1.d-9/ 
	
      Integer I,KIA1,KIA2,KF1,KF2,Kin1,Kin2
      Double Precision QF2,alpha_s,qg2,qg1,sh,th,uh,wt1,sum2
      Double Precision wtg,wtq1,wtq2
      Double Precision cab,cabB,SH_NEW,pcm_fin,cos_fin,phi,phio,test
      Double Precision sh_old,shf,check,dphi,dtheta
      Double Precision E1,E2,yph,yph_bar,xsi_hard(2),epsilon
      Double Precision k1q,k1q1,k1p,kq,pq,pq1,plus,minus
      Double Precision X1,X2
      Double Precision alph_em,pcm_new,nq
      Double Precision SumAab,SumAnab,SumaL
      Double Precision SumB1ab,SumB1nab,SumB2ab,SumB2nab
      Double Precision SumC1ab,SumC1nab,SumC2ab,SumC2nab
      Double Precision SUMA,SUMB,SUMC,SUMC1,SUMC2,SUMCL
      Double Precision A1,A2,A3,Asum
      double Precision XPQC1,XPQC2,QSUM
	Double Precision rntest,Rleg,xsecnorm
      Integer NFLT
      Integer ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON /CASCOL/ ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      double Precision Aas
      Integer IRPAT,IRPBT,IRPCT,IRPDT,IRPET,IRPFT,IRPGT

      Double Precision dot,dot1,dcasrn
      Double Precision sL,rt,ru
      Double Precision Del,Nc
      Double Precision rap1,rap2,rapp,ptstm1,ptstm2,zz1,zz2,QT2
      Integer npfin
      Integer IQBAR
      Integer  ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
	Double Precision MEQG
      
	Double Precision QQg,lambda
	Double Precision Qscal 
	Integer Iqqbar
	Common/GLUDAT2/Qscal,Iqqbar
c      Data IQBAR/1/
c      Data IQBAR/2/
      Double Precision PT2
      Integer Nacc
      
      Double Precision alphasca

      Integer IFLflag
      Common/FlagFL/ IFLflag

      Data Nc/3.d0/
      Logical ltest
c      data ltest/.true./
      data ltest/.false./
      Logical first
      data first/.true./
c      SAVE STMAX1,STMAX2, 
      Save XPQ2, phio   

c     call pylist(1)
      IQBAR = Iqqbar
      IRESPRO = 0
      IRPAT = IRPA
      IRPBT = IRPB
      IRPCT = IRPC
      IRPDT = IRPD
      IRPET = IRPE
      IRPGT = IRPG
c      IQbar = 2
      Xsecnorm = 1.
      IF(INTER.EQ.2) THEN
         write(6,*) ' INTER = 2 not implemented for meoffsh'
         STOP
      ENDIF
      
c      write(6,*) IPRO,IRPA,IRPB,IRPC
      IF(IPRO.EQ.11) THEN
c	   write(6,*) kpa
         QF2 = DFLOAT(PYCHGE(KPA))**2/9.D0
         NQ = 1d0
         IRPBT = 0
         IRPCT = 0
         IRPDT = 0
         IRPET = 0
      ELSE
         QF2=2.D0/3.D0
         NQ = 3d0
      ENDIF
      NPFIN = N
c	call pylist(1)
      SH = 2.D0 * AM(1)**2 + 2.D0 * DOT1(NF1,NF2)
      WT1 = 0.d0
c      write(6,*) ' meoffsh ',am(1)
c	write(6,*) ' meoffsh: ',IRPAT,IRPBT,IRPCT,IRPDT,IRPET,IRPFT,IRPGT
      IF(SH.LE.0.) RETURN
c boost  cms
      CALL PYROBO(0,N,0.D0,0.D0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
      KIA1=N+1
      KIA2=N+2
      KF1=N+3
      KF2=N+4
      Kin1=N+5
      Kin2=N+6
      N=N+6
      DO I=1,4
         K(KIA1,I) = K(NIA1,I)
         K(KIA2,I) = K(NIA2,I)
         K(KF1,I) = K(NF1,I)
         K(KF2,I) = K(NF2,I)
         P(KIA1,I) = P(NIA1,I)
         P(KIA2,I) = P(NIA2,I)
         P(KF1,I) = P(NF1,I)
         P(KF2,I) = P(NF2,I)
         DSH(I)= P(NIA1,I) + P(NIA2,I)
      ENDDO
      K(KIA1,1)=13
      K(KIA2,1)=13
      K(KF1,1)=13
      K(KF2,1)=13

      K(Kin1,1)=13
      K(Kin2,1)=13
c define incoming particles as massless, to ensure proper cancellation in
c matrix-element
      P(Kin1,4)=Dabs(P(1,3))
      P(Kin1,5)=0
      P(Kin2,4)=Dabs(P(2,3))
      P(Kin2,5)=0
      Do I=1,3
         P(Kin1,I)=P(1,I)
         P(Kin2,I)=P(2,I)
      enddo
c
      SH_old = DOT(DSH,DSH)

      if(ltest) then
         SH_NEW = DOT(DSH,DSH)
         goto 10
      endif

      K12 = DOT1(KIA1,KIA1)
      K22 = DOT1(KIA2,KIA2)
      IF(ICCFM.EQ.0) Then
	   If(ires(1).eq.1) then
            P(KIA1,1) = 0D0
            P(KIA1,2) = 0D0
            P(KIA1,5) = 0
            K12 = 0
         endif
	   If(ires(2).eq.1) then
            P(KIA2,1) = 0D0
            P(KIA2,2) = 0D0
            P(KIA2,5) = 0
            K22 = 0
         endif

      Else
         CALL PYROBO(0,0,0.D0,0.D0,-DSH(1)/DSH(4), -DSH(2)/DSH(4), -DSH(3)
     +   /DSH(4))
         DPHI = PYANGL(P(KIA1,1),P(KIA1,2))
         call PYROBO(0,0,0.D0,-dphi,0.d0,0.d0,0.d0)
         DTHETA = PYANGL(P(KIA1,3),P(KIA1,1))
         call PYROBO(0,0,-DTHETA,0.D0,0.d0,0.d0,0.d0)
c      call pylist(1)
         P(KIA1,3) = dsqrt(((sh_old-K12-K22)**2 - 4d0*K12*K22)/4d0/sh_old)
         P(KIA1,4) = 0.5d0*(sh_old+(-k22+k12))/sqrt(sh_old)
         P(KIA2,3) =-dsqrt(((sh_old-K12-K22)**2 - 4d0*K12*K22)/4d0/sh_old)
         P(KIA2,4) = 0.5d0*(sh_old-(-k22+k12))/sqrt(sh_old) 

c      write(6,*) ' Kia1,2 = ',KIA1,KIA2,' KF1,2 = ',KF1,KF2
c      call pylist(1)
            
         do i=1,4
            dshnf(i) = P(kf1,i)+P(kf2,i)
         enddo
         shf = dot(dshnf,dshnf)

         check = abs(P(kf1,1)+P(kf2,1)+
     &     P(kf1,2)+P(kf2,2) +
     &     P(kf1,3)+P(kf2,3) +
     &     P(kf1,4)+P(kf2,4)-P(KIA2,4)-P(KIA1,4))
         if(check.gt.1e-5) then 
            write(6,*) 'shat ',sqrt(sh_old),sqrt(sh_new),sqrt(shf),check
            write(6,*) ' meoffsh momentum checks ',P(kf1,1)+P(kf2,1),
     &      P(kf1,2)+P(kf2,2), 
     &      P(kf1,3)+P(kf2,3), 
     &      P(kf1,4)+P(kf2,4)-P(KIA2,4)-P(KIA1,4)
         endif
ctest
         PCM_FIN=SQRT(P(NF1,4)**2 - AM(1)**2)
         COS_FIN=P(NF1,3)/PCM_FIN
cend     
         call PYROBO(0,0,DTHETA,0.D0,0.d0,0.d0,0.d0)
         call PYROBO(0,0,0.D0,dphi,0.d0,0.d0,0.d0)

         CALL PYROBO(0,0,0.D0,0.D0,DSH(1)/DSH(4) ,DSH(2)/DSH(4), DSH(3)/
     +   DSH(4))
      Endif
c      call pylist(1)

ctest
      P(KIA1,3) = 0.5d0*(P(KIA1,3)-P(KIA1,4))
      P(KIA1,4) = DABS(P(KIA1,3))
      P(KIA2,3) = 0.5d0*(P(KIA2,3)+P(KIA2,4))
      P(KIA2,4) = DABS(P(KIA2,3))
       DO I=1,4
         DSHN(I)= P(KIA1,I) + P(KIA2,I)
      ENDDO
      SH_NEW = DOT(DSHN,DSHN)
      if(sh_new.lt.0) goto 20
      IF(P(NF1,5)-P(NF2,5).lt.0.01) then
         CALL PYROBO(0,0,0.D0,0.D0,-DSHN(1)/DSHN(4), -DSHN(2)/DSHN(4),
     +   -DSHN(3)/DSHN(4))

         P(KF1,4) = 0.5d0*DSQRT(SH_NEW)
         if(P(KF1,4)**2.le. AM(1)**2) goto 20
         PCM_NEW=SQRT(P(KF1,4)**2 - AM(1)**2)
         P(KF1,3) =PCM_NEW*COS_FIN
         IF(IMEPHI.EQ.0) THEN
            phi=2.d0*PI*dcasrn()
            phio = phi
         ELSEIF(IMEPHI.EQ.1) THEN
            phi=phio
         ENDIF
         P(KF1,1) = DSQRT(PCM_NEW**2-P(KF1,3)**2)*cos(phi)
         P(KF1,2) = DSQRT(PCM_NEW**2-P(KF1,3)**2)*sin(phi) 
         DO I=1,4
            P(KF2,I) = P(KIA1,I) + P(KIA2,I) - P(KF1,I)
         ENDDO
c	Call pylist(1)
         CALL PYROBO(0,0,0.D0,0.D0,DSHN(1)/DSHN(4) ,DSHN(2)/DSHN(4),
     +   DSHN(3)/DSHN(4))
      Endif
     
ctesting
      K22 = DOT1(KIA2,KIA2)
     
	IF(Ires(1).ne.0.and.Ires(2).ne.0.and.IPRO.eq.10.and.
     &	(IRPB.eq.1.or.IRPC.eq.1).and.iccfm.ne.0) then

c         write(6,*) ' enter Rleg '
         do i=1,4
           pfin(i)=P(KF1,i)+P(KF2,i)
         end do
         plus=(pfin(4)-pfin(3))/2D0
         minus=(pfin(4)+pfin(3))/2D0
         P(KIA2,4)=minus
         P(KIA2,3)=minus
         P(KIA1,4)=plus
         P(KIA1,3)=-plus
c  select which leg should be on-shell
         Rleg=dcasrn()
ccc         Rleg = 0.6
cccc         xsecnorm = 2.
         if(Rleg.gt.0.5) then
           P(KIA2,2)=0D0
           P(KIA2,1)=0D0
           P(KIA1,2)=pfin(2)
           P(KIA1,1)=pfin(1)
         else
           P(KIA1,2)=0D0
           P(KIA1,1)=0D0
           P(KIA2,2)=pfin(2)
           P(KIA2,1)=pfin(1)
         endif
         
	else
c      write(6,*) ' should not enter ',IRES(1),IRPBT,IRPCT
	ENDIF
Cend testing

   10 continue

      K1K2 = DOT1(KIA1,KIA2)
      K1P2 = DOT1(KIA1,Kin2)
      K1P3 = DOT1(KIA1,KF2)
      K1P4 = DOT1(KIA1,KF1) 
      K2P1 = DOT1(KIA2,Kin1)
      K2P3 = DOT1(KIA2,KF2)
      K2P4 = DOT1(KIA2,KF1)
      K12 = DOT1(KIA1,KIA1)
      K22 = DOT1(KIA2,KIA2)
      P1P2 = DOT1(Kin1,Kin2)
      P1P3 = DOT1(Kin1,KF2)
      P1P4 = DOT1(Kin1,KF1)
      P2P3 = DOT1(Kin2,KF2)
      P2P4 = DOT1(Kin2,KF1)
      
      M2 = AM(1)**2

      TH = AM(1)**2 - 2.D0 * DOT1(KIA1,KF1) + K12
      UH = AM(1)**2 - 2.D0 * DOT1(KIA2,KF1) + K22
      SH = K12 + 2.D0 * DOT1(KIA1,KIA2) + K22
      


cc      test=Sh_new+th+uh-2*AM(1)**2-K12-K22
      test=Sh+th+uh-2*AM(1)**2-K12-K22
	
	
      if(abs(test).gt.1e-4) then
         write(6,*) ' sh+th+uh ,m**2 ',Sh_new+th+uh,2*AM(1)**2+K12+
     +   K22,AM(1)
         write(6,*) ' sh+th+uh-m**2 ',Sh_new+th+uh-2*AM(1)**2-K12-
     +   K22
      endif
C Catani npb 366 (1991) 135
      IF(M2.gt.0.001) then
        E1 = -K12/M2
        E2 = -K22/M2
c        write(6,*) ' mass ',m2
        else
        E1 = -K12*CUT
        E2 = -K22*CUT
c        write(6,*) ' no mass ',m2
      Endif
      SumAL=0d0
      SumAab = 0d0
      SumAnab = 0d0
      SumC1ab = 0d0
      SumC1nab = 0d0
      SumB1ab = 0d0
      SumB1nab = 0d0
      SumB2ab = 0d0
      SumB2nab = 0d0
      SumC2ab = 0d0
      SumC2nab = 0d0
      SumC2 = 0
      Sumcl = 0
      epsilon = 0D0
c      write(6,*) ' K12 ',K12,' K22 ',K22
c      IF(E2.GT.CUT) write(6,*) ' meoffsh ',E1,E2,CUT,M2
      IF(E1.LT.CUT.AND.E2.LT.CUT) THEN
C all incoming partons are on mass shell
C parton 1 on shell, parton 2 on shell
C eq B23,B24
         SumAab = K1K2/K1P4/K1P3*((K1P3**2 + K1P4**2)/K1K2 + 2D0*M2*
     +   (1D0 - M2*K1K2/2D0/K1P4/K1P3))
         SumAnab = -2D0/K1K2*((K1P3**2 + K1P4**2)/K1K2 + 2D0*M2*(1D0 -
     +   M2*K1K2/2D0/K1P4/K1P3))
         p3p4 = DOT1(KF1,KF2)
c	for g*g -> gg	   
	   SumB1ab= 0
	   SumB2ab= 0
C labels as in paper JHEP 092009121 (note: p3 <-> p4 in code) (Formula is not in paper)
         SumB1nab=((k1k2**4+k2p3**4+k2p4**4)*
     &   (k1p4*k2p3+k1p3*k2p4+k1k2*p3p4))/(k1k2*k1p3*k1p4*k2p3*k2p4*p3p4)
c exchange momenta
         SumB2nab=((k1k2**4+k1p4**4+k1p3**4)*
     &   (k2p3*k1p4+k2p4*k1p3+k1k2*p3p4))/(k1k2*k2p4*k2p3*k1p4*k1p3*p3p4)
c        divide by 2, because of identical particals in final state (see Ellis,Stirling,Webber eq. 7.22         
         SumB1nab = (SumB1nab + SumB2nab)/2.
         
c	for g*q -> gq
C labels as in paper JHEP 092009121 (note: p3 <-> p4 in code) (Formula is not in paper)
	   SumC1ab=(k1k2**2 + k2p4**2)/(k1p3*p3p4)
 	   SumC1nab= SumC1ab*(-1d0 + (k1p3*k2p4)/(k1p4*k2p3) + (k1k2*p3p4)/(k1p4*k2p3))   
c exchange momenta
	   SumC2ab=(k1k2**2 + k1p3**2)/(k2p4*p3p4)
 	   SumC2nab= SumC2ab*(-1d0 + (k2p4*k1p3)/(k2p3*k1p4) + (k1k2*p3p4)/(k2p3*k1p4))
         
         SumC1ab = SumC1ab + SumC2ab 
         SumC1nab = SumC1nab + SumC2nab 
c factor 2: according to K.Kutak & factor 2 (p1 <-> p2)             
         SumC1ab = SumC1ab * 2. *2.
         SumC1nab = SumC1nab * 2. *2.
c         write(6,*) ' both partons are on shell:  K12 ',K12,' K22 ',K22
      ELSEIF(E1.LT.CUT.AND.E2.GT.CUT) THEN
C photoproduction, gluon off mass
C parton 1 on shell, parton 2 off shell
C eq. B21, B22
         SumAab = (K1K2/K1P2)**2*((P2P3**2+P2P4**2)/K1P4/K1P3 - 2D0*M2/
     +   K22*(P2P4/K1P3-P2P3/K1P4)**2)
         SumAnab = 0D0
         p3p4 = DOT1(KF1,KF2)
c	for g*g -> gg	   
	   SumB1ab= 0
C labels as in paper JHEP 092009121 (note: p3 <-> p4 in code)
	   SumB1nab=(k1k2/k1p2)**2*(p2p3**4+k1p2**4+p2p4**4)*
     &   (p3p4*k1p2+k1p3*p2p4+p2p3*k1p4)/(p2p3*k1p3*p3p4*k1p2*p2p4*k1p4)
c factor 2 (p1 <-> p2) and divide by 2, because of identical particals in final state (see Ellis,Stirling,Webber eq. 7.22         

c	for g*q -> gq
C labels as in paper JHEP 092009121 (note: p3 <-> p4 in code)
 	   SumC1ab=(k1k2/k1p2)**2*(k1p2**2+p2p4**2)/(k1p3*p3p4)
	   SumC1nab= SumC1ab*((p3p4*k1p2+k1p3*p2p4)/(k1p4*p2p3)-1.)        
c factor 2: according to K.Kutak & factor 2 (p1 <-> p2)             
         SumC1ab = SumC1ab * 2. *2.
         SumC1nab = SumC1nab * 2. *2.
      ELSEIF(E1.GT.CUT.AND.E2.LT.CUT) THEN
C DIS gluon onshell
C parton 1 off shell, parton 2 on shell
C eq. B21, B22
         SumAab = (K1K2/K2P1)**2*((P1P3**2+P1P4**2)/K2P4/K2P3 - 2D0*M2/
     +   K12*(P1P4/K2P3-P1P3/K2P4)**2)
         SumAnab = 0D0
c NIA1 = gluon
C NIA2 = quark
C NF1 = gluon
C NF2 = quark
         p3p4 = DOT1(KF1,KF2)
c	for g*g -> gg	   
	   SumB1ab= 0
C labels as in paper JHEP 092009121 (note: p3 <-> p4 in code)
	   SumB1nab=(k1k2/k2p1)**2*(p1p4**4+k2p1**4+p1p3**4)*
     &   (p3p4*k2p1+k2p4*p1p3+p1p4*k2p3)/(p1p4*k2p4*p3p4*k2p1*p1p3*k2p3)
     
c factor 2 (p1 <-> p2) and divide by 2, because of identical particals in final state (see Ellis,Stirling,Webber eq. 7.22         

c	for g*q -> gq
C labels as in paper JHEP 092009121 (note: p3 <-> p4 in code)
	   SumC1ab=(k1k2/k2p1)**2*(k2p1**2+p1p3**2)/(k2p4*p3p4)
	   SumC1nab= SumC1ab*((p3p4*k2p1+k2p4*p1p3)/(k2p3*p1p4)-1.)    
c factor 2: according to K.Kutak & factor 2 (p1 <-> p2)             
         SumC1ab = SumC1ab * 2. *2.
         SumC1nab = SumC1nab * 2. *2.
c         write(6,*) ' K1 = ', sqrt(-K12),' K2 ',sqrt(-K22)      
cccc         write(6,*) 'cut = ',cut,'e1,e2 ',e1,e2,sqrt(M2), ' K12 = ', K12,' K22 ',K22    

      ELSEIF(E1.GT.CUT.AND.E2.GT.CUT) THEN
C DIS, gluon off mass
C eq B13,B14
C parton 1 off shell, parton 2 off shell
         SumAab = 4D0*(K1P2*K2P1/P1P2**2)**2*(P1P2**2/(TH-M2)/(UH-M2) -
     +   1D0/K12/K22*(P1P2+2D0*P1P4*P2P3/(TH-M2)+ 2D0*P1P3*P2P4/(UH-
     +   M2))**2)
         Del = 2D0*P1P2*(2D0*P1P3*P2P4/P1P2 - 2D0*P1P4*P2P3/P1P2 - 
     +   K12*P1P4/K2P1 + K22*P2P4/K1P2 + K1P4-K2P4)
         SumAnab = (K1P2*K2P1/P1P2**2)**2*(4D0*P1P2**2*( -1D0/(TH-M2)/
     +   (UH-M2)-1D0/SH*(1D0/(TH-M2)-1D0/(UH-M2))* (P1P3/K2P1-P2P3/
     +   K1P2)+1D0/SH*P1P2/K2P1/K1P2)+2D0/K12/K22* (P1P2+4D0*P1P4
     +   *P2P3/(TH-M2)-Del/SH)* (P1P2+4D0*P1P3*P2P4/(UH-M2)+Del/SH))
         SumaL = SumAab*
     +    (1.d0 -1d0/3d0 - 8d0/3d0/SumAab*k2p1**2/p1p2**2*(
     +   (p2p4**2+p2p3**2)/(TH-M2)/(UH-M2) - (2d0*M2 + k12)/k22*
     +   (p2p3/(TH-M2) - p2p4/(UH-M2))**2))

         p3p4 = DOT1(KF1,KF2)

C for DIS gamma*q -qg
	   SumC2ab=(k1k2/k2p1)**2*(k2p1**2+p1p3**2)/(k2p4*p3p4)
	   SumC2nab= SumC2ab*((p3p4*k2p1+k2p4*p1p3)/(k2p3*p1p4)-1.)   
c         write(6,*) ' SumC2ab ',SumC2ab,SumC2nab 

      ENDIF


      IF(FIRST) THEN
	   Xdum = 0.001
	   SFSdum = 2.
	   STMdum = 5.
         call cauniglu(K(2,2),Xdum,SFSdum,STMdum,XPQdum)
         IQBAR = Iqqbar
         IF(IQBAR.EQ.1) THEN
            write(6,*) ' meoffsh - scale: qq-pair '
         ELSE
            write(6,*) ' meoffsh - scale: q or q_bar '
         ENDIF
	      write(6,*) ' meoffsh - renorm. scale ',Qscal
	      write(6,*) ' meoffsh - fact. scale   ',scalfaf
         FIRST=.FALSE.
      ENDIF

c	write(6,*) K12,k22,K(nia1,2),nia1,K(nia1,2),nia2
c	call pylist(1)
      QG1 = P(NIA1,1)**2+ P(NIA1,2)**2
      QG2 = P(NIA2,1)**2+ P(NIA2,2)**2
      SFSCAL1 = DABS(QG1)
      SFSCAL2 = DABS(QG2)
c      write(6,*) 'k12 ',k12,' k22 ',k22
c      write(6,*) 'qg1 ',qg1,' qg2 ',qg2
c new scale
      PT2Q = (P(KF1,1)**2 +P(KF1,2)**2+P(KF2,1)**2 +P(KF2,2)**2)/2.
c      PT2Q = (P(KF1,1)**2 +P(KF1,2)**2) 
      xg_bar =  -qg2/(1.d0-xg)/SSS
      yph = dble(YY)
      yph_bar =  -qg1/(1.d0 - yph)/SSS

      xg_bar = 0.  ! hannes hack
      yph_bar = 0. ! hannes hack

	if(sfscal1.le.sfcut) sfscal1=sfcut 
	if(sfscal2.le.sfcut) sfscal2=sfcut 
      xsi_hard(1) = (xg+yph_bar)/(yph+xg_bar)	
      xsi_hard(2) = (yph+xg_bar)/(xg+yph_bar)

c which of the quarks is closer to the proton:
      rap1 = rapp-PYP(NF1,17)
      rap2 = rapp-PYP(NF2,17)
      zz1 = dabs(qg1)/yph/sss/xg
      zz2 = dabs(qg2)/yph/sss/xg
c      write(6,*) ' old zz1,2 ',zz1,zz2
      zz1 = 1. - (P(NF1,4)+ P(NF1,3))/2./P(2,4)/yph
      zz2 = 1. - (P(NF2,4)+ P(NF2,3))/2./P(2,4)/xg
* new
      ptstm1=dsqrt(dabs(qg1)+P(NF1,5)**2)
      ptstm2=dsqrt(dabs(qg2)+P(NF2,5)**2)

      If(IQBAR.eq.1) then
         stmax1 = yph*sqrt(xsi_hard(1)*SSS)
         stmax2 = xg*sqrt(xsi_hard(2)*SSS)
      elseif(IQBAR.EQ.2) then
         stmax1 = ptstm1/(1.-zz1)
         stmax2 = ptstm2/(1.-zz2)
      elseif(IQBAR.EQ.3) then
         QT2 = (P(KIA1,1)+P(KIA2,1))**2+(P(KIA1,2) +P(KIA2,2))**2
         stmax1 = sqrt(sh + QT2)
         stmax2 = stmax1
      else
         write(6,*) ' meoffsh. no valid scale selected ',IQBAR
      endif
	
      CALL PYROBO(0,N,0.D0,0.D0,-CM(1)/CM(4),-CM(2)/CM(4),-CM(3)/CM(4)
     +)

c new scale
c check on scale:
      SCALE = PT2Q+M2
      IF(IQ2.EQ.1) THEN
         SCALE = 4.D0*M2
      ELSEIF(IQ2.EQ.2) THEN
         SCALE = SH_NEW
      ELSEIF(IQ2.EQ.3) THEN
         SCALE = 4D0*M2 + PT2Q
c        SCALE = M2 + PT2Q
      ELSEIF(IQ2.EQ.4) THEN
         SCALE = Q2
      ELSEIF(IQ2.EQ.5) THEN
         SCALE = Q2 + PT2Q + 4.D0*M2
      ELSEIF(IQ2.EQ.6) THEN
         SCALE = DABS(K12)+DABS(K22)
c this is for testing
c      ELSEIF(IQ2.EQ.7) THEN
c         sL=2d0*k1p2*k2p1/p1p2
c         rt=1d0/th
c         ru=1d0/uh
c         Scale=2d0*(sL*rt*ru)/((sL*rt*ru)**2+rt**2+ru**2)
      ELSEIF(IQ2.EQ.7) THEN
         Scale = max(stmax1,stmax2)**2
c test average
c         PT2Q = (P(KF1,1)**2 +P(KF1,2)**2+P(KF2,1)**2 +P(KF2,2)**2)/2.
c         SCALE = 4D0*M2 + PT2Q
c test maximum
c         PT2Q = max(P(KF1,1)**2 +P(KF1,2)**2,P(KF2,1)**2 +P(KF2,2)**2)
c         SCALE = 4D0*M2 + PT2Q
      ELSE
         WRITE(6,*) ' NO VALID Q2 SCALE. STOP'
         STOP
      ENDIF
	
c	write(6,*) ' scales stmax = ',stmax,' sqrt(scale) = ',sqrt(scale)
      Dscale = scale
c      alpha_s = ALPHASCA(sqrt(dscale))
	QQg = Qg0

c         write(6,*) ' new alphas ',alpha_s
c      alpha_s = 0.1

      ALPH_EM = ALPH
c	alph_em = 1./137.
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)
      
      do i=-6,6
         XPQ1(i) = 0.d0
         XPQ2(i) = 0.d0
      end do
      stmax1 = stmax1 *scalfaf
      stmax2 = stmax2 *scalfaf
      If(Ires(1).ne.0) then
         if(stmax1.lt.sqrt(sss)) then
            call cauniglu(K(1,2),XG1SF,SFSCAL1,STMAX1,XPQ1)
         endif
      endif
      If(Ires(2).ne.0) then
         if(stmax2.lt.sqrt(sss)) then
            call cauniglu(K(2,2),XG2SF,SFSCAL2,STMAX2,XPQ2)
         endif
      endif
c      write(6,*) ' call alphas '
      alpha_s = ALPHASCA(sqrt(dscale))
     
      
      IF(Ires(1).eq.0.AND.Ires(2).ne.0) THEN
C gamma glu -> q qbar
         Cab = 16.D0*PI*PI*ALPHA_S*ALPH_EM*QF2
         CabB = 0D0
         epsilon = yph**2/2d0/(1d0-yph+yph**2/2d0)
         
         If(SumAL.lt.0) SumAl = 0  ! hannes hack
c         SumAl = 0.
         
         if ((SumAab - epsilon*SumAL).le.0) then
cc            write(6,*) ' meoffsh: problem with long. photon:', (Aab -
cc     +      epsilon*AL)/Aab
         endif
         If(SumAL.gt.0) then
            if ((SumAab - epsilon*SumAL)/SumAab.gt.10.0.or.
     +        (SumAab - epsilon*SumAL)/SumAab.lt.0.1) then
c               write(6,*) ' meoffsh: problem with long. photon:', (SumAab -
c     +         epsilon*SumAL)/SumAab,SumAab,SumAL,epsilon,yph
            endif
c	      write(6,*) 'here entered ',K(NIA1,2),K(NIA2,2),QF2
            if(epsilon.le.0.or.epsilon.gt.1.) then
               write(6,*) ' meoffsh: problem with long. photon:', (SumAab -
     +         epsilon*SumAL)/SumAab,SumAab,SumAL,epsilon,yph
            endif
         endif
c         write(6,*) ' after corr meoffsh ', SumAab,SumAL
c         SumAab = max(0.,SumAab - epsilon*SumAL)  ! hannes hack
c         write(6,*) ' meoffsh IFLflag = ',IFLflag
         If(IFLflag.eq.0) then 
            SumAab = SumAab - epsilon*SumAL ! hannes hack
            Elseif(IFLflag.eq.1) then
            SumAab = SumAab
            Elseif(IFLFLag.eq.2) then
            SumAab = SumAl
            Else
            write(6,*) ' Flag IFLflag = ',IFLflag,' not found: STOP'
            Stop
         Endif
         WTG= dble(XPQ2(0))
         SUMA = (Cab*SumAab + CabB*SumAnab)*WTG
cc this is for gamma*q -> gq
c         XPQC2 = 0
c	   XPQC2 = dble(XPQ2(3)+xpq2(-3))
c	   Cab = 16.D0*PI*PI*ALPHA_S*ALPH_EM*QF2*8./3.
c	   CabB = 0D0
c         SUMC2 = (Cab*SumC2ab + CabB*SumC2nab)*XPQC2
c	   SUMC2 = 0
c         
c         SUMC =  SUMC2
c       write(6,*) ' meoffsh checks 4',SumC2ab,SumC,Cab,XPQC2
   	   Rleg = 0.6
         K(NIA1,2) = 22
         K(NIA2,2) = 21
cc this is for gamma*q -> gq
c	   XPQC2 = dble(XPQ2(3)+xpq2(-3))
c	   Cab = 16.D0*PI*PI*ALPHA_S*ALPH_EM*QF2*8./3.
c	   CabB = 0D0
c         SUMC2 = (Cab*SumC2ab + CabB*SumC2nab)*XPQC2
c	   SUMC2 = 0
c         SUMC =  SUMC2
c ME from Coll-Fact
         xpqc2 = 0.
c         do i=-3,3
c         XPQC2= xpqc2 + XPQ2(I)*DFLOAT(PYCHGE(I))**2/9.D0
c         write(6,*) ' xpqc2 ',i,xpqc2
c         end do
	   XPQC2 = dble(XPQ2(3)+xpq2(-3))*DFLOAT(PYCHGE(3))**2/9.D0
         Q2 = -k12
         if(abs(th).gt.10.and.abs(uh).gt.1.and.abs(sh).gt.1) then
            SUMC2 = 16.D0*8.D0/3.D0*PI*PI*ALPHA_S*ALPH_EM*(-TH/SH - SH/TH + 2.D0*UH*Q2/SH/
     +      TH - 2.D0*UH*Q2/(SH+Q2)**2)
            SUMCL = 16.D0*32.D0/3.D0*PI*PI*Q2*(Q2+SH+TH)/(SH+Q2)**2
            SUMCL = 0.
cc         SUMC2 = 16.D0*8.D0/3.D0*PI*PI*ALPHA_S*ALPH_EM*(-uH/SH - SH/TH + 2.D0*tH*Q2/SH/
c     +   uH - 2.D0*tH*Q2/(SH+Q2)**2)
         else
            SUMC2=0
            SUMCL = 0.
         endif
                   If(IFLflag.eq.0) then 
C transverse and longitudinal          
            SumC = SumC2 - epsilon*SumCL ! hannes hack
            Elseif(IFLflag.eq.1) then
C transverse         
            SumC = SumC2  
            Elseif(IFLFLag.eq.2) then
C longitudinal          
            SumC = SumCL ! hannes hack
         Endif
c ensure DGLAP treatment
	   If(sfscal2.gt.scale*scalfa) SUMC = 0 
	   If(sfscal2.gt.STMAX2**2) SUMC = 0 
         SUMC =  SUMC*XPQC2
      ELSEIF(Ires(1).ne.0.and.Ires(2).ne.0) THEN
C glu glu -> qqbar
         Cab = 16D0*PI**2*ALPHA_S**2/2D0/NC
         CabB = 16D0*PI**2*ALPHA_S**2/2D0*NC/(NC**2-1D0)
         Cab  = Cab*NQ
         CabB = CabB*NQ
         SUMA = 0
         If(IRPAT.ne.0) then
           WTG= dble(XPQ1(0)*XPQ2(0))
           SUMA = (Cab*SumAab + CabB*SumAnab)*WTG
         Endif
C glu glu -> glu glu
         SUMB = 0
         IF(IRPBT.ne.0) then
           Cab = 16D0*PI**2*ALPHA_S**2*2D0/9D0
           CabB = 16D0*PI**2*ALPHA_S**2*9D0/8D0
           WTG= dble(XPQ1(0)*XPQ2(0))
           SUMB = (Cab*SumB1ab + CabB*SumB1nab)*WTG
         Endif
C glu q -> glu q 
         SUMC = 0
         SUMC1 =0
         SUMC2 =0
         IF(IRPCT.ne.0) then        
           Cab = 16D0*PI**2*ALPHA_S**2*2D0/9D0
           CabB = 16D0*PI**2*ALPHA_S**2/4D0 
           XPQC1=XPQ2(1)+XPQ2(2)+XPQ2(3)+XPQ2(-1)+XPQ2(-2)+XPQ2(-3)
           WTQ1= dble(XPQ1(0)*XPQC1)
           XPQC2=XPQ1(1)+XPQ1(2)+XPQ1(3)+XPQ1(-1)+XPQ1(-2)+XPQ1(-3)
           WTQ2= dble(XPQ2(0)*XPQC2)
           If(Rleg.gt.0.5) then
c KIA2 is onshell          
c factor  2 because of Rleg selection, but both contributions must be added
   	        SUMC2=(Cab*SumC2ab+CabB*SumC2nab)*WTQ2*2
c              write(6,*) ' Rleg > 0.5 ',Sumc1,SumC1ab,SumC1nab,XPQ1(0),XPQC1
              else
c KIA1 is onshell           
c factor  2 because of Rleg selection, but both contributions must be added
   	        SUMC2=(Cab*SumC1ab+CabB*SumC1nab)*WTQ1*2
c              write(6,*) ' Rleg < 0.5 ',Sumc2,SumC1ab,SumC1nab,XPQ2(0),XPQC2
           endif
           SUMC = SUMC1 + SUMC2
         Endif

      ELSEIF(Ires(1).eq.0.AND.Ires(2).eq.0) THEN
C gamma gamma -> q qbar
         Cab = 2D0*(4*PI*ALPH_EM)**2*QF2*QF2*NC
         CabB = 0
         SUMA = (Cab*SumAab + CabB*SumAnab)         
         K(NIA1,2) = 22
         K(NIA2,2) = 22
         K(NF1,2) = -KPA
         K(NF2,2) =  KPA
      ENDIF
      
      SUM2 = SUMA + SUMB + SUMC
      rntest = dcasrn()
      if(sum2.le.0) goto 20
c      if(Sum2.gt.0) then
c         write(6,*) ' suma,sumb,sumc ',suma,sumb,sumc
c         write(6,*) ' K1 = ', sqrt(-K12),sqrt(qg1),
c     &    ' K2 ',sqrt(-K22),sqrt(qg2)      
c      endif
	IF((SUMA+SUMB)/SUM2.LT.rntest) THEN
         IRESPRO = 3
      ELSEIF(SUMA/SUM2.LT.rntest) THEN
         IRESPRO = 2
      ELSE
         IRESPRO = 1
      ENDIF

c      write(6,*) ' meoffsh NFLT ',NFLT,NFLAV

C new cut for on-shell DGLAP
	If(Iccfm.eq.0) then
	   If(ires(1).eq.1.and.(sfscal1.gt.scale*scalfa)) SUM2 = 0 
	   If(ires(1).eq.1.and.(sfscal1.gt.STMAX1**2)) SUM2 = 0 
	   If(ires(2).eq.1.and.(sfscal2.gt.scale*scalfa)) SUM2 = 0 
	   If(ires(2).eq.1.and.(sfscal2.gt.STMAX2**2)) SUM2 = 0 
	Endif

c      Write(6,*) ' meoffsh Cab,Aab,CabB,Anab ',Cab,Aab,CabB,Anab
c      Write(6,*) ' meoffsh Wtg ',Wtg,XG1SF,XG2SF
c      Write(6,*) ' meoffsh Cab,Aab,CabB,Anab ',Cab,Aab,CabB,Anab
      IF(SUM2.LT.0D0) THEN
c     Write(6,*) ' meoffsh Cab,Aab,CabB,Anab ',Cab,Aab,CabB,Anab
c     write(6,*) ' SUM2 ',SUM2
         if(sum2.ne.sum2) then
            write(6,*) ' meoffsh SUM2,xsi_hard ',SUM2,
     &        xsi_hard(1),xsi_hard(2)
            Write(6,*) ' meoffsh Cab,Aab,CabB,Anab ',Cab,SumAab,CabB,
     +      SumAnab
            write(6,*) ' meoffsh k12,k22 ',k12,k22
            write(6,*) ' meoffsh th,uh ',th,uh
         endif
         SUM2 = 0D0
      ENDIF
c apply pt cuts here
cc      PT2 = PYP(KF1,9)
      PT2 = PYP(NF1,9)
      CALL CUTG(PT2,NACC)
      IF(NACC.EQ.0) SUM2 = 0
cc      PT2 = PYP(KF2,9)
      PT2 = PYP(NF2,9)
      CALL CUTG(PT2,NACC)
      IF(NACC.EQ.0) SUM2 = 0

      PT2 = min(PYP(KF1,9),PYP(KF2,9)) ! hannes new
      If(PT2.le.PARU(112)**2) SUM2 = 0 ! hannes new

c cut to ensure collinear treatment
      IF(Ires(1).ne.0.and.IPRO.EQ.10.and.(IRPB.eq.1.or.IRPC.eq.1)) then
         If(Rleg.gt.0.5) then
            if(sfscal2.ge.scale) then
               sum2 = 0
            endif
         Else
            if(sfscal1.ge.scale) then
               sum2 = 0
            endif
         Endif
      endif
c apply additonal factor of 2 in case we select one on shell leg
c since we do not add the xsections but rather put it symmetric only
 
      WT1 = SUM2*xsecnorm

      if(sum2.ne.sum2) then

         write(6,*) ' meoffsh xpq(0),alpha_s,alph_em',xpq1(0),xpq2(0),
     +      alpha_s,alph_em
         write(6,*) ' meoffsh SUM2,xsi_hard ',SUM2,
     &      xsi_hard(1),xsi_hard(2)
         Write(6,*) ' meoffsh Cab,Aab,CabB,Anab ',Cab,SumAab,CabB,SumAnab
         write(6,*) ' meoffsh k12,k22 ',k12,k22
         write(6,*) ' meoffsh th,uh ',th,uh
         write(6,*) ' P(KIA1) ',(P(KIA1,I),I=1,5)
         write(6,*) ' P(KIA2) ',(P(KIA2,I),I=1,5)
         write(6,*) ' P(KF1) ',(P(KF1,I),I=1,5)
         write(6,*) ' P(KF2) ',(P(KF2,I),I=1,5)

      endif
      CALL PYEDIT(13)
      RETURN
   20 WT1 = 0.d0
      CALL PYEDIT(13)
      RETURN
      END


      FUNCTION FXN1(X,WEIGHT)

      IMPLICIT none
      Double Precision X(20)
#include "caslujets.inc"
#include "casparas.inc"
#include "caspart.inc"
#include "casluco.inc"
      DOUBLE PRECISION WEIGL,WEIGC,WEIGB
      COMMON /weiglh/weigl,weigc,weigb
      Integer IGENFL,IMEPHI
      COMMON/GENWEI/IGENFL,IMEPHI
      Real NCALL
      DATA NCALL/0./
      Integer IPROO
      Double Precision WEIG,FXN1,XSEC,WEIGHT
      Double Precision xsecval

      logical first
      double precision au 
      common/f2fit/au(50),first
      double precision xnorm 

      real qg22t,shat1t
      common/mytest/qg22t,shat1t


      IMEPHI = 0
      WEIG=0.D0
      FXN1=0.D0
      IPROO = IPRO

      Weig = 0.D0
      Weigl = 0.D0
      Weigc = 0.D0
      Weigb = 0.D0
      Xsec = 0.D0

      if(ipro.eq.2) then
c         Call Xsecpsi(X,Xsec)
         Weig = Xsec
      elseif(ipro.eq.3) then
c         Call Xsecchi(X,Xsec)
         Weig = Xsec
      elseif(ipro.eq.4) then
         Call Xsecqg(X,Xsec)
         Weig = Xsec
      elseif(ipro.eq.102) then
c         Call Xhiggs(X,Xsec)
         Weig = Xsec
      elseif(ipro.eq.500) then 
c         Call XV1(X,Xsec)
         Weig = Xsec
      elseif(ipro.eq.501) then
c         Call XVq(X,Xsec)
         Weig = Xsec
      elseif(ipro.ge.503.and.ipro.le.509) then
c         Call XVQQ(X,Xsec)
         Weig = Xsec
      elseif(ipro.ge.513.and.ipro.le.519) then
c         Call XVQQ(X,Xsec)
         Weig = Xsec

      else

         Call Xsec1(X,Xsec)

         Weigl = Xsec
         Weig = Xsec
         ncall = ncall + 1
         If(xsec.lt.0) then
         write(6,*) ' fxn1 light q < 0 '
         endif
         IF(NFLAV.GT.3.and.IPRO.EQ.10) THEN
            IMEPHI = 1
            IPRO=11
            IHFLA = 4
            Call Xsec1(X,Xsec)
         If(xsec.lt.0) then
         write(6,*) ' fxn1  charm < 0 '
         endif

            Weigc = Xsec

            WEIG = WEIGL + WEIGC
            ncall = ncall + 1
            IF(NFLAV.EQ.5) THEN
               IMEPHI = 1
               IPRO=11
               IHFLA = 5
               Call Xsec1(X,Xsec)
         If(xsec.lt.0) then
         write(6,*) ' fxn1  bottom < 0 '
         endif

               Weigb = Xsec

               WEIG = WEIGL + WEIGC + WEIGB
               ncall = ncall + 1
            ENDIF
         ENDIF
c add here also the valence quarks
c         IPRO=4
c         Call Xsecqg(X,Xsecval)
c         xnorm = au(7)
c         if(xnorm.lt.0) xnorm = 0
c         write(6,*) ' fxn1 xnorm ',xnorm
c         Weig = Weig + xnorm*Xsecval
c         if(weig.lt.0) then
c         write(6,*) ' valence q le 0 ',Weig,Xsecval,xnorm
c         endif
         weig=abs(weig)
c        
      endif
      weig=abs(weig)
      if(weig.lt.0) then 
         write(6,*) ' weig < 0', weig,ncall
c         weig = 1E-43         
         weig = 0d0    
      endif
      IMEPHI = 0
      IPRO = IPROO
      FXN1=WEIG
      WEIGHT = WEIG
c      write(6,*)' fxn1 ',weig,ncall
      
      RETURN
      END



      SUBROUTINE EVENT
      IMPLICIT NONE
      Integer NDIM
      COMMON/INTDIM/ NDIM
#include "caseffic.inc"
#include "caspara.inc"
#include "caspara1.inc"
#include "casluco.inc"
#include "caspart.inc"
#include "casludat1.inc"
#include "casgki.inc"
#include "caslujets.inc"
#include "casvalues.inc"
#include "caspara5.inc"
#include "casmult.inc"
#include "casherup.inc"
#include "caspar.inc"

	Integer MSTP,MSTI
	Double Precision PARP,PARI
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)

      Double Precision X(20)
      Double Precision  XGR(20)
      Integer IMIX
      COMMON /OALPHAS/ IMIX

      Double Precision XV
      Integer NDIMEN
      COMMON /XVAL/ XV(20),NDIMEN
      Double Precision dcasrn
      Integer IGENFL,IMEPHI
      COMMON/GENWEI/IGENFL,IMEPHI

      Integer Nqqb,Nccb,Nbbb
      Common/NEVOUT/Nqqb,Nccb,Nbbb

      Integer LERR
      COMMON/ERR/LERR(100)
      Integer IDEBUG
      COMMON/ INTERN/IDEBUG
      Double Precision remPARJ32
      LOGICAL FIRST,first2

      Integer NPRIN,IDEBUGO,IPROO,KPAO,KINT22,MSTP128
      Integer IHFO
      Integer I,J
      Double Precision WEIGL,WEIGC,WEIGB
      COMMON /weiglh/weigl,weigc,weigb

      Double Precision WTOT,Xsec
      Integer IQUARK
      Double Precision p1,p2,p3,p4,pc1,pc2,pc3,pc4,chec
      EXTERNAL dcasrn
      DATA FIRST/.TRUE./
      DATA FIRST2/.TRUE./
      DATA NPRIN/0/
      SAVE remPARJ32, IHFO, IPROO, KPAO, KINT22   



      RETURN
      END


      SUBROUTINE caend(IK)
      IMPLICIT None
      CHARACTER* 16 CNAM
      CHARACTER* 6 CBEAM1,CBEAM2
      Integer IMIX,IK,LERR
      COMMON /OALPHAS/IMIX
      COMMON/ERR/LERR(100)
      Integer Nqqb,Nccb,Nbbb
      Common/NEVOUT/Nqqb,Nccb,Nbbb
      REAL ACC1,ACC2
      Integer IINT,NCB
      COMMON /INTEGR/ ACC1,ACC2,IINT,NCB
	Integer n1min,n1max,n2min,n2max,n3min,n3max,ncall
	common/caerrstf/n1min,n1max,n2min,n2max,n3min,n3max,ncall

#include "caspart.inc"
#include "caseffic.inc"
#include "casparas.inc"
#include "caslujets.inc"
#include "casluco.inc"
	Integer haprint
	common/haprint/haprint
      Double Precision xs_l,xs_c,xs_b
	if(haprint.eq.1) goto 1000
      IF(IK.EQ.10) WRITE(6,*) 'TIME LIMIT REACHED........'
      CALL PYNAME(KBEAM(1,2),CNAM)
      CBEAM1 = CNAM(1:6)
      CALL PYNAME(KBEAM(2,2),CNAM)
      CBEAM2 = CNAM(1:6)
      CALL PYNAME(KINT(2,2),CNAM)
      WRITE(6,*) ' x - section for ',CBEAM1,CBEAM2,' --> Q Q_bar '
     +//'X"'
      WRITE(6,10000) AVGI,SD
10000 FORMAT('  sigma  = ',G16.7,' nb   +/- ',G16.7)
      IF(IK.EQ.20) THEN
         WRITE(6,10100) AVGI,SD
10100 FORMAT('  sigma visible = ',G16.7,' nb   +/- ',G16.7)
         xs_l = dble(Nqqb)/dble(nout)*avgi
         xs_c = dble(Nccb)/dble(nout)*avgi
         xs_b = dble(Nbbb)/dble(nout)*avgi
         write(6,*) ' light  quarks : xsec = ',xs_l,' Nqqb = ',Nqqb
         write(6,*) ' charm  quarks : xsec = ',xs_c,' Nccb = ',Nccb
         write(6,*) ' bottom quarks : xsec = ',xs_b,' Nbbb = ',Nbbb

         WRITE(6,*) ' Nr of events generated : ',NIN
         WRITE(6,*) ' Nr of events written: ',NOUT
         WRITE(6,*) ' Error summary on event generation '
         WRITE(6,*) ' Errors and their meaning meaning:'
         WRITE(6,*) ' CAPS x > 0.999...........: ',LERR(45)
         WRITE(6,*) ' CAPS boost PS error......: ',LERR(46)
         WRITE(6,*) ' CAPS energy not conserved: ',LERR(100)
         WRITE(6,*) ' CAREMN frag. cuts.........: ',LERR(48),LERR(49)
         WRITE(6,*) ' PYPREP error..............: ',LERR(50)
         WRITE(6,*) ' cascps check 1st..........: ',LERR(55)
         WRITE(6,*) ' cascps no more memory.....: ',LERR(51),LERR(52)
         WRITE(6,*) ' cascps boost error........: ',LERR(53)
         WRITE(6,*) ' cascps xfb(iflb)=0........: ',LERR(54)
         WRITE(6,*) ' cascps xsi < 0 ...........: ',LERR(58)
         WRITE(6,*) ' cascps NTRY > 15000.......: ',LERR(56)
         WRITE(6,*) ' cabran stopping < 0 ......: ',LERR(59)
         WRITE(6,*) ' cabran Wt > 1.............: ',LERR(60)
	   
	   write(6,*) ' calls to ccfmstfu:  ',ncall
	   write(6,*) ' ccfmstfu: x out of range (min/max) ',n2min,n2max
	   write(6,*) ' ccfmstfu: k2 out of range (min/max) ',n1min,n1max
	   write(6,*) ' ccfmstfu: p out of range (min/max) ',n3min,n3max

      ENDIF
1000      RETURN
      END
      

      SUBROUTINE calcglu 
      Implicit None
#include "casludat1.inc"      
	character *72 TXT
      Double Precision  XI,KT2I,TI
      Real X,KT2,T,scale
	integer kfl
	common /gluon_int/kfl,x,scale,t
      Logical FIRSTmy
      Data FIRSTmy/.TRUE./
      Logical start
      data start/.true./
      LOGICAL ex
      Integer IRR
	Integer Iccfmset
      
      Integer ikincut,Ipgg,ns_sel,Iqqbar
      Real Qg,Qg0,Qscal
      Double Precision x1minL,x1maxL,x1dif,x2dif,x3dif
      Double Precision x2minL,x2maxL,q2min,q2max,p3max,p3min 
	Double Precision x3lmin,x3lmax,x3ldif
      Double Precision x1min,x1max,x2min,x2max,x3min,x3max
	Double Precision Qgmin
      
      logical first
      double precision au 
      common/f2fit/au(50),first
	
	Real A0,B0,EPS,SUM,XA_INT,Sum1,XA_INTq
      Real result,err
	External XA_INT,XA_INTq
      CHARACTER   filename1*132, filname*132
	character adum
      data filename1/'ccfm-test.dat'/

      character filename2*132
      Common/updfout/filename2
      
      Integer Noevol
      Common/Evoltest/Noevol
      Integer Nbx,Nbk,Nbp
      Parameter (Nbx=50)
      Parameter (Nbk=50)
      Parameter (Nbp=50)
      Integer i,j,k
      Real X1(0:Nbx+1),X2(0:Nbk+1),X3(0:Nbp+1)
      Real XPQS(-6:6)
      Integer NbxL,Nbxlin
      Double Precision x1minlin,x1maxlin
      real sigm,xnorm,xmean,xnor
      
      Real abstol 
      Integer nseg
      Real QCDlam

c Evolution or only starting distribution
      Noevol = 0

      write(6,*) ' here in calglu1  ',filname,filename2,filename1
      if(filename2.eq.' ' ) then
        filname = filename1
        write(6,*) ' calcglu old naming: ',filname,filename1
        else
        filname = filename2
        write(6,*) ' calcglu new name: ',filname,filename2
      endif
      IF(FIRSTmy) THEN
         FIRSTmy=.FALSE.
         Write(6,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
         Write(6,*) '++++  calculate gluon density    ++++++++++++++'
         Write(6,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
         write(6,*) 'Cascade/updfGrids/reading ccfm-grid.dat '
         open(30,FILE='Cascade/updfGrids/ccfm-grid.dat', FORM='formatted',
     +   STATUS= 'OLD',IOSTAT=IRR,ERR=50 )
c         write(6,*) ' ccfm-kernel.dat '
c         open(30,FILE='ccfm-kernel.dat', FORM='formatted',
c     +   STATUS= 'OLD',IOSTAT=IRR,ERR=50 )
200      Read(30,101) TXT
  101    Format(A72)
C         WRITE(6,101) '  line ',TXT
         If(TXT(1:4).EQ.'  Qg') then 
	      read(txt,1000) adum,Qg0,adum,ikincut
1000        format(A7,f12.8,A10,I6)
c	      read(txt,1000) Qg0,ikincut
c1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Qg0') then 
c	      read(txt,1000) Qg0,ikincut
c         WRITE(6,101) ' 1st line ',TXT
	      read(txt,1000) adum,Qg0,adum,ikincut
C1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
c1000        format(A7,f12.8,A10,I6)
            goto 200
         Endif
         If(TXT(1:4).EQ.' Ipg') then 
c	      read(txt,1001) Ipgg,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
	      read(txt,1001) adum,Ipgg,adum,ns_sel
1001        format(A8,I4,A10,I4)
c         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:6).EQ.' Qscal') then 
	      read(txt,1002)  adum,Qscal, adum,Iqqbar
1002        format(A9,f7.3,A10,I4)
c	      read(txt,1002) Qscal,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
C         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:7).EQ.' QCDlam') then 
c         WRITE(6,101) '3rd line ',TXT
            read(txt,1003) adum,QCDLam
1003        format(A9,f12.8)
c            PARU(112)=QCDLam
            goto 200
         Endif
         If(TXT(1:4).EQ.' ln(') then 
c         WRITE(6,101) '2 or 3rd line',TXT
         Endif
         close(30)
      ENDIF
      If(QCDlam.le.0) then
         write(6,*) ' No value for QCDlam found on file. Set to ',PARU(112)
         QCDlam = PARU(112)
      Endif
      inquire(FILE=filname,EXIST=ex)
      if(ex) then
         open(31,FILE=filname, FORM='formatted',STATUS=
     +   'old',IOSTAT=IRR,ERR=50 )
          write(6,*) '  old file found: overwritten ',filname
          else
          open(31,FILE=filname, FORM='formatted',STATUS=
     +    'new',IOSTAT=IRR,ERR=50 )
          write(6,*) ' new file created ',filname
	Endif
	write(31,10000) Qg0,ikincut
10000 format(' Qg0 = ',f12.8,' ikincut= ',I6)
      write(31,10001) Ipgg,ns_sel
10001 format(' Ipgg = ',I4,' ns_sel = ',I4)
      write(31,10002) QCDlam
10002 format(' QCDlam = ',f12.8)
      write(31,10100)
10100 format(' ln(xg),  ln(k_t^2), ln(p),  xgx')

	   x1min = 1e-7
c         x1max =.9999999
         x1max = 0.1d0
         x1minL=log(x1min)
         x1maxL=log(x1max)
         NbxL = Nbx-10
         x1dif=(x1maxL-x1minL)/Real(NbxL)
         Do I=0,NbxL+1
            x1(I) = x1minL + x1dif*Real(I)
         Enddo
c lin grid in x from 0.1 - 1
         x1minlin=0.1d0
         x1maxlin=1.0d0
         Nbxlin = 10
         x1dif=(x1maxlin-x1minlin)/Real(Nbxlin)
         Do I=NbxL,NbxL+Nbxlin+1
            x1(I) = log(x1minlin + x1dif*Real(I-NbxL))
         Enddo

         q2min = 0.001 ! new         
         q2max = 40000.
	   x2minL=log(q2min)
	   x2maxL=log(q2max)
	   x2dif=(x2maxL-x2minL)/real(nbk)
         Do I=0,Nbk+1
            x2(I) = x2minL + x2dif*Real(I)
	   Enddo

         p3min=Qg0
 	   p3max= 14000.
	   x3lmin=log(p3min)
	   x3lmax=log(p3max)
	   x3ldif = (x3lmax-x3lmin)/Real(Nbp)	
         Do I=0,Nbp+1
            x3(I) = x3lmin + x3ldif*Real(I)
         Enddo
      
      if(first) then
        write(6,*) ' new parameters in calcglu :', AU(1),AU(2),
     &        AU(3),AU(4),Au(5),Au(6),Au(7),AU(8),AU(9)
      endif

      do j=0,Nbk
         do i=0,Nbx
            do k=0,Nbp
	         A0 = 0.000000001
	         B0 = 1.
	         eps = 0.01 
		   kfl = 2212
               x = exp(x1(i))
               t=exp(x3(k))
               scale=exp(x2(j))
c               write(6,*) ' gadap sum: ',sum
               abstol = 0.0  ! abstol = 0, value of abstol will be ignored
               nseg = 1000
               call readgrid(kfl,Real(x),scale,t,xpqs)
               CALL RADAPT(xa_int,A0,B0,NSEG,eps,abstol,SUM,ERR)
               CALL RADAPT(xa_intq,A0,B0,NSEG,eps,abstol,SUM1,ERR)
c               write(6,*) ' radapt sum: ',sum
cccc               Sum = result
               xpqs(0) = sum
               xpqs(3) = sum1

c               write(6,*) sum
               if(xpqs(0).ne.xpqs(0)) then
                 write(6,*) ' error in calcglu ',x,t,scale,xpqs(0)
                 xpqs(0) = 0
               endif
               if(xpqs(0).ge.1.E10) then
c                 write(6,*) ' error in calcglu xpqs too large ',
c     &                 x,t,scale,xpqs(0)
                 xpqs(0) = 1.E10
               endif
               write(31,*) x1(i),x2(j),x3(k),xpqs(-4),xpqs(-3),xpqs(-2),xpqs(-1),xpqs(0),xpqs(1),xpqs(2),xpqs(3),xpqs(4)
c                  write(6,*) x1(i),x2(j),x3(k),sum,i,j,k
            enddo
         enddo
c         write(6,*) ' calcglu: loop i =',i,j,k
c         write(6,*) ' last loop i =',i,j,k
      enddo
      close(31)
      Write(6,*) ' end of calcglu '
      RETURN
   50 write(6,*) ' calcglu: error opening file'
      Return
      END


      Function xa_int(xt)
	Implicit None
	Real xa_int,x0,z,test
      Real XPQS(-6:6)

      Real X,T,scale
	integer kfl
	common /gluon_int/kfl,x,scale,t
      DOUBLE PRECISION XXX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU,
     +VAL(20),PP2,qq2
      CHARACTER*20 PARM(20)
      LOGICAL PDFFIRST	
	Real xmin,xmax,xt,wt
	Real A1,A2,A3,A4,A5,A6
      logical first
	double precision au 
      common/f2fit/au(50),first
      real xsat
	Integer jj

	Integer Noevol
      Common/Evoltest/Noevol  
	
	xmin=x
	xmax=1.
	wt = log(xmax/xmin)
	x0= xmin*(xmax/xmin)**xt
      
c       write(6,*) ' here in xa_int '
	z = x/x0
	if(x0.lt.x) then
	   xa_int = 0
	   return
	   endif
c	write(6,*) ' x,x0,z,sclae,t ',x,x0,z,scale,t
	
        If(noevol.eq.0) then
           call readgrid(kfl,Real(z),scale,t,xpqs)
c           write(6,*) 'kfl,z,t,xpqs',kfl,z,scale,xpqs(0)
        else
           xpqs(0) = 1.
           z=0.98
           x0 = xt
           wt = 1./scale
        endif


	A1 = au(1)
	A2 = au(2)
      A3 = au(3)
	A4 = au(4)
      A5 = au(5)
      A6 = au(6)
c       write(6,*) a1,a2,a3,a4
       if(z.gt.0.99999) then 
	  xpqs(0) = 1.
	  else
	  endif
       if(z.lt.0.) then
	   test = 0.
	   else
        test =A1*5.*(1./x0**A2)*(1.-x0)**A3*(1.-A4*x0+A5*sqrt(x0)+A6*x0**2)*xpqs(0)/x0
         test = max(0.,test)	   
      endif 
c     f(x)/g(x)*int(g(x)dx) = f(x)*x*wt
	xa_int = test * wt * x0	
	Return 
	
	end
      
      Function xa_intq(xt)
	Implicit None
	Real xa_intq,x0,z,test
      Real XPQS(-6:6)

      Real X,T,scale
	integer kfl
	common /gluon_int/kfl,x,scale,t
      DOUBLE PRECISION XXX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU,
     +VAL(20),PP2,qq2
      CHARACTER*20 PARM(20)
      LOGICAL PDFFIRST	
	Real xmin,xmax,xt,wt
	Real A1,A2,A3,A4,A5,A6
      logical first
	double precision au 
      common/f2fit/au(50),first
      real xsat
	Integer jj

	Integer Noevol
      Common/Evoltest/Noevol        
	
	xmin=x
	xmax=1.
	wt = log(xmax/xmin)
	x0= xmin*(xmax/xmin)**xt
      
c       write(6,*) ' here in xa_int '
	z = x/x0
	if(x0.lt.x) then
	   xa_intq = 0
	   return
	   endif
c	write(6,*) ' x,x0,z,sclae,t ',x,x0,z,scale,t
	
        If(noevol.eq.0) then
           call readgrid(kfl,Real(z),scale,t,xpqs)
c           write(6,*) 'kfl,z,t,xpqs',kfl,z,scale,xpqs(0)
	   xpqs(3) = xpqs(3)
        else
           xpqs(3) = 1.
           z=0.98
           x0 = xt
           wt = 1./scale
        endif


	A1 = au(1)
	A2 = au(2)
      A3 = au(3)
	A4 = au(4)
      A5 = au(5)
      A6 = au(6)
c       write(6,*) a1,a2,a3,a4
       if(z.gt.0.99999) then 
	  xpqs(0) = 1.
	  else
	  endif
       if(z.lt.0.) then
	   test = 0.
	   else
        test =A1*5.*(1./x0**A2)*(1.-x0)**A3*(1.-A4*x0+A5*sqrt(x0)+A6*x0**2)*xpqs(3)/x0
         test = max(0.,test)	   
      endif 
c     f(x)/g(x)*int(g(x)dx) = f(x)*x*wt
	xa_intq = test * wt * x0	
	Return 
	
	end
      SUBROUTINE readgrid(KF,X,Q2,P,XPQ)
      Implicit None
      Integer n1,n2,n3
      REAL XPQ(-6:6),X,Q2,P
*! new
      Parameter (n1=51,n2=51,n3=51)
      Real Q2x,xx,px,xpd0,xpd1,xpd2,xpd3
      DIMENSION Q2X(0:n1),XX(0:n2),PX(0:n3),XPD0(0:n2,0:n1,0:n3),XPD1(0:n2,0:n1,0:n3),XPD2(0:n2,0:n1,0:n3),XPD3(0:n2,0:n1,0:n3)
      Real amax
      Dimension amax(0:n2)
      Real kl,kt,xpqm
c      Double Precision XA(3),f_grid(n1,n2,n3),A(N1+N2+N3)
c      Double Precision f_gridq(n1,n2,n3)
      Double Precision XA(3),f_grid0(n1,n2,n3),f_grid1(n1,n2,n3),f_grid2(n1,n2,n3),f_grid3(n1,n2,n3),A(N1+N2+N3)
      INTEGER NA(3)
      DATA NA/n1,n2,n3/
      Double Precision DHFINT
	Integer  ikincut,Ipgg,ns_sel
	Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Integer IGLU
      Common/CAGLUON/Iglu
	Double Precision Qscal
	Integer Iqqbar
	Common/GLUDAT2/Qscal,Iqqbar
      Integer i,j,k,in,irr,igrid,ip,n2test
      Integer ipoint,kf,it
	Integer Iset,Iccfmset
      Double Precision scal,rx,rq2,rp,glu,uval,dval,ssea,xpqmax,xm
      Double Precision RKMS4B,RKMS3B,RKMS2B,RKMS1B,RKMS0,RKMS1,RKMS2,RKMS3,RKMS4
      Real xl,pl
	Integer n1min,n1max,n2min,n2max,n3min,n3max,ncall
	common/caerrstf/n1min,n1max,n2min,n2max,n3min,n3max,ncall
	character *72 TXT
	character adum
      Real QCDlam
      LOGICAL FIRST,afirst
      DATA FIRST/.TRUE./
      Data Afirst/.TRUE./
	if(iabs(kf).eq.2212) then 
	else
         write(6,*) ' readgrid: no valid kf = ',kf
         write(6,*) ' readgrid: only applicable for proton' 	   
	   stop
	endif
      IF(FIRST) THEN
	   n1min = 0
	   n1max = 0
	   n2min = 0
	   n2max = 0
	   n3min = 0
	   n3max = 0
	   ncall = 0
         i=0
         scal = 1.0
         Ipgg = 0
	   ns_sel = -1
	   If(IGLU.GT.1000) then
	   ISET = Int(IGLU/1000)
	   ELSE
         ISET = IGLU
	   ENDIF
         write(6,*) ' start reading Cascade/updfGrids/ccfm-grid.dat '
c         write(6,*) ' start reading ccfm-kernel.dat '
         if(iset.eq.1) then
               open(30,FILE='Cascade/updfGrids/ccfm-grid.dat', FORM='formatted',
c               open(30,FILE='ccfm-kernel.dat', FORM='formatted',
     +         STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
c               write(6,*) ' read unintegrated gluon from file: '
c     +         //'Cascade/updfGrids/ccfm-grid.dat'
         else
            write(6,*) ' No valid IGLU selected for readgrid: IGLU = ',
     +      Iset
            STOP
         endif
c         read(30,10000) Qg0,ikincut
c10000    format(' Qg0 = ',f12.8,' ikincut= ',I6)
c         read(30,10100)
c10100    format('xg,  kt, p  xgx') 
200      Read(30,101) TXT
  101    Format(A72)
C         WRITE(6,101) '  line ',TXT
         If(TXT(1:4).EQ.'  Qg') then 
c	      read(txt,1000) Qg0,ikincut
c1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
	      read(txt,1000) adum,Qg0,adum,ikincut
C1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
1000        format(A7,f12.8,A10,I6)
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Qg0') then 
c	      read(txt,1000) Qg0,ikincut
	      read(txt,1000) adum,Qg0,adum,ikincut
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Ipg') then 
c	      read(txt,1001) Ipgg,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
cc         WRITE(6,101) '2nd line ',TXT
	      read(txt,1001) adum,Ipgg,adum,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
1001        format(A8,I4,A10,I4)
            goto 200
         Endif
         If(TXT(1:7).EQ.' QCDlam') then 
c         WRITE(6,101) '3rd line ',TXT
            read(txt,1003) adum,QCDLam
1003        format(A9,f12.8)
c            PARU(112)=QCDLam
            goto 200
         Endif
         If(TXT(1:6).EQ.' Qscal') then 
c	      read(txt,1002) Qscal,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
C         WRITE(6,101) '2nd line ',TXT
	      read(txt,1002)  adum,Qscal, adum,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
1002        format(A9,f7.3,A10,I4)
            goto 200
         Endif
         If(TXT(1:4).EQ.' ln(') then 
c         WRITE(6,101) '2 or 3rd line',TXT
         Endif
	   
c         write(6,*) ' soft cut Q0 ',Qg0,' scale factor = ',scal
c         write(6,*) ' kin cut ',ikincut,' Ipgg = ',Ipgg,
c     &	   ' ns_sel = ',ns_sel
c         write(6,*) ' type of splitting function and non-sudakov: '
	   If(Ipgg.eq.0) then
c            write(6,*) ' Pgg = asb(k)/z*Delta_ns + asb(q)/(1-z) '
		elseif(Ipgg.eq.1) then
c		write(6,*) 
c     &	' Pgg = asb(k)*B*(1-z)/z*Delta_ns + asb(q)*B*z/(1-z)'
c		write(6,*) ' with BB = ',BB
            elseif(Ipgg.eq.2) then
c	     write(6,*) ' Pgg = asb(q)[1/z*Delta_ns + 1/(1-z)]'
            elseif(Ipgg.eq.3) then
c		write(6,*) 
c     &	' Pgg = asb(q)*[B*(1-z)/z*Delta_ns + B*z/(1-z)]'
		else
c     Ipgg = 0 ! Pgg=(as(k)*1/z +as(q)*1/(1-z))   
c     Ipgg = 1 ! Pgg=(as(k)*B(1-z)/z +as(q)*B*z/(1-z)), Gavins full splitting
c     Ipgg = 2 ! Pgg=(as(q)*1/z +as(q)*1/(1-z)), as(q) eveywhere
c     Ipgg = 3 ! Pgg=(as(q)*B(1-z)/z +as(q)*B*z/(1-z)), as(q) eveywhere + full split
 		write(6,*) ' non - valid splitting fct selected ',Ipgg
		stop
	   Endif
c         write(6,*) ' scal factor = ',Qscal,' fact. scale = ',Iqqbar
		
         do j=1,n1
            do i=1,n2
               do k=1,n3
c                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS
cccc                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS,RKMS2
                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,RKMS1,RKMS2,RKMS3,RKMS4
                  xx(i) = rx
                  q2x(j) = rq2
                  px(k) = rp
cc                  xpd(i,j,k) = max(0.,rkms)*scal
cc                  xpdq(i,j,k) = max(0.,rkms2)*scal
                  xpd0(i,j,k) = max(0.,rkms0)*scal
                  xpd1(i,j,k) = max(0.,rkms1)*scal
                  xpd2(i,j,k) = max(0.,rkms2)*scal
                  xpd3(i,j,k) = max(0.,rkms3)*scal
               enddo
            enddo
         enddo
         write(6,*) ' after reading Cascade/updfGrids/ccfm-grid.dat ',j,i,k
c         write(6,*) ' after reading ccfm-kernel.dat ',j,i,k
c we have only gluons
         DO 10 J=1,N2
            DO 10 I=1,N1
               DO 10 k=1,N3

   10    CONTINUE
         IN=0
         DO I=1,n1
            IN=IN+1
            A(IN) = q2x(I)
         ENDDO
         DO I=1,N2
            IN=IN+1
            A(IN) = xx(I)
         ENDDO
         DO I=1,N3
            IN=IN+1
            A(IN) = px(I)
         ENDDO
         DO I=1,N2
            DO k=1,N3
               DO J=1,N1
c                  f_grid(j,i,k) = xpd(i,j,k)
c                  f_gridq(j,i,k) = xpdq(i,j,k)
                  f_grid0(j,i,k) = xpd0(i,j,k)
                  f_grid1(j,i,k) = xpd1(i,j,k)
                  f_grid2(j,i,k) = xpd2(i,j,k)
                  f_grid3(j,i,k) = xpd3(i,j,k)
c                  write(6,*)f_grid(j,i,k) 
               ENDDO
            ENDDO
         ENDDO
         FIRST=.FALSE.
c         write(6,*) '  parton densities read from file unit 30 '
         Close(30)
      ENDIF

      ncall = ncall +1 


c      XA(1) = Q2
c      XA(2) = X
c      XA(3) = P
      XA(1) = log(Q2)
      XA(2) = log(X)
      XA(3) = log(P)
      if(xa(3).lt.px(1)) then
	   n3min = n3min + 1
	  
	   if(n3min.lt.10) then 
c            write(6,*) ' p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' p out of range ',xa(3),px(1)
	   elseif(n3min.eq.10) then 
c            write(6,*) ' p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' last message printed: min p'
	   endif
         xa(3)=px(1)
      endif
      if(xa(3).gt.px(n3)) then
	   n3max = n3max + 1
	   if(n3max.lt.10) then 
c            write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
	   elseif(n3max.eq.10) then 
c            write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
c            write(6,*) ' last message printed: max p'
	   endif
         xa(3)=px(n3)
      endif
      if(xa(2).ge.xx(n2)) xa(2)=xx(n2)-0.0001
      if(xa(2).lt.xx(1)) then
	   n2min = n2min + 1
	   if(n2min.lt.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
	   elseif(n2min.eq.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
c            write(6,*) ' last message printed: min x'
	   endif
	   xa(2)=xx(1)
	   endif 
      if(xa(1).lt.q2x(1)) then
	   n1min = n1min + 1
	   if(n1min.lt.10) then 
c            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
	   elseif(n1min.eq.10) then 
c            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
c            write(6,*) ' last message printed: min k2'
	   endif
         xa(1)=q2x(1)
      endif
      if(xa(1).ge.q2x(n1)) then
	   n1max = n1max + 1
         if(n1max.lt.10) then 
c            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
         elseif(n1max.eq.10) then 
c            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
c            write(6,*) ' last message printed: max k2'
	   endif
         xa(1)=q2x(n1)-0.1
      endif
c check if interpolation or grid wanted
      igrid = 0
      if(igrid.eq.1) then
         if(xa(1).lt.q2x(1)) then
            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
            xa(1)=q2x(1)
         endif
         if(xa(1).ge.q2x(n1)) then
            write(6,*) '  k2 out of range: x = ', x,n1,q2x(n1)
            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
            xa(1)=q2x(n1)-0.1
         endif
         i=0
   20    i=i+1
         if(xa(1).gt.a(na(1))) then
c       write(6,*) ' q2  ',xa(1),a(na(1))
            write(6,*) ' q2 not found ',q2,a(na(1)),q2x(n1),xa(1)
            i=na(1)
         else
            if(xa(1).ge.A(i).and.xa(1).lt.a(i+1)) Then
            else
               if(i.le.na(1)) then
                  goto 20
               else
                  write(6,*) ' at 20: q2 not found ',i,q2
               endif
            endif
         endif
         j=0
   30    j=j+1
         if(xa(2).ge.A(na(1)+j).and.xa(2).lt.a(na(1)+j+1)) Then
         else
            if(j.le.na(2)) then
               goto 30
            else
               write(6,*) ' at 30: x not found ',x,xa(2),j
            endif
         endif
         k=0
   40    k=k+1
         if(xa(3).ge.a(na(1)+na(2)+na(3))) then
            k=na(3)
c       write(6,*) ' p  ',xa(3),a(na(1)+na(2)+na(3))
         else
            if(xa(3).ge.A(na(1)+na(2)+k).and. xa(3).lt.a(na(1)+na(2)+k+1)) Then
            else
               if(k.le.na(3)) then
                  goto 40
               else
                  write(6,*) ' at 40: p not found ',k,p
               endif
            endif
         endif
c do not divide by q2 here, because it is done in ccfmstfu
c         glu = f_grid(i,j,k)/q2
c         glu = f_grid(i,j,k)
c         quark = f_gridq(i,j,k)
         glu = DHFINT(3,XA,NA,A,f_grid0)
         uval = DHFINT(3,XA,NA,A,f_grid1)
         dval = DHFINT(3,XA,NA,A,f_grid2)
         ssea = DHFINT(3,XA,NA,A,f_grid3)
         
      else
c do not divide by q2 here, because it is done in ccfmstfu
c         glu = DHFINT(3,XA,NA,A,f_grid)/q2
c         glu = DHFINT(3,XA,NA,A,f_grid)  
c         quark = DHFINT(3,XA,NA,A,f_gridq)  
         glu = DHFINT(3,XA,NA,A,f_grid0)
         uval = DHFINT(3,XA,NA,A,f_grid1)
         dval = DHFINT(3,XA,NA,A,f_grid2)
         ssea = DHFINT(3,XA,NA,A,f_grid3)
      endif

c           write(6,*) ' new readgrid x,q2,p,xpq',x,q2,p,glu
    
      DO  IP=-6,6
         XPQ(IP)=0.0
      ENDDO
      if(glu.gt.0.) xpq(0) = glu
      if(uval.gt.0.) xpq(1) = uval
      if(dval.gt.0.) xpq(2) = dval
      if(ssea.gt.0.) xpq(3) = ssea
      if(ssea.gt.0.) xpq(-3) = ssea
      if(glu.gt.5E6) then
         write(6,*) ' new readgrid gluon x,q2,p,xpq',x,q2,p,glu
c         write(6,*) ' i,j,k ',i,j,k
      endif
      if(ssea.gt.5E6) then
         write(6,*) ' new readgrid sea x,q2,p,xpq',x,q2,p,ssea
c         write(6,*) ' i,j,k ',i,j,k
      endif
      return


   80 write(6,*) ' readgrid: error in opening file '
      stop
   90 continue
      write(6,*) ' readgrid: end of file at ',i,j,k,RX,RQ2,RP,RKMS1,RKMS2,RKMS3
      stop
      END
      SUBROUTINE CCFMSTFQ(KF,X,Q2,P,XPQ)
      Implicit None
      Integer n1,n2,n3
      REAL XPQ(-6:6),X,Q2,P
*! new
      Parameter (n1=51,n2=51,n3=51)
      Real Q2x,xx,px,xpd1,xpd2
      DIMENSION Q2X(0:n1),XX(0:n2),PX(0:n3)
      DIMENSION XPD1(0:n2,0:n1,0:n3),XPD2(0:n2,0:n1,0:n3)
      Real amax
      Dimension amax(0:n2)
      Real kl,kt,xpqm
      double Precision f_grid1(n1,n2,n3)
      double Precision f_grid2(n1,n2,n3)
      double Precision XA(3),A(N1+N2+N3)
      INTEGER NA(3)
      DATA NA/n1,n2,n3/
      Double Precision DHFINT
#include "casshwr.inc"
#include "casludat1.inc"
	Integer  ikincut,Ipgg,ns_sel
	Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Integer IGLU
      Common/CAGLUON/Iglu
	Double Precision Qscal
	Integer Iqqbar
	Common/GLUDAT2/Qscal,Iqqbar
      Integer i,j,k,in,irr,igrid,ip,n2test
      Integer ipoint,kf,it
	Integer Iset,Iccfmset
      Real scal,rx,rq2,rp,rkms1,rkms2,quark1,quark2,xpqmax,xm
      Real xl,pl
	Integer n1min,n1max,n2min,n2max,n3min,n3max,ncall
	common/caerrstf/n1min,n1max,n2min,n2max,n3min,n3max,ncall
      Real QCDLam
	character *72 TXT
	character adum
      character*100 pdfname
      character*512 PDFPATH
      common/CASPDF/PDFPATH
      INTEGER PDFPATHLEN

      Real Qg0q
      
      LOGICAL FIRST,afirst
      DATA FIRST/.TRUE./
      Data Afirst/.TRUE./
            
	if(iabs(kf).eq.2212) then 
	else
         write(6,*) ' ccfmstfq: no valid kf = ',kf
         write(6,*) ' ccfmstfq: only applicable for proton' 	   
	   stop
	endif
      IF(FIRST) THEN
	   n1min = 0
	   n1max = 0
	   n2min = 0
	   n2max = 0
	   n3min = 0
	   n3max = 0
	   ncall = 0
         i=0
         scal = 1.0
         Ipgg = 0
	   ns_sel = -1
c         write(6,*) ' ccfmtfq pdfpath:',pdfpath
c         PDFPATHLEN=INDEX(PDFPATH,' ')-1
         write(6,*)'read grid files from Cascade/updfGrids/ccfm-quark.dat'
         open(30,FILE="Cascade/updfGrids/ccfm-quark.dat", 
     +   FORM='formatted',STATUS= 'OLD',
     +   IOSTAT=IRR,ERR=80 )
c         read(30,10000) Qg0,ikincut
c10000    format(' Qg0 = ',f12.8,' ikincut= ',I6)
c         read(30,10100)
c10100    format('xg,  kt, p  xgx') 
200      Read(30,101) TXT
  101    Format(A72)
c         WRITE(6,101) '  line ',TXT
         If(TXT(1:4).EQ.'  Qg') then 
	      read(txt,1000) adum,Qg0q,adum,ikincut
C1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
1000        format(A7,f12.8,A10,I6)
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Qg0') then 
	      read(txt,1000) adum,Qg0q,adum,ikincut
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Ipg') then 
	      read(txt,1001) adum,Ipgg,adum,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
1001        format(A8,I4,A10,I4)
c         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:6).EQ.' Qscal') then 
	      read(txt,1002)  adum,Qscal, adum,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
1002        format(A9,f7.3,A10,I4)
C         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:7).EQ.' QCDlam') then 
c         WRITE(6,101) '3rd line ',TXT
            read(txt,1003) adum,QCDLam
1003        format(A9,f12.8)
            PARU(112)=QCDLam
            goto 200
         Endif
         If(TXT(1:4).EQ.' ln(') then 
c         WRITE(6,101) '2 or 3rd line',TXT
         Endif
	   
         write(6,*) ' soft cut Q0 ',Qg0q,' scale factor = ',scal
         write(6,*) ' kin cut ',ikincut,' Ipgg = ',Ipgg,
     &	   ' ns_sel = ',ns_sel
         write(6,*) ' QCD_lam used in uPDF: ',PARU(112)
         write(6,*) ' type of splitting function and non-sudakov: '
	   If(Ipgg.eq.0) then
            write(6,*) ' Pgg = asb(k)/z*Delta_ns + asb(q)/(1-z) '
		elseif(Ipgg.eq.1) then
		write(6,*) 
     &	' Pgg = asb(k)*B*(1-z)/z*Delta_ns + asb(q)*B*z/(1-z)'
c		write(6,*) ' with BB = ',BB
            elseif(Ipgg.eq.2) then
	     write(6,*) ' Pgg = asb(q)[1/z*Delta_ns + 1/(1-z)]'
            elseif(Ipgg.eq.3) then
		write(6,*) 
     &	' Pgg = asb(q)*[B*(1-z)/z*Delta_ns + B*z/(1-z)]'
		else
c     Ipgg = 0 ! Pgg=(as(k)*1/z +as(q)*1/(1-z))   
c     Ipgg = 1 ! Pgg=(as(k)*B(1-z)/z +as(q)*B*z/(1-z)), Gavins full splitting
c     Ipgg = 2 ! Pgg=(as(q)*1/z +as(q)*1/(1-z)), as(q) eveywhere
c     Ipgg = 3 ! Pgg=(as(q)*B(1-z)/z +as(q)*B*z/(1-z)), as(q) eveywhere + full split
 		write(6,*) ' non - valid splitting fct selected ',Ipgg
		stop
	   Endif
         write(6,*) ' scal factor = ',Qscal,' fact. scale = ',Iqqbar
		
         do j=1,n1
            do i=1,n2
               do k=1,n3
                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS1,RKMS2
                  xx(i) = rx
                  q2x(j) = rq2
                  px(k) = rp
                  xpd2(i,j,k) = max(0.,rkms1)*scal  ! up quarks
                  xpd1(i,j,k) = max(0.,rkms2)*scal  ! down quarks
               enddo
            enddo
         enddo
         write(6,*) ' end of file at ',j,i,k
c we have only valence quarks
         IN=0
         DO I=1,n1
            IN=IN+1
            A(IN) = q2x(I)
         ENDDO
         DO I=1,N2
            IN=IN+1
            A(IN) = xx(I)
         ENDDO
         DO I=1,N3
            IN=IN+1
            A(IN) = px(I)
         ENDDO
         DO I=1,N2
            DO k=1,N3
               DO J=1,N1
                  f_grid1(j,i,k) = xpd1(i,j,k)
                  f_grid2(j,i,k) = xpd2(i,j,k)
               ENDDO
            ENDDO
         ENDDO
         FIRST=.FALSE.
         write(6,*) '  quark densities read from file unit 30 '
         Close(30)
      ENDIF

      ncall = ncall +1 


c      XA(1) = Q2
c      XA(2) = X
c      XA(3) = P
      XA(1) = log(Q2)
      XA(2) = log(X)
      XA(3) = log(P)
      if(xa(3).lt.px(1)) then
	   n3min = n3min + 1
	  
	   if(n3min.lt.10) then 
c            write(6,*) ' p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' p out of range ',xa(3),px(1)
	   elseif(n3min.eq.10) then 
c            write(6,*) ' p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' last message printed: min p'
	   endif
         xa(3)=px(1)
      endif
      if(xa(3).gt.px(n3)) then
	   n3max = n3max + 1
	   if(n3max.lt.10) then 
c            write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
	   elseif(n3max.eq.10) then 
c            write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
c            write(6,*) ' last message printed: max p'
	   endif
         xa(3)=px(n3)
      endif
      if(xa(2).ge.xx(n2)) xa(2)=xx(n2)-0.0001
      if(xa(2).lt.xx(1)) then
	   n2min = n2min + 1
	   if(n2min.lt.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
	   elseif(n2min.eq.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
c            write(6,*) ' last message printed: min x'
	   endif
	   xa(2)=xx(1)
	   endif 
      if(xa(1).lt.q2x(1)) then
	   n1min = n1min + 1
	   if(n1min.lt.10) then 
c            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
	   elseif(n1min.eq.10) then 
c            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
c            write(6,*) ' last message printed: min k2'
	   endif
         xa(1)=q2x(1)
      endif
      if(xa(1).ge.q2x(n1)) then
	   n1max = n1max + 1
         if(n1max.lt.10) then 
c            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
         elseif(n1max.eq.10) then 
c            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
c            write(6,*) ' last message printed: max k2'
	   endif
         xa(1)=q2x(n1)-0.1
      endif
c check if interpolation or grid wanted
      igrid = 0
      if(igrid.eq.1) then
         if(xa(1).lt.q2x(1)) then
            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
            xa(1)=q2x(1)
         endif
         if(xa(1).ge.q2x(n1)) then
            write(6,*) '  k2 out of range: x = ', x,n1,q2x(n1)
            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
            xa(1)=q2x(n1)-0.1
         endif
         i=0
   20    i=i+1
         if(xa(1).gt.a(na(1))) then
c       write(6,*) ' q2  ',xa(1),a(na(1))
            write(6,*) ' q2 not found ',q2,a(na(1)),q2x(n1),xa(1)
            i=na(1)
         else
            if(xa(1).ge.A(i).and.xa(1).lt.a(i+1)) Then
            else
               if(i.le.na(1)) then
                  goto 20
               else
                  write(6,*) ' at 20: q2 not found ',i,q2
               endif
            endif
         endif
         j=0
   30    j=j+1
         if(xa(2).ge.A(na(1)+j).and.xa(2).lt.a(na(1)+j+1)) Then
         else
            if(j.le.na(2)) then
               goto 30
            else
               write(6,*) ' at 30: x not found ',x,xa(2),j
            endif
         endif
         k=0
   40    k=k+1
         if(xa(3).ge.a(na(1)+na(2)+na(3))) then
            k=na(3)
c       write(6,*) ' p  ',xa(3),a(na(1)+na(2)+na(3))
         else
            if(xa(3).ge.A(na(1)+na(2)+k).and. xa(3).lt.a(na(1)+na(2)+k+
     +      1)) Then
            else
               if(k.le.na(3)) then
                  goto 40
               else
                  write(6,*) ' at 40: p not found ',k,p
               endif
            endif
         endif
         quark1 = f_grid1(i,j,k)/q2
         quark2 = f_grid2(i,j,k)/q2
      else
         quark1 = DHFINT(3,XA,NA,A,f_grid1)/q2
         quark2 = DHFINT(3,XA,NA,A,f_grid2)/q2
      endif
      DO  IP=-6,6
         XPQ(IP)=0.0
      ENDDO
      if(quark1.gt.0.) then 
         xpq(1) = quark1
         endif 
      if(quark2.gt.0.) then 
         xpq(2) = quark2
         endif 
      if(quark1.gt.1000000) then
         write(6,*) ' new ccfmstfq x,q2,p,xpq',x,q2,p,quark1,quark2
         write(6,*) ' i,j,k ',i,j,k
      endif
      return


   80 write(6,*) ' error in opening file ccfmstfQ'
   90 write(6,*) ' end of file at ',i,j,k,RX,RQ2,RP,RKMS1,RKMS2
      stop
      END

      SUBROUTINE CCFMSTFU(KF,X,Q2,P,XPQ)
      Implicit None
      Integer n1,n2,n3
      REAL XPQ(-6:6),X,Q2,P
*! new
      Parameter (n1=51,n2=51,n3=51)
c      Real Q2x,xx,px,xpd
c      DIMENSION Q2X(0:n1),XX(0:n2),PX(0:n3),XPD(0:n2,0:n1,0:n3)
      Real Q2x,xx,px,xpd0,xpd1,xpd2,xpd3
      DIMENSION Q2X(0:n1),XX(0:n2),PX(0:n3),XPD0(0:n2,0:n1,0:n3),XPD1(0:n2,0:n1,0:n3),XPD2(0:n2,0:n1,0:n3),XPD3(0:n2,0:n1,0:n3)
      Real amax
      Dimension amax(0:n2)
      Real kl,kt,xpqm
c      Double Precision XA(3),f_grid(n1,n2,n3),A(N1+N2+N3)
      Double Precision XA(3),f_grid0(n1,n2,n3),f_grid1(n1,n2,n3),f_grid2(n1,n2,n3),f_grid3(n1,n2,n3),A(N1+N2+N3)
      INTEGER NA(3)
      DATA NA/n1,n2,n3/
      Double Precision DHFINT
#include "casshwr.inc"
#include "casludat1.inc"
	Integer  ikincut,Ipgg,ns_sel
	Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Integer IGLU
      Common/CAGLUON/Iglu
	Double Precision Qscal 
	Integer Iqqbar
	Common/GLUDAT2/Qscal,Iqqbar 
      Integer i,j,k,in,irr,igrid,ip,n2test
      Integer ipoint,kf,it
	Integer Iset,Iccfmset
c      Real scal,rx,rq2,rp,rkms,glu,xpqmax,xm
      Double Precision scal,rx,rq2,rp,glu,uval,dval,ssea,xpqmax,xm
      Double Precision RKMS4B,RKMS3B,RKMS2B,RKMS1B,RKMS0,RKMS1,RKMS2,RKMS3,RKMS4
      Real xl,pl
	Integer n1min,n1max,n2min,n2max,n3min,n3max,ncall
	common/caerrstf/n1min,n1max,n2min,n2max,n3min,n3max,ncall
      Real QCDLam
	character *72 TXT
	character adum
      LOGICAL FIRST,afirst
      DATA FIRST/.TRUE./
      Data Afirst/.TRUE./
      Logical First_dum,Fccfm1,Fccfm2
      Common/ myfirst/First_dum,Fccfm1,Fccfm2

      Integer Itheory
      Common/theory/Itheory
      
      First=Fccfm1
      Afirst=Fccfm2
      
	if(iabs(kf).eq.2212) then 
	else
         write(6,*) ' ccfmstfu: no valid kf = ',kf
         write(6,*) ' ccfmstfu: only applicable for proton' 	   
	   stop
	endif
      IF(FIRST) THEN
	   n1min = 0
	   n1max = 0
	   n2min = 0
	   n2max = 0
	   n3min = 0
	   n3max = 0
	   ncall = 0
         i=0
         scal = 1.0
         Ipgg = 0
	   ns_sel = -1
	   If(IGLU.GT.1000) then
	   ISET = Int(IGLU/1000)
	   ELSE
         ISET = IGLU
	   ENDIF
         if(iset.eq.1) then
            IF(ICCFM.EQ.0) THEN
cc               open(30,FILE='dglap.dat', FORM='formatted',STATUS=
cc     +         'OLD', IOSTAT=IRR,ERR=80 )
cc               write(6,*) ' read unintegrated gluon from file: '
cc     +         //'dglap.dat'
		   if(iglu.lt.1000) then
                    write(6,*) ' reading: Cascade/updfGrids/dglap.dat'
                    open(30,FILE='Cascade/updfGrids/dglap.dat', FORM='formatted',STATUS=
     +               'OLD', IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: dglap.dat'
               else
	           Iccfmset =  MOD(IGLU,1000)
	           write(6,*) 'Iccfmset',Iccfmset
	            if(iccfmset.eq.111) then
                    write(6,*) ' reading: ccfm-test.dat'
                    open(30,FILE='ccfm-test.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-test.dat'
	            elseif(iccfmset.eq.113) then
                    write(6,*) ' reading: ICCFM=0 for itheory=103 ccfm-grid.dat'
                    open(30,FILE='Cascade/updfGrids/ccfm-grid.dat', FORM='formatted',
     +              STATUS= 'OLD',IOSTAT=IRR,ERR=50 )
                    write(6,*) ' read: for itheory=103 ccfm-grid.dat'
                  else
                    write(6,*) ' reading: Cascade/updfGrids/dglap.dat'
                    open(30,FILE='Cascade/updfGrids/dglap.dat', FORM='formatted',STATUS=
     +               'OLD', IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: dglap.dat'
                  endif
               endif
            ELSEIF(ICCFM.EQ.1) THEN
		   if(iglu.lt.1000) then
               open(30,FILE='Cascade/updfGrids/ccfm.dat', FORM='formatted',STATUS= 'OLD',
     +         IOSTAT=IRR,ERR=80 )
               write(6,*) ' read unintegrated gluon from file: '
     +         //'ccfm.dat'
		   else
	         Iccfmset =  MOD(IGLU,1000)
	         write(6,*) 'Iccfmset',Iccfmset
	            if(iccfmset.eq.1) then
                    open(30,FILE='Cascade/updfGrids/ccfm-ipgg0ns-1-0-all-lhc.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-ipgg0ns-1-0-all-lhc.dat'
	            elseif(iccfmset.eq.2) then
                    open(30,FILE='Cascade/updfGrids/ccfm-ipgg1ns1-0-all-lhc.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-ipgg1ns1-0-all-lhc.dat'
	            elseif(iccfmset.eq.3) then
                    open(30,FILE='Cascade/updfGrids/ccfm-ipgg2ns2-0-all-lhc.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-ipgg2ns2-0-all-lhc.dat'
		      elseif(iccfmset.eq.10) then
                    open(30,FILE='Cascade/updfGrids/ccfm-final-setA0.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-final-setA0.dat'
	            elseif(iccfmset.eq.11) then
                    open(30,FILE='Cascade/updfGrids/ccfm-final-setA0+.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-final-setA0+.dat'
	            elseif(iccfmset.eq.12) then
                    open(30,FILE='ccfm-final-setA0-.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-final-setA0-.dat'
	            elseif(iccfmset.eq.13) then
                    open(30,FILE='Cascade/updfGrids/ccfm-final-setA1.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-final-setA1.dat'
		      elseif(iccfmset.eq.20) then
                    open(30,FILE='Cascade/updfGrids/ccfm-final-setB0.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-final-setB0.dat'
	            elseif(iccfmset.eq.21) then
                    open(30,FILE='ccfm-final-setB0+.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-final-setB0+.dat'
 	            elseif(iccfmset.eq.22) then
                    open(30,FILE='ccfm-final-setB0-.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-final-setB0-.dat'
	            elseif(iccfmset.eq.23) then
                    open(30,FILE='Cascade/updfGrids/ccfm-final-setA1.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-final-setB1.dat'
	            elseif(iccfmset.eq.101) then
                    open(30,FILE='ccfm-setC.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-setC.dat'
	            elseif(iccfmset.eq.111) then
                    open(30,FILE='ccfm-test.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-test.dat'
	            elseif(iccfmset.eq.113) then		     
                    open(30,FILE='Cascade/updfGrids/ccfm-grid.dat', FORM='formatted',
     +              STATUS= 'OLD',IOSTAT=IRR,ERR=50 )
                    write(6,*) ' read: for itheory=103 ccfm-grid.dat'                 
                  endif		   
               endif		   
            ENDIF
         elseif(iset.eq.6) then
            open(30,FILE='Cascade/updfGrids/kmr.dat', FORM='formatted',STATUS= 'OLD',
     +      IOSTAT=IRR,ERR=80 )
            write(6,*) ' read unintegrated gluon from file: '
     +      //'kmr.dat'
         else
            write(6,*) ' No valid IGLU selected for CCFMSTFU: IGLU = ',
     +      Iset
            STOP
         endif
c         read(30,10000) Qg0,ikincut
c10000    format(' Qg0 = ',f12.8,' ikincut= ',I6)
c         read(30,10100)
c10100    format('xg,  kt, p  xgx') 
200      Read(30,101) TXT
  101    Format(A72)
c         WRITE(6,101) '  line ',TXT
         If(TXT(1:4).EQ.'  Qg') then 
	      read(txt,1000) adum,Qg0,adum,ikincut
C1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
1000        format(A7,f12.8,A10,I6)
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Qg0') then 
	      read(txt,1000) adum,Qg0,adum,ikincut
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Ipg') then 
	      read(txt,1001) adum,Ipgg,adum,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
1001        format(A8,I4,A10,I4)
c         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:6).EQ.' Qscal') then 
	      read(txt,1002)  adum,Qscal, adum,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
1002        format(A9,f7.3,A10,I4)
C         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:7).EQ.' QCDlam') then 
c         WRITE(6,101) '3rd line ',TXT
            read(txt,1003) adum,QCDLam
1003        format(A9,f12.8)
            PARU(112)=QCDLam
            goto 200
         Endif
         If(TXT(1:4).EQ.' ln(') then 
c         WRITE(6,101) '2 or 3rd line',TXT
         Endif
	   
         write(6,*) ' soft cut Q0 ',Qg0,' scale factor = ',scal
         write(6,*) ' kin cut ',ikincut,' Ipgg = ',Ipgg,
     &	   ' ns_sel = ',ns_sel
         write(6,*) ' QCD_lam used in uPDF: ',PARU(112)
         write(6,*) ' type of splitting function and non-sudakov: '
	   If(Ipgg.eq.0) then
            write(6,*) ' Pgg = asb(k)/z*Delta_ns + asb(q)/(1-z) '
		elseif(Ipgg.eq.1) then
		write(6,*) 
     &	' Pgg = asb(k)*B*(1-z)/z*Delta_ns + asb(q)*B*z/(1-z)'
c		write(6,*) ' with BB = ',BB
            elseif(Ipgg.eq.2) then
	     write(6,*) ' Pgg = asb(q)[1/z*Delta_ns + 1/(1-z)]'
            elseif(Ipgg.eq.3) then
		write(6,*) 
     &	' Pgg = asb(q)*[B*(1-z)/z*Delta_ns + B*z/(1-z)]'
		else
c     Ipgg = 0 ! Pgg=(as(k)*1/z +as(q)*1/(1-z))   
c     Ipgg = 1 ! Pgg=(as(k)*B(1-z)/z +as(q)*B*z/(1-z)), Gavins full splitting
c     Ipgg = 2 ! Pgg=(as(q)*1/z +as(q)*1/(1-z)), as(q) eveywhere
c     Ipgg = 3 ! Pgg=(as(q)*B(1-z)/z +as(q)*B*z/(1-z)), as(q) eveywhere + full split
 		write(6,*) ' non - valid splitting fct selected ',Ipgg
		stop
	   Endif
         write(6,*) ' scal factor = ',Qscal,' fact. scale = ',Iqqbar
		
         do j=1,n1
            do i=1,n2
               do k=1,n3
cc                  write(6,*) j,i,k,RX,RQ2,RP,RKMS
cc                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS
cc                  write(6,*) ' in ccfmstfu ',iccfmset
                  if(iccfmset.lt.111) then
                     READ(30,*,Err=90 ) RX,RQ2,RP,RKMS0
                  else
                    READ(30,*,Err=90 ) RX,RQ2,RP,RKMS4b,RKMS3b,RKMS2b,RKMS1b,RKMS0,RKMS1,RKMS2,RKMS3,RKMS4
                  endif
                  xx(i) = rx
                  q2x(j) = rq2
                  px(k) = rp
c                  xpd(i,j,k) = max(0.,rkms)*scal
                  xpd0(i,j,k) = max(0.,rkms0)*scal
                  xpd1(i,j,k) = max(0.,rkms1)*scal
                  xpd2(i,j,k) = max(0.,rkms2)*scal
                  xpd3(i,j,k) = max(0.,rkms3)*scal
               enddo
            enddo
         enddo
         write(6,*) ' end of file at ',j,i,k
c we have only gluons
         DO 10 J=1,N2
            DO 10 I=1,N1
               DO 10 k=1,N3

   10    CONTINUE
         IN=0
         DO I=1,n1
            IN=IN+1
            A(IN) = q2x(I)
         ENDDO
         DO I=1,N2
            IN=IN+1
            A(IN) = xx(I)
         ENDDO
         DO I=1,N3
            IN=IN+1
            A(IN) = px(I)
         ENDDO
         DO I=1,N2
            DO k=1,N3
               DO J=1,N1
c                  f_grid(j,i,k) = xpd(i,j,k)
                  f_grid0(j,i,k) = xpd0(i,j,k)
                  f_grid1(j,i,k) = xpd1(i,j,k)
                  f_grid2(j,i,k) = xpd2(i,j,k)
                  f_grid3(j,i,k) = xpd3(i,j,k)
               ENDDO
            ENDDO
         ENDDO
         FIRST=.FALSE.
         Fccfm1=.FALSE.
         write(6,*) '  parton densities read from file unit 30 '
         Close(30)
      ENDIF

      ncall = ncall +1 


c      XA(1) = Q2
c      XA(2) = X
c      XA(3) = P
      XA(1) = log(Q2)
      XA(2) = log(X)
      XA(3) = log(P)
      if(xa(3).lt.px(1)) then
	   n3min = n3min + 1
	  
	   if(n3min.lt.10) then 
c            write(6,*) ' p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' p out of range ',xa(3),px(1)
	   elseif(n3min.eq.10) then 
c            write(6,*) ' p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' last message printed: min p'
	   endif
         xa(3)=px(1)
      endif
      if(xa(3).gt.px(n3)) then
	   n3max = n3max + 1
	   if(n3max.lt.10) then 
c            write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
	   elseif(n3max.eq.10) then 
c            write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
c            write(6,*) ' last message printed: max p'
	   endif
         xa(3)=px(n3)
      endif
      if(xa(2).ge.xx(n2)) xa(2)=xx(n2)-0.0001
      if(xa(2).lt.xx(1)) then
	   n2min = n2min + 1
	   if(n2min.lt.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
	   elseif(n2min.eq.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
c            write(6,*) ' last message printed: min x'
	   endif
	   xa(2)=xx(1)
	   endif 
      if(xa(1).lt.q2x(1)) then
	   n1min = n1min + 1
	   if(n1min.lt.10) then 
c            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
	   elseif(n1min.eq.10) then 
c            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
c            write(6,*) ' last message printed: min k2'
	   endif
         xa(1)=q2x(1)
      endif
      if(xa(1).ge.q2x(n1)) then
	   n1max = n1max + 1
         if(n1max.lt.10) then 
c            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
         elseif(n1max.eq.10) then 
c            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
c            write(6,*) ' last message printed: max k2'
	   endif
         xa(1)=q2x(n1)-0.1
      endif
c check if interpolation or grid wanted
      igrid = 0
      if(igrid.eq.1) then
         if(xa(1).lt.q2x(1)) then
            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
            xa(1)=q2x(1)
         endif
         if(xa(1).ge.q2x(n1)) then
            write(6,*) '  k2 out of range: x = ', x,n1,q2x(n1)
            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
            xa(1)=q2x(n1)-0.1
         endif
         i=0
   20    i=i+1
         if(xa(1).gt.a(na(1))) then
c       write(6,*) ' q2  ',xa(1),a(na(1))
            write(6,*) ' q2 not found ',q2,a(na(1)),q2x(n1),xa(1)
            i=na(1)
         else
            if(xa(1).ge.A(i).and.xa(1).lt.a(i+1)) Then
            else
               if(i.le.na(1)) then
                  goto 20
               else
                  write(6,*) ' at 20: q2 not found ',i,q2
               endif
            endif
         endif
         j=0
   30    j=j+1
         if(xa(2).ge.A(na(1)+j).and.xa(2).lt.a(na(1)+j+1)) Then
         else
            if(j.le.na(2)) then
               goto 30
            else
               write(6,*) ' at 30: x not found ',x,xa(2),j
            endif
         endif
         k=0
   40    k=k+1
         if(xa(3).ge.a(na(1)+na(2)+na(3))) then
            k=na(3)
c       write(6,*) ' p  ',xa(3),a(na(1)+na(2)+na(3))
         else
            if(xa(3).ge.A(na(1)+na(2)+k).and. xa(3).lt.a(na(1)+na(2)+k+
     +      1)) Then
            else
               if(k.le.na(3)) then
                  goto 40
               else
                  write(6,*) ' at 40: p not found ',k,p
               endif
            endif
         endif
         glu = f_grid0(i,j,k)/q2
      else
c         glu = DHFINT(3,XA,NA,A,f_grid)/q2
         glu = DHFINT(3,XA,NA,A,f_grid0)/q2
         uval = DHFINT(3,XA,NA,A,f_grid1)/q2
         dval = DHFINT(3,XA,NA,A,f_grid2)/q2
         ssea = DHFINT(3,XA,NA,A,f_grid3)/q2
      endif
      DO  IP=-6,6
         XPQ(IP)=0.0
      ENDDO
      if(glu.gt.0.) xpq(0) = glu
      if(uval.gt.0.) xpq(1) = uval
      if(dval.gt.0.) xpq(2) = dval
      if(ssea.gt.0.) xpq(3) = ssea
      if(ssea.gt.0.) xpq(-3) = ssea
      if(glu.gt.5000000) then
         
ccc         write(6,*) ' new ccfmstfu x,q2,p,xpq',x,q2,p,glu
ccc         write(6,*) ' i,j,k ',i,j,k
      endif
      return
      ENTRY ccfmmax(x,xpqmax)

      xpqmax = 0.
c first calculate maxima and store it in vector amax
      if(Afirst) then
         write(6,*) ' calculate now maxima '
         DO I=1,N2
            n2test = I
            xpqmax = 0.
            do it=n2test,n2
               DO k=1,N3
                  DO J=1,N1
                     xm=f_grid0(j,it,k)/exp(q2x(j))
c              write(6,*) f_grid(j,it,k),exp(q2x(j))
                     if(xm.gt.xpqmax) xpqmax=xm
                  Enddo
               Enddo
            Enddo
            amax(i)=xpqmax
         Enddo
         Afirst=.false.
         Fccfm2=Afirst
         write(6,*) ' finished with maxima '
      endif
      xl = log(x)
      DO I=1,N2
         If(xx(I).gt.xl) then
            ipoint=i
            goto 50
         Endif
      Enddo
   50 continue
      xpqmax = amax(ipoint)

      return

      ENTRY ccfmmin(x,kt,p,xpqm)
      xpqm = 100000.
      xl=log(x)
      kl=log(kt)
      pl=log(p)

      if(kl.lt.q2x(1)) then
c         write(6,*) '  k2 out of range ',kt,' min ',exp(q2x(1))
         kl=q2x(1)
      endif
      if(kl.ge.q2x(n1)) then
c         write(6,*) ' ccfmmin: k2 out of range: p =',p,' x= ',x,n1
c         write(6,*) ' ccfmmin: k2 out of range ',kt,' max ',exp(q2x(n1))
         kl=q2x(n1)-0.1
      endif
      i=0
   60 i=i+1
      if(kl.gt.a(na(1))) then
c       write(6,*) ' q2  ',xa(1),a(na(1))
         write(6,*) ' q2 not found ',q2,a(na(1)),q2x(n1),xa(1)
         i=na(1)
      else
         if(kl.ge.A(i).and.kl.lt.a(i+1)) Then
         else
            if(i.le.na(1)) then
               goto 60
            else
               write(6,*) ' ccfmmin q2 not found ',i,kt
            endif
         endif
      endif
      j=0
   70 j=j+1
      if(xl.ge.A(na(1)+j).and.xl.lt.a(na(1)+j+1)) Then
      else
         if(j.le.na(2)) then
            goto 70
         else
            write(6,*) ' ccfmmin x not found: set j=na(2) ',x,xl,j
		j = na(2)
         endif
      endif
      DO K=1,N3
         if(pl.ge.px(k)) then
            xm=max(f_grid0(i,j,k)/kt,0.005)
            if(xm.lt.xpqm) xpqm=xm
         endif
      Enddo
c        write(6,*) ' ccfm x,p,xpqm ',x,p,xpqm
      return


   80 write(6,*) ' CCFMSTFU: error in opening file ',iccfmset
   90 write(6,*) ' end of file at ',i,j,k,RX,RQ2,RP,RKMS0
      stop
      END

      FUNCTION PYR(IDUMMY)

      double precision pyr,dcasrn
	external dcasrn
      Logical first
	data first/.true./
	if(first) then
	  first=.false.
	  write(6,*) ' cascade version of pyr used...'
	  endif
      PYR= dcasrn(idummy)

      RETURN
      END

      DOUBLE PRECISION FUNCTION dcasrn(i)

      DOUBLE PRECISION DR(1)
      DOUBLE PRECISION TWOM24
      PARAMETER (TWOM24 = 2.D0**(-24))
      Real Rvec(2)
	Integer i
c      CALL RANLUX(RVEC,2)
c      DR(1) = DBLE (RVEC(1)) + DBLE (RVEC(2)) * TWOM24
      CALL RM48 (DR,1)
      dcasrn = DR(1)


      RETURN
      END


      SUBROUTINE dcasrnv(RVEC,LENV)

      Integer lenv,I
      Double Precision RVEC(LENV)
      Double Precision dcasrn,t
      Do I=1,Lenv
         t=dcasrn(i)
         if(t.ne.t) then
            write(6,*) ' DCASRNv: serious problem with Rn = ',t
         elseif(t.gt.1d0.or.t.lt.0.d0) then
            write(6,*) ' DCASRNv: out of range Rn = ',t
         endif
         Rvec(i)=t
      enddo

      RETURN
      END

C call cascade to calc total xsection

	Subroutine Sigcalc(xn,q2n,f2,f2e)
      Implicit None
	Integer N1,N2
	Integer K1,K2,LUX
#include "caslujets.inc"
#include "casluco.inc"
#include "casparas.inc"
#include "casvalues.inc"
#include "caspart.inc"
#include "caseffic.inc"
#include "casparam.inc"
#include "casshwr.inc"
#include "casludat1.inc"
      Integer Iglu
      Common/CAGLUON/Iglu

      Integer I
      Integer Iseed
      Common/random/ISEED
	
	Integer Nevent
	Common/steer1/Nevent

	Logical First,Fccfm1,Fccfm2
	Common/ myfirst/First,Fccfm1,Fccfm2
      logical lcalc
	Common/ myfirst2/lcalc
      Integer Itheory
      Common/theory/Itheory
      Integer IcasHF
      Common/CasHF/IcasHF
	Integer IPROO
	real timeleft
	Integer Minuts
	External Minuts
      Double Precision DOT
      Double Precision xn,q2n,f2,f2e,alph_em
	Double Precision GEV2NB
      DATA GEV2NB/.3893857E+6/
      
      Real xm_test
      
	Integer haprint
	common/haprint/haprint
      Double Precision PARU112
      
      PARU112 = PARU(112)


      f2 = 0
      f2e = 0
      
      LUX = 3
      K1=0
      K2=0
cc      CALL RLUXGO(LUX,ISEED,K1,K2)

      haprint = 1

      n1=0
      n2=0
      CALL RM48IN(ISEED,N1,N2)
c      write(6,*) ' cascade rm48 initialised from ',iseed


C initialise PYTHIA 6 parameters
ccc      CALL GPYINI
C initialise CASCADE parameters
      CALL CASINI

      CALL STEER_new
	
C-- change standard parameters of CASCADE
   
	Call cascha
c
C-- change standard parameters of JETSET/PYTHIA
      Call pytcha
c      write(6,*) ' paru(112) ',PARU(112),PARU112
      IF(PARU112.ge.0) PARU(112) = PARU112
c
c Start integration of x-section
c
c      CALL PTIME(' cascade  ',1,0)
ccc	Fccfm1 = .true. Now set xfitter
	Fccfm2 = .true.
	First = .true.
      
c iglu is set by steer-ep  for itheory=102

      if(itheory.eq.103.or.itheory.eq.104) Iglu = 1113
      if(itheory.eq.101) Iglu = 1111
      if(itheory.eq.105) Iglu = 1111
     
c      write(6,*) ' sigcalc itheory = ',itheory,' iglu = ',iglu

c     here now calcualte y,q2 from x,q2
      P(1,1) = 0.0D0
      P(1,2) = 0.0D0
      P(1,3) = DBLE(PLEPIN)
      P(1,5) = PYMASS(KE)
      P(1,4) = DSQRT(P(1,1)**2+P(1,2)**2+P(1,3)**2+P(1,5)**2)
      P(2,1) = 0.0D0
      P(2,2) = 0.0D0
      P(2,3) = DBLE(PPIN)
      P(2,5) = PYMASS(KP)
      P(2,4) = DSQRT(P(2,1)**2+P(2,2)**2+P(2,3)**2+P(2,5)**2)
      DO 20 I=1,4
         CM(I) =P(1,I)+P(2,I)
   20 CONTINUE
      SSS = DOT(CM,CM)
      qmi = q2n
      qMA = qmi 
      
      ymi = q2n/xn/sss
      yma=ymi  
 
            
      if(ymi.gt.1.or.q2n.gt.sss) goto 345
        
c      write(6,*) ' after x,y calc ',xn,yma,q2n,sqrt(sss)
c     check on xg_min = smim/y/sss
      xm_test = 1./ymi/sss
c      write(6,*) xm_test
      if(xm_test.gt.0.99) goto 345
      if(ymi.le.0) goto 345
      if(ymi.gt.0.99) goto 345
c      write(6,*) ' using IGLU = ',IGLU,
c     &  ' for xsection calc in cascade-uitl.F',ymi,yma,qmi,qma
c      write(6,*) ' sigcalc x,y, q2, sqrts ',xn,yma,q2n,sqrt(sss)
      IPROO = IPRO
      if(IcasHF.eq.4) then
         IPRO =11
         IHFLA = 4
         elseif(IcasHF.eq.5) then
         IPRO = 11
         IHFLA = 5 
         else
         IPRO = 10
         NFLAV = 3
      endif
c      write(6,*) ' sigcalc ipro , ihfla ', ipro,nflav,ihfla,icasHF
      CALL CASCADE
c
c Print out result of integration of x-section
c
      CALL CAEND(1)
      
 	ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2n)      
c note in p_semih Fgam=Fwei=Fwei1=1
      f2  = avgi*q2n/GEV2NB/4d0/pi**2/alph_em
	f2e = sd*q2n/GEV2NB/4d0/pi**2/alph_em
cc      write(6,*) ' sigcalc Iccfm ',ICCFM
cc      write(6,*) ' sigcalc Ipro, xn,q2n,ymi,f2: ',ipro,IcasHF,NFlav,xn,qmi,ymi,f2
      IPRO = IPROO
123   Continue
	
      Return
345   continue
      f2 = 0
      f2e =0 
      return      
c	stop
      END
	Subroutine Sigqpm(xn,q2n,f2,f2e)
      Implicit None
	Integer N1,N2
	Integer K1,K2,LUX
#include "caslujets.inc"
#include "casluco.inc"
#include "casparas.inc"
#include "casvalues.inc"
#include "caspart.inc"
#include "caseffic.inc"
#include "casparam.inc"
#include "casludat1.inc"
      Integer Iglu
      Common/CAGLUON/Iglu

      Integer I
      Integer Iseed
      Common/random/ISEED
	
	Integer Nevent
	Common/steer1/Nevent

	Logical First,Fccfm1,Fccfm2
	Common/ myfirst/First,Fccfm1,Fccfm2
      logical lcalc
	Common/ myfirst2/lcalc
      Integer Itheory
      Common/theory/Itheory
	
      Integer IPROOLD
      
	real timeleft
	Integer Minuts
	External Minuts
      Double Precision DOT
      Double Precision xn,q2n,f2,f2e,alph_em
	Double Precision GEV2NB
      DATA GEV2NB/.3893857E+6/
      
      Real xm_test
      
	Integer haprint
	common/haprint/haprint
      Double Precision PARU112
      
      PARU112 = PARU(112)

      LUX = 3
      K1=0
      K2=0
c      CALL RLUXGO(LUX,ISEED,K1,K2)

      haprint = 1

      n1=0
      n2=0
      CALL RM48IN(ISEED,N1,N2)
c      write(6,*) ' cascade rm48 initialised from ',iseed

C initialise PYTHIA 6 parameters
ccc      CALL GPYINI
C initialise CASCADE parameters
      CALL CASINI

      CALL STEER_new
	
C-- change standard parameters of CASCADE
   
	Call cascha
c
C-- change standard parameters of JETSET/PYTHIA
      Call pytcha
      IF(PARU112.ge.0) PARU(112) = PARU112
c
c Start integration of x-section
c
c      CALL PTIME(' cascade  ',1,0)
ccc	Fccfm1 = .true. Now set xfitter
	Fccfm2 = .true.
	First = .true.
c iglu is set by steer-ep  for itheory=102

      if(itheory.eq.103.or.itheory.eq.104) Iglu = 1113
      if(itheory.eq.101) Iglu = 1111
      if(itheory.eq.105) Iglu = 1111
      
c     here now calcualte y,q2 from x,q2
      P(1,1) = 0.0D0
      P(1,2) = 0.0D0
      P(1,3) = DBLE(PLEPIN)
      P(1,5) = PYMASS(KE)
      P(1,4) = DSQRT(P(1,1)**2+P(1,2)**2+P(1,3)**2+P(1,5)**2)
      P(2,1) = 0.0D0
      P(2,2) = 0.0D0
      P(2,3) = DBLE(PPIN)
      P(2,5) = PYMASS(KP)
      P(2,4) = DSQRT(P(2,1)**2+P(2,2)**2+P(2,3)**2+P(2,5)**2)
      DO 20 I=1,4
         CM(I) =P(1,I)+P(2,I)
   20 CONTINUE
      SSS = DOT(CM,CM)
      qmi = q2n
      qMA = qmi 
      
      ymi = q2n/xn/sss
      yma=ymi  
 
            
      if(ymi.gt.1.or.q2n.gt.sss) goto 345
        
c      write(6,*) ' after x,y calc ',xn,yma,q2n,sqrt(sss)
c     check on xg_min = smim/y/sss
      xm_test = 1./ymi/sss
c      write(6,*) xm_test
      if(xm_test.gt.0.99) goto 345
      if(ymi.le.0) goto 345
      if(ymi.gt.0.99) goto 345
      IPROOLD = IPRO
      IPRO = 4  ! switch for qpm
      
c      write(6,*) ' sigqpm using IGLU = ',IGLU,ipro
      CALL CASCADE
c
c Print out result of integration of x-section
c
      CALL CAEND(1)
      
c      write(6,*) ' end of cascade ',avgi,sd
 	ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2n)      
c note in p_semih Fgam=Fwei=Fwei1=1
      f2  = avgi*q2n/GEV2NB/4d0/pi**2/alph_em
	f2e = sd*q2n/GEV2NB/4d0/pi**2/alph_em
c      write(6,*) ' end of cascade ',f2,f2e,pi,q2n,PYALEM(Q2n)
123   Continue
	IPRO = IPROOLD
      Return
345   continue
      f2 = 0
      f2e =0 
      return      
c	stop
      END
 	
      Subroutine xsec1(X,XSEC)
      IMPLICIT None
	Double Precision X(20)
#include "caslujets.inc"
#include "casparas.inc"
#include "caspart.inc"
#include "casluco.inc"
      Integer IDEBUG
      COMMON/ INTERN/IDEBUG
#include "caspara5.inc"
	Integer NCALL,I,J
	Double Precision GEV2NB,XSEC,WEIG,WPART,WT1
      DATA NCALL/0/
      DATA GEV2NB/.3893857E+6/
      DO 10  I=1,20
         DO 10 J=1,5
            K(I,J) = 0
   10 P(I,J) = 0.0D0
      N=0
      XSEC=0.D0
      WEIG=0.D0
      WT1 = 0D0
      WPART=0
	if(x(1).ne.x(1)) then
        write(6,*) 'fatal in xsec ',(x(i),i=1,20)
      endif
c	IDEBUG = 1
c      write(6,*) KBEAM(1,2),KBEAM(2,2)
c      write(6,*) ' calling meoffsh in xsec1 '
      IF(IPRO.ge.30) then
         IF(IPRO.eq.30) call p_gbw(x,wpart)
         if(IPRO.eq.31) call p_new(x,wpart)
         if(IPRO.eq.32) call p_full(x,wpart)
         WEIG=WPART*GEV2NB
      else
      IF((KBEAM(1,2).eq.22.or.Iabs(KBEAM(1,2)).eq.11).and.
     & 	(KBEAM(2,2).eq.22.or.Iabs(KBEAM(2,2)).eq.11)) Then
c         CALL P_gamgam(X,WPART)
c	   write(6,*) ' gamgam called'
	Else
c	   write(6,*) ' before p_semih'
         CALL P_SEMIH(X,WPART)
c	   write(6,*) ' standard called', wpart,ipro
	Endif
      IF(WPART.GT.0.) THEN
         CALL MEOFFSH(WT1) 
      ENDIF
c      write(6,*) ' wpart,wt1', wpart,wt1,ncall
      WEIG=WPART*WT1*GEV2NB
      endif
C BOOST BACK TO OVERALL CMS
c      write(6,*) ' wpart,wt1', IPRO,wpart,wt1,ncall
      ncall = ncall + 1
c      call PYLIST(1)
      XSEC=WEIG
      RETURN
      END

	
	subroutine steer_new
	
	Implicit None
	
      character *132 TXT,TXTS
	character *6 char
	character *10 char2
	Integer i1,i2,i3
	Integer inte,i,j
	Integer Irr
	real r3
	
	
	Integer nloop,nmax
	Parameter (nmax=1000)
      character *4   para
	Real Rval 
	Integer Ival,Ld,Le
	Common/steering/Nloop,Ld(nmax),Le(nmax),Ival(nmax),
     &	para(nmax),Rval(nmax)
	Integer Nevent
	Common/steer1/Nevent
	character *132 filnam
	Integer lunhbk
	Common/steer2/lunhbk,filnam
      LOGICAL         HBKOUT
      COMMON /QHBKLO/ HBKOUT
	
      open(41,FILE='steer-ep', FORM='formatted',STATUS=
     +         'OLD', IOSTAT=IRR,ERR=80 )
	
	nloop = 0
   20 Read(41,101,END=100) TXT
      If(TXT(1:1).EQ.'*') then
c         WRITE(6,103) '*',TXT
         GOTO 20
      Endif
c check for special cards	
      If(TXT(1:6).EQ.'NEVENT') then
         read(txt,*) char,i1
         NEVENT = i1
c	   write(6,*) ' Nevents = ',i1
         GOTO 20
      Endif
      If(TXT(1:6).EQ.'HBKOUT') then
         read(txt,*) char,lunhbk,filnam
c	   write(6,*) char,lunhbk,filnam
	   HBKOUT=.TRUE.
         GOTO 20
      Endif
      
	If(TXT(1:4).EQ.'END$') goto 90
	inte = 0 
      do j=2,132
	i=j-1
	if(TXT(i:i).EQ."'".and.
     &(TXT(J:J).EQ.'H'.OR.TXT(J:J).EQ.'I'.OR.TXT(J:J).EQ.'J'.
     &OR.TXT(J:J).EQ.'K'.OR.TXT(J:J).EQ.'L'.OR.TXT(J:J).EQ.'M'.
     &OR.TXT(J:J).EQ.'N')) then
      inte = 1
	endif
	enddo
	nloop = nloop + 1
      if(nloop.ge.nmax) then
        write(6,*) ' fatal : steer nllop > nmax ',nloop,nmax
        stop
        endif
	Ival(nloop) = -9999
	Rval(nloop) = -9999.
	if(inte.eq.1) then
      read(txt,*) char,i1,i2,i3,char2
c	write(6,*) ' output Integer :',char,i1,i2,i3 
	Ival(nloop) = i3
      else
      read(txt,*) char,i1,i2,r3,char2
c	write(6,*) ' output real    :',char,i1,i2,r3
	Rval(nloop) = r3
	endif
	Ld(nloop) = i1
	Le(nloop) = i2
	para(nloop) = char
	
c	write(6,*) ' output :',nloop,para(nloop),Ival(nloop),rval(nloop)
	
      goto 20
	
   90 Continue
  100 Continue
      close(41)
	rewind(41)
      Return
   80 write(6,*) ' Steer-new: error in opening file steer-ep'
	stop
  101 Format(A132)
      End


C#######################################################################
C
C   One- and two-dimensional adaptive Gaussian integration routines.
C
C **********************************************************************

      SUBROUTINE GADAP(A0,B0,F,EPS,SUM)
	IMPLICIT NONE
C
C   PURPOSE           - INTEGRATE A FUNCTION F(X)
C   METHOD            - ADAPTIVE GAUSSIAN
C   USAGE             - CALL GADAP(A0,B0,F,EPS,SUM)
C   PARAMETERS  A0    - LOWER LIMIT (INPUT,REAL)
C               B0    - UPPER LIMIT (INPUT,REAL)
C               F     - FUNCTION F(X) TO BE INTEGRATED. MUST BE
C                       SUPPLIED BY THE USER. (INPUT,REAL FUNCTION)
C               EPS   - DESIRED RELATIVE ACCURACY. IF SUM IS SMALL EPS
C                       WILL BE ABSOLUTE ACCURACY INSTEAD. (INPUT,REAL)
C               SUM   - CALCULATED VALUE FOR THE INTEGRAL (OUTPUT,REAL)
C   PRECISION         - SINGLE
C   REQ'D PROG'S      - F
C   AUTHOR            - T. JOHANSSON, LUND UNIV. COMPUTER CENTER, 1973
C   REFERENCE(S)      - THE AUSTRALIAN COMPUTER JOURNAL,3 P.126 AUG. -71
C
	Integer Num,Ifu,N,L,I
	Real F,A,B,F1,F2,F3,S,F1F,F2F,F3F,AA,BB,EPS,RED,SUM,C
	Real A0,B0,DSUM,W1,U2,SOLD,SS
      COMMON/GADAP1/ NUM,IFU
      EXTERNAL F
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)

      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)
      IF(EPS.LT.1.0E-8) EPS=1.0E-8
      RED=1.3
      L=1
      I=1
      SUM=0.
      C=SQRT(15.)/5.
      A(1)=A0
      B(1)=B0
      F1(1)=F(0.5*(1+C)*A0+0.5*(1-C)*B0)
      F2(1)=F(0.5*(A0+B0))
      F3(1)=F(0.5*(1-C)*A0+0.5*(1+C)*B0)
      IFU=3
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)
   10 CONTINUE
      L=L+1
      N(L)=3
      EPS=EPS*RED
      A(I+1)=A(I)+C*(B(I)-A(I))
      B(I+1)=B(I)
      A(I+2)=A(I)+B(I)-A(I+1)
      B(I+2)=A(I+1)
      A(I+3)=A(I)
      B(I+3)=A(I+2)
      W1=A(I)+(B(I)-A(I))/5.
      U2=2.*W1-(A(I)+A(I+2))/2.
      F1(I+1)=F(A(I)+B(I)-W1)
      F2(I+1)=F3(I)
      F3(I+1)=F(B(I)-A(I+2)+W1)
      F1(I+2)=F(U2)
      F2(I+2)=F2(I)
      F3(I+2)=F(B(I+2)+A(I+2)-U2)
      F1(I+3)=F(A(I)+A(I+2)-W1)
      F2(I+3)=F1(I)
      F3(I+3)=F(W1)
      IFU=IFU+6
      IF(IFU.GT.5000) GOTO 40
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))
      SS=S(I+1)+S(I+2)+S(I+3)
      I=I+3
      IF(I.GT.300)GOTO 30
      SOLD=S(I-3)
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 10
      SUM=SUM+SS
      I=I-4
      N(L)=0
      L=L-1
   20 CONTINUE
      IF(L.EQ.1) GOTO 40
      N(L)=N(L)-1
      EPS=EPS/RED
      IF(N(L).NE.0) GOTO 10
      I=I-1
      L=L-1
      GOTO 20
   30 WRITE(6,10000)
10000 FORMAT(' GADAP:I TOO BIG')
   40 RETURN
      END
      SUBROUTINE cascha
*#**********************************************************************
*#
*#    SUBROUTINE cascha
*#
*#**********************************************************************
      character *4   param 
*
	Integer nloop,nmax
	Parameter (nmax=1000)
      character *4   para 
	Real Rval 
	Integer Ival,Ld,Le
	Common/steering/Nloop,Ld(nmax),Le(nmax),Ival(nmax),
     &	para(nmax),Rval(nmax)

#include "casparas.inc"
#include "casluco.inc"
#include "caspara5.inc"
#include "caspart.inc"
#include "casvalues.inc"
#include "casshwr.inc"
#include "casprkt.inc"
#include "casscal.inc"
#include "casglu.inc"
#include "casmult.inc"
#include "casherup.inc"

      REAL ACC1,ACC2
      COMMON /INTEGR/ ACC1,ACC2,IINT,NCB
      Integer IGLU
	Common/CAGLUON/Iglu
	Double Precision Qscal
	Integer Iqqbar
	Common/GLUDAT2/Qscal,Iqqbar

C---Determine the machine accuacy EPSMAC (D= 2.E-06)
      EPSMAC = 0.5
      DO 10 I= 1, 100
         EPSMAC = EPSMAC * 0.5
         IF ((1.0+EPSMAC) .EQ. 1.0) GO TO 20
   10 CONTINUE
      EPSMAC = 1.0E-6
   20 EPSMAC = 2.0 * EPSMAC
c      WRITE(6,10000) EPSMAC
10000 FORMAT(1X,/,' computed machine accuracy = ',G10.2,/)

C---Interpret the text containing the parameters to be changed
      If(Nloop.ge.1) Then
	      Do 30 I=1,Nloop

            PARAM = Para(I)
            INT = Ival(I)
            REA = Rval(I)

            IF (PARAM.EQ.'KBE1') THEN
               IDEF = KE
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  KE = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRE1') THEN
               IDEF = Ires(1)
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Ires(1) = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'KBE2') THEN
               IDEF = KP
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  KP = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRE2') THEN
               IDEF = Ires(2)
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Ires(2) = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'NFRA') THEN
               IDEF = NFRAG
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  NFRAG = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IPST') THEN
               IDEF = IPST
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IPST = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IFPS') THEN
               IDEF = IFPS
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IFPS = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ITIM') THEN
               IDEF = ITIMSHR
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ITIMSHR = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRAS') THEN
               IDEF = IRUNA
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRUNA = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRAM') THEN
               IDEF = IRUNAEM
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRUNAEM = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IQ2S') THEN
               IDEF = IQ2
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IQ2 = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IPRO') THEN
               IDEF = IPRO
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IPRO = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'NFLA') THEN
               IDEF = NFLAV
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  NFLAV = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'INTE') THEN
               IDEF = INTER
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  INTER = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'IHFL') THEN
               IDEF = IHFLA
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IHFLA = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF


            IF (PARAM.EQ.'ISEH') THEN
               IDEF = ISEMIH
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ISEMIH = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ICCF') THEN
               IDEF = ICCFM
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ICCFM = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IGLU') THEN
               IDEF = IGLU
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IGLU = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ICGL') THEN
               IDEF = ICOLGLU
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ICOLGLU = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IINT') THEN
               IDEF= IINT
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IINT = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'NCAL') THEN
               IDEF= NCB
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  NCB = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IREM') THEN
               IDEF= IRspl
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRspl = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'MULT') THEN
               IDEF= Mult
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Mult = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'ILHA') THEN
               IDEF = Ilha
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Ilha = INT
c                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
c                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
           IF (PARAM.EQ.'PBE2') THEN
               RDEF = PPIN
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PPIN = REA
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'PBE1') THEN
               RDEF = PLEPIN
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PLEPIN = REA
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'QMIN') THEN
               RDEF = SNGL(QMI)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  QMI = DBLE(REA)
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'QMAX') THEN
               RDEF = SNGL(QMA)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  QMA = DBLE(REA)
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'YMIN') THEN
               RDEF = SNGL(YMI)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  YMI = DBLE(REA)
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'YMAX') THEN
               RDEF = SNGL(YMA)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  YMA = DBLE(REA)
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'THEM') THEN
               RDEF = SNGL(THEMA)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  THEMA = DBLE(REA)
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'THEI') THEN
               RDEF = SNGL(THEMI)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  THEMI = DBLE(REA)
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'PTCU') THEN
               RDEF = SNGL(PT2CUT(IPRO))
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PT2CUT(IPRO) = DBLE(REA)
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF


            IF (PARAM.EQ.'ACC1') THEN
               RDEF = ACC1
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  ACC1 = REA
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ACC2') THEN
               RDEF = ACC2
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  ACC2 = REA
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'XKT1') THEN
               RDEF = SNGL(PRKT1)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PRKT1 = DBLE(REA)
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'XKTR') THEN
               RDEF = SNGL(PRKTREM)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PRKTREM = DBLE(REA)
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'XKT2') THEN
               RDEF = SNGL(PRKT2)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PRKT2 = DBLE(REA)
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'XKTR') THEN
               RDEF = SNGL(PRKTREM)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PRKTREM = DBLE(REA)
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'SCAL') THEN
               RDEF = SNGL(SCALFA)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  SCALFA = DBLE(REA)
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'SCAF') THEN
               RDEF = SNGL(SCALFAF)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  SCALFAF = DBLE(REA)
c                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
c                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF

   30    CONTINUE
      ENDIF

10100 FORMAT(' cascha: ',A4,'(',I4,')',
     +       ' changed from the default ',I8,' TO',I8)
10200 FORMAT(' cascha (WARNING):',A4,'(',I4,')',
     +       ' not changed from the default ',I8,' TO',I8)
10300 FORMAT(' cascha: ',A4,'(',I4,')',
     +       ' changed from the default ',E10.3,' TO',E10.3)
10400 FORMAT(' cascha (WARNING): ',A4,'(',I4,')',
     +       ' not changed from the default ',E10.3,' TO',E10.3)

      RETURN
      END

      SUBROUTINE PYTCHA

*#**********************************************************************
*#
*#    SUBROUTINE PYTCHA
*#
*# PURPOSE: Interpret the seering to change the default
*#          parameters of PYTHIA 5.7.
*#          note: only those parameters can be changed which the user
*#                may change according to the PYTHIA manual.
*#
*#**********************************************************************
      character *4   param 
*
	Integer nloop,nmax
	Parameter (nmax=1000)
      character *4   para 
	Real Rval 
	Integer Ival,Ld,Le
	Common/steering/Nloop,Ld(nmax),Le(nmax),Ival(nmax),
     &	para(nmax),Rval(nmax)


      Double Precision paru,parj,pmas,parf,vckm,brat
      CHARACTER CHAF*16
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
 	COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT4/CHAF(500,2)

C---additional PYTHIA commons (not part of standard PYTHIA)
      CHARACTER*8     FRAME,BEAM,TARGET
      Double Precision        win,pbeam,ptarge
      COMMON /PBEAM1/ FRAME,BEAM,TARGET
      COMMON /PBEAM2/ WIN,PBEAM(3),PTARGE(3)
      SAVE /PBEAM1/,/PBEAM2/

      Double Precision ckin,parp,pari
C---PYTHIA commons
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      SAVE /PYSUBS/,/PYPARS/

      CHARACTER*4     CHAINT


C---Determine the machine accuacy EPSMAC (D= 2.E-06)
      EPSMAC = 0.5
      DO 33 I= 1, 100
      EPSMAC = EPSMAC * 0.5
      IF ((1.0+EPSMAC) .EQ. 1.0)  GO TO 35
   33 CONTINUE
      EPSMAC = 1.0E-6
   35 EPSMAC = 2.0 * EPSMAC
c      WRITE(6,330) EPSMAC
 330  FORMAT(1X,/,' COMPUTED MACHINE ACCURACY = ',G10.2,/)

C---Interpret the text containing the parameters to be changed
      If(Nloop.ge.1) Then
	      Do 30 I=1,Nloop

            LND = Ld(I)
            LNE = Le(I)

            PARAM = Para(I)
            INT = Ival(I)
            REA = Rval(I)

            IF (PARAM.EQ.'MSEL') THEN
                IDEF = MSEL
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    MSEL = INT
c                    PRINT 91001, PARAM,LND,IDEF,INT
                ELSE
c                    PRINT 92001, PARAM,LND,IDEF,INT
                ENDIF
            ENDIF

            IF (PARAM.EQ.'MSUB') THEN
                IDEF = MSUB(LND)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    MSUB(LND) = INT
c                    PRINT 91001, PARAM,LND,IDEF,INT
                ELSE
c                    PRINT 92001, PARAM,LND,IDEF,INT
                ENDIF
            ENDIF

            IF (PARAM.EQ.'MSTP') THEN
                IDEF = MSTP(LND)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    MSTP(LND) = INT
c                    PRINT 91001, PARAM,LND,IDEF,INT
                ELSE
c                    PRINT 92001, PARAM,LND,IDEF,INT
                ENDIF
            ENDIF

            IF (PARAM.EQ.'CKIN') THEN
                RDEF = CKIN(LND)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    CKIN(LND) = REA
c                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
c                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'WIN ') THEN
                RDEF = WIN
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    WIN = dble(REA)
c                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
c                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'PBEA') THEN
                RDEF = PBEAM(LND)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    PBEAM(LND) = REA
c                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
c                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'PTAR') THEN
                RDEF = PTARGE(LND)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    PTARGE(LND) = REA
c                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
c                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'PARP') THEN
                RDEF = PARP(LND)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    PARP(LND) = dble (REA)
c                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
c                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'KFIN') THEN
                IDEF = KFIN(LND,LNE)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    KFIN(LND,LNE) = REA
c                    PRINT 91003, PARAM,LND,LNE,IDEF,INT
                ELSE
c                    PRINT 92003, PARAM,LND,LNE,IDEF,INT
                ENDIF
            ENDIF


            IF (PARAM.EQ.'MSTU') THEN
                IDEF = MSTU(LND)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    MSTU(LND) = INT
c                    PRINT 91001, PARAM,LND,IDEF,INT
                ELSE
c                    PRINT 92001, PARAM,LND,IDEF,INT
                ENDIF
            ENDIF

            IF (PARAM.EQ.'MSTJ') THEN
                IDEF = MSTJ(LND)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    MSTJ(LND) = INT
c                    PRINT 91001, PARAM,LND,IDEF,INT
                ELSE
c                    PRINT 92001, PARAM,LND,IDEF,INT
                ENDIF
            ENDIF

            IF (PARAM.EQ.'PARU') THEN
                RDEF = PARU(LND)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    PARU(LND) = dble(REA)
c                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
c                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'PARJ') THEN
                RDEF = PARJ(LND)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    PARJ(LND)= dble(REA)
c                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
c                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'KCHG') THEN
                IDEF = KCHG(LND,LNE)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    KCHG(LND,LNE) = INT
c                    PRINT 91003, PARAM,LND,LNE,IDEF,INT
                ELSE
c                    PRINT 92003, PARAM,LND,LNE,IDEF,INT
                ENDIF
            ENDIF

            IF (PARAM.EQ.'PMAS') THEN
                RDEF = PMAS(LND,LNE)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    PMAS(LND,LNE)=dble(REA)
c                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
c                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'PARF') THEN
                RDEF = PARF(LND)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    PARF(LND) = dble(REA)
c                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
c                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'VCKM') THEN
                RDEF = VCKM(LND,LNE)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    VCKM(LND,LNE) = dble(REA)
c                    PRINT 91002, PARAM,LND,LNE,RDEF,REA
                ELSE
c                    PRINT 92002, PARAM,LND,LNE,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'MDCY') THEN
                IDEF = MDCY(LND,LNE)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    MDCY(LND,LNE) = INT
c                    PRINT 91004, PARAM,LND,LNE,IDEF,INT
                ELSE
c                    PRINT 92004, PARAM,LND,LNE,IDEF,INT
                ENDIF
            ENDIF

            IF (PARAM.EQ.'MDME') THEN
                IDEF = MDME(LND,LNE)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    MDME(LND,LNE) = INT
c                    PRINT 91004, PARAM,LND,LNE,IDEF,INT
                ELSE
c                    PRINT 92004, PARAM,LND,LNE,IDEF,INT
                ENDIF
            ENDIF
            IF (PARAM.EQ.'BRAT') THEN
                RDEF = BRAT(LND)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    BRAT(LND) = dble(REA)
cc                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
cc                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'KFDP') THEN
                IDEF = KFDP(LND,LNE)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    KFDP(LND,LNE) = INT
c                    PRINT 91004, PARAM,LND,LNE,IDEF,INT
                ELSE
c                    PRINT 92004, PARAM,LND,LNE,IDEF,INT
                ENDIF
            ENDIF


30      CONTINUE
      ENDIF

91001 FORMAT(' pytcha: ',A4,'(',I4,')',
     &       ' changed from the default ',I8,' TO',I8)
92001 FORMAT(' pytcha (warning):',A4,'(',I4,')',
     &       ' not changed from the default ',I8,' TO',I8)
91004 FORMAT(' pytcha: ',A4,'(',I4,')',
     &       ' changed from the default ',I8,I8,' TO',I8,I8)
92004 FORMAT(' pytcha (warning):',A4,'(',I4,')',
     &       ' not changed from the default ',I8,I8,' TO',I8,I8)
91002 FORMAT(' pytcha: ',A4,'(',I4,')',
     &       ' changed from the default ',E10.3,' TO',E10.3)
92002 FORMAT(' pytcha (warning):',A4,'(',I4,')',
     &       ' not changed from the default ',E10.3,' TO',E10.3)

91003 FORMAT(' pytcha:',
     &'CHAF(',I4,') changed from the default:',A4,A4)
92003 FORMAT(' pytcha:',
     &'CHAF(',I4,') not changed from the default :',A4,A4)
*91004 FORMAT('                                     TO:',A4,A4)


      RETURN
      END
      SUBROUTINE BCASIN

      IMPLICIT NONE
      CHARACTER*23  FLAG,CHECK

C...Give sensible default values to switches and parameters.
      INTEGER I
      DOUBLE PRECISION EPARL
      REAL X,Y,Q2
      DOUBLE PRECISION PARL
      Integer LST
      COMMON /CASTOU/PARL(30),X,Y,Q2,LST(40)
      Integer IPY
      Double Precision PYVAR,PYPAR
      COMMON /MYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      INTEGER KLST,KIPY
      REAL EPYPAR,EPYVAR
      DIMENSION KLST(40),EPARL(30)
      DIMENSION KIPY(80),EPYPAR(80),EPYVAR(80)
C...CASTOU: Cuts, basic switches and parameters.
      DATA KLST/0,1,5,1,3,1,1,12,5,1,0,4,5,1,1,1,0,2,3,21*0/
**      DATA EPARL/1.,1.,0.44,0.75,0.226,0.,0.,0.015,2.,0.,0.01,4.,
**     +0.1,0.44,0.01,7.29735E-03,1.16637E-05,0.044,0.03,1.,10*0./
      DATA EPARL/1.,0.7,0.44,0.75,0.226,0.,0.,0.015,2.,0.,0.01,4.,
     +0.1,0.44,0.01,7.29735E-03,1.16637E-05,0.044,0.03,1.,10*0./
      DATA KIPY/
     + 0,     0,     2,     2,     6,     1,     1,     4,     3,     1,
     + 3,     1,     1,     2,     1,     1,     4,     1,     1,     1,
     + 0,     1,     1,     1,     1,     1,     1,     0,     0,     0,
     + 1,     2,     1,     1,    30,    33,     1,     1,     7,     0,
     + 0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     + 0,     0,     0,     1,   100,     0,     0,     0,     0,     0,
     + 0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     + 0,     0,     0,     0,     0,     0,     0,     0,     0,     0/
      DATA (EPYPAR(I),I=1,40)/
     +   7.299E-03,   2.290E-01,   2.000E-01,   2.500E-01,   4.000E+00,
     +   1.000E+00,   4.400E-01,   4.400E-01,   7.500E-02,   0.000E+00,
     +   2.000E+00,   2.000E+00,   1.000E+00,   0.000E+00,   3.000E+00,
     +   1.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   1.000E+00,
     +   2.500E-01,   1.000E+00,   2.000E+00,   1.000E-03,   1.000E+00,
     +   1.000E+00,   1.000E+00,  -2.000E-02,  -1.000E-02,   0.000E+00,
     +   0.000E+00,   1.600E+00,   0.500E+00,   0.200E+00,   3.894E-01,
     +   1.000E+00,   3.300E-01,   6.600E-01,   0.000E+00,   1.000E+00/
      DATA (EPYPAR(I),I=41,80)/
     +   2.260E+00,   1.000E+04,   1.000E-04,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00/
      DATA EPYVAR/80*0./
      DATA CHECK /'CASCADE 2.3 INITIALISED'/
      DATA FLAG  /'                      '/

      IF (FLAG .EQ. CHECK) THEN
c         WRITE(6,*) ' BLOCK DATA for CASCADE 2.0 REINITIALISED.'
      ENDIF

      CALL UCOPY(KLST,LST,40)
      CALL UCOPY(EPARL,PARL,30)
      CALL UCOPY(KIPY,IPY,80)
      CALL UCOPY(EPYPAR,PYPAR,80)
      CALL UCOPY(EPYVAR,PYVAR,80)
      FLAG = CHECK
c      WRITE(6,*) ' BLOCK DATA for CASCADE  initialised'

      RETURN
      END
      SUBROUTINE casini

      IMPLICIT NONE
#include "casparas.inc"
#include "caspara5.inc"
#include "casparam.inc"
#include "casluco.inc"
#include "casludat1.inc"
#include "casvalues.inc"
#include "casshwr.inc"
#include "casprkt.inc"
#include "casscal.inc"
#include "casglu.inc"
#include "casmult.inc"
#include "casherup.inc"
      CHARACTER CHIN*100
      Double Precision PYPAR,PYVAR
	Integer IPY
      COMMON /MYPARA/ IPY(80),PYPAR(80),PYVAR(80)
	Integer IDEBUG
      COMMON/ INTERN/IDEBUG
      REAL ACC1,ACC2
	Integer IINT,NCB
      COMMON /INTEGR/ ACC1,ACC2,IINT,NCB
      Integer IGLU
	Common/CAGLUON/Iglu
      DOUBLE PRECISION SIN2W,XMW2
      COMMON/ELWEAK/SIN2W,XMW2
	DOUBLE PRECISION XMW
	Double Precision Qscal 
	Integer Iqqbar
	Common/GLUDAT2/Qscal,Iqqbar 

      CALL BCASIN
c select unintegrated gluon : 	
c IGLU = 1: CCFM old set JS2001
c IGLU =  1001: CCFM J2003 set 1 
c IGLU =  1002: CCFM J2003 set 2 
c IGLU =  1003: CCFM J2003 set 3 
c IGLU =  1010: CCFM set A0
c IGLU =  1011: CCFM set A0+
c IGLU =  1012: CCFM set A0-
c IGLU =  1013: CCFM set A1
c IGLU =  1020: CCFM set B0
c IGLU =  1021: CCFM set B0+
c IGLU =  1022: CCFM set B0-
c IGLU =  1023: CCFM set B1
c
c IGLU=2 --> derivative GRV
c IGLU=3 --> Bluemlein
c IGLU=4 --> KMS
c IGLU=5 --> GBW
c IGLU=6 --> KMR
c IGLU=6 --> RS
      IGLU = 1010
      IDEBUG = 0
c collinear gluon(for tests...)
      ICOLGLU=1005006	
c select integration procedure
      NCB = 20000
      ACC1 = 1.
      ACC2 = 0.5
C FINAL STATE PARTON SHOWER
      IFPS = 1
      ILEPTO = 1
c     scale factor for final state PS
      SCALFAF = 1.      
C LEPIN KE = 11
      KE = 11
	Ires(1) = 0
      KGL = 21
      KPH=22
      KEB = KPH
      KPA = 1
      NFLAV = 5
      IHF = 0
	
C PPIN KP = 2212
      KP = 2212
	Ires(2) = 1
C select interaction type
C INTER = 0 --> photon exchange
      INTER = 0
C select BGF in semihard approach a la Catani et al. ISEMIH = 1
      ISEMIH = 1
C RUNNING ALPHA_S IRUNA = 1
      IRUNA=1
C RUNNING ALPHA_EM IRUNAEM = 1
      IRUNAEM=0
C SELECT SCALE FOR STRUCTURE FKT AND ALPHA_S
C IQ2 =1 SCALE Q2 = MASS**2
C IQ2 =2 SCALE Q2 = SHAT
C IQ2 =3 SCALE Q2 = MASS**2 + PT**2
C IQ2 =4 SCALE Q2 = Q2
C IQ2 =5 SCALE Q2 = Q2 + pt**2
      IQ2=3
C SELECT PROCESS
      IPRO= 10
	PT2CUT(10) = 0
	PT2CUT(11) = 0
C Minimum Q^2 of electron to be generated
      QMI = 5.d0
C Maximum Q^2 of electron to be generated
      QMA = 10D8
C Minimum y of electron to be generated
      YMI=0.0d0
C Minimum y of electron to be generated
      YMA=1.0d0
C Maximium theta angle of scattered electron
      THEMA = 180.0D0
C Minimum  theta angle of scattered electron
      THEMI =   0.0D0
C ELECTRON MOMENTUM
      PLEPIN =-30.
C PROTON MOMENTUM
      PPIN   = 820.
C PERFORM FRAGMENTATION NFRAG=0/1
      NFRAG = 1
c printout of intermediate/virtual partons during initial state cascade
      IPST = 0
c safety for fragmentation
      PYPAR(12) = 1.5
      MSTU(12) = 0
      CHIN='PMAS(4,1)=1.5'
c      CALL PYGIVE(CHIN)
c      CALL PYGIVE('PMAS(1,1)=0.450')
c      CALL PYGIVE('PMAS(2,1)=0.450')
      IHFLA = 3
c

	ICCFM = 1
	
	IORDER = 1
c primordial kt
      PRKT1 = 0d0
      PRKT2 = 0d0
c width of Gaussian distribution in transverse momentum
c when a non-trivial target remnant is split into two particles
c (PARL(14) in LEPTO )
      PRKTREM = 0.35d0
c energy sharing of remnant:
      IRspl = 4	
	SCALFA = 1d0	
c multiple scattering
      Mult = 0	
C SET ELECTROWEAK PARAMETERS
      SIN2W = 0.23D0
      XMW = 80.D0
      XMW2 = XMW**2
c parameters for scales
      Qscal = 1d0
	Iqqbar = 1

c use les houches accord for frag and final state ps 
      ILHA = 0	
c      WRITE(6,*) ' casini: initialisation of CASCADE '

      RETURN
      END
C...  END INITIALIZE



*CMZ :  5.10/01 30/06/98  14.37.25  by  Hannes Jung
*CMZ :  5.10/00 02/05/98  17
*-- Author :
************************************************************************
*    ====================================================              *
      SUBROUTINE BASES( FXN, S, SIGMA, CTIME, IT1, IT2 )
*    ====================================================              *
*      Subroutine BASES for the Numerical integration.                 *
*      In terms of this program Integration can be done, furthermore   *
*      a probability distribution can be made for the event generation.*
*      The event with weight one is generated by program SPRING.       *
* ((Input))                                                            *
*    from the arguement                                                *
*      FXN    : Name of function program                               *
*    from the labeled common /BASE1/                                   *
*      XL(50) : Lower limits of the integration variabels              *
*      XU(50) : upper limits of the integration variabels              *
*      NDIM   : Dimension of the integration                           *
*      NCALL  : Number of sampling points per iteration                *
*    from the lebeled common /BASE2/                                   *
*      ITMX*  : Number of iteration                                    *
*      ACC*   : Required accuracies                                    *
* ((Output))                                                           *
*      S      : Estimate of the integral                               *
*      SIGMA  : Standard deviation of the estimate                     *
*      CTIME  : Computing time required for integration                *
*      IT1    : Number of iterations for the grid defining step        *
*      IT2    : Number of iterations for the integration step          *
C*                                                                     *
C*       Coded by S.Kawabata         April '94                         *
C*                                                                     *
C***********************************************************************
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL FXN
      PARAMETER (MXDIM = 50)
*
*     JFLAG =  0 : First trial for defining grids.
*     JFLAG =  1 : First trial for event accumulation.
*     JFLAG =  2 : Second or more trial for defining grids.
*     JFLAG =  3 : Second or more trial for accumulation.
*                                                                      *
      COMMON /BASE0/ JFLAG,IBASES
      COMMON /BASE1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,
     .               IG(MXDIM),NCALL
      COMMON /BASE2/ ACC1,ACC2,ITMX1,ITMX2
      REAL*4 STIME
      COMMON /BSRSLT/AVGI,SD,CHI2A,STIME,ITG,ITF
      CHARACTER*80 ERROR
      COMMON /BWARN1/ NERROR
      COMMON /BWARN2/ ERROR(3,3)
*        INTV = ( 0 / 1 / any ) = ( Batch / Batch(Unix) / Interactive )
*        IPNT = ( 0 / any ) = ( IBM Type / Ascii printer )
      COMMON /BSCNTL/ INTV, IPNT, NLOOP, MLOOP

       COMMON/NINFO/ NODEID, NUMNOD
       COMMON /BDATE/ IDATE(3),ITIME(2)
*            IDATE(1) : year        ITIME(1) : hour
*            IDATE(2) : month       ITIME(2) : minute
*            IDATE(3) : day
      REAL*4 TIMEBS,TIMINT,TIMESP,TIME0,RTIME,TIMEB1,TIMEB2,TIMES1
      COMMON /BTIME1/ TIME0,RTIME,TIMEB1,TIMEB2,TIMES1
      COMMON /BTIME2/ TIMEBS(0:2),TIMINT,TIMESP(0:2)

*-------------------------------------------------
*     Check the parameters defined by user
*------------------------------------------------------

      CALL BSCHCK

* ---------------------------------------------------------------
*          Initialize timer
* ---------------------------------------------------------------

       CALL BSDATE

       JFLAG  = 0
       LU     = 16
c       LU     = 6

       IF( INTV .GT. 1 ) THEN
           CALL BSPRNT( LU, 1, IDUM1, IDUM2 )
       ENDIF

C  -----------------------------------------------------
C     Defining grids
C  -----------------------------------------------------
*
       DO 100 I = 1, NWILD
          IG(I) = 1
  100  CONTINUE

       CALL BSETGU

       IF( INTV .GT. 1 ) THEN
           CALL BSPRNT( LU, 4, IDUM1, IDUM2 )
       ENDIF

       CALL BSUTIM( 0, 2 )

*     ===================
       CALL BSINTG( FXN )
*     ===================        For a parallel computer
C                                      CALL BSCAST( JFLAG, 1 )

*  ----------------------------------------------------
*     Accumulation to make probability distribution
*  ----------------------------------------------------
*     ===================
       CALL BSINTG( FXN )
*     ===================        For a parallel computer
C                                      CALL BSCAST( JFLAG, 1 )
       S     = AVGI
       SIGMA = SD
       CTIME = STIME
       IT1   = ITG
       IT2   = ITF

       CALL BSUTIM( 0, 2 )
       TIMEB2 = RTIME

       IF( NERROR .GT. 0 ) THEN
           WRITE(6,9900)
 9900      FORMAT(1X,'****************************************',
     .               '***************************************',
     .           /1X,'* (((( Warning in the integration step ',
     .               '))))                                   *',
     .           /1X,'*                                      ',
     .               '                                       *')
           DO 990 J = 1,NERROR
           DO 990 I = 1,3
              WRITE(6,9901) ERROR(I,J)
 9901         FORMAT(1X,A79)
  990      CONTINUE
           WRITE(6,9902)
 9902      FORMAT(1X,'*                                      ',
     .               '                                       *',
     .           /1X,'*(( Suggestion ))                      ',
     .               '                                       *',
     .           /1X,'* (1) Try integration again with larger ',
     .               'number of sample points than this job.*',
     .           /1X,'* or                                   ',
     .               '                                       *',
     .           /1X,'* (2) The integral variables are not sui',
     .               'ted for the function.                 *',
     .           /1X,'*     Take another integral variables !!',
     .               '                                      *',
     .           /1X,'*                                       ',
     .               '                                      *',
     .           /1X,'****************************************',
     .               '***************************************')
       ENDIF

       IF( INTV .GT. 1 ) THEN
           CALL BSPRNT( LU, 2, IDUM1, IDUM2 )
       ENDIF

       RETURN
       END
*CMZ :  5.10/00 02/05/98  17.48.43  by  Unknown
*-- Author :
***********************************************************************
*                                                                     *
*    ==========================                                       *
      SUBROUTINE BSINTG( FXN )
*    ==========================                                       *
*((Function))                                                         *
*    Subroutine performs N-dimensional Monte Carlo integration        *
*    for four vector generation of simulated events                   *
*                                                                     *
*       JFLAG = 0 ; First Trial of Defining Grid                      *
*       JFLAG = 1 ; First Trial of Data Accumulation                  *
*       JFLAG = 2 ; Second Trial of Defining Grid                     *
*       JFLAG = 3 ; Second Trial of Data Accumulation                 *
*                                                                     *
*    Coded   by S.Kawabata    July 1980 at DESY, Hamburg              *
*    Last update              March 1994                              *
*                                                                     *
***********************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)

      EXTERNAL FXN
      PARAMETER (MXDIM = 50, NDMX = 50, LENG = 32768)
      COMMON /BASE0/ JFLAG,IBASES
      COMMON /BASE1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,
     .               IG(MXDIM),NCALL
      COMMON /BASE2/ ACC1,ACC2,ITMX1,ITMX2
      COMMON /BASE3/ SCALLS,WGT,TI,TSI,TACC,IT
      COMMON /BASE4/ XI(NDMX,MXDIM),DX(MXDIM),DXD(LENG),DXP(LENG),
     .               ND,NG,NPG,MA(MXDIM)
      PARAMETER (ITM = 50)
      REAL*4 TIME, EFF, WRONG, TRSLT, TSTD, PCNT
      COMMON /BASE5/ ITRAT(ITM,0:1),TIME(ITM,0:2),EFF(ITM,0:1),
     .               WRONG(ITM,0:1),RESLT(ITM,0:1),ACSTD(ITM,0:1),
     .               TRSLT(ITM,0:1),TSTD(ITM,0:1),PCNT(ITM,0:1)
      COMMON /BASE6/ D(NDMX,MXDIM),
     .               ALPH,XSAVE(NDMX,MXDIM),XTI,XTSI,XACC,ITSX
      REAL*4 STIME
      COMMON /BSRSLT/AVGI,SD,CHI2A,STIME,ITG,ITF
      CHARACTER*80 ERROR
      COMMON /BWARN1/ NERROR
      COMMON /BWARN2/ ERROR(3,3)
*
*        INTV = ( 0 / 1 / any ) = ( Batch / Batch(Unix) / Interactive )
*        IPNT = ( 0 / any ) = ( IBM Type / Ascii printer )
      COMMON /BSCNTL/ INTV, IPNT, NLOOP, MLOOP

      REAL*8  X(MXDIM)
      INTEGER KG(MXDIM),IA(MXDIM)

      COMMON/NINFO/ NODEID, NUMNOD
      REAL*4 TIMEBS,TIMINT,TIMESP,TIME0,RTIME,TIMEB1,TIMEB2,TIMES1
      COMMON /BTIME1/ TIME0,RTIME,TIMEB1,TIMEB2,TIMES1
      COMMON /BTIME2/ TIMEBS(0:2),TIMINT,TIMESP(0:2)
C     REAL*8  TX(2)
      INTEGER NCNODE(2,512),NPNODE(2,512)
C     INTEGER NEFF(2)
*
*     Parameters for checking convergency
*
      DATA ACLMT,FC / 25.0D0, 5.0D0 /


      DATA  ONE/ 1.0D0/, ZERO/0.0D0/, LU / 16/
      DATA  N0/0/, N1/1/, HUNDRT/100.0D0/

************************************************************************
*                       Initialization Part
************************************************************************
*=======================================================================
*          Determine the number of hypercubes NSP
*=======================================================================

      XND     = ND
      NSP     = NG**NWILD
      XJAC    = 1.0D0
      DO  5 I = 1, NDIM
         XJAC = XJAC*DX(I)
    5 CONTINUE
      CALLS   = NSP*NPG
      DXG     = 1.0D0/NG
      DV2G    = DXG**(2*NWILD)/NPG/NPG/(NPG-1)
      DXG     = DXG*XND

      IF( NSP .EQ. 1 ) THEN
*=======================================================================
*           Determination of the number of sampling points
*               per node in the single hypercube case
*=======================================================================
          MEX     = MOD(NPG,NUMNOD)
          NPERCP  = NPG/NUMNOD
          NPGT    = 0
          DO  12 NODEX = 1,NUMNOD
             NPGS  = NPGT + 1
             NPGT  = NPGT + NPERCP
             IF( NODEX .LE. MEX ) NPGT = NPGT + 1
             NCNODE(1,NODEX) = 1
             NCNODE(2,NODEX) = 1
             NPNODE(1,NODEX) = NPGS
             NPNODE(2,NODEX) = NPGT
   12     CONTINUE
      ELSE
*=======================================================================
*          Determination of the number of hypercubes
*              per node in many hypercubes case
*=======================================================================
          MEX     = MOD(NSP,NUMNOD)
          NPERCP  = NSP/NUMNOD
          NSPT    = 0
          DO  15 NODEX = 1,NUMNOD
             NSPS  = NSPT + 1
             NSPT  = NSPT + NPERCP
             IF( NODEX .LE. MEX ) NSPT = NSPT + 1
             NCNODE(1,NODEX) = NSPS
             NCNODE(2,NODEX) = NSPT
             NPNODE(1,NODEX) = 1
             NPNODE(2,NODEX) = NPG
   15     CONTINUE
      ENDIF
*=======================================================================
      NEND    = N0
      ATACC   = ZERO
      NERROR  = N0
      NER1    = N0
      NER2    = N0
      NER3    = N0
      SUMTI   = ZERO
      SUMTSI  = ZERO

      IF(JFLAG .EQ. N0 .OR. JFLAG .EQ. N1 ) THEN
*-----------------------------------------------------------------------
*        JFLAG = 0  : The first trial of the grid optim. step
*        JFLAG = 1  : The first trial of the integration step
*-----------------------------------------------------------------------
         DO 10 J  = N1,NSP
           DXD(J) = ZERO
           DXP(J) = ZERO
   10    CONTINUE
*       -----------------
         ISTEP   = JFLAG
*       -----------------
         IT1   = N1
         SI    = ZERO
         SI2   = ZERO
         SWGT  = ZERO
         SCHI  = ZERO
*       =============
         CALL BHRSET
*       =============
         NSU     = N0
         SCALLS= ZERO
      ELSE
*-----------------------------------------------------------------------
*        JFLAG = 2  : The continuation of the grid optim. step
*        JFLAG = 3  : The continuation of the integration step
*-----------------------------------------------------------------------
C        IF( JFLAG .EQ. 2 ) THEN
*           -------------
C            ISTEP  = N0
*           -------------
C        ELSE
C    .   IF( JFLAG .EQ. 3 ) THEN
*           -------------
C            ISTEP  = N1
*           -------------
C        ELSE
C                *****************
C                      STOP
C                *****************
C         ENDIF
C
C         IT1   = IT + 1
      ENDIF

*------- Set the expected accuracy and the max. iteration number -------

      ITMX   = ITMX1
      ACC    = ACC1*0.01D0
      IF( ISTEP .EQ. N1 ) THEN
         ITMX = ITMX2
         ACC  = ACC2*0.01D0
      ENDIF

*-------- Print the title of the convergency behavior table -----------
*                  in the interactive mode
      IF( INTV .GT. 1 ) THEN
*         -----------------------------------
           CALL BSPRNT( LU, 5, ISTEP, IDUM2 )
*         -----------------------------------
      ENDIF
      NEGFLG     = 0

*    =====================
      CALL BSUTIM( 0, 2 )
*    =====================

*********************************************************************
*               Main Integration Loop
*********************************************************************
*    ========
      DO 500  IT = IT1,ITMX
*    ========
*=======================================================================
*                 Initialization for the iteration
*=======================================================================

         SCALLS  = SCALLS + CALLS
         NGOOD   = N0
         NEGTIV  = N0
         TI      = ZERO
         TSI     = TI

         IF( ISTEP .EQ. N0 ) THEN
             DO 200 J= N1,NDIM
             DO 200 I=1,ND
                D(I,J)= TI
  200        CONTINUE
         ENDIF

         NODEX  = NODEID
         IF( NODEID .EQ. 0 )  NODEX = NUMNOD

*---------------------------------------------------------------------
*        Distributing hyper cubes to NumNode nodes
*           NCNODE(1,NODEX)   : 1st cube number for the node NODEX
*           NCNODE(2,NODEX)   : Last cube number for the node NODEX
*                    NODEX    : node number 1 => NumNode(=0)
*                    NODEX    : node number 1 => NumNode(=0)
*---------------------------------------------------------------------

         NSP1  = NCNODE(1,NODEX)
         NSP2  = NCNODE(2,NODEX)
*                                 Dummy loopfor a parallel processor
C                                 IF( NSP1 .GT. 1 ) THEN
C                                     CALL DRLOOP( NDIM*NPG*(NSP1-1) )
C                                 ENDIF

*=====================================================================
*      Loop for hypercube from NSP1 to NSP2 in the NodeX-th node
*=====================================================================
*       ========
         DO 400 NCB = NSP1, NSP2
*       ========
            FB      = 0.0
            F2B     = 0.0
            NP      = NCB - 1
            IF( NWILD .GT. 1 ) THEN
                DO 210 J = 1,NWILD-1
                   NUM   = MOD(NP,MA(J+1))
                   KG(J) = NUM/MA(J) + 1
  210           CONTINUE
            ENDIF
            KG(NWILD)     = NP/MA(NWILD) + 1

*---------------------------------------------------------------------
*       If number of hypercubes is only one,
*        Distributing sampling points to NumNode nodes
*           NPNODE(1,NODEX)   : 1st sample point for the node NODEX
*           NPNODE(2,NODEX)   : Last sample point for the node NODEX
*                    NODEX    : node number 1 => NumNode(=0)
*---------------------------------------------------------------------

            NPG1  = NPNODE(1,NODEX)
            NPG2  = NPNODE(2,NODEX)
*                                 Dummy loop for a parallel processor
C                                 IF( NPG1 .GT. 1 ) THEN
C                                     CALL DRLOOP( NDIM*(NPG1-1) )
C                                 ENDIF

*=====================================================================
*          Loop for sampling points from NPG1 to NPG2
*                in the single hypercube case
*=====================================================================
*          ========
            DO 300 NTY = NPG1,NPG2
*          ========
*---------------------------------------------------------------------
*        Determine the integration variables by random numbers
*---------------------------------------------------------------------

               WGT   = XJAC
               DO 250 J= 1,NDIM
                  IF( J .LE. NWILD ) THEN
                      XN  = (KG(J)-DRN(IDUMY))*DXG+1.D0
                  ELSE
                      XN  = ND*DRN(IDUMY)+1.D0
                  ENDIF
                  IA(J)   = XN
                  IAJ     = IA(J)
                  IF( IAJ .EQ. 1) THEN
                      XO  = XI(IAJ,J)
                      RC  = (XN-IA(J))*XO
                  ELSE
                      XO  = XI(IAJ,J)-XI(IAJ-1,J)
                      RC  = XI(IAJ-1,J)+(XN-IAJ)*XO
                  ENDIF
                  X(J)    = XL(J)+RC*DX(J)
                  WGT     = WGT*XO*XND
  250          CONTINUE
*-----------------------------------------------------------------------
*                     =======
               FXG  =  FXN(X)*WGT
*                     =======
*-----------------------------------------------------------------------
*             Check the value of the integrand
*-----------------------------------------------------------------------

               IF( FXG .NE. 0.0 ) THEN
                   NGOOD = NGOOD + 1
                   IF( ISTEP .EQ. 1 ) THEN
                       DXD(NCB) = DXD(NCB) + FXG
                       IF( FXG .GT. DXP(NCB) ) DXP(NCB) = FXG
                   ENDIF
                   IF( FXG .LT. 0.0 ) THEN
                       NEGTIV= NEGTIV+ 1
                       IF( NEGFLG .EQ. 0 ) THEN
                          WRITE(6,9200) IT,NODEID
 9200                     FORMAT(1X,
     .                       '******* WARNING FROM BASES ********',
     .                       '***********',
     .                       /1X,'*  Negative FUNCTION at IT =',I3,1X,
     .                       ', node = ',I3,1X,'*',
     .                       /1X,'***********************************',
     .                       '***********')
                          NEGFLG  = 1
                       ENDIF
                   ENDIF
               ENDIF

*-----------------------------------------------------------------------
*              Accumulation of FXG and FXG*FXG
*-----------------------------------------------------------------------

               F2    = FXG*FXG
               FB    = FB + FXG
               F2B   = F2B + F2

               IF( ISTEP .EQ. 0 ) THEN
                   DO 260  J = 1,NDIM
                      D(IA(J),J)= D(IA(J),J)+F2
  260              CONTINUE
               ENDIF
*======
  300       CONTINUE
*======
*------------------------------------------- for a parallel processor
*                                 Dummy loop for a parallel processor
C                                 IF( NPG2 .LT. NPG ) THEN
C                                     CALL DRLOOP(NDIM*(NPG-NPG1))
C                                 ENDIF
*                                 Global sum of FB and F2B
C                                 IF( NSP .EQ. 1 ) THEN
C                                     CALL BSDSUM(  FB, 1 )
C                                     CALL BSDSUM( F2B, 1 )
C                                 ENDIF
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*         Calculate the estimate and variance in the hypercube
*-----------------------------------------------------------------------

            F2B   = DSQRT(F2B*NPG)
            F2S   = (F2B-FB)*(F2B+FB)
            TI    = TI+FB
            TSI   = TSI + F2S

*======
  400    CONTINUE
*======
*------------------------------------------- for a parallel processor
*                                 Dummy loop
C                                 IF( NSP2 .LT. NSP ) THEN
C                                     CALL DRLOOP(NDIM*NPG*(NSP-NSP2))
C                                 ENDIF

*                                 Global sum of efficiency and frequency
*                                     of negative valued function
C                                 NEFF(1) = NGOOD
C                                 NEFF(2) = NEGTIV
C                                 CALL BSISUM( NEFF, 2 )

C                                 TX(1) = TI
C                                 TX(2) = TSI
C                                 IF( NSP .EQ. 1 ) THEN
C                                     CALL BSDSUM(   TX, 2 )
C                                 ENDIF

*                                 Global sum of grid information
C                                 IF( ISTEP .EQ. 0 ) THEN
C                                     NOWORK = NDMX*NDIM
C                                     CALL BSDSUM(    D, NOWORK )
C                                 ENDIF

*=====================================================================
*           Compute Result of this Iteration
*=====================================================================
*--------------------------------------------------------------------
*           Accumulate the histogram entries
*--------------------------------------------------------------------
*       -------------
         CALL BHSAVE
*       -------------
*--------------------------------------------------------------------

C        TI     = TX(1)
C        TSI    = TX(2)
C        NGOOD  = NEFF(1)
C        NEGTIV = NEFF(2)

         TI    = TI/CALLS
         TSI   = TSI*DV2G
**
         IF( TSI .LE. 1.0D-37 ) TSI = 1.0D-37
**
         TI2   = TI*TI

         IF( NGOOD .LE. 10 ) THEN
*           --------------------------------
             CALL BSPRNT( LU, 9, IDUM1, IDUM2 )
*           --------------------------------
*            *****************
c hannes hack
ccc                    write(6,*) ' integration has not enough good points '
                    Avgi = 0.
c                   STOP
                   return
*            *****************

         ENDIF

*--------------------------------------------------------------------
*               Calculate the cumulative result
*--------------------------------------------------------------------

         WGT   = ONE/TSI
         SI    = SI+TI*WGT
         SWGT  = SWGT+WGT
         SCHI  = SCHI+TI2*WGT
         AVGI  = SI/SWGT
         CHI2A = ZERO
         IF(IT .GT. N1 ) CHI2A = (SCHI - SI*AVGI)/(IT-.999D0)
         SD    = DSQRT(ONE/SWGT)

*---------------------------------------------------------------------
*             Save the results in the buffer
*---------------------------------------------------------------------

         TSI   = DSQRT(TSI)
         ITX         = MOD( IT, ITM)
         IF( ITX .EQ. 0 ) ITX = ITM
         ITRAT(ITX,ISTEP)  = IT
         EFF  (ITX,ISTEP)  = NGOOD/CALLS*HUNDRT
         WRONG(ITX,ISTEP)  = NEGTIV/CALLS*HUNDRT
         RESLT(ITX,ISTEP)  = AVGI
         ACSTD(ITX,ISTEP)  = SD
         TRSLT(ITX,ISTEP)  = TI
         TACC              = ABS(TSI/TI*HUNDRT)
         TSTD (ITX,ISTEP)  = TACC
         PCNT (ITX,ISTEP)  = ABS(SD/AVGI*HUNDRT)

*----------------------------------------------------------------------
*                  Check cumulative accuracy
*----------------------------------------------------------------------

         IF( NODEID .EQ. 0 ) THEN

*-------------------  Check cumulative accuracy -----------------------

             SDAV  = SD/AVGI
             IF((ABS(SDAV) .LE. ACC)) NEND = N1

             IF( ISTEP .EQ. N1 ) THEN
                 IF( TACC .GT. ACLMT ) THEN
                     IF( NER1 .EQ. 0 ) THEN
                         NERROR = NERROR + 1
                         WRITE(ERROR(1,NERROR),9900) NERROR,IT,ACLMT
 9900                    FORMAT('* (',I1,') Temp. accuracy of it-#',
     .                         I3,' is too large comparing to',
     .                         F6.2,' percent.',6X,'*')
                         WRITE(ERROR(2,NERROR),9901) TACC,ACLMT
 9901                    FORMAT('*',8X,'Temp. accuracy (',
     .                         F7.4,' % )  >>   (',
     .                         F7.4,' % )',23X,'*')
                         WRITE(ERROR(3,NERROR),9902)
 9902                    FORMAT('*',77X,'*')
                         NER1  = 1
                     ENDIF
                 ENDIF
                 IF( IT .GT. 1 ) THEN
                     IF(( TI .GT. AVTI+FDEVI ) .OR.
     .                  ( TI .LT. AVTI-FDEVI )      ) THEN
                          IF( NER2 .EQ. 0 ) THEN
                              NERROR = NERROR + 1
                              WRITE(ERROR(1,NERROR),9910) NERROR,IT,FC
 9910                         FORMAT('* (',I1,') Temp. estimate of ',
     .                        'it-#',I3,' fluctuates more than ',
     .                               F4.1,'*average-sigma.',6X,'*')
                              RE = TI
*patch TI:1995/08/25
                              ARE = ABS(RE)
*old                          CALL BSORDR( RE, FX2, ORDER, IORDR )
                              CALL BSORDR( ARE, FX2, ORDER, IORDR )
*patch end
                              RE = TI/ORDER
                              RE1 = AVTI
                              AC  = FDEVI
*patch TI:1995/08/25
                              ARE1 = ABS(AVTI)
                              AAC  = ABS(FDEVI)
                              IF( ARE1 .GE. AAC ) THEN
                                  CALL BSORDR( ARE1, FX2, ORDR1, IORDR1)
                              ELSE
                                  CALL BSORDR( AAC, FX2, ORDR1, IORDR1)
                              ENDIF
*                             IF( RE1 .GE. AC ) THEN
*                                 CALL BSORDR( RE1, FX2, ORDR1, IORDR1)
*                             ELSE
*                                 CALL BSORDR( AC, FX2, ORDR1, IORDR1)
*                             ENDIF
*patch end
                              RE1 = AVTI/ORDR1
                              AC  = AC/ORDR1
                              WRITE(ERROR(2,NERROR),9911) RE,IORDR,
     .                                          RE1,AC,IORDR1
 9911                         FORMAT('*        Temp. Estimate (',
     .                         F10.6,' E',I3,')  >  (',F10.6,'+',F8.6,
     .                         ' ) E',I3,', or',1X,'*')
                              WRITE(ERROR(3,NERROR),9912) RE,IORDR,
     .                                          RE1,AC,IORDR1
 9912                         FORMAT('*        Temp. Estimate (',
     .                         F10.6,' E',I3,')  <  (',F10.6,'-',F8.6,
     .                         ' ) E',I3,5X,'*')
                              NER2 = 1
                          ENDIF
                     ENDIF
                     IF( TSI .GT. FDEVI ) THEN
                         IF( NER3 .EQ. 0 ) THEN
                             NERROR = NERROR + 1
                             WRITE(ERROR(1,NERROR),9920) NERROR,IT,FC
 9920                        FORMAT('* (',I1,') Error of it-#',
     .                              I3,' fluctuates more than',F4.1,
     .                              '*average-sigma.',16X,'*')
                             RE1 = TSI
*patch TI:1995/08/25
                             ARE1 = ABS(TSI)
*                            CALL BSORDR( RE1, FX2, ORDER, IORDR)
                             CALL BSORDR( ARE1, FX2, ORDER, IORDR)
*patch end;
                             RE1 = TSI/ORDER
                             AC  = FDEVI
*patch TI:1995/08/25
                             AAC  = ABS(FDEVI)
*                            CALL BSORDR( AC, FX2, ORDR1, IORDR1)
                             CALL BSORDR( AAC, FX2, ORDR1, IORDR1)
*patch end;
                             AC  = AC/ORDR1
                             WRITE(ERROR(2,NERROR),9921) RE1,IORDR,
     .                                         AC,IORDR1
 9921                        FORMAT('*        Temp. Error (',
     .                         F10.6,' E',I3,')  >  (',F10.6,
     .                         ' E',I3,')',18X,'*')
                             WRITE(ERROR(3,NERROR),9902)
                             NER3  = 1
                         ENDIF
                     ENDIF
                 ENDIF
                 SUMTSI = SUMTSI + TSI
                 SUMTI  = SUMTI  + TI
                 AVTSI  = SUMTSI/FLOAT(IT)
                 AVTI   = SUMTI/FLOAT(IT)
                 FDEVI  = FC*AVTSI
             ENDIF
         ENDIF

*------------------------------------------- for a parallel processor

*                                  Broadcast
C                                  CALL BSCAST( NEND, 1 )

*----------------------------------------------------------------------
*        Smoothing the Distribution D(I,J) and refine the grids
*----------------------------------------------------------------------

         IF( ISTEP .LE. N0 ) THEN
             IF( IT .EQ. ITMX ) NEND = N1
*           ---------------------
             CALL BSETGV( NEND )
*           ---------------------
         ENDIF
*       ==========================
         CALL BSUTIM( 0, ISTEP )
*       ==========================

         TIME (ITX,ISTEP)  = TIMINT
         STIME             = TIMINT

*---- Print the convergency behavior table in the interactive mode ----
         IF( INTV .GT. 1 ) THEN
*            ---------------------------------
              CALL BSPRNT ( LU, 6, ISTEP, IDUM2 )
*            ---------------------------------
         ENDIF

         IF( NEND .EQ. N1 ) GO TO 600

*       ======================
         CALL BSUTIM( 0, 2 )
*       ======================
*======
  500 CONTINUE
*======
      IT    = IT - N1
      NEND  = N1

***********************************************************************
*                   Termination of BASES
***********************************************************************
*======
  600 CONTINUE
*======
*---------------------------------------------- For a parallel computer

*                                 Global sum of histograms
C                                 CALL BHSUM
*                                 Global sum of probabilities
C                                 CALL BSDSUM(  DXD, NSP )
*                                 Global sum of the max.value in each HC
C                                 CALL BSDSUM(  DXP, NSP )


*======================= End of the step ? ============================

      IF( NEND .EQ. N1 ) THEN
          IF( INTV .GT. 1 ) THEN
*            ---------------------------------
              CALL BSPRNT ( LU, 7, IDUM1, IDUM2 )
*            ---------------------------------
          ENDIF
          IF( ISTEP .EQ. N0) THEN
              JFLAG   = N1
              ITG     = IT
          ELSE
              JFLAG   = N0
              ITF     = IT
          ENDIF
      ENDIF
*    ======================
       CALL BSUTIM( 0, 2 )
*    ======================

      RETURN
      END

*
* $Id: rm48.F,v 1.2 1996/12/12 16:32:06 cernlib Exp $
*
* $Log: rm48.F,v $
* Revision 1.2  1996/12/12 16:32:06  cernlib
* Variables ONE and ZERO added to SAVE statement, courtesy R.Veenhof
*
* Revision 1.1.1.1  1996/04/01 15:02:55  mclareni
* Mathlib gen
*
*
      SUBROUTINE RM48(RVEC,LENV)
C     Double-precision version of
C Universal random number generator proposed by Marsaglia and Zaman
C in report FSU-SCRI-87-50
C        based on RANMAR, modified by F. James, to generate vectors
C        of pseudorandom numbers RVEC of length LENV, where the numbers
C        in RVEC are numbers with at least 48-bit mantissas.
C   Input and output entry points: RM48IN, RM48UT.
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RM48:                                    ++
C!!!      CALL RM48 (RVEC, LEN)     returns a vector RVEC of LEN     ++
C!!!                   64-bit random floating point numbers between  ++
C!!!                   zero and one.                                 ++
C!!!      CALL RM48IN(I1,N1,N2)   initializes the generator from one ++
C!!!                   64-bit integer I1, and number counts N1,N2    ++
C!!!                  (for initializing, set N1=N2=0, but to restart ++
C!!!                    a previously generated sequence, use values  ++ 
C!!!                    output by RM48UT)                            ++ 
C!!!      CALL RM48UT(I1,N1,N2)   outputs the value of the original  ++
C!!!                  seed and the two number counts, to be used     ++
C!!!                  for restarting by initializing to I1 and       ++  
C!!!                  skipping N2*100000000+N1 numbers.              ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C for 32-bit machines, use IMPLICIT DOUBLE PRECISION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RVEC(*)
      COMMON/R48ST1/U(97),C,I97,J97
      PARAMETER (MODCNS=1000000000)
      SAVE CD, CM, TWOM24, NTOT, NTOT2, IJKL,TWOM49, ONE, ZERO
      DATA NTOT,NTOT2,IJKL/-1,0,0/
      Logical FIRST
      Data FIRST/.TRUE./
C
      IF (NTOT .GE. 0)  GO TO 50
C
C        Default initialization. User has called RM48 without RM48IN.
      IJKL = 54217137
      NTOT = 0
      NTOT2 = 0
      KALLED = 0
      GO TO 1
C
      ENTRY      RM48IN(IJKLIN, NTOTIN,NTOT2N)
C         Initializing routine for RM48, may be called before
C         generating pseudorandom numbers with RM48.   The input
C         values should be in the ranges:  0<=IJKLIN<=900 OOO OOO
C                                          0<=NTOTIN<=999 999 999
C                                          0<=NTOT2N<<999 999 999!
C To get the standard values in Marsaglia's paper, IJKLIN=54217137
C                                            NTOTIN,NTOT2N=0
      IJKL = IJKLIN
      NTOT = MAX(NTOTIN,0)
      NTOT2= MAX(NTOT2N,0)
      KALLED = 1
C          always come here to initialize
    1 CONTINUE
      IJ = IJKL/30082
      KL = IJKL - 30082*IJ
      I = MOD(IJ/177, 177) + 2
      J = MOD(IJ, 177)     + 2
      K = MOD(KL/169, 178) + 1
      L = MOD(KL, 169)
      If(FIRST) then
         WRITE(6,'(A,I10,2X,2I10)') 
     &	 ' My-RM48 INITIALIZED:',IJKL,NTOT,NTOT2
         FIRST=.FALSE.
      ENDIF
CCC      PRINT '(A,4I10)', '   I,J,K,L= ',I,J,K,L
      ONE = 1.
      HALF = 0.5
      ZERO = 0.
      DO 2 II= 1, 97
      S = 0.
      T = HALF
      DO 3 JJ= 1, 48
         M = MOD(MOD(I*J,179)*K, 179)
         I = J
         J = K
         K = M
         L = MOD(53*L+1, 169)
         IF (MOD(L*M,64) .GE. 32)  S = S+T
    3    T = HALF*T
    2 U(II) = S
      TWOM49 = T
      TWOM24 = ONE
      DO 4 I24= 1, 24
    4 TWOM24 = HALF*TWOM24
      C  =   362436.*TWOM24
      CD =  7654321.*TWOM24
      CM = 16777213.*TWOM24
      I97 = 97
      J97 = 33
C       Complete initialization by skipping
C            (NTOT2*MODCNS + NTOT) random numbers
      DO 45 LOOP2= 1, NTOT2+1
      NOW = MODCNS
      IF (LOOP2 .EQ. NTOT2+1)  NOW=NTOT
      IF (NOW .GT. 0)  THEN
      WRITE(6,'(A,I15)') ' RM48IN SKIPPING OVER ',NOW
          DO 40 IDUM = 1, NTOT
          UNI = U(I97)-U(J97)
          IF (UNI .LT. ZERO)  UNI=UNI+ONE
          U(I97) = UNI
          I97 = I97-1
          IF (I97 .EQ. 0)  I97=97
          J97 = J97-1
          IF (J97 .EQ. 0)  J97=97
          C = C - CD
          IF (C .LT. ZERO)  C=C+CM
   40     CONTINUE
      ENDIF
   45 CONTINUE
      IF (KALLED .EQ. 1)  RETURN
C
C          Normal entry to generate LENV random numbers
   50 CONTINUE
      DO 100 IVEC= 1, LENV
      UNI = U(I97)-U(J97)
      IF (UNI .LT. ZERO)  UNI=UNI+ONE
      U(I97) = UNI
      I97 = I97-1
      IF (I97 .EQ. 0)  I97=97
      J97 = J97-1
      IF (J97 .EQ. 0)  J97=97
      C = C - CD
      IF (C .LT. ZERO)  C=C+CM
      UNI = UNI-C
      IF (UNI .LT. ZERO) UNI=UNI+ONE
      RVEC(IVEC) = UNI
C             Replace exact zeros by 2**-49
         IF (UNI .EQ. ZERO)  THEN
            RVEC(IVEC) = TWOM49
         ENDIF
  100 CONTINUE
      NTOT = NTOT + LENV
         IF (NTOT .GE. MODCNS)  THEN
         NTOT2 = NTOT2 + 1
         NTOT = NTOT - MODCNS
         ENDIF
      RETURN
C           Entry to output current status
      ENTRY RM48UT(IJKLUT,NTOTUT,NTOT2T)
      IJKLUT = IJKL
      NTOTUT = NTOT
      NTOT2T = NTOT2
      RETURN
      END
      
*
* importing sources
*
* Revision 1.1.1.1  2008/10/12 13:56:39  jung
* importing sources
*
* Revision 1.1.1.1  1996/02/15 17:48:36  mclareni
* Kernlib
*
*
          FUNCTION DHFINT(NARG,ARG,NENT,ENT,TABLE)
          Implicit None
C
C   INTERPOLATION ROUTINE. AUTHOR C. LETERTRE.
C   MODIFIED BY B. SCHORR, 1.07.1982.
C

          INTEGER   NENT(500)
          Double Precision     ARG(500),   ENT(500),TABLE(1000000)
          INTEGER   INDEX(32)
          Double Precision       WEIGHT(32)
          Double Precision DHFINT
          Integer NARG,LMAX,ISTEP,KNOTS,N,NDIM,ISHIFT,I,K
          Integer LMIN,LOCA,LOCB,LOCC,LGFILE
          Double Precision X,H,ETA
          LOGICAL   MFLAG,    RFLAG
          DHFINT  =  0.
          IF(NARG .LT. 1  .OR.  NARG .GT. 5)  GOTO 300
          LMAX      =  0
          ISTEP     =  1
          KNOTS     =  1
          INDEX(1)  =  1
          WEIGHT(1) =  1.
          DO 100    N  =  1, NARG
             X     =  ARG(N)
             NDIM  =  NENT(N)
             LOCA  =  LMAX
             LMIN  =  LMAX + 1
             LMAX  =  LMAX + NDIM
             IF(NDIM .GT. 2)  GOTO 10
             IF(NDIM .EQ. 1)  GOTO 100
             H  =  X - ENT(LMIN)
             IF(H .EQ. 0.)  GOTO 90
             ISHIFT  =  ISTEP
             IF(X-ENT(LMIN+1) .EQ. 0.)  GOTO 21
             ISHIFT  =  0
             ETA     =  H / (ENT(LMIN+1) - ENT(LMIN))
             GOTO 30
  10         LOCB  =  LMAX + 1
  11         LOCC  =  (LOCA+LOCB) / 2
             IF(X-ENT(LOCC))  12, 20, 13
  12         LOCB  =  LOCC
             GOTO 14
  13         LOCA  =  LOCC
  14         IF(LOCB-LOCA .GT. 1)  GOTO 11
             LOCA    =  MIN0( MAX0(LOCA,LMIN), LMAX-1 )
             ISHIFT  =  (LOCA - LMIN) * ISTEP
             ETA     =  (X - ENT(LOCA)) / (ENT(LOCA+1) - ENT(LOCA))
             GOTO 30
  20         ISHIFT  =  (LOCC - LMIN) * ISTEP
  21         DO 22  K  =  1, KNOTS
                INDEX(K)  =  INDEX(K) + ISHIFT
  22            CONTINUE
             GOTO 90
  30         DO 31  K  =  1, KNOTS
                INDEX(K)         =  INDEX(K) + ISHIFT
                INDEX(K+KNOTS)   =  INDEX(K) + ISTEP
                WEIGHT(K+KNOTS)  =  WEIGHT(K) * ETA
                WEIGHT(K)        =  WEIGHT(K) - WEIGHT(K+KNOTS)
  31            CONTINUE
             KNOTS  =  2*KNOTS
  90         ISTEP  =  ISTEP * NDIM
 100         CONTINUE
          DO 200    K  =  1, KNOTS
             I  =  INDEX(K)
             DHFINT  =  DHFINT + WEIGHT(K) * TABLE(I)
 200         CONTINUE
          RETURN
 300      Continue
c          CALL KERMTR('E104.1',LGFILE,MFLAG,RFLAG)
          IF(MFLAG) THEN
             IF(LGFILE .EQ. 0) THEN
                WRITE(*,1000) NARG
             ELSE
                WRITE(LGFILE,1000) NARG
             ENDIF
          ENDIF
          IF(.NOT. RFLAG) CALL ABEND
          RETURN
1000      FORMAT('  HFUNCTION DHFINT ... NARG =',I6,
     +              '  NOT WITHIN RANGE')
          END
      
C #include "/Users/jung/jung/cvs/pythia6410/pydata.f"
      subroutine p_gbw(X,WPART)
      IMPLICIT None
#include "casgki.inc"
#include "caslujets.inc"
#include "casluco.inc"
#include "caspart.inc"
#include "caspara1.inc"
#include "caspara.inc"
#include "casvalues.inc"
#include "caspara5.inc"
#include "casshwr.inc"
#include "casmult.inc"
      Double Precision X,WPART,WMAT
      DIMENSION X(20)
      Integer ndimc
      COMMON /DIMEN/ NDIMC
      Integer NDIM
      COMMON/INTDIM/ NDIM
      Integer NDIMEN,NRN
      Double precision kt2min,kt2max,kt2,z1,z,mq,mc,weight
      Double precision charge,alpha_s,xbj,sh
      Double precision phi_L,phi_T,phic_L,phic_T
      Double Precision XGSF,SFSCAL2,STMAX2,XPQ2(-6:6),ALPH_EM
      Double Precision scale,sigm0,lam,x0,r02,glu 
	logical first
      double precision au 
      common/f2fit/au(50),first
      
      Double Precision alphasca
      
      integer imass
c treat charm in massless (imass=0) or massive approach (imass=1)    
      imass = 0

      q2 = qmi
      xbj=qmi/ymi/sss
      

      kt2min = 0.00001
      kt2max = Q2
      kt2max = 500000.
c      kt2max = 500.
      
      NDIMEN = NDIM
      
      NRN = 0

      NRN = NRN + 1      
      kt2=kt2min*(kt2max/kt2min)**x(NRN)
      weight = kt2*log(kt2max/kt2min)
      
      NRN = NRN + 1
      z = X(NRN)
      NRN = NRN + 1
      z1 = X(NRN)
c      write(6,*) x(1),x(2),x(3)
c	write(6,*) ' dummy routine for p_gbw',q2,xgsf,ymi,sss,nrn

      charge = 1. + 4. + 1.
c      mq=0.140
c      mq=0.07
      mq = au(8)
      mc = au(9)
c      mq=0.07
c      mq=0.0
c      mq=0.1
c	write(6,*) ' mq ,mc ',mq,mc
      sh = -Q2 -kt2 +Q2/z
c      if(sh.le.0) goto 345
      charge = charge/9.
      alpha_s=ALPHASCA(sqrt(kt2))
      alpha_s=0.2
	ALPH_EM = ALPH
	IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)      
	phi_L=8d0*alpha_s*ALPH_EM*charge*
     &	(kt2*z*(1d0-z)*z1**2*(1-z1)**2*Q2)/
     &      ((Q2*z1*(1d0-z1)+kt2*z*(1-z)+mq**2)*(Q2*z1*(1-z1)+mq**2))

      phi_T=alpha_s*ALPH_EM*charge*
     &      (kt2*Q2*(z1**2 +(1-z1)**2)*z1*
     &      (1d0-z1)*(z**2 +(1d0-z)**2)+
     &       kt2*(z1**2 +(1-z1)**2)*mq**2+4d0*z1*(1-z1)*z*(1-z)*mq**2)/
     &      ((Q2*z1*(1d0-z1)+kt2*z*(1d0-z)+mq**2)*
     &      (Q2*z1*(1d0-z1)+mq**2))
c      if(sh.gt.4*mc**2) then  ! including charm
      if(sh.gt.0) then  ! including charm
         charge = 4./9
c for the massless limit set mc = mq
         if(imass.eq.0) mc = mq
 	   phic_L=8d0*alpha_s*ALPH_EM*charge*
     &	(kt2*z*(1d0-z)*z1**2*(1-z1)**2*Q2)/
     &      ((Q2*z1*(1d0-z1)+kt2*z*(1-z)+mc**2)*(Q2*z1*(1-z1)+mc**2))

         phic_T=alpha_s*ALPH_EM*charge*
     &      (kt2*Q2*(z1**2 +(1-z1)**2)*z1*
     &      (1d0-z1)*(z**2 +(1d0-z)**2)+
     &       kt2*(z1**2 +(1-z1)**2)*mc**2+4d0*z1*(1-z1)*z*(1-z)*mc**2)/
     &      ((Q2*z1*(1d0-z1)+kt2*z*(1d0-z)+mc**2)*
     &      (Q2*z1*(1d0-z1)+mc**2))
     
      endif 
c include fator of pi coming from d^2k = pi dk^2     
cc      wmat = 3.1415*(phi_T +phic_T + phi_L + phic_L)/kt2
cc	phi_L = 0
      wmat = 3.1415*(phi_T  + phi_L )/kt2

      wmat = wmat * weight
c      write(6,*) ' WMAT = ',WMAT,weight,PYALEM(Q2),pi
      
      sfscal2= kt2
      stmax2 = 10. 
      
c       write(6,*) ' 1  NDIM = ',NDIM,NRN
      NDIMC = NRN
      NDIM = NRN
ccc      XGSF = xbj*(q2+4*mq**2)/Q2
      XGSF = xbj
c      XGSF = xbj*(q2+mq**2)/Q2
c      write(6,*) xbj,xgsf
      XGSF = xbj/z
      stmax2 = sqrt(Q2/z)  ! 
      call cauniglu(2212,XGSF,SFSCAL2,STMAX2,XPQ2)
cc    F_2(x,Q^2)=Q**2/(4d0*pi**2)\int_0^Q2 d^2k/k2 *(phi_L+phi_T)*A(x,k2)
c       write(6,*) ' 1  NDIM = ',NDIM,NRN
      WPART = WMAT*XPQ2(0)
c      write(6,*) ' WPART = ',WPART,WMAT,xpq2(0)
      RETURN
 345  Continue
      Wpart = 0.
      Return     
      END
      subroutine p_new(X,WPART)
      IMPLICIT None
#include "casgki.inc"
#include "caslujets.inc"
#include "casluco.inc"
#include "caspart.inc"
#include "caspara1.inc"
#include "caspara.inc"
#include "casvalues.inc"
#include "caspara5.inc"
#include "casshwr.inc"
#include "casmult.inc"
      Double Precision X,WPART,WMAT
      DIMENSION X(20)
      Integer ndimc
      COMMON /DIMEN/ NDIMC
      Integer NDIM
      COMMON/INTDIM/ NDIM
      Integer NDIMEN,NRN
      Double precision qt2min,qt2max,qt2,qtp2,z,mq,weight
      Double precision kt2min,kt2max,kt2,phi,zmin,zmax
      Double precision charge,alpha_s,xbj
      Double precision phi_L,phi_T,Pqg,Tr,xq,sig,sh,mass
      Double Precision XGSF,SFSCAL2,STMAX2,XPQ2(-6:6),ALPH_EM
      Double Precision scale,sigm0,lam,x0,r02,glu 
      Integer Icolinear
      
      Double Precision alphasca
      
      Icolinear = 0 ! use full qt dependent splitting fct (eq 4.9)
      Icolinear = 1 ! use integrated spitting fct (eq 4.6) 
      NDIMEN = NDIM
      
      charge = 1. + 4. + 1.   ! w/o charm
c      charge = 1. + 4. + 1. + 4  ! including charm
      Tr = 0.5 
      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)      
      NRN = 0
      sig = 0.
      q2 = qmi
      xbj=qmi/ymi/sss
      
      NRN = NRN + 1
      zmin = xbj
      zmax = 0.99
c      zmax = 0.9
      z = zmin + X(NRN)*zmax
      NRN = NRN + 1      
      kt2min = 0.00001
c      kt2min = 0.01
c      kt2max = 5000.
      kt2max = Q2/z
      kt2=kt2min*(kt2max/kt2min)**x(NRN)
      weight = kt2*log(kt2max/kt2min)*(zmax-zmin) 
      sh = -Q2 -kt2 +Q2/z
      mass = 1.5
c      mass = 1.3
cc      if(sh.gt.4*mass**2) charge = 1. + 4. + 1. + 4  ! including charm
      if((Q2 -z*kt2).lt.0.) then
         write(6,*) ' q2 -zkt2 < 0 ',q2,z*kt2
         goto 123 
      endif
      If(Icolinear.eq.1) then
c         alpha_s=ALPHAS(sqrt(kt2))
         alpha_s=ALPHASCA(sqrt(q2))
         alpha_s= 0.2
         Pqg = ((1.-z)**2+z**2)*log(Q2/z/kt2) 
         Pqg = Pqg -(1.-6.*z*(1.-z))*(1.-z*kt2/Q2)
         xq =alpha_s/2./pi*Pqg*Tr * pi 
c         write(6,*) ' in icoll =1 ',alpha_s,kt2
      else 
         qt2min = 0.000001
         qt2max = (1.-z)*(Q2 -z*kt2)
         if(qt2max.lt.qt2min) goto 123 
         NRN = NRN + 1      
         qt2=qt2min*(qt2max/qt2min)**x(NRN)
         weight = qt2*log(qt2max/qt2min)*weight 
         alpha_s=ALPHASCA(sqrt(kt2))
         alpha_s=ALPHASCA(sqrt(qt2))
         alpha_s=ALPHASCA(sqrt(qt2)/(1-z))
c         alpha_s=ALPHAS(sqrt(q2))
        alpha_s=0.2
c        qtp2 = qt2 + z**2*kt2 - 2*z*sqrt(qt2)*sqrt(kt2)*cos(phi)
         qtp2 = qt2 
         Pqg = (1.-z)**2+z**2+4.*z**2*(1.-z)**2*kt2/qtp2
         Pqg = (qtp2/(qtp2+z*(1.-z)*kt2))**2*Pqg*Tr
c factor of pi coming from kt integration ?
         xq =alpha_s/2./pi*Pqg * pi 
         xq = xq/qtp2
      
      endif

      if(xq.lt.0) then
         write(6,*) ' p_new: xq < 0 ',xq,qtp2
      endif
      charge = charge/9.
      xq = xq * weight 
     
      sfscal2= kt2
      stmax2 = 10. 
c      stmax2 = kt2*2 
c     qbar = sqrt (xg**2* \ksi *s) = sqrt(xg * y *s) = sqrt(Q2/z) 
      stmax2 = sqrt(Q2/z)  ! 
     
c       write(6,*) ' 1  NDIM = ',NDIM,NRN
      NDIMC = NRN
      NDIM = NRN
      XGSF = xbj/z
      call cauniglu(2212,XGSF,SFSCAL2,STMAX2,XPQ2)
cc    F_2(x,Q^2)=Q**2/(4d0*pi**2)\int_0^Q2 d^2k/k2 *(phi_L+phi_T)*A(x,k2)
c       write(6,*) ' 1  NDIM = ',NDIM,NRN
      xq = xq * XPQ2(0)
      sig = xq * charge * 4. * pi**2 * ALPH_EM/q2
 123  continue      
      WPART = sig  
      if(wpart.lt.0) then 
        write(6,*) ' p_new: WPART = ',WPART,WMAT,xpq2(0),xgsf,
     &     SFSCAL2,STMAX2
      endif
      RETURN
      END
      subroutine p_full(X,WPART)
      IMPLICIT None
#include "casgki.inc"
#include "caslujets.inc"
#include "casluco.inc"
#include "caspart.inc"
#include "caspara1.inc"
#include "caspara.inc"
#include "casvalues.inc"
#include "caspara5.inc"
#include "casshwr.inc"
#include "casmult.inc"
      Double Precision X,WPART,WMAT
      DIMENSION X(20)
      Integer ndimc
      COMMON /DIMEN/ NDIMC
      Integer NDIM
      COMMON/INTDIM/ NDIM
      Integer NDIMEN,NRN
      
      Double precision charge
      Double precision sig,kt2
      Double Precision XGSF,SFSCAL2,STMAX2,XPQ2(-6:6)
	Double Precision ECM,WT,Wt2
	Double Precision PZ
	Double Precision KT_1,KT_2,X1,X2,x1bar,x2bar
	Double Precision CMS(4)
      Double Precision phi1,phi2
      Double Precision ktmin,ktmax, weight
      
      DOUBLE PRECISION K1K2,K1P2,K1P3,K1P4,K2P1,K2P3,K2P4,K12,K22
      DOUBLE PRECISION P1P2,P1P3,P1P4,P2P3,P2P4,M2,PT2Q,QF2
      Double Precision Scale,alpha_s,alph_em,sh,th,uh,xsi_hard
	Double Precision Aab,Anab,al,cab,epsilon,mq,mc,wt1,xbj
	Double precision DOT,Dot1
      Double Precision XY
      DOUBLE PRECISION STHETA,SPHI
      Double Precision sL,rt,ru
      
      Integer KIA1,KIA2,KF1,KF2
      Double Precision phi,cos_fin,pcm_fin,pcm_new,sh_new
      DOUBLE PRECISION DSHN(4)
      Double Precision dcasrn
      double Precision cphi,test,E1,E2,cut
      DATA CUT /1.d-6/ ! as smallx
      Logical corr
c      Data corr/.true./ ! std	
      Data corr/.false./
      Integer ibar
c      Data ibar/1/  ! std
      Data ibar/0/
      
      COMMON/XVAR/ XY(10)
	Integer I,J,KPF,NP,ncount
	logical first
      double precision au 
      common/f2fit/au(50),first
      
      Double Precision alphasca
      
      mq = au(8)
      mc = au(9)
      
      ktmin = 0.000001d0
      ktmax = sqrt(sss)/2.
      
      xbj = qmi/ymi/sss
      xmin = 0.0000001
c      xmin = dmax1(4*mq**2,qmi)/ymi/sss
      xmin = 4*mq**2/ymi/sss
      xmax = 0.99
      
      NRN = 0
      sig = 0.

      kt_1 = qmi*(1.-ymi)
      if(ibar.eq.0) kt_1 = qmi
      x1=ymi
      
      
      if(x(1).ne.x(1)) then
        write(6,*) ' fatal ', (x(i),i=1,10)
      endif
c incoming particle momenta
	PZ = Sqrt(sss)/2.
c      write(6,*) ' s old ',sss
      K(1,2)=11
      K(2,2)=2212
      K(1,1)=21
      K(2,1)=21
      P(1,1) = 0d0
	P(1,2) = 0d0
	P(1,3) = -PZ
	P(1,4) = DABS(PZ)
	P(1,5) = 0d0
	P(2,1) = 0d0
	P(2,2) = 0.d0
	P(2,3) = PZ
	P(2,5) = 0d0
	P(2,4) =DABS(PZ)	
      CMS(1)=  P(1,1) + P(2,1)
      CMS(2)=  P(1,2) + P(2,2)
      CMS(3)=  P(1,3) + P(2,3)
      CMS(4)=  P(1,4) + P(2,4)
      SSS=DOT(CMS,CMS)	
c      write(6,*) ' s new ',sss
c call Phase-space 
      NIA1=3
	NIA2=4	
      K(Nia1,1)=21
      K(Nia2,1)=21
	KPF=5
	KPA = KPF
      
      IPRO = 11
      IQ2 = 3      
     
      
      AM(1) = 1.5d0
      AM(1) = 0.25d0
      AM(2) = AM(1)
c	if(kt_1.le.0.000001) kt_1=0.000001
	NRN = NRN + 1 
      phi1=2.d0*PI*X(NRN)	
	NRN = NRN + 1
      phi2=2.d0*PI*X(NRN)
      NRN = NRN + 1
      kt_2 = ktmin*(ktmax/ktmin)**X(NRN)	
      Weight = kt_2 * log(ktmax/ktmin)
      NRN = NRN + 1
      x2 = xmin * (xmax/xmin)**X(NRN)
      Weight = Weight * log(xmax/xmin)
      if(ibar.eq.0) then
c         x1bar = -kt_1/(1.-x1)/sss
         x1bar = 0.
         x2bar = 0.
         else
         x1bar = -kt_1/(1.-x1)/sss
         x2bar = -kt_2/(1.-x2)/sss
      endif
c      write(6,*) ' x1bar ',x1bar,x2bar
	P(Nia1,1) = Dsqrt(KT_1)*cos(phi1)
	P(Nia1,2) = Dsqrt(KT_1)*sin(phi1)
	P(Nia1,3) = x1*P(1,3) + x1bar*P(2,3)
	P(Nia1,4) = x1*P(1,4) + x1bar*P(2,4)
	
	P(Nia1,5) = -Dsqrt(KT_1)
      K(Nia1,2)=22
c	if(kt_2.le.0.000001) kt_2=0.000001
c 2pi from each phi integration      
ccc      Weight = Weight *  pi**2 
	P(Nia2,1) = Dsqrt(KT_2)*cos(phi2)
	P(Nia2,2) = Dsqrt(KT_2)*sin(phi2)
	P(Nia2,3) = x2*P(2,3) + x2bar*P(1,3)
	P(Nia2,4) = x2*P(2,4) + x2bar*P(1,4)
	P(Nia2,5) = -Dsqrt(KT_2)
      K(Nia2,2)=21
c      write(6,*) ' 1st kt ',P(NIA1,1)**2+P(NIA1,2)**2,
c     &  P(NIA2,1)**2+P(NIA2,2)**2
c      write(6,*) ' test ',(-dot1(nia1,nia1)-dot1(nia2,nia2)+4*mq**2)/x1/sss,x2
      ncount = ncount + 1
	if(ncount.lt.10) then
c         write(6,*) ' 1st kt^2: ',kt_1,kt_2,ncount
      endif
	NRN = NRN + 1
      XY(1) = X(NRN)
      NRN = NRN + 1
      XY(2) = X(NRN)
         
       
      NF1=NIA1+2
      NF2=NIA1+3
      K(NF1,1)=2
      K(NF1,2)=KPF
      K(NF2,1)=1
      K(NF2,2)=-KPF
      K(NF1,3)=NIA1
      K(NF2,1)=1
      K(NF2,3)=NIA1
	
      DBCMS(1)=  P(NIA1,1) + P(NIA2,1)
      DBCMS(2)=  P(NIA1,2) + P(NIA2,2)
      DBCMS(3)=  P(NIA1,3) + P(NIA2,3)
      DBCMS(4)=  P(NIA1,4) + P(NIA2,4)
      SHAT=DOT(DBCMS,DBCMS)	      
      if(shat.ne.shat) then
         write(6,*) ' shat ',shat,' ecm ',ecm, ' pi ',pi
         write(6,*) ' P(nia1) ', (P(NIA1,I),i=1,5)
         write(6,*) ' P(nia2) ', (P(NIA2,I),i=1,5)
         write(6,*) ' x1,kt_1 ',x1,kt_1,' phi ',phi1,phi2,' nrn ',nrn
         write(6,*) ' X() ',(x(i),i=1,6)
      endif
      NP=2
	N=6
      wt1 = 0.
      do j = 1,2
         wt = 0.
         if(j.eq.1) then 
            AM(1) = mq
            QF2 = 2./3.
            elseif(j.eq.2) then
            AM(1) = mc
            QF2 = 4./9.
         endif
         AM(2) = AM(1)
         P(NF1,5)=AM(1)
         P(NF2,5)=AM(1)
  
	   IF(SHAT.LE.0) goto 170
	   IF(SHAT.LE.4*AM(1)**2) goto 170
         ECM =DSQRT(SHAT)
         sh_new = shat
c         if(x2.le.xbj) then
c      write(6,*) ' test ',(-dot1(nia1,nia1)-dot1(nia2,nia2)+4*mq**2)/x1/sss,x2
c         write(6,*) ' xmin,x2, xbj, sqrt(shat) ',xmin,x2,xbj,sqrt(shat)
c         endif

         CALL PYROBO(0,0,0.D0,0.D0,-DBCMS(1)/DBCMS(4), -DBCMS(2)/DBCMS(4), 
     +     -DBCMS(3)/DBCMS(4))
         SPHI = PYANGL(P(NIA1,1),P(NIA1,2))
         call PYROBO(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
         STHETA = PYANGL(P(NIA1,3),P(NIA1,1))
         call PYROBO(0,0,-STHETA,0.D0,0.d0,0.d0,0.d0)
         CALL PHASE(NP,ECM,AM,PCM,WT)
c         write(6,*) ' phase wt = ',wt
         IF(WT.LE.0.D0) GOTO 170
         DO 80  I=1,4
            P(NF1,I)=PCM(I,1)
            P(NF2,I)=PCM(I,2)
   80    CONTINUE
c         PCM_FIN=SQRT(P(NF1,4)**2 - AM(1)**2)
c         COS_FIN=P(NF1,3)/PCM_FIN
         call PYROBO(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
         call PYROBO(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
         PCM_FIN=SQRT(P(NF1,4)**2 - AM(1)**2)
         COS_FIN=P(NF1,3)/PCM_FIN
         CALL PYROBO(0,0,0.D0,0.D0,DBCMS(1)/DBCMS(4), DBCMS(2)/DBCMS(4), 
     +     DBCMS(3)/DBCMS(4))
         if(shat.ne.shat) then
           write(6,*) ' shat ',shat,' ecm ',ecm
           write(6,*) dbcms
           stop
         endif
         PT2Q = P(NF2,1)**2 + P(NF2,2)**2
c         call pylist(1)
c      write(6,*) ' 2nd kt ',P(NIA1,1)**2+P(NIA1,2)**2,
c     &  P(NIA2,1)**2+P(NIA2,2)**2
         if(corr) then
            KIA1=N+1
            KIA2=N+2
            KF1=N+3
            KF2=N+4
            N=N+4
            K(KIA1,1)=13
            K(KIA2,1)=13
            K(KF1,1)=13
            K(KF2,1)=13
            DO I=1,4
               P(KIA1,I) = P(NIA1,I)
               P(KIA2,I) = P(NIA2,I)
               P(KF1,I) = P(NF1,I)
               P(KF2,I) = P(NF2,I)
            ENDDO
            P(KIA1,3) = 0.5d0*(P(KIA1,3)-P(KIA1,4))
            P(KIA1,4) = DABS(P(KIA1,3))
            P(KIA2,3) = 0.5d0*(P(KIA2,3)+P(KIA2,4))
            P(KIA2,4) = DABS(P(KIA2,3))
            DO I=1,4
              DSHN(I)= P(KIA1,I) + P(KIA2,I)
            ENDDO
            SH_NEW = DOT(DSHN,DSHN)
            CALL PYROBO(0,0,0.D0,0.D0,-DSHN(1)/DSHN(4), -DSHN(2)/DSHN(4),
     +      -DSHN(3)/DSHN(4))
            P(KF1,4) = 0.5d0*DSQRT(SH_NEW)
            PCM_NEW=SQRT(P(KF1,4)**2 - AM(1)**2)
            P(KF1,3) =PCM_NEW*COS_FIN
            phi=2.d0*PI*dcasrn()         
            P(KF1,1) = DSQRT(PCM_NEW**2-P(KF1,3)**2)*cos(phi)
            P(KF1,2) = DSQRT(PCM_NEW**2-P(KF1,3)**2)*sin(phi)

            DO I=1,4
               P(KF2,I) = P(KIA1,I) + P(KIA2,I) - P(KF1,I)
            ENDDO
            CALL PYROBO(0,0,0.D0,0.D0,DSHN(1)/DSHN(4) ,DSHN(2)/DSHN(4),
     +      DSHN(3)/DSHN(4))
c            call pylist(1)
            DO I=1,4
               P(NIA1,I) = P(KIA1,I)
               P(NIA2,I) = P(KIA2,I)
               P(NF1,I) = P(KF1,I)
               P(NF2,I) = P(KF2,I)
            ENDDO
c            write(6,*) N
	      N=6
         endif
         NDIM = NRN
         NDIMC = NRN
         K1K2 = DOT1(NIA1,NIA2)
         K1P2 = DOT1(NIA1,2)
         K1P3 = DOT1(NIA1,NF2)
         K1P4 = DOT1(NIA1,NF1)
         K2P1 = DOT1(NIA2,1)
         K2P3 = DOT1(NIA2,NF2)
         K2P4 = DOT1(NIA2,NF1)
         K12 = DOT1(NIA1,NIA1)
         K22 = DOT1(NIA2,NIA2)
         P1P2 = DOT1(1,2)
         P1P3 = DOT1(1,NF2)
         P1P4 = DOT1(1,NF1)
         P2P3 = DOT1(2,NF2)
         P2P4 = DOT1(2,NF1)
         M2 = AM(1)**2
         TH = AM(1)**2 - 2.D0 * DOT1(NIA1,NF1) + K12
         UH = AM(1)**2 - 2.D0 * DOT1(NIA2,NF1) + K22
         SH = K12 + 2.D0 * DOT1(NIA1,NIA2) + K22
         test=Sh_new+th+uh-2*AM(1)**2-K12-K22
c         write(6,*) ' kt ',kt_1,k12,kt_2,k22
	
	
      if(abs(test).gt.1e-5) then
         write(6,*) ' sh+th+uh ,m**2 ',Sh_new+th+uh,2*AM(1)**2+K12+
     +   K22,AM(1)
         write(6,*) ' sh+th+uh-m**2 ',Sh_new+th+uh-2*AM(1)**2-K12-
     +   K22
         write(6,*) ' sh,th,uh ',Sh_new,th,uh
      endif
c      write(6,*) ' kt ',P(NIA1,1)**2+P(NIA1,2)**2,k12,
c     &  P(NIA2,1)**2+P(NIA2,2)**2,k22

c         PT2Q = P(NF2,1)**2 + P(NF2,2)**2
         SCALE = PT2Q + 4.*M2
         IF(IQ2.EQ.1) THEN
            SCALE = Real(4.D0*M2)
         ELSEIF(IQ2.EQ.2) THEN
            SCALE = Real(SH_NEW)
         ELSEIF(IQ2.EQ.3) THEN
            SCALE = Real(4D0*M2 + PT2Q)
c            SCALE = Real(M2 + PT2Q)
         ELSEIF(IQ2.EQ.4) THEN
            SCALE = Real(Q2)
         ELSEIF(IQ2.EQ.5) THEN
            SCALE = Real(Q2 + PT2Q + 4.D0*M2)
         ELSEIF(IQ2.EQ.6) THEN
            SCALE = Real(DABS(K22))
c this is for testing
         ELSEIF(IQ2.EQ.7) THEN
            write(6,*) ' iq2 = 7 incorrect for p_full'
            stop
            sL=2d0*k1p2*k2p1/p1p2
            rt=1d0/th
            ru=1d0/uh
            Scale=real(2d0*(sL*rt*ru)/((sL*rt*ru)**2+rt**2+ru**2))
         ELSE
            WRITE(6,*) ' NO VALID Q2 SCALE. STOP'
            STOP
         ENDIF
         alpha_s = ALPHASCA(sqrt(scale))
         ALPH_EM = ALPH
	   IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(KT_1)
         E1 = -K12/M2
         E2 = -K22/M2
         AL = 0
         IF(E1.GT.CUT.AND.E2.LT.CUT) THEN
c           write(6,*) 'on mass e1,e2 ',e1,e2
C DIS, gluon on mass
C eq. B21, B22
           Aab = (K1K2/K2P1)**2*((P1P3**2+P1P4**2)/K2P4/K2P3 - 2D0*M2/
     +     K12*(P1P4/K2P3-P1P3/K2P4)**2)
           Anab = 0D0
         if(Aab.le.0) then
           write(6,*) ' on mass Aab = ',Aab,' am = ',am(1),e2,m2,k22
c           call pylist(1)
c           stop 
         endif
         ELSEIF(E1.GT.CUT.AND.E2.GT.CUT) THEN
c            write(6,*) 'off mass e1,e2 ',e1,e2
          Aab = 4D0*(K1P2*K2P1/P1P2**2)**2*(P1P2**2/(TH-M2)/(UH-M2) -
     +     1D0/K12/K22*(P1P2+2D0*P1P4*P2P3/(TH-M2)+ 2D0*P1P3*P2P4/(UH-
     +     M2))**2)
           aL = Aab*(1.d0 -1d0/3d0 - 8d0/3d0/Aab*k2p1**2/p1p2**2*(
     +     (p2p4**2+p2p3**2)/(TH-M2)/(UH-M2) - (2d0*M2 + k12)/k22*
     +     (p2p3/(TH-M2) - p2p4/(UH-M2))**2))
           Cab = 16.D0*PI*PI*ALPHA_S*ALPH_EM*QF2
           epsilon = ymi**2/2d0/(1d0-ymi+ymi**2/2d0)
           if(Aab.lt.0)  write(6,*) ' Aab = ',Aab,' AL = ',al,epsilon
           Aab = Aab - epsilon*AL
         endif 
         if(Aab.le.0) then
           write(6,*) ' Aab = ',Aab,' am = ',am(1),e2,m2,k22
c           call pylist(1)
c           stop 
         endif
c         Aab=1
c         Cab = 1
	   wt1 =  wt1 + Cab*Aab*wt
170	   Continue      
c         call pylist(1)
c	   IF(SHAT.LE.4*AM(1)**2) then
c         	write(6,*) ' shat < mass ', wt1,wt,am(1),np
c         endif
         sfscal2= Dabs(P(NIA2,1)**2+ P(NIA2,2)**2)
	   if(sfscal2.le.1e-3) sfscal2=1.d-3 
         xsi_hard = (x1)/(x2)
         stmax2 = x2*sqrt(xsi_hard*SSS)
      end do
      weight = weight * wt1 
c       write(6,*) ' 1  NDIM = ',NDIM,NRN
         
      XGSF = x2
      call cauniglu(2212,XGSF,SFSCAL2,STMAX2,XPQ2)
      sig = weight * XPQ2(0)/2./x1/x2/sss/4./pi**2
      WPART = sig  
      if(wpart.lt.0) then 
        write(6,*) ' p_full: WPART = ',WPART,xpq2(0),xgsf,
     &     SFSCAL2,STMAX2
      endif
      RETURN
 180  Continue
      Wpart = 0.
	Return
      END
      
      
      SUBROUTINE cauniglu(KF,XI,KT2I,TI,XPQ)
      Implicit None
      Integer kf,i,Iccfmset
      Double Precision  XI,KT2I,TI
      Double Precision  XPQ(-6:6)
      Real XPQS(-6:6),XPQQ(-6:6)
      Real X,KT2,T,scale
      Double Precision DGRV_NLN
      DOUBLE PRECISION DGRV_NL,BLUEML,d_xgx,XGLU
      Double Precision SIGM0,LAM,X0,R02,alpha_s,glu
      Integer IGLU
      Common/CAGLUON/Iglu
      Integer ISEL
      Double Precision UPV,DNV,USEA,DSEA
      Common /collpdf/UPV,DNV,USEA,DSEA 
	Double Precision xminpdf,xmaxpdf,q2minpdf,q2maxpdf
	common/W50513/xminpdf,xmaxpdf,q2minpdf,q2maxpdf
      EXTERNAL DGRV_NLN,DGRV_NL,BLUEML,D_XGX,XGLU
      Double Precision d_dnv,d_upv
      External d_dnv,d_upv
      LOGICAL FIRST,PDFFIRST
      Data PDFFIRST/.TRUE./
      Data FIRST/.TRUE./
      Integer Iquark
      Data Iquark/-1/
c      Data Iquark/0/
      Integer ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      SAVE ISEL  

c ISEL = 1 CCFM
c      = 2 derivative of collinear gluon
c      = 3 Bluemlein
c      = 4 KMS
c      = 5 Golec-Biernat Wuesthoff
c      = 6 KMR
	qg0=max(0.3,qg0)
      IF(FIRST) THEN
	   If(IGLU.GT.1000) then
	   ISEL = Int(IGLU/1000)
c	   write(6,*) ' isel ',ISEL
	   ELSE
         ISEL = IGLU
	   ENDIF
         FIRST=.FALSE.
         Write(6,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
         Write(6,*) '++++  gluon density selection    ++++++++++++++'
         Write(6,*) '++++  ISEL= ',ISEL,
     +      '                  ++++'
         If(Isel.eq.1) then 
	     if(iglu.lt.1000) then
           Write(6,*) '++++  ISEL=1 --> CCFM   JS2001             ++++'
	     elseif(iglu.gt.1000) then
	        Iccfmset =  MOD(IGLU,1000)
     	        write(6,*) '++++  Iccfmset = ',Iccfmset
	        if(iccfmset.eq.1) then
           Write(6,*) '++++  CCFM  J2003 set 1 (Ipgg0ns-1)        ++++'
	        elseif(iccfmset.eq.2) then
           Write(6,*) '++++  CCFM J2003 set 2 (Ipgg1ns1)          ++++'
	        elseif(iccfmset.eq.3) then
           Write(6,*) '++++  CCFM J2003 set 3 (Ipgg2ns2)          ++++'
	        elseif(iccfmset.eq.10) then 
           Write(6,*) '++++  CCFM set A0                          ++++'
	        elseif(iccfmset.eq.11) then
           Write(6,*) '++++  CCFM set A0+                         ++++'
	        elseif(iccfmset.eq.12) then
           Write(6,*) '++++  CCFM set A0-                         ++++'
	        elseif(iccfmset.eq.13) then
           Write(6,*) '++++  CCFM set A1                          ++++'
	        elseif(iccfmset.eq.20) then
           Write(6,*) '++++  CCFM set B0                          ++++'
	        elseif(iccfmset.eq.21) then
           Write(6,*) '++++  CCFM set B0+                         ++++'
	        elseif(iccfmset.eq.22) then
           Write(6,*) '++++  CCFM set B0-                         ++++'
	        elseif(iccfmset.eq.23) then
           Write(6,*) '++++  CCFM set B1                          ++++'
	        elseif(iccfmset.eq.101) then
           Write(6,*) '++++  CCFM set C                           ++++'
	        elseif(iccfmset.eq.111) then
           Write(6,*) '++++  CCFM test set                        ++++'
	        elseif(iccfmset.eq.113) then
           Write(6,*) '++++  CCFM new test set                        ++++'
	        else
	     Write(6,*) '++++ no valid set specified:',iglu,'++++'
	        endif
	     endif
	     elseif(isel.eq.2) then
         Write(6,*) '++++  ISEL=2 --> derivative of coll. gluon ++++'
	     elseif(isel.eq.3) then
         Write(6,*) '++++  ISEL=3 --> Bluemlein                 ++++'
	     elseif(isel.eq.4) then
         Write(6,*) '++++  ISEL=4 --> Kwiecinski,Martin,Stasto  ++++'
	     elseif(isel.eq.5) then
         Write(6,*) '++++  ISEL=5 --> Wuesthoff Golec-Biernat   ++++'
	     elseif(isel.eq.6) then
         Write(6,*) '++++  ISEL=6 --> Kimber,Martin,Ryskin      ++++'
	     elseif(isel.eq.7) then
         Write(6,*) '++++  ISEL=7 --> Ryskin,Shabelski          ++++'
	   endif
         Write(6,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
      ENDIF

      X = real(XI)
      KT2 = real(KT2I)
c for new p_max
      T=max(0.0001,real(TI))
      DO 10 I=-6,6
         XPQ(i)=0.d0
   10 XPQS(I)=0.0
      SCALE = KT2
      If(ISEL.EQ.1) THEN
	   if(iabs(kf).eq.2212) then 
            call ccfmstfu(kf,Real(x),scale,t,xpqs)
	   elseif(kf.eq.22.or.iabs(kf).eq.11) then
            call ccfmstfg(kf,Real(x),scale,t,xpqs)
	   else
            write(6,*) ' cauniglu: no valid kf = ',kf	      
	   endif 
         if(iquark.eq.0) then
            call ccfmstfq(kf,Real(x),scale,t,xpqq)
            xpqs(1) = xpqq(1)
            xpqs(2) = xpqq(2)
            
         else
c	   XPQS(1) = dnv/q2minpdf
c	   XPQS(2) = upv/q2minpdf
c	   XPQS(1) = dnv/scale 
c	   XPQS(2) = upv/scale
c         SCALE = MAX(q2minpdf,KT2)
cccc         call xpdf(dble(x),dble(scale))
c         call xpdf(dble(x),dble(t**2))
c  	      XPQS(1) = (dnv)/scale
c	      XPQS(2) = (dnv)/scale
         endif
c         write(6,*) ' cauniglu kf,x,scale,t ',kf,x,scale,t
c         write(6,*) ' cauniglu ',xpqs(0),xpqs(1),xpqs(2)
      ELSE
         write(6,*) ' ISEL = ',ISEL,' not implemented '
         STOP
      ENDIF
      DO I=-6,6
         XPQ(I)=DBLE(XPQS(I))
      ENDDO
      RETURN
      END

      FUNCTION ALPHASCA(RQ)
      IMPLICIT NONE
#include "caspara.inc"
#include "casludat1.inc"
#include "casvalues.inc"
#include "casscal.inc"
	Double Precision ALPHASCA,RQ,RQ1,PYALPS,ALAM
      Integer ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel

      Double Precision SCALE,Qgcut
	Double Precision Qscal
      Double Precision test_scale,test_as
	Integer Iqqbar
	Common/GLUDAT2/Qscal,Iqqbar
	Logical alpfirst
	Data alpfirst/.true./
      RQ1 = RQ
      IF(IRUNA.EQ.0) THEN
         ALPHASCA=ALPHS
      ELSE

         ALAM=PARU(112)
	   QGcut=max(ALAM,Dble(Qg0))
         SCALE = RQ1*RQ1
	   SCALE = SCALE * SCALFA * Qscal
c          write(6,*) 'lambda ',alam
         IF(SCALE.LT.Qgcut**2) SCALE=Qgcut**2
c         if(sqrt(scale).lt.1.5)  	
c     &      write(6,*) 'lambda ',alam,qgcut,sqrt(scale),PYALPS(SCALE)
         ALPHASCA =PYALPS(SCALE)
	   If(alpfirst) then
	     write(6,*) ' alpha_s for ME: PARU(112) = ',PARU(112)
	     write(6,*) ' alpha_s for ME: MSTU(111) = ',MSTU(111)
	     write(6,*) ' alpha_s for ME: Scalfa = ',Scalfa,' Qscal ',Qscal
           test_scale = 91.188d0*91.188d0
           test_as = PYALPS(test_scale)
           write(6,*) ' alpha_s at M_z: ',test_as
	     alpfirst=.false.
	   endif
	
      ENDIF
      RETURN
      END
      
	Subroutine RAINT(FUNC,EPS,RESULT,ERR)
      Implicit None
      Real FUNC,RESULT,ERR,EPS,x,y
      Double Precision pyr
      External Func
      double Precision fsum,ffsum,gsum,ggsum,hsum,hhsum,sigma2
      Double Precision error,error2,f
      Integer n1,n2,np1,np2,npoints
      
      np1 = 200
      np2 = 20
      npoints = np1*np2
      
      fsum = 0.
      ffsum = 0.
      hsum = 0.
      hhsum = 0.
      
      do n1=0, np1-1
	    gsum = 0
          ggsum = 0
          do n2 = 0, np2-1
             x = pyr(1)
             y = x
             f = func(x)             
             fsum = fsum + f
             ffsum = ffsum + f*f
             x = (y+real(n1))/real(np1)
		 f = func(x)
             gsum = gsum + f
             ggsum = ggsum + f*f
           end do
           gsum = gsum/real(np2)
           ggsum = ggsum /real(np2)
	     sigma2 = ggsum - gsum*gsum
	     error2 = sigma2/real(np2-1)
	     hhsum = hhsum+error2
	     hsum = hsum+gsum/real(np1)
	end do

	fsum = fsum /real(npoints)
	ffsum = ffsum /real(npoints)
	sigma2 = ffsum - fsum*fsum
	error = sqrt(sigma2/real(npoints))
      
c      write(6,*) ' For one range of ',npoints,' points:', fsum,error
	error = sqrt(hhsum)/real(np1);
c      write(6,*) 'For ',np1,'ranges of each ',np2,' points:', hsum,error
      Result = hsum
      Err = error/max(1d0,hsum)
      if(err.gt.eps) then
         write(6,*) ' RAINT: desired precision ',EPS,' not reached:',err
      endif 
	return
      end
 
      SUBROUTINE meqg(WT1)
      IMPLICIT None
C
C   q g* ----> q
C
C         P1-----*
C                  *
C                    *-------- q
C                  *
C         P2-----*
C
C

#include "caslujets.inc"
#include "caspara1.inc"
#include "caspara.inc"
#include "casluco.inc"
#include "caspara5.inc"
#include "caspart.inc"
#include "casgki.inc"
#include "casshwr.inc"
#include "casmult.inc"
#include "casscal.inc"
      Integer IGENFL,IMEPHI
      COMMON/GENWEI/IGENFL,IMEPHI
      Double Precision XPQ1(-6:6),XPQ2(-6:6),SFSCAL1,SFSCAL2
	Double Precision STMAX1,STMAX2
	common/myscale/STMAX1,STMAX2
      DOUBLE PRECISION DSH(4),DSHN(4)
      Double Precision DScale
      Double Precision SH_NEW,M2,xh,sig,PT2Q
	Double Precision K1K2,K12,K22
      Double Precision scale
      Double Precision Sfcut
      DATA Sfcut /1.d-3/ 
c      DATA Sfcut /1.d-9/ 
	Double Precision WT1,SH,qg2,qg1,wtg,sum2
      Double Precision alph_em,pcm_new,alpha_s,alpha_em
      Double Precision yph,yph_bar,xsi_hard(2),gf
      Double Precision dot,dot1,dcasrn
      Double Precision rap1,rap2,rapp,ptstm,zz1,zz2,XPQC2
      Integer npfin
      Integer IQBAR
	Integer KIA1,KIA2,KF1,KF2,Kin1,Kin2
	Integer I
      Integer ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Double Precision E1,E2
	Double Precision CUT
      Double Precision QFT
      DATA CUT /1.d-6/ ! as smallx
c      DATA CUT /1.d-2/ ! 
	
      DOUBLE PRECISION SIN2W,XMW2
      COMMON/ELWEAK/SIN2W,XMW2
	Double Precision QQg,lambda
      Integer NFL
      Data IQBAR/1/
c      Data IQBAR/2/

      Logical ltest
c      data ltest/.true./
      data ltest/.false./
      Logical first
      data first/.true./
      
      Double Precision alphasca
      
      ALPH_EM = ALPH
c     call pylist(1)
      NPFIN = N
c	call pylist(1)
      SH = DOT1(NF1,NF1) + 2.D0 * DOT1(NF1,NF2) * DOT1(NF2,NF2)
      WT1 = 0.d0
      IF(SH.LE.0.) RETURN
c boost  cms
c      write(6,*) ' before boost'
c      call pylist(1)
      CALL PYROBO(0,N,0.D0,0.D0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
c      call pylist(1)
      KIA1=N+1
      KIA2=N+2
      KF1=N+3
      KF2=N+4
      Kin1=N+5
      Kin2=N+6
      N=N+6
            
c select incoming flavor
      IF(Iabs(KE).eq.11) then
        K(NIA2,2) = K(NF1,2)       
      else
        K(NIA1,2) = K(NF1,2)
      endif
      DO I=1,4
         P(KIA1,I) = P(NIA1,I)
         P(KIA2,I) = P(NIA2,I)
         P(KF1,I) = P(NF1,I)
         P(KF2,I) = P(NF2,I)
         DSH(I)= P(NIA1,I) + P(NIA2,I)
      ENDDO
      K(KIA1,1)=13
      K(KIA2,1)=13
      K(KF1,1)=13
      K(KF2,1)=13

      K(Kin1,1)=13
      K(Kin2,1)=13
c define incoming particles as massless, to ensure proper cancellation in
c matrix-element
      P(Kin1,4)=Dabs(P(1,3))
      P(Kin1,5)=0
      P(Kin2,4)=Dabs(P(2,3))
      P(Kin2,5)=0
      Do I=1,3
         P(Kin1,I)=P(1,I)
         P(Kin2,I)=P(2,I)
      enddo
c
      if(ltest) then
         SH_NEW = DOT(DSH,DSH)
         goto 10
      endif
      P(KIA1,3) = 0.5d0*(P(KIA1,3)-P(KIA1,4))
      P(KIA1,4) = DABS(P(KIA1,3))
      P(KIA2,3) = 0.5d0*(P(KIA2,3)+P(KIA2,4))
      P(KIA2,4) = DABS(P(KIA2,3))

      IF(ICCFM.EQ.0) Then
	   If(ires(1).eq.1) then
            P(KIA1,1) = 0D0
            P(KIA1,2) = 0D0
         endif
	   If(ires(2).eq.1) then
            P(KIA2,1) = 0D0
            P(KIA2,2) = 0D0
         endif
      ENDIF

      DO I=1,4
         DSHN(I)= P(KIA1,I) + P(KIA2,I)
      ENDDO
      SH_NEW = DOT(DSHN,DSHN)
c      IF(SH_NEW.LE.AM(1)**2) GOTO 20


   10 continue


c      K1K2 = DOT1(KIA1,KIA2)
c      K12 = DOT1(KIA1,KIA1)
c      K22 = DOT1(KIA2,KIA2)
      k12 = P(KIA1,1)**2 + P(KIA1,2)**2
      k22 = P(KIA2,1)**2 + P(KIA2,2)**2
	K1K2 = P(KIA1,1)*P(KIA2,1) + P(KIA1,2)*P(KIA2,2)
      M2 = AM(1)**2
      E1 = K12/M2
      E2 = K22/M2
	sig=0
c	write(6,*) k12,e1,k22,e2
	  
 	sig = 0.5


c	write(6,*) K12,k22,K(nia1,2),nia1,K(nia1,2),nia2
c	call pylist(1)
      QG1 = P(NIA1,1)**2+ P(NIA1,2)**2
      QG2 = P(NIA2,1)**2+ P(NIA2,2)**2
      SFSCAL1 = DABS(QG1)
      SFSCAL2 = DABS(QG2)
c      write(6,*) 'k12 ',k12,' k22 ',k22
c      write(6,*) 'qg1 ',qg1,' qg2 ',qg2
c new scale
      PT2Q = P(KF1,1)**2 + P(KF1,2)**2  
      xg_bar =  -qg2/(1.d0-xg)/SSS
      yph = dble(YY)
      yph_bar =  -qg1/(1.d0 - yph)/SSS

	if(sfscal1.le.sfcut) sfscal1=sfcut 
	if(sfscal2.le.sfcut) sfscal2=sfcut 
      xsi_hard(1) = (xg+yph_bar)/(yph+xg_bar)	
      xsi_hard(2) = (yph+xg_bar)/(xg+yph_bar)
 
c which of the quarks is closer to the proton:
      rapp = PYP(NPFIN,17)
      zz1 = 1. - (P(NF1,4)+ P(NF1,3))/2./P(2,4)/yph
      zz2 = 1. - (P(NF1,4)+ P(NF1,3))/2./P(2,4)/xg
      ptstm=dsqrt(dabs(qg1)+P(NF1,5)**2)

      If(IQBAR.eq.1) then
         if(xsi_hard(1).le.0) xsi_hard(1)=0
         if(xsi_hard(2).le.0) xsi_hard(2)=0
         stmax1 = yph*sqrt(xsi_hard(1)*SSS)
         stmax2 = xg*sqrt(xsi_hard(2)*SSS)
      elseif(IQBAR.EQ.2) then
         stmax1 = ptstm/(1.-zz1)
         stmax2 = ptstm/(1.-zz2)
      else
         write(6,*) ' meqg. no valid scale selected ',IQBAR
      endif
		
	if(stmax1.ne.stmax1) then
c	   write(6,*) ' meqg: stmax problem ',xg,xsi_hard,SSS
	   stmax1=1d20
	endif
      IF(FIRST) THEN
         IF(IQBAR.EQ.1) THEN
            write(6,*) ' meqg - scale: qq-pair '
         ELSE
            write(6,*) ' meqg - scale: q or q_bar '
         ENDIF
         FIRST=.FALSE.
      ENDIF
c            write(6,*) 'meqg xsi ',xsi_hard,' stmax ',stmax
c            write(6,*) ' new xg,xg_bar ',xg,xg_bar
c            write(6,*) ' new yph,yph_bar ',yph,yph_bar
      CALL PYROBO(0,N,0.D0,0.D0,-CM(1)/CM(4),-CM(2)/CM(4),-CM(3)/CM(4)
     +)

c new scale
c check on scale:
      SCALE = PT2Q+M2
      IF(IQ2.EQ.1) THEN
         SCALE = M2
      ELSEIF(IQ2.EQ.2) THEN
         SCALE = SH_NEW
      ELSEIF(IQ2.EQ.3) THEN
         SCALE = M2 + PT2Q
      ELSEIF(IQ2.EQ.5) THEN
         SCALE = Q2 + PT2Q
      ELSEIF(IQ2.EQ.6) THEN
         SCALE = DABS(K12)+DABS(K22)
c this is for testing
      ELSEIF(IQ2.EQ.7) THEN
         Scale = max(stmax1,stmax2)**2
      ELSE
         WRITE(6,*) ' NO VALID Q2 SCALE. STOP'
         STOP
      ENDIF
	
      Dscale = scale
c      write(6,*) ' in meqg, before alphas'
c      alpha_s = ALPHASCA(sqrt(dscale))
c	QQg = Qg0
c	Lambda = 0.2	  
c	alpha_s=12.*3.1415927/(33.-8)
c     &	   /Log(dmax1(Dscale,QQg**2)/Lambda**2)

c      write(6,*) ' new alphas ',alpha_s

c      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)
     
      IF(K(NIA1,2).eq.22.and.K(NIA2,2).NE.21) then
         Q2 = -Dot1(nia1,nia1)
         IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)
c         stmax2 = sqrt(q2)
         if(stmax2.ge.sqrt(sss)) then
            XPQ2(0) = 0.d0
         else
            call cauniglu(K(2,2),XG2SF,SFSCAL2,STMAX2,XPQ2)
         endif
c	   write(6,*) ' meqg: stmax,xpq1,xpq2 ',
c     &	   STMAX,XPQ1(0),XPQ2(0)
         xpqc2 = 0
         NFL=2
         do i=-NFL,NFL
         XPQC2= xpqc2 + XPQ2(I)*DFLOAT(PYCHGE(I))**2/9.D0
c         write(6,*) 'meqg: xpqc2 ',i,xpq2(i)
         end do
cccc         xpqc2 = 0.1
c select the flavour         
         QFT = - dcasrn()*XPQC2
         KPA=-NFL-1
   90    KPA=KPA+1
         QFT = QFT + DBLE(PYCHGE(KPA))**2/9.D0*DBLE(XPQ2(KPA))
         IF(QFT.LT.0.0D0) GOTO 90
         if(abs(kpa).gt.3) then
           write(6,*) ' fatal error in meqg: KPA = ',KPA
         endif
         K(NIA2,2) = KPA
         K(NF1,2) = KPA
         WTG= dble(XPQC2)
c         If(SFSCAL2.gt.STMAX2**2) WTG = 0 
c use photon flux in p_semih
         Sig = 4D0* ALPH_EM*PI*PI/Q2
c         write(6,*) ' meqg: yy,q2 ',yy,q2
c         write(6,*) 'meqg: ',sig,wtg,kpa
      else IF(K(NIA1,2).NE.21.AND.K(NIA2,2).EQ.21) THEN
         if(stmax1.ge.sqrt(sss)) then
            XPQ1(0) = 0.d0
         elseif(stmax2.ge.sqrt(sss)) then
            XPQ2(0) = 0.d0
         else
            call cauniglu(K(1,2),XG1SF,SFSCAL1,STMAX1,XPQ1)
            call cauniglu(K(2,2),XG2SF,SFSCAL2,STMAX2,XPQ2)
         endif
c	   write(6,*) ' meqg: stmax,xpq1,xpq2 ',
c     &	   STMAX,XPQ1(0),XPQ2(0)
         XPQC2=XPQ1(1)+XPQ1(2)+XPQ1(3)+XPQ1(-1)+XPQ1(-2)+XPQ1(-3)
         WTG= dble(XPQ2(0)*XPQC2)
      ELSE
	   write(6,*) ' no valid process for ipro = 4 ',
     &   K(NIA1,2),K(NIA2,2)
         WTG= 0.d0
      ENDIF
      
      
      
      SUM2 = Sig 
c	write(6,*) ' meqg :',sig,wtg
C new cut for on-shell DGLAP
	If(Iccfm.eq.0) then
	   If(ires(1).eq.1.and.(Qg1.gt.scale*scalfa)) SUM2 = 0 
	   If(ires(1).eq.1.and.(Qg1.gt.STMAX1**2)) SUM2 = 0 
	   If(ires(2).eq.1.and.(Qg2.gt.scale*scalfa)) SUM2 = 0 
	   If(ires(2).eq.1.and.(Qg2.gt.STMAX2**2)) SUM2 = 0 
	Endif

c      Write(6,*) ' meqg Wtg ',Wtg,XG1SF,XG2SF
      IF(SUM2.LT.0D0) THEN
c     write(6,*) ' SUM2 ',SUM2
         if(sum2.ne.sum2) then
            write(6,*) ' meqg SUM2,xsi_hard ',SUM2,xsi_hard
            write(6,*) ' meqg k12,k22 ',k12,k22
         endif
         SUM2 = 0D0
      ENDIF
 
c      write(6,*) ' SUM2 ',SUM2,wtg
      SUM2 = SUM2*WTG
c      write(6,*) ' meqg:  SUM2 ',SUM2,wtg

      WT1 = SUM2
      if(sum2.ne.sum2) then

         write(6,*) ' meqg xpq(0),alpha_s',xpq1(1),xpq2(0),alpha_s
         write(6,*) ' meqg SUM2,xsi_hard ',SUM2,xsi_hard
         write(6,*) ' meqg k12,k22 ',k12,k22
         write(6,*) ' P(KIA1) ',(P(KIA1,I),I=1,5)
         write(6,*) ' P(KIA2) ',(P(KIA2,I),I=1,5)
         write(6,*) ' P(KF1) ',(P(KF1,I),I=1,5)

      endif
      CALL PYEDIT(13)
      
      RETURN
   20 WT1 = 0.d0
      CALL PYEDIT(13)
      RETURN
      END
C From Leonard J. Moss of SLAC:

C Here's a hybrid QuickSort I wrote a number of years ago.  It's
C based on suggestions in Knuth, Volume 3, and performs much better
C than a pure QuickSort on short or partially ordered input arrays.  

      SUBROUTINE SORTRX(N,DATA,INDEX)
C===================================================================
C
C     SORTRX -- SORT, Real input, indeX output
C
C
C     Input:  N     INTEGER
C             DATA  REAL
C
C     Output: INDEX INTEGER (DIMENSION N)
C
C This routine performs an in-memory sort of the first N elements of
C array DATA, returning into array INDEX the indices of elements of
C DATA arranged in ascending order.  Thus,
C
C    DATA(INDEX(1)) will be the smallest number in array DATA;
C    DATA(INDEX(N)) will be the largest number in DATA.
C
C The original data is not physically rearranged.  The original order
C of equal input values is not necessarily preserved.
C
C===================================================================
C
C SORTRX uses a hybrid QuickSort algorithm, based on several
C suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
C "pivot key" [my term] for dividing each subsequence is chosen to be
C the median of the first, last, and middle values of the subsequence;
C and the QuickSort is cut off when a subsequence has 9 or fewer
C elements, and a straight insertion sort of the entire array is done
C at the end.  The result is comparable to a pure insertion sort for
C very short arrays, and very fast for very large arrays (of order 12
C micro-sec/element on the 3081K for arrays of 10K elements).  It is
C also not subject to the poor performance of the pure QuickSort on
C partially ordered data.
C
C Created:  15 Jul 1986  Len Moss
C
C===================================================================
 
      INTEGER   N,INDEX(N)
      REAL      DATA(N)
 
      INTEGER   LSTK(31),RSTK(31),ISTK
      INTEGER   L,R,I,J,P,INDEXP,INDEXT
      REAL      DATAP
 
C     QuickSort Cutoff
C
C     Quit QuickSort-ing when a subsequence contains M or fewer
C     elements and finish off at end with straight insertion sort.
C     According to Knuth, V.3, the optimum value of M is around 9.
 
      INTEGER   M
      PARAMETER (M=9)
 
C===================================================================
C
C     Make initial guess for INDEX
 
      DO 50 I=1,N
         INDEX(I)=I
   50    CONTINUE
 
C     If array is short, skip QuickSort and go directly to
C     the straight insertion sort.
 
      IF (N.LE.M) GOTO 900
 
C===================================================================
C
C     QuickSort
C
C     The "Qn:"s correspond roughly to steps in Algorithm Q,
C     Knuth, V.3, PP.116-117, modified to select the median
C     of the first, last, and middle elements as the "pivot
C     key" (in Knuth's notation, "K").  Also modified to leave
C     data in place and produce an INDEX array.  To simplify
C     comments, let DATA[I]=DATA(INDEX(I)).
 
C Q1: Initialize
      ISTK=0
      L=1
      R=N
 
  200 CONTINUE
 
C Q2: Sort the subsequence DATA[L]..DATA[R].
C
C     At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
C     r > R, and L <= m <= R.  (First time through, there is no
C     DATA for l < L or r > R.)
 
      I=L
      J=R
 
C Q2.5: Select pivot key
C
C     Let the pivot, P, be the midpoint of this subsequence,
C     P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
C     so the corresponding DATA values are in increasing order.
C     The pivot key, DATAP, is then DATA[P].
 
      P=(L+R)/2
      INDEXP=INDEX(P)
      DATAP=DATA(INDEXP)
 
      IF (DATA(INDEX(L)) .GT. DATAP) THEN
         INDEX(P)=INDEX(L)
         INDEX(L)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
      IF (DATAP .GT. DATA(INDEX(R))) THEN
         IF (DATA(INDEX(L)) .GT. DATA(INDEX(R))) THEN
            INDEX(P)=INDEX(L)
            INDEX(L)=INDEX(R)
         ELSE
            INDEX(P)=INDEX(R)
         ENDIF
         INDEX(R)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
C     Now we swap values between the right and left sides and/or
C     move DATAP until all smaller values are on the left and all
C     larger values are on the right.  Neither the left or right
C     side will be internally ordered yet; however, DATAP will be
C     in its final position.
 
  300 CONTINUE
 
C Q3: Search for datum on left >= DATAP
C
C     At this point, DATA[L] <= DATAP.  We can therefore start scanning
C     up from L, looking for a value >= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         I=I+1
         IF (DATA(INDEX(I)).LT.DATAP) GOTO 300
 
  400 CONTINUE
 
C Q4: Search for datum on right <= DATAP
C
C     At this point, DATA[R] >= DATAP.  We can therefore start scanning
C     down from R, looking for a value <= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         J=J-1
         IF (DATA(INDEX(J)).GT.DATAP) GOTO 400
 
C Q5: Have the two scans collided?
 
      IF (I.LT.J) THEN
 
C Q6: No, interchange DATA[I] <--> DATA[J] and continue
 
         INDEXT=INDEX(I)
         INDEX(I)=INDEX(J)
         INDEX(J)=INDEXT
         GOTO 300
      ELSE
 
C Q7: Yes, select next subsequence to sort
C
C     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
C     for all L <= l < I and J < r <= R.  If both subsequences are
C     more than M elements long, push the longer one on the stack and
C     go back to QuickSort the shorter; if only one is more than M
C     elements long, go back and QuickSort it; otherwise, pop a
C     subsequence off the stack and QuickSort it.
 
         IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=J+1
            RSTK(ISTK)=R
            R=I-1
         ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=L
            RSTK(ISTK)=I-1
            L=J+1
         ELSE IF (R-J .GT. M) THEN
            L=J+1
         ELSE IF (I-L .GT. M) THEN
            R=I-1
         ELSE
C Q8: Pop the stack, or terminate QuickSort if empty
            IF (ISTK.LT.1) GOTO 900
            L=LSTK(ISTK)
            R=RSTK(ISTK)
            ISTK=ISTK-1
         ENDIF
         GOTO 200
      ENDIF
 
  900 CONTINUE
 
C===================================================================
C
C Q9: Straight Insertion sort
 
      DO 950 I=2,N
         IF (DATA(INDEX(I-1)) .GT. DATA(INDEX(I))) THEN
            INDEXP=INDEX(I)
            DATAP=DATA(INDEXP)
            P=I-1
  920       CONTINUE
               INDEX(P+1) = INDEX(P)
               P=P-1
               IF (P.GT.0) THEN
                  IF (DATA(INDEX(P)).GT.DATAP) GOTO 920
               ENDIF
            INDEX(P+1) = INDEXP
         ENDIF
  950    CONTINUE
 
C===================================================================
C
C     All done
 
      END

	Subroutine Siggrid(xn,q2n,sf2,sfl,sf2c,sflc,sf2b,sflb)
      Implicit None
	Integer N1,N2
	Integer K1,K2,LUX
#include "caslujets.inc"
#include "casluco.inc"
#include "casparas.inc"
#include "casvalues.inc"
#include "caspart.inc"
#include "caseffic.inc"
#include "casparam.inc"
      Integer Iglu
      Common/CAGLUON/Iglu

      Integer I
      Integer Iseed
      Common/random/ISEED
	
	Integer Nevent
	Common/steer1/Nevent

	Logical First,Fccfm1,Fccfm2
	Common/ myfirst/First,Fccfm1,Fccfm2
      logical lcalc
	Common/ myfirst2/lcalc
	
      Double Precision q2test,xtest,x1min,x1max,x1minL,x1maxL,x1dif,q2ini
      Double Precision x1minlin,x1maxlin
      Integer NbxL,Nbxlin
      Integer j
      
      Double Precision xn,q2n,f2e,f2,fl,f2c,flc,f2b,flb
      Double Precision sf2,sfl,sf2c,sflc,sf2b,sflb
                 
      Logical firstgrid
      data firstgrid/.true./

      Real sq2,xmymin
      Integer nq2point
      Common /myq2grid/ nq2point,sq2(200),xmymin
#include "cascade-util_nbx.inc"
      INTEGER NA(1)
      double precision fgrid(Nbx+1),xa(1),a(Nbx+1)
      Common /mykernel/fgrid,xa,a,na
      Double precision x,q2
      Common /myvalues/x,q2
      
      Double Precision fgridxq2(Nbx+1,200),flgridxq2(Nbx+1,200)
      Double Precision fcgridxq2(Nbx+1,200),flcgridxq2(Nbx+1,200)
      Double Precision fbgridxq2(Nbx+1,200),flbgridxq2(Nbx+1,200)

      Integer IFLflag
      Common/FlagFL/ IFLflag
      Integer IcasHF
      Common/CasHF/IcasHF
      
      Double precision small 
      Parameter (small=1.d-3)
      
      Integer Iread

      INTEGER IN
c	Real A0,B0,EPS,SUM
c      Real ERR,ABSTOL
c      Real Tol
c      Real f2_int
	Double Precision A0,B0,EPS,SUM
      Double Precision ERR,ABSTOL
      Double Precision Tol
      Integer Nseg
      
      Double Precision f2_int
      Logical ex
      Integer Irr       
      External f2_int
      
      double precision au 
      Logical firstt
      common/f2fit/au(50),firstt
      
      Integer nq2test,Nbxtest
      Real q2grid(400)
      
      Integer Itheory
      Common/theory/Itheory
      
      Integer jfound
      
      Iread = 0
      
      
      x = xn
      q2 = q2n
c      write(6,*) ' siggrid ',firstgrid,x,q2,nbx,nq2point
      If(Firstgrid) then
         inquire(FILE='Cascade/updfGrids/f2tilde-grid.dat',EXIST=ex)
c create grid
         if(ex) then
            open(3,FILE='Cascade/updfGrids/f2tilde-grid.dat', FORM='formatted',
     +      STATUS='OLD',IOSTAT=IRR,ERR=220)
            write(6,*) ' Cascade/updfGrids/f2tilde-grid.dat existing '
            if(Iread.eq.0) write(6,*) ' overwriting existing file '
         else 
            write(6,*) '  creating f2tilde-grid.dat ' 
            open(3,FILE='Cascade/updfGrids/f2tilde-grid.dat', FORM='formatted',STATUS='NEW',
     +      IOSTAT=IRR,ERR=220)
         Endif
         
ccc         if(Itheory.eq.103) Iread = 1
         
         NA(1)=Nbx+1

	   x1min = 1e-5
	   x1min = xmymin
c         x1max =0.1d0
         x1max =0.9999d0
         x1minL=log(x1min)
         x1maxL=log(x1max)
c         NbxL = Nbx-10
         NbxL = Nbx
         x1dif=(x1maxL-x1minL)/Real(NbxL)         
         In=0
c         Do I=0,NbxL-1
         Do I=0,NbxL
           in=in+1
           A(IN)= exp(x1minL + x1dif*Real(I))
         enddo
c lin grid in x from 0.1 - 1
c         x1minlin=0.1d0
c         x1maxlin=1.0d0
c         Nbxlin = 10
c         x1dif=(x1maxlin-x1minlin)/Real(Nbxlin)
c         Do I=NbxL+1,NbxL+Nbxlin+1
c           in=in+1
c           A(in) = x1minlin + x1dif*Real(I-NbxL-1)
c         Enddo

c lin grid in x from xymin - 1
c         x1minlin=xmymin
c         x1maxlin=1.0d0
c         Nbxlin = Nbx
c         x1dif=(x1maxlin-x1minlin)/Real(Nbx)
c         in=0
c         Do I=0,Nbx
c           in=in+1
c           A(in) = x1minlin + x1dif*Real(I)
c         Enddo
         Nbxtest = Nbx
         nq2test = nq2point
         if(nq2test.eq.0) then 
            write(6,*) ' siggrid check: nq2test = 0 - STOP',nq2test
         endif
         write(6,*) ' siggrid check',Nbxtest,nq2test,iread
c         if(Iread.eq.1) then
c            read(3,*) Nbxtest,nq2test
c            write(6,*) ' siggrid check',Nbx,Nbxtest,nq2point,nq2test
c         endif
         do j = 1, nq2test
            q2test = sq2(j)
            q2grid(j) = sq2(j)
            Do I=0,Nbxtest 
c               xtest  = exp(x1minL + x1dif*Real(I))
               xtest  = A(I+1)
               if(iread.eq.0) then
               f2 = 0
               fl = 0
               f2c = 0
               flc = 0
               f2b = 0
               flb = 0
c do the light quarks here
               IcasHF=0
c calcualte SumAab
               IFLflag = 1
               call sigcalc(xtest,q2test,f2,f2e)
               fgridxq2(i+1,j)=f2
c calcualte SumAL
               IFLflag = 2
               call sigcalc(xtest,q2test,fl,f2e)
               flgridxq2(i+1,j)=fl
c do charm here
               IcasHF=4
c calcualte SumAab
               IFLflag = 1
               call sigcalc(xtest,q2test,f2c,f2e)
               fcgridxq2(i+1,j)=f2c
c calcualte SumAL
               IFLflag = 2
               call sigcalc(xtest,q2test,flc,f2e)
               flcgridxq2(i+1,j)=flc
c do bottom here
               IcasHF=5
c calcualte SumAab
               IFLflag = 1
               call sigcalc(xtest,q2test,f2b,f2e)
               fbgridxq2(i+1,j)=f2b
c calcualte SumAL
               IFLflag = 2
               call sigcalc(xtest,q2test,flb,f2e)
               flbgridxq2(i+1,j)=flb
               write(3,*) xtest,q2test,f2,fl,f2c,flc,f2b,flb
c               write(6,*) j,I,' q2test = ',q2test,' x = ',xtest,' f2 ',f2,' f2e ',f2e
               endif
            End do
         end do
         if(Iread.eq.1) then
               I=0
               J=0
               q2ini = 0.
330            read(3,*,END=331,ERR=331) xtest,q2test,f2,fl,f2c,flc,f2b,flb
               if(q2test.ne.q2ini) then
                  q2ini=q2test
                  J = J + 1
                  I=0
               endif
               I=I+1
               fgridxq2(i,j)=f2
               flgridxq2(i,j)=fl
               fcgridxq2(i,j)=f2c
               flcgridxq2(i,j)=flc
               fbgridxq2(i,j)=f2b
               flbgridxq2(i,j)=flb
               A(I) = xtest    
               q2grid(j) = q2test   
c               write(6,*) ' siggrid ',i,j,xtest,q2test, f2,fl,f2c,flc,f2b,flb     
               goto 330
 331           continue
         Nbxtest = I-1
         nq2test = j
         NA(1)=Nbx+1
         write(6,*) ' Nbxtest ',Nbxtest,nq2test,na(1),nbx
         endif
         close(3)
      endif
      firstgrid = .false.
      
      
      do i=0,Nbxtest
         fgrid(i+1) = 0.
      end do 
      jfound = 0
      do j =1,nq2test
      
      if(q2.ge.q2grid(j)-small.and.q2.le.q2grid(j)+small) then
c      write(6,*) ' q2 found : ',q2,small,j,q2grid(j)
      jfound = j 
cc      q2 point found
        do i=0,Nbxtest
          fgrid(i+1) = fgridxq2(i+1,j)
        end do 
      goto 345
      else 
c      write(6,*) ' q2 = ',q2,q2grid(j),small
      endif
      enddo        
345   continue
      If(jfound.eq.0) then
         write(6,*) ' siggrid: problem finding q2 point ',q2
      endif
c      integrate with starting distribution
	A0 = 0.0001
	A0 = 0.0000001
	B0 = 1.0
	eps = 0.001 
      abstol = 0  ! abstol = 0, value of abstol will be ignored
      nseg = 50
c      nseg = 1 ! fully automatic, adapt until tolerance attained, 
      CALL DADAPT(f2_int,A0,B0,NSEG,eps,abstol,SUM,ERR)
      sf2 = SUM
c now calculate Fl
      do i=0,Nbxtest
          fgrid(i+1) = flgridxq2(i+1,jfound)
      end do 
      abstol = 0  ! abstol = 0, value of abstol will be ignored
      CALL DADAPT(f2_int,A0,B0,NSEG,eps,abstol,SUM,ERR)      
      sfl = SUM
c now calculate F2c
      do i=0,Nbxtest
          fgrid(i+1) = fcgridxq2(i+1,jfound)
c          write(6,*) ' A(I) ',A(i),fgrid(i+1)
      end do 
      abstol = 0  ! abstol = 0, value of abstol will be ignored
      CALL DADAPT(f2_int,A0,B0,NSEG,eps,abstol,SUM,ERR)      
      sf2c = SUM
c now calculate Flc
      do i=0,Nbxtest
          fgrid(i+1) = flcgridxq2(i+1,jfound)
      end do 
      abstol = 0  ! abstol = 0, value of abstol will be ignored
      CALL DADAPT(f2_int,A0,B0,NSEG,eps,abstol,SUM,ERR)      
      sflc = SUM
c now calculate F2b
      do i=0,Nbxtest
          fgrid(i+1) = fbgridxq2(i+1,jfound)
      end do 
      abstol = 0  ! abstol = 0, value of abstol will be ignored
      CALL DADAPT(f2_int,A0,B0,NSEG,eps,abstol,SUM,ERR)      
      sf2b = SUM
c now calculate Flb
      do i=0,Nbxtest
          fgrid(i+1) = flbgridxq2(i+1,jfound)
      end do 
      abstol = 0  ! abstol = 0, value of abstol will be ignored
      CALL DADAPT(f2_int,A0,B0,NSEG,eps,abstol,SUM,ERR)      
      sflb = SUM
    
      
      tol = 0.0
      f2e = sqrt(ERR**2 + (tol*sum)**2) 
c      write(6,*) ' siggrid: ',xn,q2n,f2,fl
      return  
220   write(6,*) ' f2tilde_grid.dat  on unit 3 cannot be opened '          
c	stop
      END
      
      Function f2_int(xt)
	Implicit None
c	Real f2_int,x0,z
c	Real xmin,xmax,xt,wt
	Double Precision f2_int,x0,z
	Double Precision xmin,xmax,xt,wt
      Double precision test
      Double Precision f2tilde
      Real sq2,xmymin
      Integer nq2point
      Common /myq2grid/ nq2point,sq2(200),xmymin
#include "cascade-util_nbx.inc"
      INTEGER NA(1)
      double precision fgrid(Nbx+1),xa(1),a(Nbx+1)
      Common /mykernel/fgrid,xa,a,na
      Double precision x,q2
      Common /myvalues/x,q2
      Double Precision DHFINT
      Double Precision DHPOLINT

	Double Precision A1,A2,A3,A4,A5,A6
      logical first
	double precision au 
      common/f2fit/au(50),first
      real xsat
	Integer jj

	Integer Noevol
      Common/Evoltest/Noevol        
	
	xmin=x
	xmax=0.99999
	wt = log(xmax/xmin)
	x0= xmin*(xmax/xmin)**xt
      
	if(x0.lt.x) then
	   f2_int = 0
	   return
	   endif
	XA(1) = x/x0
      if(xa(1).gt.1. ) then
	   f2_int = 0
	   return
	   endif
        
c      write(6,*) 'f2_int',(AU(jj),jj=1,4)
      f2tilde = DHFINT(1,XA,NA,A,fgrid)    


	A1 = au(1)
	A2 = au(2)
      A3 = au(3)
	A4 = au(4)
      A5 = au(5)
      A6 = au(6)
c       write(6,*) a1,a2,a3,a4
	 test = 0.
       if(x0.lt.0.) then
	   test = 0.
	   else
c         test =A1*5.*(1./x0**A2)*(1.-x0)**A3*(1.-A4*x0)/x0
c updf style
         test =A1*5.*(1./x0**A2)*(1.-x0)**A3*(1.-A4*x0+A5*sqrt(x0)+A6*x0**2)/x0
c herapdf style 
c         test =A1*5.*(1./x0**A2)*(1.-x0)**A3*(1.+ A4*x0 + A5*x0**2 + A6*x0**3)/x0
c cteq style 
c        test =A1*5.*(1./x0**A2)*(1.-x0)**A3*exp(x0*A4)*(1.+exp(A5*x0) + exp(A6*x0**2))/x0
c	  write(6,*) 'xa_int ',f2tilde,test,x 
        if(test.le.0) test=0
         test=test*f2tilde
         test = max(0.,test)	   
        endif 
c	write(6,*) 'xa_int ',xpqs(0),test,x 
c     f(x)/g(x)*int(g(x)dx) = f(x)*x*wt
	f2_int = test * wt * x0	
	Return 
	
	end
      
	Subroutine Siggrid1(xn,q2n,f2,f2e)
      Implicit None
	Integer N1,N2
	Integer K1,K2,LUX
#include "caslujets.inc"
#include "casluco.inc"
#include "casparas.inc"
#include "casvalues.inc"
#include "caspart.inc"
#include "caseffic.inc"
#include "casparam.inc"
      Integer Iglu
      Common/CAGLUON/Iglu

      Integer I
      Integer Iseed
      Common/random/ISEED
	
	Integer Nevent
	Common/steer1/Nevent

	Logical First,Fccfm1,Fccfm2
	Common/ myfirst/First,Fccfm1,Fccfm2
      logical lcalc
	Common/ myfirst2/lcalc
	
      Double Precision q2test,xtest,x1min,x1max,x1minL,x1maxL,x1dif
      Integer j
      
      Double Precision xn,q2n,f2,f2e
                  
      Logical firstgrid
      data firstgrid/.true./

      Real sq2,xmymin
      Integer nq2point
      Common /myq2grid/ nq2point,sq2(200),xmymin
#include "cascade-util_nbx.inc"
      INTEGER NA(1)
      double precision fgrid(Nbx+1),xa(1),a(Nbx+1)
      Common /mykernel/fgrid,xa,a,na
      Double precision x,q2
      Common /myvalues/x,q2
      
      Double Precision fgridxq2(Nbx+1,200)
      
      Double precision small 
      Parameter (small=1.d-3)
      
      Integer Iread
      Parameter (Iread=0)

      INTEGER IN
	Real A0,B0,EPS,SUM
      Real ERR,ABSTOL
      Real Tol
      Integer Nseg
      
      Real f2_int
      Logical ex
      Integer Irr       
      External f2_int

      x = xn
      q2 = q2n
c      write(6,*) ' siggrid ',firstgrid
      If(Firstgrid) then
            q2test = q2
               xtest  = x
               write(6,*) ' q2test = ',q2test,' x = ',xtest
               call sigcalc(xtest,q2test,f2,f2e)
               write(6,*) ' q2test = ',q2test,' x = ',xtest,
     &           ' f2 ',f2,' f2e ',f2e
      endif
      return  
220   write(6,*) ' f2tilde_grid.dat  on unit 3 cannot be opened '          
c	stop
      END
      
	function dhpolint(NARG,ARG,NENT,ENT,TABLE)
      Implicit None
C
C
      Integer NARG
      INTEGER   NENT(500)
      Double Precision     ARG(500),   ENT(500)
      Double Precision TABLE(NENT(1))
      Double Precision DHPOLINT
      Integer N,NDIM,LMAX,LMIN
      Double Precision X
      Integer N1,N2,N3
      Integer jj(narg)
      Parameter (N1=3, N2=4, N3=4)
c      Double precision x1a(N1),x2a(N2),x3a(N3),ya(N1,N2),yb(N1,N2,N3)
      Double precision x1a(N1),ya(N1),yb(N1)
      Integer i,j,k,i1,j1,k1,nn
      Integer istart,iend
      Double Precision y,dy
      
     
      dhpolint = 0.
      LMAX = 0
      
      DO N  =  1, NARG
         X     =  ARG(N)
         NDIM  = NENT(N)
         LMIN  = LMAX + 1
         LMAX  = LMAX + NDIM
         call locate1(ENT(LMIN),NDIM,X,jj(N))
c         write(6,*) ' after locate ', x,jj(n)
c         stop
         if(n.eq.1) then 
           nn = n1
           else
           write(6,*) ' N > 1: stop '
           endif 
         i1=0
         if(nn.le.2) then
           istart = jj(n)
           iend   = jj(n)+nn-1
         elseif(nn.eq.3) then
           istart = jj(n)-nn+2
           iend   = jj(n)+nn-2
         elseif(nn.eq.4) then
           istart = jj(n)-nn+3
           iend   = jj(n)+nn-2
         elseif(nn.eq.5) then
           istart = jj(n)-nn+3
           iend   = jj(n)+nn-3
         elseif(nn.eq.6) then
           istart = jj(n)-nn+4
           iend   = jj(n)+nn-3
         elseif(nn.eq.7) then
           istart = jj(n)-nn+4
           iend   = jj(n)+nn-4
         elseif(nn.eq.8) then
           istart = jj(n)-nn+5
           iend   = jj(n)+nn-4
         elseif(nn.eq.9) then
           istart = jj(n)-nn+5
           iend   = jj(n)+nn-5
         endif  
         
      End Do
         i1 = 0
         do i=istart,iend
           i1=i1+1
c           write(6,*) 'test1 ', i1,n1,i,istart,iend
           if(i1.gt.n1) write(6,*) 'test1 ', i1,n1,i,istart,iend
c           write(6,*) ' index i,i1,j = ',i,i1,j
           x1a(i1) = ent(i)
           if(i.gt.51.or.i.le.0) then
c            write(6,*) X,ent(i),i,istart,iend
c            x1a(i1) = 0
           endif
         end do 
      i1=0
      do i=istart,iend
        i1=i1+1
           if(i1.gt.n1) write(6,*) ' test2 ', i1,n1
            ya(i1) = table(i)
      end do
c      write(6,*) ' x1a ',(x1a(i),I=1,n1)
c      write(6,*) ' ya ',(ya(i),I=1,n1)
      call polint(x1a,ya,n1,arg(1),y,dy)
ccc      call polin2(x1a,x2a,ya,n1,n2,arg(1),arg(2),y,dy)
ccc      call polin3(x1a,x2a,x3a,yb,n1,n2,n3,arg(1),arg(2),arg(3),y,dy)
c      write(6,*) ' result: y = ',y,' dy = ',dy
      dhpolint = y
c      dhpolint = ya(1,1)
      RETURN
1000  FORMAT('  HFUNCTION DPOLINT ... NARG =',I6,
     +       '  NOT WITHIN RANGE')
      END
      SUBROUTINE locate1(xx,n,x,j)

c     From Numerical recipes:
c     -----------------------
c     Given an array XX of length N, and a given value of X, returns a
c     value of J such that X is between XX(J) and XX(J+1).  XX must be
c     monotonic, either increasing or decreasing. J=0 or J=N is
c     returned to indicate that X is out of range. 
      Implicit None     
      Integer j,n
      Double Precision x, xx(n)
      Integer jl,jm,ju
      jl=0
      ju=n+1
 10   IF(ju-jl.gt.1)THEN
          jm=(ju+jl)/2
          IF((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))THEN
              jl=jm
          ELSE
              ju=jm
          ENDIF
          go to 10
      ENDIF
      if(x.eq.xx(1)) then
        j=1
      else if(x.eq.xx(n)) then
        j=n-1
      else
        j=jl
      endif
      RETURN
      END
      SUBROUTINE POLINT (XA,YA,N,X,Y,DY)
C     Adapted from "Numerical Recipes" 
      IMPLICIT NONE
      DOUBLE PRECISION XA(N),YA(N),X,Y,DY
      Integer NMAX
      Parameter (NMAX = 10)
      DOUBLE PRECISION C(NMAX),DIF,DIFT,HO,HP,W,D(NMAX),DEN
      INTEGER NS,I,M,N

      NS=1
      DIF=DABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
cz          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
      
*CMZ :          04/10/2001  22.14.38  by  Hannes Jung
*CMZ :  5.10/00 05/08/2000  08.56.38
*-- Author :
***********************************************************************
*    =======================================                          *
      SUBROUTINE BSPRNT( LU, ID, IP1, IP2 )
*    =======================================                          *
* ((purpose))                                                         *
*     Print out routine of BASES.                                     *
*  (Argument)                                                         *
*     ID  : Identity number of printouts.                             *
*     IP1... IP2 : Integer                                            *
*  (Author)                                                           *
*     S. Kawabata   May 1992                                          *
*     Last update   March 1994                                        *
***********************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MXDIM = 50, NDMX = 50, LENG = 32768)
      COMMON /BASE0/ JFLAG,IBASES
      COMMON /BASE1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,
     .               IG(MXDIM),NCALL
      COMMON /BASE2/ ACC1,ACC2,ITMX1,ITMX2
      COMMON /BASE3/ SCALLS,WGT,TI,TSI,TACC,IT
      COMMON /BASE4/ XI(NDMX,MXDIM),DX(MXDIM),DXD(LENG),DXP(LENG),
     .               ND,NG,NPG,MA(MXDIM)
      PARAMETER (ITM = 50)
      REAL*4 TIME, EFF, WRONG, TRSLT, TSTD, PCNT
      COMMON /BASE5/ ITRAT(ITM,0:1),TIME(ITM,0:2),EFF(ITM,0:1),
     .               WRONG(ITM,0:1),RESLT(ITM,0:1),ACSTD(ITM,0:1),
     .               TRSLT(ITM,0:1),TSTD(ITM,0:1),PCNT(ITM,0:1)
      REAL*4 STIME
      COMMON /BSRSLT/AVGI,SD,CHI2A,STIME,IT1,ITF
      CHARACTER*51 ICH(0:1)
      CHARACTER*1 CN
*        INTV = ( 0 / 1 / any ) = ( Batch / Batch(Unix) / Interactive )
*        IPNT = ( 0 / any ) = ( IBM Type / Ascii printer )
      COMMON /BSCNTL/ INTV, IPNT, NLOOP, MLOOP
*
       COMMON /BDATE/ IDATE(3),ITIME(2)
*            IDATE(1) : year        ITIME(1) : hour
*            IDATE(2) : month       ITIME(2) : minute
*            IDATE(3) : day
      REAL*4 TIMEBS,TIMINT,TIMESP,TIME0,RTIME,TIMEB1,TIMEB2,TIMES1
      COMMON /BTIME1/ TIME0,RTIME,TIMEB1,TIMEB2,TIMES1
      COMMON /BTIME2/ TIMEBS(0:2),TIMINT,TIMESP(0:2)
      REAL*4 XTIME
*
       COMMON/NINFO/ NODEID, NUMNOD
*
      DATA  ICH / 'Convergency Behavior for the Grid Optimization Step',
     .            'Convergency Behavior for the Integration Step      '/

      IF( NODEID .NE. 0 ) RETURN
      CN = CHAR(12)

      GO TO ( 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 ), ID
C----------------------------------------------------------- BSMAIN

  100 IF( IPNT .EQ. 0 ) THEN
          WRITE(LU,9600)
 9600     FORMAT(/1H1,/1H )
      ELSE
          WRITE(LU,9610) CN
 9610     FORMAT(A1)
      ENDIF
      WRITE(LU,9620) (IDATE(I),I=1,3),(ITIME(J),J=1,2)
 9620 FORMAT(53X,'Date: ',I4,'/',I2,'/',I2,2X,I2.2,':',I2.2)
      WRITE(LU,9050)
 9050 FORMAT(
     . 8X,'**********************************************************',
     ./8X,'*                                                        *',
     ./8X,'*     BBBBBBB     AAAA     SSSSSS   EEEEEE   SSSSSS      *',
     ./8X,'*     BB    BB   AA  AA   SS    SS  EE      SS    SS     *',
     ./8X,'*     BB    BB  AA    AA  SS        EE      SS           *',
     ./8X,'*     BBBBBBB   AAAAAAAA   SSSSSS   EEEEEE   SSSSSS      *',
     ./8X,'*     BB    BB  AA    AA        SS  EE            SS     *',
     ./8X,'*     BB    BB  AA    AA  SS    SS  EE      SS    SS     *',
     ./8X,'*     BBBB BB   AA    AA   SSSSSS   EEEEEE   SSSSSS      *',
     ./8X,'*                                                        *',
     ./8X,'*                   BASES Version 5.1                    *',
     ./8X,'*           coded by S.Kawabata KEK, March 1994          *',
     ./8X,'**********************************************************')

          RETURN
C----------------------------------------------------------- BSMAIN

  200     IF( IPNT .EQ. 0 ) THEN
              WRITE(LU,9600)
          ELSE
              WRITE(LU,9610) CN
          ENDIF
          WRITE(LU,9300)
 9300     FORMAT(20X,
     .         '****** END OF BASES *********')

C----------------------------------------------------------- BSMAIN

  300 CONTINUE
      WRITE(LU,9305)
 9305 FORMAT(
     .//5X,'<<   Computing Time Information   >>')

*     WRITE(LU,9310) (IDATE(I),I=1,3),(ITIME(J),J=1,2)
*9310 FORMAT(/15X,'Start at: ',I2,'/',I2,'/',I2,2X,I2.2,':',I2.2)
*     CALL BSDATE
*     WRITE(LU,9320) (IDATE(I),I=1,3),(ITIME(J),J=1,2)
*9320 FORMAT(15X,'End   at: ',I2,'/',I2,'/',I2,2X,I2.2,':',I2.2)
      WRITE(LU,9330)
 9330 FORMAT(/15X,'(1) For BASES              H: M:  Sec')
      CALL BSTCNV(TIMEBS(2),IH,MN,IS1,IS2)
      WRITE(LU,9340) IH, MN, IS1, IS2
 9340 FORMAT(19X,'Overhead           : ',I3,':',I2,':',I2,'.',I2.2)
      CALL BSTCNV(TIMEBS(0),IH,MN,IS1,IS2)
      WRITE(LU,9350) IH, MN, IS1, IS2
 9350 FORMAT(19X,'Grid Optim. Step   : ',I3,':',I2,':',I2,'.',I2.2)
      CALL BSTCNV(TIMEBS(1),IH,MN,IS1,IS2)
      WRITE(LU,9360) IH, MN, IS1, IS2
 9360 FORMAT(19X,'Integration Step   : ',I3,':',I2,':',I2,'.',I2.2)
      XTIME = TIMEB2 - TIMEB1
      CALL BSTCNV(XTIME,IH,MN,IS1,IS2)
      WRITE(LU,9365) IH, MN, IS1, IS2
 9365 FORMAT(19X,'Go time for all    : ',I3,':',I2,':',I2,'.',I2.2)
      EXTIM  = TIMEBS(1)*1000.0/SCALLS/0.7
      WRITE(LU,9375)
 9375 FORMAT(/15X,'(2) Expected event generation time')
      WRITE(LU,9376) EXTIM
 9376 FORMAT(19X,'Expected time for 1000 events :',F10.2,' Sec')
      RETURN

C----------------------------------------------------------- BASES

  400 NSP   = NG**NWILD
      MCALL = NSP*NPG
      WRITE(LU,9400) NDIM,NWILD,MCALL,NCALL,ND,NG,NSP
 9400 FORMAT(
     .//5X,'<<   Parameters for BASES    >>',
     .//5X,' (1) Dimensions of integration etc.',
     . /5X,'     # of dimensions :    Ndim    =',I9,3X,'( 50 at max.)',
     . /5X,'     # of Wilds      :    Nwild   =',I9,3X,'( 15 at max.)',
     . /5X,'     # of sample points : Ncall   =',I9,'(real)',
     .                                         I9,'(given)',
     . /5X,'     # of subregions    : Ng      =',I9,' / variable',
     . /5X,'     # of regions       : Nregion =',I9,' / variable',
     . /5X,'     # of Hypercubes    : Ncube   =',I9,
     .//5X,' (2) About the integration variables')
      WRITE(LU,9405)
 9405 FORMAT(10X,'------',2('+---------------'),'+-------+-------')
      WRITE(LU,9410)
 9410 FORMAT(10X,'    i       XL(i)           XU(i)     ',
     .           '  IG(i)   Wild')
      WRITE(LU,9405)
       DO 450 I = 1,NDIM
          IF( I .LE. NWILD ) THEN
          WRITE(LU,9420) I,XL(I),XU(I),IG(I)
 9420     FORMAT(10X,I5,1P,2('  ',E14.6),'  ',3X,0P,I1,3X,
     .                       '   yes')
          ELSE
          WRITE(LU,9421) I,XL(I),XU(I),IG(I)
 9421     FORMAT(10X,I5,1P,2('  ',E14.6),'  ',3X,0P,I1,3X,
     .                        '    no')
          ENDIF
  450  CONTINUE
       WRITE(LU,9405)
       WRITE(LU,9450) ITMX1,ACC1,ITMX2,ACC2
 9450  FORMAT(
     . /5X,' (3) Parameters for the grid optimization step',
     . /5X,'     Max.# of iterations: ITMX1 =',I9,
     . /5X,'     Expected accuracy  : Acc1  =',F9.4,' %',
     .//5X,' (4) Parameters for the integration step',
     . /5X,'     Max.# of iterations: ITMX2 =',I9,
     . /5X,'     Expected accuracy  : Acc2  =',F9.4,' %')

          RETURN
C----------------------------------------------------------- BASES

  500    IF( INTV .LE. 1 )    RETURN
         ISTEP  = IP1
         IF( IPNT .EQ. 0 ) THEN
             WRITE(LU,9600)
         ELSE
             WRITE(LU,9610) CN
         ENDIF
         WRITE(LU,9620) (IDATE(I),I=1,3),(ITIME(J),J=1,2)
         WRITE(LU,9500) ICH(ISTEP)
 9500    FORMAT(15X,A)
         WRITE(LU,9570)
         WRITE(LU,9550)
 9550    FORMAT(1X,'<- Result of  each iteration ->',
     .          2X,'<-     Cumulative Result     ->',
     .          1X,'< CPU  time >',
     .         /1X,' IT Eff R_Neg   Estimate  Acc %',
     .          2X,'Estimate(+- Error )order  Acc %',
     .          1X,'( H: M: Sec )')
         WRITE(LU,9570)
 9570    FORMAT(1X,7('----------'),'--------')
         RETURN

C----------------------------------------------------------- BASES

  600    IF( INTV .LE. 1 ) RETURN
         ISTEP  = IP1
         ITX = MOD( IT, ITM)
         IF( ITX .EQ. 0 ) ITX = ITM

         CALL BSLIST( LU, ITX, ISTEP )

         RETURN

  700    IF( INTV .LE. 1 ) RETURN
         WRITE(LU,9570)

         RETURN
C----------------------------------------------------------- BASES

  800    ITJ    = IP1
         ISTEP  = IP2
         ITX  = MOD( ITJ, ITM )
         IF( ITX .EQ. 0 ) ITX = ITM

         IF( ITRAT(1,ISTEP) .EQ. 1 ) THEN
             NDEV   = 1
         ELSE
             NDEV   = 2
             ITFN   = ITM
             ITMN   = 10000
             DO 610 I = 1,ITM
                IF( ITRAT(I,ISTEP) .LT. ITMN ) THEN
                    ITST = I
                    ITMN = ITRAT(I,ISTEP)
                ENDIF
  610        CONTINUE
             IF( ITST .EQ. 1 ) NDEV = 1
         ENDIF

         IF( IPNT .EQ. 0 ) THEN
             WRITE(LU,9600)
         ELSE
             WRITE(LU,9610) CN
         ENDIF
         WRITE(LU,9620) (IDATE(I),I=1,3),(ITIME(J),J=1,2)
         WRITE(LU,9500) ICH(ISTEP)
         WRITE(LU,9570)
         WRITE(LU,9550)
         WRITE(LU,9570)

  625    IF( NDEV .EQ. 1 ) THEN
             ITST = 1
             ITFN = ITX
         ENDIF

         DO 650 I = ITST, ITFN

            CALL BSLIST( LU, I, ISTEP )

  650    CONTINUE
         NDEV  = NDEV - 1
         IF( NDEV .GT. 0 ) GO TO 625
         WRITE(LU,9570)

      RETURN

C----------------------------------------------------------- BASES

  900 WRITE(LU,9950)
 9950 FORMAT(1X,'******** FATAL ERROR IN BASES **************',
     .      /1X,'There are no enough good points in this iteration.',
     .      /1X,'Process was terminated due to this error.')

      RETURN

C-----------------------------------------------------------------
 1000 LOOP = IP1
      IF( IP2 .NE. 0 ) THEN
          IF( IPNT .EQ. 0 ) THEN
              WRITE(LU,9600)
           ELSE
              WRITE(LU,9610) CN
           ENDIF
           WRITE(LU,9620) (IDATE(I),I=1,3),(ITIME(J),J=1,2)
           WRITE(LU,9650)
 9650      FORMAT(
     .      20X,'Results of Integration',
     .     /10X,5('----------'),'------',
     .     /10X,' Loop#  Estimate(+- Error )order',
     .                     '  It1  It2 ( H: M: Sec )',
     .     /10X,5('----------'),'------')
      ENDIF

      RE  = AVGI
      AC  = ABS(SD)
      ARE = ABS(RE)
      IF( ARE .GE. AC) THEN
          CALL BSORDR( ARE, F2, ORDER, IORDR)
      ELSE
          CALL BSORDR(  AC, F2, ORDER, IORDR )
      ENDIF
      RE  = RE/ORDER
      AC  = AC/ORDER
      CALL BSTCNV( STIME, IH, MN, IS1, IS2)
      WRITE(LU,9660) LOOP,RE,AC,IORDR,IT1,IT,IH,MN,IS1,IS2
 9660 FORMAT(10X,I6,F10.6,'(+-',F8.6,')E',I3.2,2I5,
     .        1X,I3,':',I2,':',I2,'.',I2.2,
     .      /10X,5('----------'),'------')

      RETURN
      END
***********************************************************************
*    ===================================                              *
      SUBROUTINE BSLIST( LU, I, ISTEP )
*    ===================================                              *
* ((purpose))                                                         *
*     Print out results of each iteration and cumulative result       *
* ((Argument))                                                        *
*  (Input)                                                            *
*     LU      : Logical unit number for the printer                   *
*     I       : Address in the arrays of common /BASE5/               *
*     ISTEP   : The Set-Identifier                                    *
*               ISTEP = ( 0 / 1 ) = ( Grid opt. / Integration step )  *
*                                                                     *
*     S. Kawabata   March '94                                         *
***********************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ITM = 50)
      REAL*4 TIME, EFF, WRONG, TRSLT, TSTD, PCNT
      COMMON /BASE5/ ITRAT(ITM,0:1),TIME(ITM,0:2),EFF(ITM,0:1),
     .               WRONG(ITM,0:1),RESLT(ITM,0:1),ACSTD(ITM,0:1),
     .               TRSLT(ITM,0:1),TSTD(ITM,0:1),PCNT(ITM,0:1)

      CALL BSTCNV( TIME(I,ISTEP), IH, MN, IS1, IS2 )

      RE  = RESLT(I,ISTEP)
      AC  = ABS(ACSTD(I,ISTEP))
      ARE = ABS(RE)
      IF( ARE .GE. AC) THEN
          CALL BSORDR( ARE, F2, ORDER, IORDR)
      ELSE
          CALL BSORDR(  AC, F2, ORDER, IORDR )
      ENDIF
      RE  = RE/ORDER
      AC  = AC/ORDER
      IEFF = EFF(I,ISTEP)
      WRITE(LU,9631) ITRAT(I,ISTEP),IEFF,WRONG(I,ISTEP),
     .              TRSLT(I,ISTEP),TSTD(I,ISTEP),
     .              RE,AC,IORDR,PCNT(I,ISTEP),IH,MN,IS1,IS2
 9631 FORMAT(I4,I4,F6.2,1P,E11.3, 0P,1X,F6.3,
     .              F10.6,'(+-',F8.6,')E',I3.2,1X,F6.3,
     .          1X,I3,':',I2,':',I2,'.',I2.2)


      RETURN
      END
 
