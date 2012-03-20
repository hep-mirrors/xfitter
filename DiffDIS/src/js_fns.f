c function FF integrated by DAIND1 to give the structure function
c=================================================
      real*8 function ffm(z)
      implicit none
      real*8 m2, mq2, z
      real*8 xbj, scale, pi, alphas, nlf, eh2
      real*8 fborn, fgcorr, fqcorr, xz
      integer nborn, ngcorr, nqcorr, ifch, ncount
      parameter (pi = 3.1415926535897932384626433832795D0)
      common/mass/m2, mq2, xbj
      common/include/nborn, ngcorr, nqcorr, ifch
      common/coupling/scale, alphas, nlf, eh2
      ffm = 0.d0
c here we get alpha_s and the parton density values from FNEW
      call fnew(xbj,z)
      xz = xbj/z
      if (nborn .ne. 0) ffm = ffm + eh2*xz*fborn(z)/z
      if (ngcorr .ne. 0) ffm = ffm + xz*fgcorr(z)/z
      if (nqcorr .ne. 0) ffm = ffm + xz*fqcorr(z)/z
      return
      end





c=================================================
      DOUBLE PRECISION FUNCTION DAIND1(A,B,FUN,EPS,KEY,MAX,KOUNT,EST)
C     ---------------------------------------------------------------
C  INPUTPARAMETERS
C  A,B      LIMITS OF THE INTEGRATION INTERVAL
C  FUN      FUNCTION TO BE INTEGRATED (TO BE DECLARED EXTERNAL IN THE MAIN PR.)
C  EPS      ABSOLUTE OR RELATIVE TOLERANCE,DEPENDING OF THE VALUE OF 'KEY'
C  KEY      =1 THEN 'EPS' DENOTES AN ABSOLUTE, =2 THEN A RELATIVE TOLERANCE
C  MAX      UPPER BOUND ON THE NUMBERS OF INTEGRAND EVALUATIONS (MAX.LE.10000)
C
C  OUTPUTPARAMETERS
C  KOUNT    NUMBER OF INTEGRAND EVALUATIONS
C  EST      ESTIMATION OF THE ABSOLUTE ERROR OF THE APPROXIMATION
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION MAXIM,MINIM,MODUL1,MODUL2
      INTEGER RANG(130)
      DIMENSION
     &AINIT(250),END(250),EPSIL(250),PART(250),W1(5),W2(5),W3(6),
     &     X1(5),X2(5)
      DATA X1/0.973906528517D+0,0.865063366689D+0,0.679409568299D+0,
     *         0.433395394129D+0,0.148874338981D+0/
      DATA X2/0.995657163026D+0,0.930157491356D+0,0.780817726586D+0,
     *         0.562757134669D+0,0.294392862701D+0/
      DATA W1/0.666713443087D-1,0.149451349151D+0,0.219086362516D+0,
     *         0.269266719310D+0,0.295524224715D+0/
      DATA W2/0.325581623080D-1,0.750396748109D-1,0.109387158802D+0,
     *         0.134709217311D+0,0.147739104901D+0/
      DATA W3/0.116946388674D-1,0.547558965744D-1,0.931254545837D-1,
     *         0.123491976262D+0,0.142775938577D+0,0.149445554003D+0/
      DATA TOL/0.23D-15/
      MAX1 = (MAX+21)/42+1
      MAX2 = MAX1/2+2
      ALFA = A
      BETA = B
      MAAT = 1
C EVALUATION OF GAUSSIAN AND KRONROD FORMULAS
   10 S = 0.5D+0*(BETA-ALFA)
      U = 0.5D+0*(BETA+ALFA)
      RES1 = 0.0D+0
      RES2 = W3(6)*FUN(U)
      DO 20 K = 1,5
        C = S*X1(K)
        C = FUN(C+U)+FUN(U-C)
        RES1 = RES1+W1(K)*C
        RES2 = RES2+W2(K)*C
        C = S*X2(K)
   20   RES2 =RES2+W3(K)*(FUN(C+U)+FUN(U-C))
      PAT = RES2*S
      MODUL2 = ABS(PAT-RES1*S)
      IF(MAAT.GT.1) GOTO 50
      EST = MODUL2
      BINT = PAT
      KOUNT =21
      PART(1) = BINT
      GOTO 90
   30 RANG(1) = 1
      AINIT(1) = A
      END(1) = B
      EPSIL(1) = EST
   40 NR = RANG(1)
      BINT = BINT-PART(NR)
      EST =EST-EPSIL(NR)
C THE SUBINTERVAL WITH LARGEST ERROR IS SPLIT UP INTO TWO EQUAL PARTS
      ALFA = AINIT(NR)
      BETA = (AINIT(NR)+END(NR))*0.5
      JJ = 1
      MAAT = MAAT+1
      GOTO 10
   50 EST = EST+MODUL2
      BINT = BINT+PAT
      IF(JJ.EQ.0) GOTO 60
      MODUL1 = MODUL2
      PAT1 = PAT
      ALFA = BETA
      BETA = END(NR)
      JJ = 0
      GOTO 10
   60 MA = MAAT
      IF(MAAT.GT.MAX2) MA = MAX1+3-MAAT
      IF(MODUL1.GT.MODUL2) GOTO 70
      EPSIL(NR) = MODUL2
      EPSIL(MAAT) = MODUL1
      AINIT(MAAT) = AINIT(NR)
      AINIT(NR) = ALFA
      END(MAAT) = ALFA
      MAXIM = MODUL2
      MINIM = MODUL1
      PART(NR) = PAT
      PART(MAAT) = PAT1
      GOTO 80
   70 EPSIL(NR) = MODUL1
      EPSIL(MAAT) = MODUL2
      END(MAAT) = BETA
      END(NR) = ALFA
      AINIT(MAAT) = ALFA
      MAXIM = MODUL1
      MINIM = MODUL2
      PART(NR) = PAT1
      PART(MAAT) = PAT
   80 KOUNT = KOUNT+42
C TEST ON THE NUMBER OF FUNCTION EVALUATIONS
      IF(KOUNT.GE.MAX) GOTO 190
   90 GOTO (100,110),KEY
C TEST ON ABSOLUTE ACCURACY
  100 IF(EST.LE.EPS) GOTO 190
      GOTO 120
C TEST ON RELATIVE ACCURACY
  110 IF(ABS(EPS*BINT).LE.TOL) GOTO 100
      IF(EST.LE.ABS(EPS*BINT)) GOTO 190
  120 IF(MAAT.EQ.1) GOTO 30
      IF(MAAT.GT.2) GOTO 130
      RANG(2) = 2
      GOTO 40
  130 MB = MA-1
C SEARCH FOR THE SUBINTERVAL WITH LARGEST ERROR
      DO 140 I = 2,MB
        IR = RANG(I)
        IF(MAXIM.GE.EPSIL(IR)) GOTO 150
  140   RANG(I-1) = RANG(I)
      RANG(MB) = NR
      RANG(MA) = MAAT
      GOTO 40
  150 RANG(I-1) = NR
      DO 160 K = I,MB
        IR = RANG(K)
        IF(MINIM.GE.EPSIL(IR)) GOTO 170
  160   CONTINUE
      RANG(MA) = MAAT
      GOTO 40
  170 DO 180 I = K,MB
        KK = MB-I+K
  180   RANG(KK+1) = RANG(KK)
      RANG(K) = MAAT
      GOTO 40
C CALCULATION OF THE INTEGRAL
  190 AIND1 = 0.0D+0
      DO 200 K = 1,MAAT
  200   AIND1 = AIND1+PART(K)
      IF(AIND1.EQ.0.0D+0)
     &        WRITE(6,*) '**** AIND=0.**** EST NOT CALCULATED'
      IF(AIND1.NE.0.0D+0) EST=EST/AIND1
      DAIND1=AIND1
      RETURN
      END





c equation (23) in PLB347 (1995) 143 - 151
c=================================================
      real*8 function fjjm(xi)
      implicit none
      real*8 pi, xi, term1, term2
      parameter (pi = 3.1415926535897932384626433832795D0)
      term1 = dsqrt(xi)
      term2 = dsqrt(4.d0 + xi)
      fjjm = 4.d0/term1/term2*dlog((term2 + term1)/(term2 - term1))
      return
      end

c  This gives the transverse Born coefficient as shown in fig 6a of
c  NPB392 (1993) 162 - 229.  For QCD take tf = 1d0/2d0, for QED take
c tf = 1d0 
c eta = (s - 4d0*m2)/4d0/m2, s is the gamma* gluon (gamma) CM Energy
c xi = Q^2/m2
c=================================================
      double precision function born_t(eta,xi)
      implicit none
      double precision eta, xi, pi, ca, cf, tf
      parameter (ca = 3d0, cf = 4d0/3d0, tf = 0.5d0)
      parameter(pi = 3.1415926535897932384626433832795D0)
      born_t = 0.5d0*pi*tf*(1.d0 + eta + 0.25d0*xi)**(-3)*
     #         (-2.d0*((1.d0 + eta - 0.25d0*xi)**2 + eta + 1.d0)*
     #         dsqrt(eta/(1.d0 + eta)) + (2.d0*(1.d0 + eta)**2 +
     #         0.125d0*xi**2 + 2.d0*eta + 1.d0)*
     #         dlog((dsqrt(1.d0 + eta) + dsqrt(eta))/
     #              (dsqrt(1.d0 + eta) - dsqrt(eta))))
      return
      end

c Longitudinal coefficient function, see above for additional comments
c Fig 6b of NPB392(1993) 162 - 229.
c=================================================
      double precision function born_l(eta,xi)
      implicit none
      double precision eta, xi, pi, tf, ca, cf
      parameter (ca = 3d0, cf = 4d0/3d0, tf = 0.5d0)
      parameter(pi = 3.1415926535897932384626433832795D0)
      born_l = 0.5d0*pi*tf*xi*(1.d0 + eta + 0.25d0*xi)**(-3.d0)*
     #         (2.d0*dsqrt(eta*(1.d0 + eta)) -
     #         dlog((dsqrt(1.d0 + eta) + dsqrt(eta))/
     #               (dsqrt(1.d0 + eta) - dsqrt(eta))))
      return
      end

c=================================================
       double precision function dilogjs(x)
       implicit double precision  (a-z)
       dimension b(8)
       integer ncall
       data ncall/0/,pi6/1.644934066848226d+00/,een,vier/1.d+00,.25d+00/
       ncall = 0
       if(ncall.eq.0)go to 2
1      if(x.lt.0)go to 3
       if(x.gt.0.5)go to 4
       z=-dlog(1.-x)
7      z2=z*z
       dilogjs=z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een)-z2*vier
       if(x.gt.een)dilogjs=-dilogjs-.5*u*u+2.*pi6
       return
2      b(1)=een
       b(2)=een/36.
       b(3)=-een/3600.
       b(4)=een/211680.
       b(5)=-een/(30.*362880.d+00)
       b(6)=5./(66.*39916800.d+00)
       b(7)=-691./(2730.*39916800.d+00*156.)
       b(8)=een/(39916800.d+00*28080.)
       ncall=1
       go to 1
3      if(x.gt.-een)go to 5
       y=een/(een-x)
       z=-dlog(een-y)
       z2=z*z
       u=dlog(y)
       dilogjs=z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een)-z2*vier-u*(z+.5*u)-pi6
       return
4      if(x.ge.een)go to 10
       y=een-x
       z=-dlog(x)
6      u=dlog(y)
       z2=z*z
       dilogjs=-z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een-u)+z2*vier+pi6
       if(x.gt.een)dilogjs=-dilogjs-.5*z*z+pi6*2.
       return
5      y=een/(een-x)
       z=-dlog(y)
       z2=z*z
       dilogjs=-z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een)-z2*vier
       return
10     if(x.eq.een)go to 20
       xx=1./x
       if(x.gt.2.)go to 11
       z=dlog(x)
       y=1.-xx
       go to 6
11     u=dlog(x)
       z=-dlog(1.-xx)
       go to 7
20     dilogjs=pi6
       return
       end


c  convolutes the pdf with the born level partonic cross section
c=================================================
      real*8 function fborn(z)
      implicit real*8 (a-h, o-z)
      real*8 z, m2, mq2, eta, xi, xbj, gluexz
      real*8 fbornt, fbornl
      integer ifch, nborn, ngcorr, nqcorr
      common/include/nborn, ngcorr, nqcorr, ifch
      common/glue/gluexz
      common/mass/m2, mq2, xbj
      fborn = 0.d0
      xi = mq2/m2
      eta = mq2*(1.d0 - z)/z/4.d0/m2 - 1.d0
      if (ifch .eq. 0) then
         fbornl = gluexz*born_l(eta,xi)
         fbornt = 0.d0
      elseif (ifch .eq. 1) then
         fbornt = gluexz*born_t(eta,xi)
         fbornl = 0.d0
      elseif (ifch .eq. 2) then
         fbornt = gluexz*born_t(eta,xi)
         fbornl = gluexz*born_l(eta,xi)
      endif
      fborn = fbornl + fbornt
      return
      end

c  convolutes the pdf with the O(alpha_s) gluon-initiated partonic 
cross section to calculate the appropriate structure function
c=================================================
      real*8 function fgcorr(z)
      implicit none
      real*8 z, m2, mq2, eta, xi, xbj, gluexz
      real*8 beta, rho
      real*8 fgtcorr, fglcorr
      real*8 alphas, pi
      real*8 scale, nlf, eh2, mu2
      real*8 xclca, xclcf, xctca, xctcf, xclbar, xctbar
      real*8 xsclca, xsclcf, xsctca, xsctcf, xsclba, xsctba
      real*8 thresha_t, threshf_t, threshbar_t
      real*8 thresha_l, threshf_l, threshbar_l
      real*8 asymp_t, asymp_l, asympbar_t, asympbar_l
      real*8 ca, cf, tf
      parameter (ca = 3d0, cf = 4d0/3d0, tf = 0.5d0)
      parameter (pi = 3.1415265359d0)
      integer ifch, nborn, ngcorr, nqcorr
      common/include/nborn, ngcorr, nqcorr, ifch
      common/glue/gluexz
      common/mass/m2, mq2, xbj
      common/coupling/scale, alphas, nlf, eh2
      xi = mq2/m2
      eta = mq2*(1.d0 - z)/z/4.d0/m2 - 1.d0
      beta = dsqrt(eta/(1d0 + eta))
      rho = 1d0/(1d0 + eta)
      mu2 = scale*scale
      if (ifch .eq. 0) then
         call sclca(eta,xi,xsclca)
         xclca = ca*tf*(xsclca + beta*asymp_l(xi) +
     #        rho*thresha_l(eta,xi))
         call sclcf(eta,xi,xsclcf)
         xclcf = cf*tf*(xsclcf + rho*threshf_l(eta,xi))
         call sclbar(eta,xi,xsclba)
         xclbar = ca*tf*(xsclba + beta*asympbar_l(xi) +
     #        rho*threshbar_l(eta,xi))
         fglcorr = xclca + xclcf + xclbar*dlog(mu2/m2)
         fgtcorr = 0.d0
      elseif (ifch .eq. 1) then
         call sctca(eta,xi,xsctca)
         xctca = ca*tf*(xsctca + beta*asymp_t(xi) +
     #        rho*thresha_t(eta,xi))
         call sctcf(eta,xi,xsctcf)
         xctcf = cf*tf*(xsctcf + rho*threshf_t(eta,xi))
         call sctbar(eta,xi,xsctba)
         xctbar = ca*tf*(xsctba + beta*asympbar_t(xi) +
     #        rho*threshbar_t(eta,xi))
         fgtcorr = xctca + xctcf + xctbar*dlog(mu2/m2)
         fglcorr = 0.d0
      elseif (ifch .eq. 2) then
         call sclca(eta,xi,xsclca)
         xclca = ca*tf*(xsclca + beta*asymp_l(xi) +
     #        rho*thresha_l(eta,xi))
         call sclcf(eta,xi,xsclcf)
         xclcf = cf*tf*(xsclcf + rho*threshf_l(eta,xi))
         call sclbar(eta,xi,xsclba)
         xclbar = ca*tf*(xsclba + beta*asympbar_l(xi) +
     #        rho*threshbar_l(eta,xi))
         fglcorr = xclca + xclcf + xclbar*dlog(mu2/m2)
         call sctca(eta,xi,xsctca)
         xctca = ca*tf*(xsctca + beta*asymp_t(xi) +
     #        rho*thresha_t(eta,xi))
         call sctcf(eta,xi,xsctcf)
         xctcf = cf*tf*(xsctcf + rho*threshf_t(eta,xi))
         call sctbar(eta,xi,xsctba)
         xctbar = ca*tf*(xsctba + beta*asympbar_t(xi) +
     #        rho*threshbar_t(eta,xi))
         fgtcorr = xctca + xctcf + xctbar*dlog(mu2/m2)
      endif
      fgcorr = 4d0*pi*alphas*eh2*gluexz*(fglcorr + fgtcorr)
      return
      end

c  assembles the appropriate partonic cross section for the structure function
c=================================================
      real*8 function fqcorr(z)
      implicit none
      real*8 z, m2, mq2, eta, xi, xbj
      real*8 beta
      real*8 fqtcorr, fqlcorr
      real*8 pi, alphas, scale, nlf, eh2
      real*8 mu2
      real*8 xchql, xclql, xqlbar, xchqt, xclqt, xqtbar
      real*8 xschql, xsclql, xsqlba, xschqt, xsclqt, xsqtba
      real*8 sumq3, sumelq3
      real*8 sumq4, sumelq4
      real*8 sumq5, sumelq5
      real*8 thresha_t, threshf_t, threshbar_t
      real*8 thresha_l, threshf_l, threshbar_l
      real*8 asymp_t, asymp_l, asympbar_t, asympbar_l
      real*8 ca, cf, tf
      parameter (ca = 3d0, cf = 4d0/3d0, tf = 0.5d0)
      character*10 q
      parameter (pi = 3.1415926535897932384626433832795D0)
      integer ifch, nborn, ngcorr, nqcorr
      common/coupling/scale, alphas, nlf, eh2
      common/include/nborn, ngcorr, nqcorr, ifch
      common/mass/m2, mq2, xbj
      common/lq3/sumq3,sumelq3
      common/lq4/sumq4,sumelq4
      common/lq5/sumq5,sumelq5
      common/quark/q
      mu2 = scale*scale
      xi = mq2/m2
      eta = mq2*(1.d0 - z)/z/4.d0/m2 - 1.d0
      beta = dsqrt(eta/(1d0 + eta))
      if (ifch .eq. 0) then
         call schql(eta,xi,xschql)
         call sclql(eta,xi,xsclql)
         call sqlbar(eta,xi,xsqlba)
         xchql = cf*tf*(xschql + beta**3*asymp_l(xi))
         xqlbar = cf*tf*(xsqlba + beta**3*asympbar_l(xi))
         xclql = cf*tf*xsclql
         if (q .eq. 'charm') then
            fqlcorr = sumq3*eh2*(xchql + xqlbar*dlog(mu2/m2))
     #           + sumelq3*xclql
         elseif (q .eq. 'bottom') then
            fqlcorr = sumq4*eh2*(xchql + xqlbar*dlog(mu2/m2))
     #           + sumelq4*xclql
         elseif (q .eq. 'top') then
            fqlcorr = sumq5*eh2*(xchql + xqlbar*dlog(mu2/m2))
     #           + sumelq5*xclql
         endif
         fqtcorr = 0.d0
      elseif (ifch .eq. 1) then
         call schqt(eta,xi,xschqt)
         call sclqt(eta,xi,xsclqt)
         call sqtbar(eta,xi,xsqtba)
         xchqt = cf*tf*(xschqt + beta**3*asymp_t(xi))
         xqtbar = cf*tf*(xsqtba + beta**3*asympbar_t(xi))
         xclqt = cf*tf*xsclqt
         if (q .eq. 'charm') then
            fqtcorr = sumq3*eh2*(xchqt + xqtbar*dlog(mu2/m2))
     #           + sumelq3*xclqt
         elseif (q .eq. 'bottom') then
            fqtcorr = sumq4*eh2*(xchqt + xqtbar*dlog(mu2/m2))
     #           + sumelq4*xclqt
         elseif (q .eq. 'top') then
            fqtcorr = sumq5*eh2*(xchqt + xqtbar*dlog(mu2/m2))
     #           + sumelq5*xclqt
         endif
         fqlcorr = 0.d0
      elseif (ifch .eq. 2) then
         call schql(eta,xi,xschql)
         call sclql(eta,xi,xsclql)
         call sqlbar(eta,xi,xsqlba)
         xchql = cf*tf*(xschql + beta**3*asymp_l(xi))
         xqlbar = cf*tf*(xsqlba + beta**3*asympbar_l(xi))
         xclql = cf*tf*xsclql
         if (q .eq. 'charm') then
            fqlcorr = sumq3*eh2*(xchql + xqlbar*dlog(mu2/m2))
     #           + sumelq3*xclql
         elseif (q .eq. 'bottom') then
            fqlcorr = sumq4*eh2*(xchql + xqlbar*dlog(mu2/m2))
     #           + sumelq4*xclql
         elseif (q .eq. 'top') then
            fqlcorr = sumq5*eh2*(xchql + xqlbar*dlog(mu2/m2))
     #           + sumelq5*xclql
         endif
         call schqt(eta,xi,xschqt)
         call sclqt(eta,xi,xsclqt)
         call sqtbar(eta,xi,xsqtba)
         xchqt = cf*tf*(xschqt + beta**3*asymp_t(xi))
         xqtbar = cf*tf*(xsqtba + beta**3*asympbar_t(xi))
         xclqt = cf*tf*xsclqt
         if (q .eq. 'charm') then
            fqtcorr = sumq3*eh2*(xchqt + xqtbar*dlog(mu2/m2))
     #           + sumelq3*xclqt
         elseif (q .eq. 'bottom') then
            fqtcorr = sumq4*eh2*(xchqt + xqtbar*dlog(mu2/m2))
     #           + sumelq4*xclqt
         elseif (q .eq. 'top') then
            fqtcorr = sumq5*eh2*(xchqt + xqtbar*dlog(mu2/m2))
     #           + sumelq5*xclqt
         endif
      endif
      fqcorr = 4d0*pi*alphas*(fqlcorr + fqtcorr)
      return
      end

c=================================================
      subroutine mrschm(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      q2=q*q
      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99
      if(x.lt.xmin.or.x.gt.xmax)       print 98,x
          if(mode.eq.1) then
        call mrsr1(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      endif 
  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ')
  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ','x= ',e10.5)
      return
      end

c=================================================
      subroutine mrsr1(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
c        open(unit=33,file='mrst.dat',status='old')
        open(unit=33,file='/userdisk/zeusqcd/develop/cor01.dat',
     +status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      

c For gamma* gluon group structure tf = 1d0/2d0,ca = 3d0, cf = 4d0/3d0
c For gamma* gamma group structure tf = 1d0 ,ca = 0d0, cf = 4d0/3d0
c eta = (W^2 - 4d0*m2)/4d0/m2, W is cm energy of gamma* gluon (gamma) system
c xi = mq2/m2 
c mq2 = Q^2 (positive)
c calculates the complete NLO longitudinal partonic cross section
c for gamma* gluon (photon) -> Q + Qbar + gluon modulo factors put in
c in the main program.
c=================================================
      double precision function gcorr_l(eta,xi)
      implicit none
      double precision ca, cf, tf, gcorra_l, gcorrf_l, gcorrbar_l
      double precision eta, xi, pi, beta, rho
      double precision scale, alphas, nlf, eh2, m2, mq2, xbj, mu2
      double precision clca, clcf, clbar
      double precision S
      double precision asymp_l, thresha_l, asympbar_l
      double precision threshbar_l, threshf_l
      parameter (pi = 3.1415926535897932384626433832795D0)
      parameter (ca = 3d0 , cf = 4d0/3d0, tf = 0.5d0)
      common/coupling/scale, alphas, nlf, eh2
c      common/mass/m2, mq2, xbj, S
      common/mass/m2, mq2, xbj
      mu2 = scale*scale
      beta = dsqrt(eta/(1.d0 + eta))
      rho = 1.d0/(1.d0 + eta)
      gcorra_l = 0.d0
      gcorrbar_l = 0.d0
      gcorrf_l = 0.d0
c We have separated the contributions by group structure
c we have also separated the threshold and asymptotic dependence from the
c full result, and what is left is contained in these subroutines
c clca corresponds to (33) of PLB347 (1995) 143 - 151.
c clcf corresponds to (34) of PLB347 (1995) 143 - 151.
c clbar corresponds to (35) of PLB347 (1995) 143 - 151.
      call sclca(eta,xi,clca)
      call sclbar(eta,xi,clbar)
      call sclcf(eta,xi,clcf)
      gcorra_l = ca*tf*(beta*asymp_l(xi) + rho*thresha_l(eta,xi)
     #     + clca)
      gcorrf_l = cf*tf*(rho*threshf_l(eta,xi) + clcf)
      gcorrbar_l = ca*tf*(beta*asympbar_l(xi) + 
     #     rho*threshbar_l(eta,xi) + clbar)*dlog(mu2/m2)
      gcorr_l = 4.d0*pi*alphas*eh2*(gcorra_l + gcorrf_l + gcorrbar_l)
      return
      end

c calculates the complete NLO transverse partonic cross section
c for gamma* gluon (photon) -> Q + Qbar + gluon modulo factors put in
c in the main program.
c for additional information see top of file 
c=================================================
      double precision function gcorr_t(eta,xi)
      implicit none
      double precision ca, cf, tf, gcorra_t, gcorrf_t, gcorrbar_t
      double precision eta, xi, pi, beta, rho
      double precision scale, alphas, nlf, eh2, m2, mq2, xbj, mu2
      double precision ctca, ctcf, ctbar
      double precision S
      double precision asymp_t, thresha_t, asympbar_t
      double precision threshbar_t, threshf_t
      parameter (pi = 3.1415926535897932384626433832795D0)
      parameter (ca = 3d0 , cf = 4d0/3d0, tf = 0.5d0)
      common/coupling/scale, alphas, nlf, eh2
c      common/mass/m2, mq2, xbj, S
      common/mass/m2, mq2, xbj
      mu2 = scale*scale
      beta = dsqrt(eta/(1.d0 + eta))
      rho = 1.d0/(1.d0 + eta)
      gcorra_t = 0.d0
      gcorrf_t = 0.d0
      gcorrbar_t = 0.d0
      call sctca(eta,xi,ctca)
      call sctbar(eta,xi,ctbar)
      call sctcf(eta,xi,ctcf)
      gcorra_t = ca*tf*(beta*asymp_t(xi) + rho*thresha_t(eta,xi)
     #     + ctca)
      gcorrf_t = cf*tf*(rho*threshf_t(eta,xi) + ctcf)
      gcorrbar_t = ca*tf*(beta*asympbar_t(xi)
     #     + rho*threshbar_t(eta,xi) + ctbar)*dlog(mu2/m2)
      gcorr_t = 4.d0*pi*alphas*eh2*(gcorra_t + gcorrf_t + gcorrbar_t)
      return
      end


c routine taken out of Numerical Recipes
c=================================================
      Subroutine Locat(xx,n,x,j)
      Integer j,n
      Double Precision x,xx(n)
      Integer jl,ju,jm
      jl = 0
      ju = n+1
 10   If (ju - jl .gt. 1) then
         jm = (ju + jl)/2
         If ((xx(n) .gt. xx(1)) .eqv. (x .gt. xx(jm))) then
            jl = jm
         else
            ju = jm
         endif
         goto 10
      endif
      j = jl
      return
      End

c eta = (W^2 - 4d0*m2)/4d0/m2, W is cm energy of gamma* light quark system
c xi = mq2/m2 
c mq2 = Q^2 (positive)
c calculates the complete NLO longitudinal partonic cross section
c modulo factors put in the main program for the reaction
c  gamma* light quark -> Q + Qbar + light quark
c  convolutes the partonic cross section with the parton density
c=================================================
      real*8 function qcorr_l(eta,xi)
      implicit none
      real*8 ca, cf, tf, qcorrh_l, qcorrl_l, qcorrbar_l
      real*8 eta, xi, beta, pi
      real*8 scale, alphas, nlf, eh2, m2, mq2, xbj, mu2
      real*8 dqll, cqhl, cqhlbar
      real*8 sumq3, sumelq3
      real*8 sumq4, sumelq4
      real*8 sumq5, sumelq5
      real*8 s
      real*8 asymp_l, asympbar_l
      character*10 quark
      parameter (pi = 3.1415926535897932384626433832795D0)
      parameter (ca = 3d0, cf = 4d0/3d0, tf = 0.5d0)
      common/coupling/scale, alphas, nlf, eh2
c      common/mass/m2, mq2, xbj, S
      common/mass/m2, mq2, xbj
c common blocks for the parton densities depending on the number of flavors
      common/lq3/sumq3,sumelq3
      common/lq4/sumq4,sumelq4
      common/lq5/sumq5,sumelq5
      common/quark/quark
      mu2 = scale*scale
      beta = dsqrt(eta/(1.d0 + eta))
      call schql(eta, xi, cqhl)
      call sqlbar(eta, xi, cqhlbar)
      call sclql(eta, xi, dqll)
c To get charm results
      if (quark .eq. 'charm') then
         qcorrh_l = cf*tf*sumq3*eh2*(
     #        (beta**3)*asymp_l(xi) + cqhl)
         qcorrbar_l = cf*tf*sumq3*eh2*((beta**3)*
     #        asympbar_l(xi) + cqhlbar)*dlog(mu2/m2)
         qcorrl_l = cf*tf*sumelq3*dqll
c To get bottom results
      elseif (quark .eq. 'bottom') then
         qcorrh_l = cf*tf*sumq4*eh2*(
     #        (beta**3)*asymp_l(xi) + cqhl)
         qcorrbar_l = cf*tf*sumq4*eh2*((beta**3)*
     #        asympbar_l(xi) + cqhlbar)*dlog(mu2/m2)
         qcorrl_l = cf*tf*sumelq4*dqll
c To get top results
      elseif (quark .eq. 'top') then
         qcorrh_l = cf*tf*sumq5*eh2*(
     #        (beta**3)*asymp_l(xi) + cqhl)
         qcorrbar_l = cf*tf*sumq5*eh2*((beta**3)*
     #        asympbar_l(xi) + cqhlbar)*dlog(mu2/m2)
         qcorrl_l = cf*tf*sumelq5*dqll
      endif
      qcorr_l = 4.d0*pi*alphas*(qcorrh_l + qcorrbar_l + qcorrl_l)
      return
      end

c calculates the complete NLO transverse partonic cross section
c for gamma* light quark -> Q + Qbar + light quark modulo factors
c put in the main program.
c for additional information see top of file
c  convolutes the partonic cross section with the parton density
c=================================================
      real*8 function qcorr_t(eta,xi)
      implicit none
      real*8 ca, cf, tf, qcorrh_t, qcorrl_t, qcorrbar_t
      real*8 eta, xi, beta, pi
      real*8 scale, alphas, nlf, eh2, m2, mq2, xbj, mu2
      real*8 d1tq, cqht, cqhtbar, d1qbar, dqtb
      real*8 sumq3, sumelq3
      real*8 sumq4, sumelq4
      real*8 sumq5, sumelq5
      real*8 S
      real*8 asymp_t, asympbar_t
      character*10 quark
      parameter (pi = 3.1415926535897932384626433832795D0)
      parameter (ca = 3d0, cf = 4d0/3d0, tf = 0.5d0)
      common/coupling/scale, alphas, nlf, eh2
c      common/mass/m2, mq2, xbj, S
       common/mass/m2,mq2,xbj
c common blocks for the parton densities depending on the number of flavors
      common/lq3/sumq3,sumelq3
      common/lq4/sumq4,sumelq4
      common/lq5/sumq5,sumelq5
      common/quark/quark
      mu2 = scale*scale
      beta = dsqrt(eta/(1.d0 + eta))
      call schqt(eta, xi, cqht)
      call sqtbar(eta, xi, cqhtbar)
      call sclqt(eta, xi, d1tq)
      call sdqtb(eta, xi, dqtb)
c To get charm results
      if (quark .eq. 'charm') then
         qcorrh_t = cf*tf*sumq3*eh2*(
     #        (beta**3)*asymp_t(xi) + cqht)
         qcorrbar_t = cf*tf*sumq3*eh2*((beta**3)*
     #        asympbar_t(xi) + cqhtbar)*dlog(mu2/m2)
         qcorrl_t = cf*tf*sumelq3*(d1tq + dqtb)
c To get bottom results
      elseif (quark .eq. 'bottom') then
         qcorrh_t = cf*tf*sumq4*eh2*(
     #        (beta**3)*asymp_t(xi) + cqht)
         qcorrbar_t = cf*tf*sumq4*eh2*((beta**3)*
     #        asympbar_t(xi) + cqhtbar)*dlog(mu2/m2)
         qcorrl_t = cf*tf*sumelq4*(d1tq + dqtb)
c To get top results
      elseif (quark .eq. 'top') then
         qcorrh_t = cf*tf*sumq5*eh2*(
     #        (beta**3)*asymp_t(xi) + cqht)
         qcorrbar_t = cf*tf*sumq5*eh2*((beta**3)*
     #        asympbar_t(xi) + cqhtbar)*dlog(mu2/m2)
         qcorrl_t = cf*tf*sumelq5*(d1tq + dqtb)
      endif
      qcorr_t = 4.d0*pi*alphas*(qcorrh_t + qcorrbar_t + qcorrl_t)
      return
      end
         
