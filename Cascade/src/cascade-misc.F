      FUNCTION DOT(A,B)
      IMPLICIT NONE
C+++++++++++
C    DOT PRODUCT OF FOUR VECTOR IN MINKOWSKI METRIK
C+++++++++++
      Double Precision A(4),B(4),DOT
      DOT = A(4)*B(4)-A(1)*B(1)-A(2)*B(2)-A(3)*B(3)
      RETURN
      END
      FUNCTION DOT1(I,J)
	IMPLICIT NONE
#include "caslujets.inc"
      Double Precision DOT1
	INTEGER I,J
C+++++++++++
C    DOT PRODUCT OF FOUR VECTOR IN MINKOWSKI METRIK
C    WITH VECTORS FROM LUJETS
C+++++++++++
      DOT1= P(I,4)*P(J,4)-(P(I,1)*P(J,1))-(P(I,2)*P(J,2))-
     .     (P(I,3)*P(J,3))
      RETURN
      END
      SUBROUTINE CCFMSTFG(KF,X,Q2,P,XPQ)
      Implicit None
#include "casshwr.inc"
      Integer n1,n2,n3
      Double Precision XPQ(-6:6),X,Q2,P
*! new
      Parameter (n1=51,n2=51,n3=51)
      Double Precision Q2x,xx,px,xpd
      DIMENSION Q2X(0:n1),XX(0:n2),PX(0:n3),XPD(0:n2,0:n1,0:n3)
      Double Precision amax
      Dimension amax(0:n2)
      Double Precision  kl,kt,xpqm
      double Precision XA(3),f_grid(n1,n2,n3),A(N1+N2+N3)
      INTEGER NA(3)
      DATA NA/n1,n2,n3/
      Double Precision DHFINT
	Integer  ikincut,Ipgg,ns_sel
	Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Integer IGLU
      Common/CAGLUON/Iglu
      Integer i,j,k,in,irr,igrid,ip,n2test
      Integer ipoint,kf,it
	Integer Iset
      Double Precision scal,rx,rq2,rp,rkms,glu,xpqmax,xm
      Double Precision xl,pl
	character *72 TXT
	character adum
      character*512 pdfname
      character*512 PDFPATH
      common/CASPDF/PDFPATH
      INTEGER PDFPATHLEN
      LOGICAL FIRST,afirst
      DATA FIRST/.TRUE./
      Data Afirst/.TRUE./
      SAVE   
      
	if(kf.eq.22.or.iabs(kf).eq.11) then
	else
         write(6,*) ' ccfmstfg: no valid kf = ',kf
         write(6,*) ' ccfmstfg: only applicable for gamma' 	   
	   stop
	endif
      IF(FIRST) THEN
         i=0
         scal = 1.0
c         scal = 0.318
	   Iset = Iglu
         write(6,*) ' ccfmtfq pdfpath:',pdfpath
         PDFPATHLEN=INDEX(PDFPATH,' ')-1
         write(6,*) ' read grid files from ',PDFPATH(1:PDFPATHLEN)
         PDFNAME=PDFPATH(1:PDFPATHLEN)//'ccfm-quark.dat'
         if(iset.eq.1) then
            IF(ICCFM.EQ.0) THEN
               PDFNAME=PDFPATH(1:PDFPATHLEN)//'dglap-gam.dat'
               open(30,FILE=PDFNAME, 
     +          FORM='formatted',STATUS=
     +         'OLD', IOSTAT=IRR,ERR=80 )
               write(6,*) ' read unintegrated gluon from file: '
     +         //'dglap-gam.dat'
            ELSEIF(ICCFM.EQ.1) THEN
               PDFNAME=PDFPATH(1:PDFPATHLEN)//'ccfm-gam.dat'
               open(30,FILE=PDFNAME, FORM='formatted',
     +         STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
               write(6,*) ' read unintegrated gluon from file: '
     +         //'ccfm-gam.dat'
            ENDIF
         else
            write(6,*) ' No valid IGLU selected for ccfmstfg: IGLU = ',
     +      Iset
            STOP
         endif
c         read(30,10000) Qg0,ikincut
c10000    format(' Qg0 = ',f12.8,' ikincut= ',I6)
c         read(30,10100)
c10100    format('xg,  kt, p  xgx')
c         write(6,*) ' soft cut Q0 ',Qg0,' scale factor = ',scal
c         write(6,*) ' kin cut ',ikincut
200      Read(30,101) TXT
  101    Format(A72)
         If(TXT(1:4).EQ.'  Qg') then 
	      read(txt,1000) Qg0,ikincut
1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Qg0') then 
	      read(txt,1000) Qg0,ikincut
c         WRITE(6,101) ' 1st line ',TXT
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
         If(TXT(1:4).EQ.' ln(') then 
c         WRITE(6,101) '2 or 3rd line',TXT
         Endif
	   
         write(6,*) ' soft cut Q0 ',Qg0,' scale factor = ',scal
         write(6,*) ' kin cut ',ikincut,' Ipgg = ',Ipgg,
     &	   ' ns_sel = ',ns_sel
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


         do j=1,n1
            do i=1,n2
               do k=1,n3
                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS
                  xx(i) = rx
                  q2x(j) = rq2
                  px(k) = rp
                  xpd(i,j,k) = max(0.,rkms)*scal
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
                  f_grid(j,i,k) = xpd(i,j,k)
               ENDDO
            ENDDO
         ENDDO
         FIRST=.FALSE.
         write(6,*) '  parton densities read from file unit 30 '
         Close(30)
	ENDIF



c      XA(1) = Q2
c      XA(2) = X
c      XA(3) = P
      XA(1) = log(Q2)
      XA(2) = log(X)
      XA(3) = log(P)
      if(xa(3).lt.px(1)) then
c         write(6,*) ' p out of range ',p,' min p ',exp(px(1))
         xa(3)=px(1)
      endif
      if(xa(3).gt.px(n3)) then
         write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
         xa(3)=px(n3)
      endif
      if(xa(2).ge.xx(n2)) xa(2)=xx(n2)-0.0001
      if(xa(1).lt.q2x(1)) then
c            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
         xa(1)=q2x(1)
      endif
      if(xa(1).ge.q2x(n1)) then
c            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
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
                  write(6,*) ' q2 not found ',i,q2
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
               write(6,*) ' x not found ',x,xa(2),j
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
                  write(6,*) ' p not found ',k,p
               endif
            endif
         endif
         glu = f_grid(i,j,k)/q2
      else
         glu = DHFINT(3,XA,NA,A,f_grid)/q2
      endif
      DO  IP=-6,6
         XPQ(IP)=0.0
      ENDDO
      if(glu.gt.0.) xpq(0) = glu
      if(glu.gt.1000000) then
         write(6,*) ' new ccfmstfg x,q2,p,xpq',x,q2,p,glu
         write(6,*) ' i,j,k ',i,j,k
      endif
      return
      ENTRY ccfmmaxg(x,xpqmax)

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
                     xm=f_grid(j,it,k)/exp(q2x(j))
c              write(6,*) f_grid(j,it,k),exp(q2x(j))
                     if(xm.gt.xpqmax) xpqmax=xm
                  Enddo
               Enddo
            Enddo
            amax(i)=xpqmax
         Enddo
         Afirst=.false.
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

      ENTRY ccfmming(x,kt,p,xpqm)
      xpqm = 100000.
      xl=log(x)
      kl=log(kt)
      pl=log(p)

      if(kl.lt.q2x(1)) then
c         write(6,*) '  k2 out of range ',kt,' min ',exp(q2x(1))
         kl=q2x(1)
      endif
      if(kl.ge.q2x(n1)) then
         write(6,*) ' ccfmmin: k2 out of range: p =',p,' x= ',x,n1
         write(6,*) ' ccfmmin: k2 out of range ',kt,' max ',exp(q2x(n1))
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
            write(6,*) ' x not found ',x,xl,j
         endif
      endif
      DO K=1,N3
         if(pl.ge.px(k)) then
            xm=max(f_grid(i,j,k)/kt,0.005)
            if(xm.lt.xpqm) xpqm=xm
         endif
      Enddo
c        write(6,*) ' ccfm x,p,xpqm ',x,p,xpqm
      return


   80 write(6,*) ' CCFMSTFG: error in opening file '
   90 write(6,*) ' end of file at ',i,j,k,RX,RQ2,RP,RKMS
      stop
      END
	subroutine mcint(fxnb,avgi,sd)
      Implicit None
      double precision fxnb,avgi,sd
      External fxnb
      Integer Mxdim
      Parameter (Mxdim = 50 )
      Double precision xlb,xub
      Integer ndim,ncall
      Common /mcint1/xlb(mxdim),xub(mxdim),ndim,ncall 
      Double precision fmax
      Common /mcint2/fmax 
      double precision ff,fsum,fsum2,sigma2,err
      Double precision xrn(mxdim),x(mxdim)
      Integer i,icall
      
      fmax=-1D20
      
      do icall = 1,ncall 
         call dcasrnv(xrn,ndim)
         do i=1,ndim
           x(i) = xlb(i) + xrn(i)*(xub(i)-xlb(I))
         end do
         ff = fxnb(x)
         if(ff.ge.fmax) fmax = ff
         fsum = fsum+ ff
         fsum2 = fsum2 + ff*ff
      end do
      
      fsum = fsum/ncall
      fsum2 = fsum2/ncall
      sigma2 = fsum2 - fsum*fsum
      err = sqrt(sigma2/ncall)
      write(6,*) ' ****************************************************'
      write(6,*) ' *                                                  *'
      write(6,*) ' *       Simple MC integration    (MCINT)           *'
      write(6,*) ' *                                                  *'
      write(6,*) ' *                                                  *'
      write(6,100) NCALL,NDIM
 100  FORMAT('  *        NCALL = ',I6,' NDIM = ',I5, '               *')
      write(6,*) ' ****************************************************'
      
      write(6,*) ' fsum = ',fsum,' err = ',err
      avgi = fsum
      sd = err
      return
      end
      FUNCTION FXNB(X)
      IMPLICIT None
	Double Precision FXNB,FXNB1,X,XG,XGF,FXN1,WEIGHT
	Integer MXDIM,I
      PARAMETER (MXDIM = 50)
      DIMENSION X(MXDIM),XGF(20)
      COMMON/XFXNB/XG(20)
      EXTERNAL FXN1
      FXNB = 0.0D0
	Weight = 0.D0
      DO 10 I=1,20
         XGF(I) = X(I)
         XG(I) = X(I)
   10 CONTINUE
      FXNB1 = FXN1(XGF,WEIGHT)
      FXNB = FXNB1
      RETURN
      END
      SUBROUTINE CUTG(PT2,NACC)
      IMPLICIT NONE
	Double Precision PT2
	INTEGER NACC
#include "casparas.inc"
      NACC=1
      IF(PT2.LT.PT2CUT(IPRO)) NACC=0
      RETURN
      END
      SUBROUTINE PHASE(NP,ET,XM,PCM,WT)
c to generate n-body phase space
c taken from: Collider Physics
c             V.D. Barger and R.J.N Phillips
c             Adison-Wesley Publishing Company, Inc.
c             1987
	Implicit none
	Integer NP
	Double Precision ET,XM,PCM,WT
	Integer N
	Double Precision WPS,PIE,M,MT,MS,P,K
      COMMON/KIN/ WPS,PIE,M(10),MT(10),MS(10),P(10,4),K(10,4),N
      DIMENSION XM(18),PCM(4,18)
	Integer I,J,Jmax
      N=NP
      PIE=4.D0*DATAN(1.D0)
      WPS=0.D0
      DO 10 I=1,N
         M(I)=0.D0
         MT(I)=0.D0
         MS(I)=0.D0
         DO 10 J=1,4
            K(I,J)=0.D0
            P(I,J)=0.D0
   10 CONTINUE
      DO 20 I=1,N
         M(I)=XM(I)
   20 CONTINUE
      K(1,1)=ET
      MT(1) =ET
      DO 30 I=1,N
   30 MS(1) = MS(1) + M(I)
      CALL ORDER
      JMAX = N - 1
      WPS = JMAX*((MT(1)-MS(1))/(4.D0*PIE**2))**(N-2)/
     +      (4.D0*PIE*MT(1))
      DO 40 J =1,JMAX
         WPS = WPS/J
         MS(J+1) = MS(J) - M(J)
         CALL GEN(J)
   40 CALL BOO(J)
      DO 50 J=1,N
         PCM(1,J) = P(J,2)
         PCM(2,J) = P(J,3)
         PCM(3,J) = P(J,4)
         PCM(4,J) = P(J,1)
   50 CONTINUE
      if(p(1,2).ne.p(1,2)) write(6,*) ' error in phase'
      WT=WPS*(2.D0*PIE)**(3*N-4)
c      write(6,*) ' phase ',wt
      RETURN
      END
      SUBROUTINE GEN(I)
      IMPLICIT None
	Integer N
	Double Precision WPS,PIE,M,MT,MS,P,K
      COMMON/KIN/ WPS,PIE,M(10),MT(10),MS(10),P(10,4),K(10,4),N
	Double Precision X,AALAM
      COMMON/XVAR/ X(10)
#include "casparas.inc"
      Double Precision PSQ,EI,PI,COST,COSTM,SINT,PHI
      Integer I
      IF(I.EQ.N-1) THEN
         MT(I+1) = M(N)
      ELSE
         MT(I+1) = MS(I+1) + (MT(1)-MS(1))*X(3*I)
      ENDIF
      PSQ = AALAM(MT(I)**2,M(I)**2,MT(I+1)**2)/(4.D0*MT(I)**2)
      IF(PSQ.LT.0.0D0) WRITE(6,*) 'PSQ.LT.0 ',PSQ
      IF(PSQ.LT.0.0D0) PSQ=0.D0
      EI = DSQRT(PSQ+M(I)**2)
      PI = DSQRT(PSQ)
      COST = 1.D0-2.D0*X(3*I-1)
      SINT = DSQRT(1.D0-COST**2)
      PHI = 2.D0*PIE*X(3*I-2)
      P(I,1) = EI
      P(I,2) = PI*SINT*DSIN(PHI)
      P(I,3) = PI*SINT*DCOS(PHI)
c      write(6,*) ' GEN: PT**2 ',P(I,2)**2+P(I,3)**2,I
      P(I,4) = PI*COST
      if(p(I,2).ne.p(i,2)) then
         write(6,*) '  error in gen '
         write(6,*) ' pi,sint,phi ',pi,sint,phi,X(3*I-2),X(3*I-1)
      endif
      WPS = WPS*PI
      RETURN
      END
      SUBROUTINE ORDER
      IMPLICIT None
	Integer N
	Double Precision WPS,PIE,M,MT,MS,P,K	
      COMMON/KIN/ WPS,PIE,M(10),MT(10),MS(10),P(10,4),K(10,4),N
      IF(N.LE.3) RETURN
      WRITE(6,*) ' N .GT. 3; N = ',N
      STOP
      END
      SUBROUTINE BOO(I)
      IMPLICIT None
	Integer N
	Double Precision WPS,PIE,M,MT,MS,P,K	
      COMMON/KIN/ WPS,PIE,M(10),MT(10),MS(10),P(10,4),K(10,4),N
	Integer I,J
	Double Precision BG1,G1,BG2,G2,BG3,G3,E,X,Z,Y
      BG1 = K(I,2)/MT(I)
      G1 = DSQRT(1.D0 + BG1**2)
      BG2 = K(I,3)/(MT(I)*G1)
      G2 = DSQRT(1.D0 + BG2**2)
      BG3 = K(I,4)/(MT(I)*G1*G2)
      G3 = DSQRT(1.D0 + BG3**2)
      E = G1*G2*G3*P(I,1) + BG1*G2*G3*P(I,2)
     +    + BG2*G3*P(I,3) + BG3*P(I,4)
      X = BG1*P(I,1) +G1*P(I,2)
      Y = BG2*G1*P(I,1) + BG2*BG1*P(I,2) +G2*P(I,3)
      Z = BG3*G2*G1*P(I,1) + BG3*G2*BG1*P(I,2)
     +    + BG3*BG2*P(I,3) +G3*P(I,4)
      P(I,1) = E
      P(I,2) = X
      P(I,3) = Y
      P(I,4) = Z
      IF(I.EQ.N-1) GOTO 20
      DO 10 J=1,4
   10 K(I+1,J) = K(I,J) - P(I,J)
C     WRITE(6,*) 'BOO ',(P(I,IK),IK=1,4)
      RETURN
   20 DO 30 J =1,4
   30 P(N,J) = K(I,J) - P(I,J)
      RETURN
      END
c get dcasrn for bases/spring
      FUNCTION DRN(IRN)
	Implicit None
	Integer IRN
      double precision drn1,DRN
      Double Precision dcasrn
	EXTERNAL dcasrn
      DRN1= dcasrn(IRN)
      DRN = DRN1
	RETURN
	END
	FUNCTION DFUN(N,X)
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION DFUN,FXN1,WEIGHT,X(20)
      EXTERNAL FXN1
      DFUN=0
      WEIGHT=0.D0
      DO 10 I=1,N
10    IF(X(I).EQ.0..OR.X(I).EQ.1.) RETURN
      DFUN = FXN1(X,WEIGHT)
      RETURN
      END
      FUNCTION AALAM(A,B,C)
      IMPLICIT NONE
	double precision A,B,C,AALAM
      AALAM = A**2 + B**2 + C**2 - 2.D0*A*B - 2.D0*B*C - 2.D0*C*A
      RETURN
      END
      Subroutine xsecqg(X,XSEC)
      IMPLICIT None
      Double Precision X(20)
#include "caslujets.inc"
#include "casparas.inc"
#include "caspart.inc"
#include "casluco.inc"
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
      WPART=0
	
      CALL P_SEMIH(X,WPART)
      IF(WPART.GT.0.) THEN
         CALL meqg(WT1)
      ENDIF
c        write(6,*) WPART,WT1,IPRO
      WEIG=WPART*WT1*GEV2NB
C BOOST BACK TO OVERALL CMS
      IF(WEIG.EQ.0.0) Return
      ncall = ncall + 1
c      call PYLIST(1)
      XSEC=WEIG
      RETURN
      END


