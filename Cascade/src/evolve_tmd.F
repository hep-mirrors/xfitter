      Subroutine evolve_tmd
      Include 'SMallx.inc'
       Integer Iev 
C---Standard LEP event common block 
      Integer NMXHEP,NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP 
      Double Precision PHEP,VHEP,EVWGT 
      PARAMETER (NMXHEP=4000) 
      COMMON/HEPEVTP/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP), 
     & JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP) 
C---Event weight 
C   N.B. weights must be renormalized to momentum sum 
      COMMON/HEPWGT/EVWGT 
      Integer nobran,ikincut
      common/myvar/nobran,ikincut

	Double Precision Qbarmy,Qbar_min,Qbar_max
      Common/mglubran/Qbarmy,Qbar_min,Qbar_max
      Integer neve
	common/myevt/neve
	Integer nloop
	common/myloop/nloop
      Double Precision x3lmin,x3lmax,x3ldif
	Integer Nbp
c      Parameter (Nbp=20)
      Parameter (Nbp=50)
c      Parameter (Nbp=1)
      Double Precision X3(0:Nbp+1)
      Double Precision X3M(0:Nbp)
      Double Precision x3b(0:Nbp)
      common/gridtt/x3m
c      save /gridtt/
	Integer ng_max
	Integer nrglu
	Common/mynrglu/nrglu
	Integer nmax,i,nx3,kev,ic
	Integer Ipgg,ns_sel
	Double Precision scal 
	common/Pggsel/Ipgg,ns_sel,scal
	Double Precision Qgmin
      character *72 TXT
      CHARACTER   FILNAME*132,testNAME*132
      Common/gludatf/filname
	Double Precision Xnorm
      common/smnorm/ Xnorm
      Logical pdflib,quark,gluon,photon,saturation
      Common /SMbran2/pdflib,quark,gluon,photon,saturation
      Integer Ioneloop,Ipdf,Itarget,Iglu,Isaturation
      Integer iparton
      Common /SMquark/iparton
      character *15 char
	Double Precision BB
	Common /splitting/ BB
      Double precision ininorm(-6:6)
      common/smininorm/ininorm
     
      Integer IRR
      
      open(50,FILE='steer_gluon-evolv', FORM='formatted',
     +         STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
      
*
*      read(50,*) filname
*      read(50,*) testname
c      write(6,*) ' output file ',filname
      Read(50,101) TXT
c      write(6,*) txt
      read(txt,1004) char,Xnorm
1004  format(a8,F16.8)
      write(6,*) txt,char,xnorm
      Read(50,101) TXT
      read(txt,1005) char,Ioneloop
      write(6,*) txt,char,Ioneloop
1005  format(a10,I8)
      Read(50,101) TXT
      read(txt,1010) char,Isaturation
1010  format(a14,I8)
      write(6,*) txt,char,Isaturation
      Read(50,101) TXT
      read(txt,1006) char,Ipdf
      write(6,*) txt,char,Ipdf
1006  format(a7,I8)
      Read(50,101) TXT
      read(txt,1007) char,Itarget
      write(6,*) txt,char,Itarget
1007  format(a10,I8)
      Read(50,101) TXT
      read(txt,1008) char,Iglu
      write(6,*) txt,char,Iglu
1008  format(a7,I8)
      Read(50,101) TXT
  101 Format(A72)
      read(txt,1000) char,Ipgg
      write(6,*) txt,char,Ipgg
1000  format(a7,I8)
      Read(50,101) TXT
      read(txt,1001) char,ns_sel
      write(6,*) txt,char,ns_sel
1001  format(a5,I8)
c      write(6,*) ' selection Ipgg = ',Ipgg,' ns_sel = ',ns_sel
      Read(50,101) TXT
c      write(6,*) txt
      read(txt,1011) char,ikincut
      write(6,*) txt,char,ikincut
1011  format(a10,I8)
      Read(50,101) TXT
      read(txt,1002) char,Qg
      write(6,*) txt,char,Qg
1002  format(a5,F16.8)
      Read(50,101) TXT
      read(txt,1003) char,Qs
      write(6,*) txt,char,Qs
1003  format(a5,F16.8)
      Read(50,101) TXT
      read(txt,1018) char,QCDlam
      write(6,*) txt,char,QCDlam
1018  format(a9,F16.8)
c      write(6,*) ' qcdlma ',qcdlam,' ikincut ',ikincut
      if(iglu.ne.0) then
          gluon = .true.
          else
          gluon = .false.
        endif 
      If(Ioneloop.eq.1) then
	  onel = .true.
	  else
	  onel=.false.
	Endif
      If(Isaturation.eq.1) then
	  saturation = .true.
	  else
	  saturation=.false.
	Endif
      if(ipdf.eq.0) then
	  pdflib=.false.
	  else
	  pdflib=.true.
	Endif
      if(itarget.eq.22) then
	  photon = .true.
	  else
	  photon = .false.
	Endif
      if(iglu.eq.0) then
        Read(50,101) TXT  
        read(txt,1009) char,Iparton
        write(6,*) txt,char,Iparton
        write(6,*) txt
1009    format(a9,I8)
      endif
	close(50)

C---Initialize run 
      Call SMinit 
	neve = 0
	nmax =nev
	Xini = 0.
      write(6,*) ' output file ',filname
      write(6,*) ' selection Ipgg = ',Ipgg,' ns_sel = ',ns_sel
      write(6,*) ' Qg = ',Qg,' Qs = ',Qs,' Xnorm = ',Xnorm
	write(6,*) ' Ipdf = ',Ipdf,' Itarget = ',Itarget
      write(6,*) ' BB = ',BB
c      if(iglu.eq.0) then
c        write(6,*) ' Iquark = ',Iparton
c      endif

         Qgmin = max(Qg-0.5d0,QCDlam)
         Qgmin = max(Qg,QCDlam)
	   x3lmin = log(Qgmin)  
	   x3lmax = log(qmax)
	   x3ldif = (x3lmax-x3lmin)/Real(Nbp)
ccc         write(6,*) ' qgmin = ',qgmin,x3lmin,x3ldif	
         Do I=0,Nbp+1
            x3(I) = exp(x3lmin + x3ldif*Real(I))
         Enddo
	   Do I=0,Nbp
	      x3m(i) = (x3(i) + x3(i+1))/2.
		x3b(i) = x3(i+1) - x3(i)
		write(6,*) ' gluon p_max ',x3m(i),i
         Enddo
	Nx3 = -1
C---Initialize analysis 
      Xini=0.
      Xfin=0.
      Call updfgrid(1)
      
c      write(6,*) ' after SMabmy BB = ',BB
c	Call updfgrid(1)


      Nx3 = 49
600   Nx3 = Nx3 + 1
      write(6,*) ' ng_max = ',ng_max,' at nx3-1 ',nx3-1
      IF(Nx3.gt.Nbp) Then
        write(6,*) ' evolve_tmd: Nx3 gt Npb -> Program stopped'
        stop
      Endif
	Qbarmy = x3m(Nx3)
	nloop = Nx3
	ng_max = 0

	write(6,*) ' Qbar_my  ',Qbarmy 
      
      do Iparton=0,2
      Jev = 0 
      write(6,*) ' evolving parton ID: ',Iparton 
        if(iparton.eq.0) then 
           gluon=.true.
           else
           gluon=.false.
        endif
        Xini = 0
        do I=1,nev
C---Initialize event 
          Call SMinfn 

C---Small-x gluon branching process 
          Call smbran 

	    Xini = Xini + x0
	    if(ng_max.le.nrglu) ng_max=nrglu
          

          If(wt.NE.0.0) then
             kev=i 
             if (kev.gt.0) ic=100000 
             if (mod(kev,ic).eq.0) write(6,*) ' event ',kev,nev,' loop P_max ',nloop,Qbarmy
	    endif
c          write(6,*) ' in loop2: Iparton = ',Iparton, I
         End do
      End do
C---Terminate analysis 

      call updfgrid(3)
      return
  80  Write(6,*) ' steering file for evolve_tmd not found '
      stop    
      End 
C------------------------------------------------------------------------
      Subroutine SMinit
C--Sets input parameters
      Include 'SMallx.inc'
      Integer nobran,ikincut
      common/myvar/nobran,ikincut
      Double Precision pbeam1,pbeam2
	Integer Iseed,minuts
	External minuts
	Double Precision Rdq,SMrgen
	Double Precision Ca
	Integer nf
      Double Precision Lambda,QQg	
	Common/hardme/Lambda,QQg,nf
	logical my_glu	
	common/glumy/my_glu
C for Sudakov (smqtem.f) B is fixed: B=0.5	
C choose BB: 0, 0.5 , 1.	
	Double Precision BB
	Common /splitting/ BB
	Integer Ipgg,ns_sel
	Double Precision scal 
	common/Pggsel/Ipgg,ns_sel,scal
	Double Precision Xnorm
      common/smnorm/ Xnorm
      Logical pdflib,quark,gluon,photon,saturation
      Common /SMbran2/pdflib,quark,gluon,photon,saturation
	Double Precision mz
	Data mz/91.188d0/
      Double Precision asbmy,test,pi
      external asbmy
	Logical First
	Common /SMstart/ First
      CHARACTER   FILENAME*132
      Common/gludatf/filename
	First = .true.


      write(6,*) ' 1st gluon = ',gluon
      if(photon) then
           gluon=.true.
        endif
        if(gluon) then
           quark=.false.
           else
           quark=.true.
        endif
        
      write(6,*) ' gluon = ',gluon
      powin=.true.
      BB=0.5
      scal = 1.	
      my_glu = .true.
	
c nobran = 0 --> standard initial state branching
c nobran = 1 --> no initial state branching
c nobran = 2 --> read ccfm.dat file
      nobran = 0
      write(6,*) ' smini: ikincut = ',ikincut,' QCDlam = ',QCDlam
      write(6,*) ' smini: Qg = ',Qg, ' Qs = ',Qs

      pow1=4.   ! F2 fit 
      pow0=0.0

c here include a factor of pi coming from d^2 k --> dkt^2	
      Xnor = Xnorm * 3.1415927

      Qx=Qs/Sqrt(2.)

      Maxpr=5
      Nmongo=27
C---Number of events
c      Nev=100000000 ! very long
      Nev=10000000 ! std
      Nev=5000000 ! old 
c      Nev=2000000
c      Nev=20000  
      Print *
      Print *,' Number of events requested:',Nev
C  LHC
C
      pbeam1=7000.
      pbeam2=7000.
C--Light quarks+ charm
      Qmas=0.25
	Ichg=4
      Idq=2

C--Set up beams: 1=electron, 2=proton
      p1(1)=0.
      p1(2)=0.
      p1(3)=-pbeam1
c      p1(5)=0.00051
      p1(5)=0.
      p1(4)=Sqrt(p1(3)**2+p1(5)**2)
      p2(1)=0.
      p2(2)=0.
      p2(3)=pbeam2
c      p2(5)=0.93827
      p2(5)=0.
      p2(4)=Sqrt(p2(3)**2+p2(5)**2)
      espec(5)=p1(5)
      pspec(5)=p2(5)
c      write(6,*) Sbig,p1(5)**2+p2(5)**2+2.*(p1(4)*p2(4)-p1(3)*p2(3))
      Sbig=p1(5)**2+p2(5)**2+2.*(p1(4)*p2(4)-p1(3)*p2(3))
C--Limit for branching
      Imax=1000
ccc	fixed xmin
	xymin = 1e-7
      xmin=xymin
      xmax=.9999999
C--y range (y=momentum fraction of virtual photon)
      ymin=xymin
      ymax=.9999999
C--Q range
      me=0.511E-3
      Qmin=ymin*me
      Qmax=sqrt(ymax*Sbig)
C  Min Q value
      Qmin=sqrt(1.)

C - now include also nf=4
      nf = 4
	QQg=Qg
	Lambda=QCDlam
	Ca = 3.
	Capb = 12.*Ca/(11*Ca - 2.*nf)
      Tg=Log(Qg/QCDLam)
      Tgg=Tg+Log(2.)
      Ts=Log(Qs/QCDLam)
      pi = 4.*atan(1.)

      Print *
      If (onel) then
        Print *,' Running alpha, Q-ordering => 1-loop'
      Else
        Print *,' Running alpha, Angle-ordering => all loops'
      Endif
      Print *
      Print *,' Qmax,Qg,Qs,QCDLam =',Qmax,Qg,Qs,QCDLam
      test=asbmy(mz)*pi/Ca
	write(6,*) ' alpha_s(m_z) = ',test
      Print *
      Print *,' Qmas,xmin =',Qmas,xmin
      Print *
      Jev=0
	Return
      End
	 
	 
C **********************************************************************

      SUBROUTINE GADAP2(A0,B0,FL,FU,F,EPS,SUM)
C
C   PURPOSE           - INTEGRATE A FUNCTION F(X,Y) OF TWO VARIABLES
C   METHOD            - ADAPTIVE GAUSSIAN IN BOTH DIRECTIONS
C   USAGE             - CALL GADAP2(A0,B0,FL,FU,F,EPS,SUM)
C   PARAMETERS  A0    - LOWER X-LIMIT (INPUT,REAL)
C               B0    - UPPER X-LIMIT (INPUT,REAL)
C               FL    - USER SUPPLIED FUNCTION FL(X) GIVING THE LOWER
C                       Y-LIMIT FOR A GIVEN X-VALUE
C                       (INPUT,REAL FUNCTION)
C               FU    - USER SUPPLIED FUNCTION FU(X) GIVING THE UPPER
C                       Y-LIMIT FOR A GIVEN X-VALUE
C                       (INPUT,REAL FUNCTION)
C               F     - USER SUPPLIED FUNCTION F(X,Y) TO BE INTEGRATED
C                       (INPUT,REAL FUNCTION)
C               EPS   - DESIRED ACCURACY (INPUT,REAL)
C               SUM   - CALCULATED VALUE FOR THE INTEGRAL (OUTPUT,REAL)
C   PRECISION         - SINGLE
C   REQ'D PROG'S      - FL,FU,F,GADAPF
C   AUTHOR            - THOMAS JOHANSSON, LDC,1973
C
	IMPLICIT NONE
	Integer NUM,IFU,I,L,N
      COMMON/GADAP1/ NUM,IFU
	Real F,FL,FU,GADAPF,A,B,F1,F2,F3,S,A0,B0,SS,SOLD
	Real DSUM,F1F,F2F,F3F,AA,BB,RED,SUM,EPS,AY,BY,W1,U2,C,X
      EXTERNAL F,FL,FU
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)

      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)
      IF(EPS.LT.1.0E-8) EPS=1.0E-8
      RED=1.4
      L=1
      I=1
      SUM=0.
      C=SQRT(15.)/5.
      A(1)=A0
      B(1)=B0
      X=0.5*(1+C)*A0+0.5*(1-C)*B0
      AY=FL(X)
      BY=FU(X)
      F1(1)=GADAPF(X,AY,BY,F,EPS)
      X=0.5*(A0+B0)
      AY=FL(X)
      BY=FU(X)
      F2(1)=GADAPF(X,AY,BY,F,EPS)
      X=0.5*(1-C)*A0+0.5*(1+C)*B0
      AY=FL(X)
      BY=FU(X)
      F3(1)=GADAPF(X,AY,BY,F,EPS)
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
      X=A(I)+B(I)-W1
      AY=FL(X)
      BY=FU(X)
      F1(I+1)=GADAPF(X,AY,BY,F,EPS)
      F2(I+1)=F3(I)
      X=B(I)-A(I+2)+W1
      AY=FL(X)
      BY=FU(X)
      F3(I+1)=GADAPF(X,AY,BY,F,EPS)
      X=U2
      AY=FL(X)
      BY=FU(X)
      F1(I+2)=GADAPF(X,AY,BY,F,EPS)
      F2(I+2)=F2(I)
      X=B(I+2)+A(I+2)-U2
      AY=FL(X)
      BY=FU(X)
      F3(I+2)=GADAPF(X,AY,BY,F,EPS)
      X=A(I)+A(I+2)-W1
      AY=FL(X)
      BY=FU(X)
      F1(I+3)=GADAPF(X,AY,BY,F,EPS)
      F2(I+3)=F1(I)
      X=W1
      AY=FL(X)
      BY=FU(X)
      F3(I+3)=GADAPF(X,AY,BY,F,EPS)
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

      FUNCTION GADAPF(X,A0,B0,F,EPS)
	IMPLICIT NONE
	Integer NUM,IFU,I,L,N
	Real F,A,B,F1,F2,F3,S,DSUM,F1F,F2F,F3F,AA,BB,EPS,RED
	Real SUM,C,A0,B0,X,W1,U2,SS,SOLD,GADAPF
      COMMON/GADAP1/ NUM,IFU
      EXTERNAL F
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)

      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)
      IF(EPS.LT.1.0E-8) EPS=1.0E-8
      RED=1.4
      L=1
      I=1
      SUM=0.
      C=SQRT(15.)/5.
      A(1)=A0
      B(1)=B0
      F1(1)=F(X,0.5*(1+C)*A0+0.5*(1-C)*B0)
      F2(1)=F(X,0.5*(A0+B0))
      F3(1)=F(X,0.5*(1-C)*A0+0.5*(1+C)*B0)
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
      F1(I+1)=F(X,A(I)+B(I)-W1)
      F2(I+1)=F3(I)
      F3(I+1)=F(X,B(I)-A(I+2)+W1)
      F1(I+2)=F(X,U2)
      F2(I+2)=F2(I)
      F3(I+2)=F(X,B(I+2)+A(I+2)-U2)
      F1(I+3)=F(X,A(I)+A(I+2)-W1)
      F2(I+3)=F1(I)
      F3(I+3)=F(X,W1)
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
   40 GADAPF=SUM
      EPS=EPS/RED
      RETURN
      END
	 
      FUNCTION SMRGEN(I)
	Implicit None
      Integer Len
      Parameter (Len=1)
      Real Rvec(Len)
      INTEGER I,j,ISEED,JSEED,K1,K2,LUX
      Double Precision Smrgen,SMrget,SMrset
	Real rndm
      Integer NTOTIN,NTO2IN,NTOTUT,NTO2UT
c
      CALL RANLUX(RVEC,LEN)
c        call ranmar(rvec,len)
	if(rvec(1).ne.rvec(1)) then
	  write(6,*) ' smrgen 1 ', (rvec(j),j=1,len),I
c          write(6,*) ' smrgen 2 ',rvec(2)       

c	  goto 123
	endif
      SMrgen=RVEC(1)
      RETURN
      END
       
C------------------------------------------------------------------------
      SUBROUTINE SMRAZM(PT,PX,PY)
C ... RANDOMLY ROTATED 2-VECTOR (PX,PY) OF LENGTH PT
      Double Precision Smrgen,SMrget,SMrset
      Double Precision RGEN,PT,PX,PY,C,S,CS,QT
   10 C=2.*SMRGEN(1)-1.
      S=2.*SMRGEN(2)-1.
      CS=C*C+S*S
      IF (CS.GT.1. .OR. CS.EQ.0) GO TO 10
      QT=PT/CS
      PX=(C*C-S*S)*QT
      PY=2.*C*S*QT
      END
C-----------------------------------------------------------------------
      SUBROUTINE SMUMAS(P)
C---PUTS MASS IN P(5) (-|MASS| IF SPACELIKE)
      Double Precision P(5),EM
      EM=(P(4)+P(3))*(P(4)-P(3))-P(1)**2-P(2)**2
      IF (EM.GE.0.) THEN
         P(5)=SQRT(EM)
      ELSE
         P(5)=-SQRT(-EM)
      ENDIF
      END
C------------------------------------------------------------------------
      SUBROUTINE SMVEQU(N,P,Q)
C---VECTOR EQUALITY
      INTEGER N,I
      Double Precision P(N),Q(N)
      DO 10 I=1,N
   10 Q(I)=P(I)
      END
       
C------------------------------------------------------------------------
      SUBROUTINE SMVDIF(N,P,Q,R)
C---VECTOR DIFFERENCE
      INTEGER N,I
      Double Precision P(N),Q(N),R(N)
      DO 10 I=1,N
   10 R(I)=P(I)-Q(I)
      END
        
      Subroutine SMinfn
C--Choose starting conditions for branching

      Include 'SMallx.inc'
      Integer Kid
      Common /SMbran1/Kid(0:1000)

      Logical pdflib,quark,gluon,photon,saturation
      Common /SMbran2/pdflib,quark,gluon,photon,saturation
      Integer iparton
      Common /SMquark/iparton
	Double Precision Xnorm
      common/smnorm/ Xnorm

      Double Precision SMrgen,SMrgau,u,yinf,frat,Qinf,Qsup
      Double Precision Rdq

	Double Precision Scale,test,xgmax,xnormn

      Double Precision  XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU,
     +VAL(20),QQ2,PP2
      Integer Ip2
      CHARACTER*20 PARM(20)
	
	Real slo,shi,xint,eps
	Integer i
	
	Double Precision xtmin,xtmax,dx,xt
	Integer nstep
	Integer nx
	Double Precision dlxmin,dlxmax,deltax,xmi,xma,sum
	Logical first,first1
	Data first/.true./,first1/.true./
	External xint
      Integer IpartonO
      Double Precision qgmin,qgmax,xmymin,WQG2,Wx0,xi,kt2i,ti,Qx2
      Double Precision xpq(-6:6),qg2max
      Double Precision kt,phi_k
      Integer nobran,ikincut
      common/myvar/nobran,ikincut
      integer iktscheme	
	Double precision wt0
      common/myweight/wt0
      Double Precision sigm, xmean
	Double precision alphas,sigm0,lam,x01,r02,xgx0
      
      
      Integer Itheory
      Common/theory/Itheory
	Real A1,A2,A3,A4,A5,A6
 
	double precision au 
      common/f2fit/au(50),first
      Double precision ininorm(-6:6)
      common/smininorm/ininorm
      
      Double Precision Ncall, sumkt, Wtest
      Double Precision pi
      pi = 4.*atan(1.)
     
     
      iktscheme =1

      Jev=Jev+1
c      Call SMrget(Js)
c      Write(6,*) ' sminit ininorm(0) ',ininorm(0),Xnor
      if(nobran.eq.2) then
         qgmin = 0.000001d0
c         qgmin = 1.
         qgmax = Sbig
c               QG2MAX = SSS*YX - Q2 - 4.D0 * AM(1)**2
         Qx2 = QGMIN*((QGMAX/QGMIN)**SMrgen(1))
         WQG2 = Qx2*LOG(QGMAX/QGMIN)
c     write(6,*) qx2,sbig
C         xmymin=max(0.95e-5,xmin) !old
         xmymin=max(1e-7,xmin)

         x0 = xmymin*(xmax/xmymin)**SMrgen(2)
         Wx0 = x0*log(xmax/xmymin)
c--choose starting Qt
         PHI_K = 2*pi*Smrgen(3)
         P(1,0)= sqrt(qx2)*cos(phi_k)
         P(2,0)= sqrt(qx2)*sin(phi_k)
c      P(1,0)= SMrgau(1,0.0d0,Qx)
c      P(2,0)= SMrgau(2,0.0d0,Qx)
        Wt=1.
      else

c--choose starting Qt
        if(iktscheme.eq.2) then
          qgmin = 0.0001
          qgmax = Sbig
          qx2 = qgmin*((qgmax/qgmin)**SMrgen(1))
c--choose starting Qt
          PHI_K = 2*pi*Smrgen(3)
          P(1,0)= sqrt(qx2)*cos(phi_k)
          P(2,0)= sqrt(qx2)*sin(phi_k)
c          write(6,*) ' sminfn_gluon ',iktscheme,qx2,P(1,0),P(2,0)
        else 
          qx2 = qx
          sigm = Qx
          xmean = 0.0
          P(1,0)= SMrgau(1,xmean,sigm)
          P(2,0)= SMrgau(2,xmean,sigm)
        endif
	  x0=0.99
        Wt=1.
        if(itheory.eq.105) then
          if(saturation) then
             qgmin = 0.000001d0
             qgmin = 0.00010001d0
             qgmax = Sbig
             qgmax = 100.
             Qx2 = QGMIN*((QGMAX/QGMIN)**SMrgen(1))
             WQG2 = Qx2*LOG(QGMAX/QGMIN)
c--choose starting Qt
             PHI_K = 2*pi*Smrgen(3)
             P(1,0)= sqrt(qx2)*cos(phi_k)
             P(2,0)= sqrt(qx2)*sin(phi_k)
	       A1 = au(1)
	       Sigm0 = au(2)
             LAM = au(3)
	       X01 = au(4)
             ALPHAS = au(5)
	       if(first) then
c               SIGM0 = 29.12*2.568
c	         LAM = 0.277
c               X01 = 0.41/10000.
c               ALPHAS = 0.2
	         write(6,*) ' sminfn: saturation selected ',saturation
	         write(6,*) ' sminfn: GBW: xgx0 = 3/4./PI^2/ALPHAS*SIGM0*R02*SCALE*EXP(-R02*SCALE) '
               write(6,*) ' sminfn: sigm0 = ',sigm0, ' alphas = ',alphas
               write(6,*) ' sminfn; r02 = (X0/X01)**LAM; with x01 = ',x01,' lam = ',lam
	         first=.false.
	       endif

             xmymin=max(1e-7,xmin)  
              
             x0 = xmymin*(xmax/xmymin)**SMrgen(4)
             Wx0 = x0*log(xmax/xmymin)
             R02 = (X0/X01)**LAM
   	       Scale=P(1,0)**2+P(2,0)**2
             glu = 3.*SIGM0*R02*SCALE*EXP(-R02*SCALE)
             xgx0 = glu/4./9.8696/ALPHAS
c	       xgmax=12.
c	       xgmax=1200.
c	       xgmax=10000.
c             test = xgx0/x0*WQG2*Wx0/xgmax
             Wt0 = xgx0/x0*Wx0*WQG2
	       Wt0 = Wt0*(1-x0)**4.2
c	       If (test.gt.1)then
c                write(6,*) ' fatal in saturation: xgx0 > xgmax :',
c     &        wt0,xgmax
c             endif
ccc	     Wt = 1.
             Wt = Wt0
c             write(6,*) sigm0,R02,scale,alphas,x0,x01,lam
cc             write(6,*)  wt,wt0,x0,scale,xgx0,Wx0,WQG2
ccc              SMtest=SMrgen(5)
c	      write(6,*) ' xmaxmin ',x0,xmax,xgx0,sigm0,r02,scale,alphas
c	      write(6,*) test,SMtest
ccc             If (test.lt.SMtest) go to 9
c             write(6,*) ' event acepted ',x0             
          else    
             xmymin=max(1e-7,xmin)   
             x0 = xmymin*(xmax/xmymin)**SMrgen(2)
             Wx0 = x0*log(xmax/xmymin)
	       A1 = au(1)
	       A2 = au(2)
             A3 = au(3)
	       A4 = au(4)
             A5 = au(5)
             A6 = au(6)
	       if(first) then
	         write(6,*) ' sminfn: initial paramterisation selected '
	         write(6,*) ' sminfn: A1*5.*(1./x0**A2)*(1.-x0)**A3*(1.-A4*x0+A5*sqrt(x0)+A6*x0**2)/x0'
               write(6,*) ' sminfn: A1 = ',A1,' A2 = ',A2,' A3 = ',A3
               write(6,*) ' sminfn: A4 = ',A4,' A5 = ',A5,' A6 = ',A6
	         first=.false.
             endif 
             Wt0 =A1*5.*(1./x0**A2)*(1.-x0)**A3*(1.-A4*x0+A5*sqrt(x0)+A6*x0**2)/x0 * Wx0
	       Wt = Wt0
c             write(6,*) wt,x0,xmax,xmymin
          endif
        endif
	  Scale=P(1,0)**2+P(2,0)**2
c        write(6,*)' scale ',scale
c	  write(6,*) quark,gluon,photon

c	write(6,*) gluon,quark,photon

	  If(iparton.eq.0) then
	    Kid(0) = 21
          if(itheory.eq.105) then
            Wtest = Wtest + Wt*x0
            ininorm(0) = Xnor*Wtest/Jev
            ininorm(3) = Xnor*Wtest/Jev
          else
            ininorm(0) = Xnor ! old
            ininorm(3) = Xnor ! old
          endif
	  else
	    Kid(0) = iparton
          if(ipartonO.ne.iparton) first1=.true.
   10     continue
          x0=xmin**SMrgen(3)
c here do the starting distribution for the valence quarks
          PARM(1) = 'DEFAULT' 
ccc          VAL(1) = 60701 ! HERAPDF15(NLO)  is not anymore in latest lhapdf
cc          VAL(1) = 60500 ! HERAPDF10(NLO) 
          VAL(1) = 10551 ! CTEQ66
	    If(first1) then
              ipartonO=iparton
              xgmax = 0.
              CALL PDFSET(PARM,VAL)
c calculate the integral over xgx	

              QQ=Qs
	        QQ2=QQ*QQ
	        xnormn = 0.
     	        nx=100
	        xmi=0.0001
	        xma=0.95
	        dlxmin=log10(xmi)
	        dlxmax=log10(xma)
	        dx = (dlxmax-dlxmin)/real(nx)
              do i=0,nx
	          xt = 10**(dlxmin+dx*real(i))
                XX=DBLE(Xt)
                CALL STRUCTM(XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
c integrate quarks
                if(iparton.eq.1) then 
                  sum = DNV
                  elseif(iparton.eq.2) then
                  sum = UPV
                  else
                  write(6,*) ' sminfn: no valid iparton specified',iparton
                  write(6,*) 'sminfn: PROGRAM STOPPED '
                  stop
                endif
                if(sum.eq.sum) then  
	            if(i.lt.nx) then
	              deltax =10**(dlxmin+dx*real(i+1))-10**(dlxmin+dx*real(i))
                    xnormn = xnormn + real(sum)*deltax
	            endif
c	            write(6,*) ' deltax,glu ',deltax,glu,QQ2,PP2
                    test=sum*1.5
c                    write(6,*) 'xx, sum ',xx,sum
                    if(test.gt.xgmax) xgmax = test
	          endif 
              enddo
              write(6,*) '  sminfn Xnorm  new :',xnormn
		  Xnor =  xnormn  ! really respect mom-sum rule for valence quarks 
              ininorm(iparton)=Xnor
c              write(6,*) '  sminfn Xnorm new :',xnormn
	        first1 = .false.
	    endif  
          QQ=Qs
	    QQ2=QQ*QQ
	    
	    xx=x0	  

          CALL STRUCTM(XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
          if(iparton.eq.1) then 
              sum = DNV
              elseif(iparton.eq.2) then
              sum = UPV
              else
              write(6,*) ' sminfn: no valid iparton specified'
     &              ,iparton
              write(6,*) 'sminfn: PROGRAM STOPPED '
              stop
          endif
	    if(sum.ne.sum) goto 10  

  	    test = sum/xgmax
	    IF(test.ge.1) 
     &      write(6,*) ' sminfn: test > 1 ',test,xx,qq,sum,xgmax,iparton
c wt = f(x)/g(x) with f(x) = q(x) and g(x) = 1/x
c    = xq(x)  	   
          If (test.lt.SMrgen(4)) go to 10
          
	  Endif
	Endif  
c--compute xgbar
      xg=x0
      xgbar=-(P(1,0)**2+P(2,0)**2+p2(5)**2*x0*(2.-x0))
     +        /(Sbig*(1.-x0))
	xgbar =0
      P(3,0)=x0*p2(3)+xgbar*p1(3)
      P(4,0)=x0*p2(4)+xgbar*p1(4)
      Call SMumas(P(1,0))
c--proton spectator
      Call SMvdif(4,p2,P(1,0),pspec)
ctest............	
	yphton = 1
	yphbar = 0.
	Wini = 1.
cend.......... test.........	
	
      q1(3)=yphton*p1(3)+yphbar*p2(3)
      q1(4)=yphton*p1(4)+yphbar*p2(4)
c--electron spectator
      Call SMvdif(4,p1,q1,espec)
      q1(5)=-Qval
      P(6,0)=Wt
      P(7,0)=0.
      P(8,0)=0.
c	write(6,*) ' sminfn_gluon :, x0 = ',x0
      Return
      End
C------------------------------------------------------------------------
      FUNCTION SMRGAU(J,A,B)
C ... Gaussian random number, mean A, standard deviation B.
      INTEGER J,I
      Double Precision A,B,X,R,SMRGAU
      Double Precision Smrgen,SMrget,SMrset
      X=0.
      DO 1 I=1,12
         R=SMRGEN(I)
    1 X=X+R
      SMRGAU=A+B*(X-6.)
      END

	
      Subroutine updfgrid(iflag)
c obtain the unitegrated gluon density as a function of xg,kt**2
c in a grid

      Implicit none
      Integer       Nmongo,Nbin,Is,Js,Imax,Nsoft,Ifin,Jfin,Jev,
     & Nev,Maxpr,Ichg,Idq
      Logical       onel,tevol,hspot,powin,afix
      Double Precision xmin,xmax,ymin,ymax,Qmin,Qmax,Xini,Xnor,
     & Xfin,Xfi2,Sigtot,Sigto2,Sigto3,Mulsum,Mulsu2,Mulsu3,xgbar,
     & x0,xg,Wt,Wsig,Qval,yphton,yphbar,Sbig,Gluch,QQt,Qtem,Pmax,
     & Qg,Qs,Qx,Ts,QCDLam,Ttem,Tg,Tgg,Told,Capb,znow,wg,Excode,Wini,
     & pow0,pow1,xmp0,q0,a0,Xmjtrg,Ppjtrg,Ptjtrg,Thjtrg,umax,xymin,
     & emax,Qmas,Qbar,Qhard,me,p1,p2,q1,q2,p3,p4,pspec,espec,P
      Common/SMicom/Nmongo,Nbin,Is,Js,Imax,Nsoft,Ifin,Jfin,Jev,
     & Nev,Maxpr,Ichg,Idq
      Common/SMlcom/onel,tevol,hspot,powin,afix
      Common/SMrcom/xmin,xmax,ymin,ymax,Qmin,Qmax,Xini,Xnor,
     & Xfin,Xfi2,Sigtot,Sigto2,Sigto3,Mulsum,Mulsu2,Mulsu3,xgbar,
     & x0,xg,Wt,Wsig,Qval,yphton,yphbar,Sbig,Gluch,QQt,Qtem,Pmax,
     & Qg,Qs,Qx,Ts,QCDLam,Ttem,Tg,Tgg,Told,Capb,znow,wg,Excode,Wini,
     & pow0,pow1,xmp0,q0,a0,Xmjtrg,Ppjtrg,Ptjtrg,Thjtrg,umax,xymin,
     & emax,Qmas,Qbar,Qhard,me,p1(5),p2(5),q1(5),q2(5),p3(5),p4(5),
     & pspec(5),espec(5),P(8,0:1000)
C--P(*,I) is for the Ith branching:
C  P(1,I)=Qtx_i
C  P(2,I)=Qty_i
C  P(3,I)=Ql_i
C  P(4,I)=E_i
C  P(5,I)=M_i (+/- for time/spacelike)
C
C  P(6,I)=current full weight
C  P(7,I)=0. for spacelike, +/- parent for int./ext. timelike
C  P(8,I)=emitted transverse momentum qt_i


      Integer NMXHEP,NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
      Double Precision PHEP,VHEP,EVWGT
      PARAMETER (NMXHEP=4000)
      COMMON/HEPEVTP/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     + JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
C---Event weight
C   N.B. weights must be renormalized to momentum sum
      COMMON/HEPWGT/EVWGT
	Double Precision BB
	Common /splitting/ BB
	Integer Ipgg,ns_sel
	Double Precision scal
	common/Pggsel/Ipgg,ns_sel,scal
      Double Precision tmax,p_max
      Common/mygrid/tmax,p_max
      Integer nobran,ikincut
      common/myvar/nobran,ikincut
	Double Precision Qbarmy,Qbar_min,Qbar_max
      Common/mglubran/Qbarmy,Qbar_min,Qbar_max
      Integer neve
	common/myevt/neve
      Double Precision sum
c	Integer nloop
c	common/myloop/nloop
	Double precision wt0
	common/myweight/wt0
c my variables
      Double Precision wt0my
      Double Precision xp2,xkt2,x_g,wtxg
c      common/cwnco1/xp2,xkt2,x_g,wtxg

      Integer Nbx,Nbk,Nbp
      Parameter (Nbx=50)
c      Parameter (Nbk=42)
      Parameter (Nbk=50)
c      Parameter (Nbk=1)
c      Parameter (Nbp=20)
      Parameter (Nbp=50)
      Double Precision X1(0:Nbx+2),X2(0:Nbk+2),X3(0:Nbp+2)
      Double Precision X1M(0:Nbx+1),X2M(0:Nbk+1),X3M(0:Nbp+1)
      Double Precision x1b(0:Nbx+1),x2b(0:Nbk+1),x3b(0:Nbp+1)
      Double Precision xgx(-4:4,0:Nbx,0:Nbk,0:Nbp)
	Double Precision xgx_f(-4:4,0:Nbx,0:Nbk,0:Nbp)
      Double Precision xm11(0:Nbx),xm21(0:Nbk),xm31(0:Nbp)
      Double Precision xm12(0:Nbx),xm22(0:Nbk),xm32(0:Nbp)
      Double Precision xm1(0:Nbx),xm2(0:Nbk),xm3(0:Nbp)
      Double Precision x3norm(-4:4,0:Nbp),x3normf(-4:4,0:Nbp),x3normN(-4:4,0:Nbp)
      Double Precision Val1,Val2,Val3
      Integer I,J,K,Iflag,Irr,Err
      Integer ix1,ix2,ix3
      double Precision nval2,Nrejx1,Nrejx2
      Double Precision x1minL,x1maxL,x1dif,x2dif,x3dif
      Double Precision x1minlin,x1maxlin,xcheck
      Integer NbxL,Nbxlin
      Double Precision x2minL,x2maxL,q2min,q2max,p3max,p3min
      Double Precision Xgnorm,Xgn,xgmy,wtmy,Xnormtest,momsum
c limits
      Double Precision x1min,x1max,x2min,x2max,x3min,x3max
      Double Precision x1dm,x2dm,x3dm
	Double Precision x3lmin,x3lmax,x3ldif
	Double Precision Qgmin,Xin
      CHARACTER   FILNAME*132
      Common/gludatf/filname
      Character*5 xxxx
      Integer Kid
      Common /SMbran1/Kid(0:1000)
	Integer ki
	Integer nold
	Double Precision kt2
      LOGICAL ex
	data nold/1/
      Data xxxx/'gluon'/
	Integer iprint
	Data iprint/0/
      Double precision ininorm(-6:6)
      common/smininorm/ininorm
      Double precision ncount0, ncount3
      
c      write(6,*) ' updfgrid ',iflag
	
      If(iflag.eq.1) Then
	   do ki=-4,4
	     do k=0,nbp
	       x3norm(ki,k) = 0.
	       x3normf(ki,k) = 0.
	       x3normN(ki,k) = 0.
	     enddo
	   enddo

c initilaze grid
cThis is the new stuff
c xg
c         x1min=0.00001
	   x1min = xmin
cccc         x1min=0.000001 ! new (not really needed for TEVATRON)
         xcheck = 0.1
         x1max=xmax
c log grid up to 0.1
         x1max = xcheck
         x1minL=log(x1min)
         x1maxL=log(x1max)
         Nbxlin = 10
         NbxL = Nbx-Nbxlin
         x1dif=(x1maxL-x1minL)/Real(NbxL)
         Do I=0,NbxL
            x1(I) = x1minL + x1dif*Real(I)
         Enddo
c	   Do I=0,NbxL
c            x1m(I) = (x1(i) + x1(i+1))/2.
c            x1b(i) = x1(i+1) - x1(i)
cc		write(6,*) ' log(x) bins x',i,exp(x1(i)),exp(x1b(i))
c         Enddo
c lin grid in x from 0.1 - 1
         x1minlin=xcheck
         x1maxlin=1.0d0
         x1dif=(x1maxlin-x1minlin)/Real(Nbxlin)
         Do I=NbxL,NbxL+Nbxlin+1
            x1(I) = log(x1minlin + x1dif*Real(I-NbxL))
         Enddo


	   Do I=0,NbxL+Nbxlin+1
            x1m(I) = (x1(i) + x1(i+1))/2.
            x1b(i) = x1(i+1) - x1(i)
c		write(6,*) ' log(x) bins x',i,exp(x1(i)),exp(x1b(i))
		write(6,*) '  bins x',i,exp(x1(i)),x1(i),x1b(i)
         Enddo


c parameters for updfgrid:
         q2min = 0.001 ! new
         q2max = qmax**2/4.
c
	   x2minL=log(q2min)
	   x2maxL=log(q2max)
	   x2dif=(x2maxL-x2minL)/real(nbk)
         Do I=0,Nbk+1
            x2(I) = x2minL + x2dif*Real(I)
	   Enddo
	   Do I=0,Nbk
            x2m(I) = (x2(i) + x2(i+1))/2.
            x2b(i) = x2(i+1) - x2(i)
c		write(6,*) ' log(kt2) bins ',i,x2m(i),exp(x2m(i))
         Enddo
c p_max
         Qgmin = max(Qg,QCDlam)
         p3min=Qgmin
	   p3max=Qmax
	   x3lmin=log(p3min)
	   x3lmax=log(p3max)
c	   write(6,*) ' p3min,p3max ',p3min,p3max,x3lmin,x3lmax
	   x3ldif = (x3lmax-x3lmin)/Real(Nbp)	
	   write(6,*) ' updfgrid-test: xl3dif ',x3ldif
         Do I=0,Nbp+1
            x3(I) = x3lmin + x3ldif*Real(I)
c	    write(6,*) ' updfgrid-test: Qbar bins = ',exp(x3(i)),I
         Enddo
	   Do I=0,Nbp
	      x3m(i) = (x3(i) + x3(i+1))/2.
		x3b(i) = x3(i+1) - x3(i)
c		write(6,*) ' log(p) bins ',i,x3b(i),x3m(i)
         Enddo
c	   Do I=0,Nbp
c	      x3(i) = (exp(x3(i)) + exp(x3(i+1)))/2.
c	    write(6,*) ' updfgrid-test2: Qbar bins = ',x3(i),I
c	      x3(i) = log(x3(i))
c         Enddo
	
         Qbar_min=x3(0)
	   Qbar_max=x3(Nbp)
c zero xg
         Do ki=1,2
         Do I=0,Nbx
            Do J=0,Nbk
               Do K=0,Nbp
                  xgx(ki,I,J,K)=0.
                  xgx_f(ki,I,J,K)=0.
               Enddo
            Enddo
         Enddo
         Enddo
         Xgnorm = 0.
         Xgn = 0.
      Elseif(iflag.eq.2) Then
c fill in values
         xgmy = xg
         wtmy = wt
         Val1 = log(xgmy)
	
         kt2=P(1,Ifin)**2+P(2,Ifin)**2
         Val2 = log(kt2)
         x_g=xgmy
         wtxg=wtmy*xgmy
c new test
         xp2 = p_max	
         Val3 = log(xp2)
         If(Val3.lt.x3(0)) Then
            write(6,*) ' val3 = ',val3,' too low: x3(0) = ',x3(0),p_max
            Val3 = x3(0)
         elseif(val3.gt.x3(Nbp+1)) Then
cc            write(6,*) ' val3 = ',val3,' too high: x3(Nbp) = ',x3(Nbp)
cc            Val3 = x3(Nbp+1)
               return
         endif
	   if(kid(Ifin).eq.21) then
	      ki=0
            ncount0 = ncount0 + 1
	      else
	      ki=kid(ifin)
            if(ki.eq.3) ncount3 = ncount3 + 1
	   Endif
         
         ix3=-1
   40    ix3 = ix3 + 1
         if(ix3.gt.Nbp) then
            write(6,*) ' ix3 outside grid ',ix3,Nbp
         else
            if(val3.gt.x3(ix3+1)) then 
c		write(6,*) ' updfgrid ',ix3,x3(ix3),val3
            goto 40
            endif
         Endif

         if(wt0.gt.0) then
            wt0my = wt0
          else
            wt0my=1.
         endif
	   x3normN(ki,ix3) = x3normN(ki,ix3) + 1.
c norm to initial condition (old style)
ccc	 x3norm(ki,ix3) = x3norm(ki,ix3) + x0*wt0my  
c norm to flavor sum rule at each scale 
	    x3normf(ki,ix3) = x3normf(ki,ix3) + wtmy
c norm to momentum sum rule at each scale 
          x3norm(ki,ix3) = x3norm(ki,ix3) + xgmy*wtmy
	      
	   Xin = Xin + x0
	   if (Excode.ne.0) then
c            write(6,*) ' updfgrid return: Excode = ',Excode 
            return
         endif
         If (Wini.eq.0.0) then
            write(6,*) ' updfgrid return: Wini = 0 ',Wini
            return
         endif
	   
	   
         If(Val1.lt.x1(0)) Then
            write(6,*) ' val1 = ',val1,' too low: x1(0) = ',x1(0)
            Nrejx1 = Nrejx1 + 1
            Val1 = x1(0)
            return
         elseif(val1.gt.x1(Nbx+1)) Then
            write(6,*) ' val1 = ',val1,' too high: x1(Nbx) = ',x1(Nbx+1)
            Val1 = x1(Nbx+1)
		return
         endif
 	   
         ix1=-1
   20    ix1 = ix1 + 1
         if(ix1.gt.Nbx) then
            write(6,*) ' ix1 outside grid ',ix1,Nbx
         else
            if(val1.gt.x1(ix1+1)) goto 20
         Endif
         If(Val2.lt.x2(0)) Then
c            write(6,*) ' val2 = ',val2,' too low: x2(0) = ',x2(0)
            Val2 = x2(0)
		return
         elseif(val2.gt.x2(Nbk+1)) Then
            nval2 = nval2 + 1
            if(nval2.lt.10) then
c               write(6,*) ' val2 = ',val2,' too high: x2(Nbk) = ',x2(Nbk+1)
            endif
            Val2 = x2(Nbk+1)
		Nrejx2 = Nrejx2 + 1
            return
         endif
         ix2=-1
   30    ix2 = ix2 + 1
         if(ix2.gt.Nbk) then
            write(6,*) ' ix2 outside grid ',ix2,Nbk
         else
            if(val2.gt.x2(ix2+1)) goto 30
         Endif
         

	   xgx(ki,ix1,ix2,ix3)= xgx(ki,ix1,ix2,ix3) + wtmy

         Xgn = Xgn + wtmy*xgmy

      Elseif(iflag.eq.3) Then
c xgx defined here is actually g(x)
c corresponds: d g(x)/dlnx = x dg(x)/dx
         if(iprint.eq.1) then
           write(6,*) Xgn,x1dif,Xgnorm,x1maxL,x1minL,Nbx
c           write(6,*) 'ncall=',ncall(ki)
	     write(6,*) ' rej x1 cut ',Nrejx1,' rej x2 cut ',Nrejx2
c	     write(6,*)  ' xgn/ncall =',xgn/real(ncall(1))
	     write(6,*) ' updfgrid: X3norm = ',(X3norm(1,i),i=0,Nbp)
	   endif
c
         Do ki = -4,4
         Xnormtest = 0.
         write(6,*) ' ininorm :',ininorm(ki),ki
         Do I = 0,Nbx
            Do J = 0,Nbk
               sum = 0
               Do K = 0, Nbp
                 Xnormtest = max(x3norm(ki,k),1.)
                 momsum = x3norm(0,k)/max(x3normN(0,k),1.d0) + x3norm(3,k)/max(x3normN(3,k),1.d0)
                 if(ki.eq.1) then
                    Xgnorm = 1./max(x3normf(ki,k),1.d0)
                 elseif(ki.eq.2) then
                    Xgnorm = 2./max(x3normf(ki,k),1.d0)                   
                 elseif(abs(ki).eq.3) then
         
c                    Xgnorm = ininorm(ki)/(0*max(x3normN(0,k),1.d0)+max(x3normN(3,k),1.d0)) 
                     Xgnorm = ininorm(ki)/Jev/momsum          
		     else
c                    Xgnorm = ininorm(ki)/Xnormtest
c                    Xgnorm = ininorm(ki)/(max(x3normN(0,k),1.d0)+0*max(x3normN(3,k),1.d0))
                     Xgnorm = ininorm(ki)/Jev/momsum          
     		     endif
c                 if(xgnorm.gt.0) then
c                   write(6,*) 'xgnorm ',xgnorm,ki,i,j,k
c                 endif
                 xgx_f(ki,I,J,K) = xgx(ki,I,J,K)/x2b(j)/x1b(i)*Xgnorm
               Enddo
            Enddo
         Enddo
         Enddo
c print out table
c         write(6,*) ' *********************************************** '
c         write(6,*) ' ******        xA(x,kt**2,p)              ****** '
c         write(6,*) ' *********************************************** '
c         write(6,*) '|   ln(xg)  |  ln(k_t^2)  |   ln(p)  |   xA     |'
c         write(6,*) ' *********************************************** '
c         write(30,*) ' Qg0 = ',Qg,' ikincut= ',ikincut
c         write(30,*) ' ln(xg),  ln(k_t^2), ln(p),  xgx'
         inquire(FILE=filname,EXIST=ex)
         write(6,*) ' updfgrid-test: filename ',filname,ex
         if(ex) then
            open(30,FILE=filname, FORM='formatted',STATUS=
     +      'old',IOSTAT=IRR,ERR=50 )
            else
            open(30,FILE=filname, FORM='formatted',STATUS=
     +      'new',IOSTAT=IRR,ERR=50 )
	   Endif
	   write(30,10000) Qg,ikincut
10000    format(' Qg0 = ',f12.8,' ikincut= ',I6)
         write(30,10001) Ipgg,ns_sel
10001    format(' Ipgg = ',I4,' ns_sel = ',I4)
         write(30,10002) QCDlam
10002    format(' QCDlam = ',f12.8)
         
         write(30,10100)
10100    format(' ln(xg),  ln(k_t^2), ln(p),  xgx')

c write output to file

         do j=0,Nbk
            do i=0,Nbx
               do k=0,Nbp
                  write(30,*) x1(i),x2(j),x3(k),(xgx_f(ki,I,J,K),ki=-4,4)
               enddo
            enddo
         enddo
         close(30)
        
         If(iprint.eq.1) then
            write(6,*) nbx,nbk,nbp
            write(6,*) ' updfgrid: after before close '
	   endif
      Endif
      Return
   50 write(6,*) ' updfgrid error opening file',filname
      Return
  100     format(' |',E12.4,' | ',F8.5,' | ',F8.5,' | ',F8.5,' |')
      End

      Subroutine SMqtem(Tlike)
C--Finds Qtem
      Include 'SMallx.inc'
      Double Precision U0,Ul,R,D,C
      Integer n
	
	
      Real sudint,suda,fl,fu,eps
      Double Precision q2hi,q2lo,qs0,Ca,nf,Lambda
      Common/sudvar/q2hi,q2lo,qs0,Ca,nf,Lambda
	External suda,fl,fu
C Parameters for Sudakov form factors
C (NMXSUD= max no of entries in lookup table)
	Real ACCUR,QEV,SUD
	Integer INTER,NQEV,NSUD
	Integer NMXSUD
      PARAMETER (NMXSUD=1024)
      COMMON/SMUSUD/ACCUR,QEV(NMXSUD,6),SUD(NMXSUD,6),INTER,NQEV,NSUD
	Real Fdiv(NMXSUD),Adiv(NMXSUD)
      Integer Isud
	Common/SMSUD1/Isud
	Double Precision Qtemgg,Ttemgg,Qtemqq,Ttemqq,Qtemqg,Ttemqg,Qtemgq,Ttemgq 
	Common/SMSUD2/ Qtemgg,Ttemgg,Qtemqq,Ttemqq,Qtemqg,Ttemqg,Qtemgq,Ttemgq 
      Logical pdflib,quark,gluon,photon,saturation
      Common /SMbran2/pdflib,quark,gluon,photon,saturation
	Real sudlo,sudup
	Double Precision power,qfac,qnow
	Real Qtest,Stest
	Real Qlim,Sudmin(4),Sudmax(4)
	 
	Real DIVDIF
	External DIVDIF
	 
	Integer IQ,ISS,I,L,L1,L2,LL
	Logical First
      Common /SMstart/ First
	Integer Iseed,minuts
	Double Precision SMrgen
	External minuts
	Common/random/Iseed
	Double Precision ana,t0,t,tmin
	Integer Ismax
	 
CCC	 
      Logical Tlike
      Data Eps/1.E-5/
      If(Tlike) then
	write(6,*) ' this version of smqtem is not valid for tlike '
	stop
      endif
c HERWIG like parameter setting	 
      NQEV=1024
      NSUD = 0
      INTER=3
      ACCUR=1.D-6
      qlim = sqrt(Sbig)
c 	 
      nf=4. 
ccc      Lambda=0.2
      Lambda=QCDlam
      Ca=3.
      eps = ACCUR
	 
      Qs0 = Qg
      Ismax = 4
	 
      IF(first) Then
c	  write(6,*) ' integrate Sudakov and produce lookup table'
         POWER=1./Real(NQEV-1)
         q2lo=Qs0**2
         QFAC=(1.1*qlim/sqrt(q2lo))**POWER
C--IS=1 FOR GLUON->GLUON+GLUON FORM FACTOR
C--IS=2 FOR quark->gluon+quark FORM FACTOR
C--IS=3 FOR gluon->quark+antiquark FORM FACTOR
C--IS=4 FOR quark->quark+GLUON FORM FACTOR
         DO 11 Is = 1,Ismax
	      ISS=IS
            SUD(1,IS)=0.
            QEV(1,IS)=q2lo
	      Sudmin(is)=99999999.
	      Sudmax(is)=-9999999.
	      Isud = Is
            DO 10 IQ=2,NQEV
               QNOW=QFAC*sqrt(QEV(IQ-1,IS))
	       if(qnow.ge.qmax) then
c	         write(6,*) ' smqtem qnow > qmax ',qnow,qmax
	       endif
	         q2hi = Qnow**2
     	         call gadap2(0.,1.,fl,fu,suda,eps,sudint)
               SUD(IQ,IS)=sudint
	         if(sudint.lt.sudmin(is)) sudmin(is)=sudint
	         if(sudint.gt.sudmax(is)) sudmax(is)=sudint
c               write(6,*) ' smqtem IQ,IS ',IQ,IS,q2hi
   10       QEV(IQ,IS)=q2hi
        write(6,*) ' new Sudakov ala Gavin '
        write(6,201) IS,NQEV
  201   FORMAT(1H1//10X,'EXTRACT FROM TABLE OF SUDAKOV FORM FACTOR NO.',
     &  I2,' (',I5,' ACTUAL ENTRIES)'//10X,'SUD IS PROBABILITY THAT',
     &  ' PARTON WITH GIVEN UPPER LIMIT ON Q WILL REACH THRESHOLD',
     &  ' WITHOUT BRANCHING'///2X,8('      Q     SUD ')/)
             L2=NQEV/8
             L1=L2/32
             IF (L1.LT.1) L1=1
             DO 40 L=L1,L2,L1
                LL=L+7*L2
                WRITE(6,30) (sqrt(QEV(I,IS)),exp(-SUD(I,IS)),I=L,LL,L2)
   30           FORMAT(2X,8(F9.2,F7.4))
   40        CONTINUE
   50         FORMAT(1H1)
	 
   11    Continue	 
	 first=.false.
      endif
 
      q2lo = QCDlam*exp(Told)
      q2lo = q2lo**2
      Qtest=q2lo
      R=SMrgen(1)
      Do Is=1,Ismax
         Qtest=q2lo
c here we ask sud(qtest)         
         sudlo=DIVDIF(SUD(1,IS),QEV(1,IS),NQEV,Qtest,INTER)	 
         If (R.eq.0.) R=Eps
            D=-Log(R)	 
	      stest = D + sudlo
            if(stest.gt.Sudmax(is)) then
c	  write(6,*) ' smqtem : stest,sudmin',stest,sudmin(is),sudmax(is)
	         stest=min(stest,sudmax(is))
	      endif 
            if(stest.lt.Sudmin(is)) then
c	  write(6,*) ' smqtem : stest,sudmin',stest,sudmin(is),sudmax(is)
	         stest=max(stest,sudmin(is))
	    endif 

	    if(sudlo.gt.stest) then
	       write(6,*) ' sudlo ',sudlo,stest,sudmin(is),sudmax(is)
	    endif
c here we ask QEV(sud)	
          Qtest = DIVDIF(QEV(1,IS),SUD(1,IS),NQEV,stest,INTER)
          
c          write(6,*) ' smqtem Qtest ',IS,sqrt(Qtest),R
	
          Qtem=sqrt(Max(0.,Qtest))
	    Qtem=Max(Qg,Qtem)
          Ttem = log(Qtem/QCDLam)

	    if(Ttem.lt.Told) then
	       write(6,*) ' smqtem_new problem: Ttem < Told ',Ttem,Told,Is
	       write(6,*) ' smqtem_new problem: Qtem < Qold ',Qtem,q2lo,Is
	    Endif
c          write(6,*) ' smqtem ',qtem,is
          If(IS.eq.1) then
             Qtem = sqrt(Max(0.,Qtest))
	       Qtem = Max(Qg,Qtem)
             Ttem = log(Qtem/QCDLam)
	
	       if(qtem.eq.0.) then
	          write(6,*) ' smqtem problem Qtem:',Qtem,Qtest,stest,D,
     &	     sudlo,q2lo,R,sudmin(is),sudmax(is)
                 Qtem=Qg+Eps
                 Ttem = log(Qtem/QCDLam)
	       endif
	       if(qtem.le.Qg) then
	          write(6,*) ' smqtem problem Qtem:',Qtem,Qg
                 Qtem=Qg+Eps
                 Ttem = log(Qtem/QCDLam)
	       endif
 
	       if(Qtem.ne.Qtem) then 
	          write(6,*) ' smqtem problem Qtem:',Qtem,Qtest,stest,
     &          D,sudlo,q2lo
	       endif
	       Qtemgg = Qtem
             Ttemgg = Ttem
          Elseif(IS.eq.2) then
             Qtemgq = sqrt(Max(0.,Qtest))
	       Qtemgq = Max(Qg,Qtemgq)
             Ttemgq = log(Qtemgq/QCDLam)
	
	       if(qtemgq.eq.0.) then
	           write(6,*) ' smqtem problem Qtem:',Qtemgq,Qtest,stest,D,
     &	           sudlo,q2lo,R,sudmin(is),sudmax(is)
                 Qtemgq=Qg+Eps
                 Ttemgq = log(Qtemgq/QCDLam)
	       endif
	       if(qtemgq.le.Qg) then
	           write(6,*) ' smqtem problem Qtem:',Qtem,Qg
                 Qtemgq=Qg+Eps
                 Ttemgq = log(Qtemgq/QCDLam)
	       endif
          Elseif(IS.eq.3) then
             Qtemqq=sqrt(Max(0.,Qtest))
	       Qtemqq=Max(Qg,Qtemqq)
             Ttemqq = log(Qtemqq/QCDLam)
	
	       if(qtemqq.eq.0.) then
	          write(6,*) ' smqtem problem Qtem:',Qtemqq,Qtest,stest,D,
     &	          sudlo,q2lo,R,sudmin(is),sudmax(is)
                Qtemqq=Qg+Eps
                Ttemqq = log(Qtemqq/QCDLam)
	       endif
	       if(qtemqq.le.Qg) then
	            write(6,*) ' smqtem problem Qtem:',Qtem,Qg
                  Qtemqq=Qg+Eps
                  Ttemqq = log(Qtemqq/QCDLam)
	       endif
          Elseif(IS.eq.4) then
             Qtemqg=sqrt(Max(0.,Qtest))
	       Qtemqg=Max(Qg,Qtemqg)
             Ttemqg = log(Qtemqg/QCDLam)
         	
	       if(qtemqg.eq.0.) then
	          write(6,*) ' smqtem problem Qtem:',Qtemqg,Qtest,stest,D,
     &	          sudlo,q2lo,R,sudmin(is),sudmax(is)
                Qtemqg=Qg+Eps
                Ttemqg = log(Qtemqg/QCDLam)
	       endif
	       if(qtemqg.le.Qg) then
	          write(6,*) ' smqtem problem Qtem:',Qtem,Qg
                Qtemqg=Qg+Eps
                Ttemqg = log(Qtemqg/QCDLam)
	       endif
	    endif
	Enddo
      End
	
      Function suda(x,y)
      Implicit None
      Real x,y
      Double Precision q2hi,q2lo,q0,Ca,nf,Lambda
      Common/sudvar/q2hi,q2lo,q0,Ca,nf,Lambda
      Real Suda
      Double Precision q,q2,z1,z1min,z1max,z
      Double Precision Capb,asb,Lambda2,Sudakov
      Double Precision Splitt,Splittgq,Splittqg,Splittqq,asbmy
      External Splitt,Splittgq,Splittqg,Splittqq,asbmy
      Double Precision BB
      Common /splitting/ BB
      Double Precision BBsave
      Integer Ipgg,ns_sel
      Double Precision scal 
      common/Pggsel/Ipgg,ns_sel,scal
      Integer Isud
      Common/SMSUD1/Isud	
      Double Precision Psplit
      Double Precision Cf
      Double Precision epsilon
      
      epsilon=1E-15
	
      Cf = 4d0/3d0
      BBsave = BB
C for the sudakov, set B always to 0.5
      BB = 0.5	 
c	 write(6,*) ' smqtem-all B =',BB,BBsave
      q2=q2lo*(q2hi/q2lo)**x
      q=sqrt(q2)
c z1= 1-z	
c z1min= = Q0/q, z1max=1
      z1min=q0/q
      z1max=1.-Epsilon
      z1= z1min*(z1max/z1min)**y
      sudakov=1.
      Capb=12.*Ca/(11.*Ca-2.*nf)
      z = 1.-z1
C--Isud=1 FOR GLUON->GLUON+GLUON FORM FACTOR
      if(Isud.eq.1) then 
	   Psplit=Splitt(z1)
      elseif(Isud.eq.2) then 
	   Psplit=Splittgq(z)*Cf/2./Ca
c	   Psplit=Splitt(z1)
      elseif(Isud.eq.3) then 
         Psplit=Splittqq(z)*Cf/2./Ca 
c	   Psplit=Splitt(z1)
      elseif(Isud.eq.4) then 
	   Psplit=Splittqg(z)*0.5/2./Ca
c	   Psplit=Splitt(z1)
      else
         write(6,*)' Isud = ',Isud,' not implemented '
      endif
      BB=BBsave
c asbmy = C_A/pi * alpha_s
      sudakov = asbmy(scal*q*z1) *Psplit/q2	 
	 
c	 
      sudakov=sudakov*q2*log(q2hi/q2lo)*z1*log(z1max/z1min)
      suda=sudakov
      return
      end
      Function Splittgq(z)
C--IS=2 FOR quark->gluon+quark FORM FACTOR
      Implicit None
      Double Precision Splittgq,z
      Splittgq = (1.+ (1.-z)**2)/z
      Return
      End 
	
      Function Splittqq(z)
C--IS=3 FOR quark->quark+GLUON FORM FACTOR
      Implicit None

      Double Precision Splittqq,z
      Splittqq = (1.+ z**2)/(1.-z)
      
      Return
      End 
	
      Function Splittqg(z)
C--IS=4 FOR gluon->quark+GLUON FORM FACTOR
      Implicit None

      Double Precision Splittqg,z
      Splittqg = z**2 + (1.-z)**2
      Return
      End 
	 
	 
      Function fl(x)
      real fl,x
      fl=0.
      return
      end
      function fu(x)
      real fu,x
      fu=1. 
      return
      end
C
C smqtem and szvalnew for full splitting function
C copy from fullsplitt.f:     
c
      Subroutine Szvalnew
C
C new version 20.5.2012
C
C--Chooses momentum fraction z for small-x branching
C  Cleaned-up version 24/10/91
C  All-loop has norm'n of small z part with no zmin
C  according to Kwiecinski Martin Sutton
C
C  changed hannes
C  All-loop has no non-sudakov normalisation
C  also alphas changed, now including cutoff for freezing alphas=0.6
C
C now also including non-singular terms according to Gavins proposal
C
C------------------------------------------------------------------------
      Include 'SMallx.inc'

      Common/gcom/g0,g1,gtot
      Double Precision asb,g0,g1,TQt,zmax,F,Ff,gtot,S2asb,arg
	Double Precision SMfreq,SMgsin,SMrgen
      Double Precision Fmin,Fmax
      Double Precision TQtp,Qlam0,xmin1
      Double Precision Pgtot,Pg1,Pg2
      Data Fmin,Fmax/0.0000001,0.9999999/
      Data Qlam0/0.71/
	Integer nsidx
	common /myns/nsidx
	Integer nf
      Double Precision Lambda,QQg	
	Common/hardme/Lambda,QQg,nf
	Real slo,shi,eps,g1p,g1f,g0p,g0f
	Double Precision g1_sudak,result
	Double Precision Splitt,z1,Splitt0,zz
	Double Precision gtotw
      Double Precision wz
      Common/smweights/wz
	Logical mytest
	Data mytest/.false./
	External g1p,Splitt
	External g0p,Splitt0
	
	Double Precision D_ns,asbold,asbold1
	Double Precision Qtemp,Capbp,QCDlamp
	Common/g1int/Qtemp,Capbp,QCDlamp,zmax
	Double Precision Qtem0,QQt0,xmin0,zmax0,Qg0
	Common/g0int/Qtem0,QQt0,xmin0,asbold,zmax0,Qg0
	Double Precision Nonsud0,Nonsud
	External Nonsud0,Nonsud
CCCCCCCCCCCCCCCCCCCCCCCCCC
C for Sudakov (smqtem.f) B is fixed: B=0.5	
C choose BB: 0, 0.5 , 1.	
	Double Precision BB
	Common /splitting/ BB
	Integer Ipgg,ns_sel
	Double Precision scal 
	common/Pggsel/Ipgg,ns_sel,scal

	Double Precision BBsave
	Double Precision asbmy
	External asbmy
	
	Double Precision alp_dns,as
	Integer Ismallz
	Common /smallz/ Ismallz,alp_dns,as
	
	Double Precision zmin
      Double Precision Dfreq,DGAUSN
      Double Precision zero
      
      Data zero/0d0/
		
	Logical First
	Data First/.True./
	Double Precision Cf,Ca
	Ca = 3d0
	Cf = 4d0/3d0

c this is alphas/(2pi) ala bryan
      TQt=Log(QQt/QCDLam)
c new try with min(QQt,Qlam0)	
	Qlam0 = Qg
      TQtp=Log((dmax1(QQt,Qlam0))**2/QCDLam**2)
      asb= Capb/TQtp

	If(First) Then
      BBsave = BB
	if(Ipgg.eq.1.or.Ipgg.eq.3) then 
	write(6,*) ' using full splitting function Pgg:',Ipgg
c	write(6,*) ' P = asb[(1-z)/z + (1-B)z(1-z)] Delta_ns +',
c     &' asb[z/(1-z) + Bz(1-z)] Delta_s'
	write(6,*) ' with B = ',BBsave
	else
	write(6,*) ' using sing part of splitting function Pgg:',Ipgg
	endif
	If(ns_sel.eq.-1) then
        write(6,*) ' Deltans is old standard one '
	elseif(ns_sel.eq.0) then
	  write(6,*) ' Deltans: as(k), sing parts, cut '
	elseif(ns_sel.eq.1) then
	  write(6,*) ' Deltans: as(k), full splitting, cut '
	elseif(ns_sel.eq.2) then
	  write(6,*) ' Deltans: as(q), sing parts, cut '
	elseif(ns_sel.eq.3) then
	  write(6,*) ' Deltans: as(q), full splitting, cut '
	endif
	write(6,*) ' new version: cutoff for asb = ',Qlam0
	write(6,*) ' asb =C_A/pi alpha_s ~ alpha_s '
	write(6,*) ' max asb = ',Capb/Log(Qlam0**2/QCDLam**2)
	First=.false.
	endif
CCCCCCCCCCCCCCCCCCCCCCCCCC
      BB = BBsave
c      write(6,*) ' szvalnew BB = ',BB
      Ismallz = 0
	alp_dns=1.	
	Wg=1.
	Wz=-99999.
	D_ns = 1.
C First get proper alphas
c asb = 2 Ca alphas/(2pi)
      asbold = asb    ! hannes new
      asbold1 = 0.5*Capb/Ttem
      if(Ipgg.ge.2) then
	  as= asbmy(scal*Qtem)
	else
	  as= asbmy(scal*QQt)
	endif
      alp_dns = as 	
	Qtemp=Qtem
	Qtem0=Qtem
	QQt0=QQt
	xmin0=xmin
	Capbp=Capb
	QCDlamp=QCDlam
      zmax=(1.-Qg/Qtem)
	zmax0 = zmax
	Qg0 = Qg
      g1 = 0.5*Capb*Log(Ttem/Tg)
      If (onel) then
	   zmin = xmin/xg
	   xmin0 = zmin
         If (zmax.gt.zmin) g0=asb*log(zmax/zmin)
         gtot=g0+g1
         If (SMrgen(1)*gtot.le.g1) then
            znow=1.-exp(Ttem*((Tg/Ttem)**SMrgen(2)-1.))
            z1=1.-znow            
	      asbold1= asbmy(Qtem*(1.-znow)) ! New 
	      wz = gtot/(asbold1/z1) ! hannes new    
	      gtotw = asbmy(scal*Qtem*(1.-znow))*Splitt(z1)*wz
         Else
            Nsoft=Nsoft+1
            znow=zmin*(zmax/zmin)**SMrgen(3)
	      wz = gtot/(asbold/znow)
	      gtotw = as*Splitt0(znow)*wz 
         Endif
      Else
         S2asb=Sqrt(2.*asb)
         arg=S2asb*(Ttem-TQt+Log(zmax))
         F=DFREQ(arg) ! changed because of numeric problems in SMFREQ
         if(QQt.ge.Qtem) then
            g0 = asb/S2asb*sqrt(2.*3.1415) *F*exp(asb*(Ttem-TQt)**2)
ctestc kt>zmax*qt
         elseif(QQt.ge.zmax*Qtem) then
            g0 = asb/S2asb*sqrt(2.*3.1415)*F
         elseif(QQt.lt.zmax*Qtem) then
            Pg1=asb*sqrt(2.*3.1415)/S2asb*DFREQ(zero)
            Pg2=asb*log(zmax*Qtem/QQt)
            g0 = Pg1+Pg2
         else
            write(6,*) ' this should never happen ',QQt,Qtem,zmax
         endif
c
         if(g0.ne.g0) then
           write(6,*) 'problem in  smzvalgg', g0
         endif
         gtot=g0+g1
c do the z weighting	
         BB = BBsave
         If (SMrgen(1)*gtot.le.g1) then
            znow=1.-exp(Ttem*((Tg/Ttem)**SMrgen(2)-1.))
            if(znow.gt.zmax) write(6,*) ' 1. smzvalgg: znow>zmax ',znow,zmax
		nsidx = 0
            z1=1.-znow            
channes	      asbold1= asbmy(Qtem*(1.-znow)) ! New 

	      wz = gtot/(asbold1/z1) ! hannes new    
	      gtotw = asbmy(scal*Qtem*(1.-znow))*Splitt(z1)*wz
         Else
            if(QQt.lt.zmax*Qtem) then
               If (SMrgen(4)*g0.gt.Pg1) then
                  xmin1 = QQt/Qtem
                  znow=xmin1*(zmax/xmin1)**SMrgen(5)
            if(znow.gt.zmax) write(6,*) ' 2. smzvalgg: znow>zmax ',znow,zmax
			nsidx = 3
               else
                  Nsoft=Nsoft+1
                  Ff=SMrgen(3)*DFREQ(zero)
                  If (Ff.lt.Fmin) then
                     znow=xmin*0.5
                     if(znow.gt.zmax) write(6,*) ' 3. smzvalgg: znow>zmax ',znow,zmax
                  Elseif (Ff.gt.Fmax) then
                     znow=zmax*Fmax
                     if(znow.gt.zmax) write(6,*) ' 4. smzvalgg: znow>zmax ',znow,zmax
                  Else
                     znow=QQt/Qtem*exp(DGAUSN(Ff)/S2asb)
                     if(znow.gt.zmax) then
                        write(6,*) ' 5. smzvalgg: znow>zmax ',znow,zmax,qqt,qtem,S2asb,asb
                        write(6,*) ' 5. smzvalgg: g0, g1 ',g0,g1,Ff,DGAUSN(Ff)        
                     endif
                  Endif
			nsidx =2
               Endif
            else
               Nsoft=Nsoft+1
               Ff=SMrgen(3)*F
               If (Ff.lt.Fmin) then
                  znow=xmin*0.5
               Elseif (Ff.gt.Fmax) then
                  znow=zmax*Fmax
               Else
                  znow=QQt/Qtem*exp(DGAUSN(Ff)/S2asb)
               Endif
		   nsidx = 1
            Endif
            if(znow.gt.zmax) then
              write(6,*) ' smzvalgg: z> zmax ',znow,' zmax ',zmax
              D_ns = 0
	  	  alp_dns = as * D_ns
              else
c Here apply the weight for the small z generation
		  Ismallz = 1
              zz = znow
              D_ns=	Nonsud(zz)
		  alp_dns = as * D_ns
	        wz = gtot/(asbold/znow*Nonsud0(zz))
	        gtotw = as*Splitt0(zz)*D_ns*wz 
            Endif 	
         Endif
      Endif	  

	if(znow.ne.znow) then
	  write(6,*) ' problem in szvalnew g0,g1:',g0,g1
        write(6,*) ' QQt,Qtem,zmax ',QQt,Qtem,zmax
	  endif
	  

	  
c calculate here g1_sudakov, for proper norm	
	slo=0.
	shi=1.
	eps=0.0001
	BB=0.5
	call gadap(slo,shi,g1p,eps,g1f)
	g1_sudak=g1f
c do the z weighting	
      BB = BBsave
      zz=znow
      if(gtotw.le.0) then
        write(6,*) ' smzvalgg: gtotw < 0 ',gtotw,g0,g1,wz
      endif
      wg = gtotw  /g1_sudak
      if(wg.le.0) then
        write(6,*) ' smzvalgg: wg < 0 ',gtotw,g1_sudak
       endif
c end new	
	if(wg.ne.wg) then
	  write(6,*) ' problem in szvalnew wg:',wg,znow
	  write(6,*) ' problem in szvalnew g0,g1:',g0,g1
        write(6,*) ' QQt,Qtem,zmax ',QQt,Qtem,zmax
	  wg=0.
	endif
      
      end
      

	function asbmy(xkt)
c  asmy = 2 Ca alphas/(2pi)      
      Include 'SMallx.inc'
      Double Precision TQtp,xkt,Qlam0
	Double Precision asbmy
	Double Precision xtest
	integer nf 
	Double Precision mass(6)
c	data mass/0.,0.,0.,1.5,5.0,175./
	data mass/0.,0.,0.,1.35,4.8,175./
      Integer Ialps
c      data Ialps/1/
      data Ialps/2/
      Integer nfmin,nfmax
      Double Precision ALAM2,Q2EFF,Q2THR,pi,B0,B1,alphas,algq
	Double Precision as
	Double Precision Cf,Ca
	Logical asbfirst
	Data asbfirst/.true./
	Ca = 3d0
	Cf = 4d0/3d0
c      Data Qlam0/0.71/
	Qlam0 = Qg	
      TQtp=Log((dmax1(xkt,Qlam0))**2/QCDLam**2)
c new approach ala shirkov
      xtest=xkt/QCDLam
      if(xtest.lt.0.9999.or.xtest.gt.1.0001) then
	else 
c	  write(6,*) ' asbmy: xtest = ',xtest,xkt,QCDLam
	Endif
c old version
c asb = Ca *alphas/pi
      asbmy = Capb/TQtp
	pi=4.*atan(1.)
      q2eff = xkt**2     
      nfmin = 3
      nfmax = 5
      nf = 4
      alam2 = qcdlam**2
c alphas from PYTHIA      
C...Evaluate first or second order alpha_strong.
  100 IF(NF.GT.MAX(2,nfmin)) THEN
        Q2THR=mass(NF)**2
        IF(Q2EFF.LT.Q2THR) THEN
          NF=NF-1
c          ALAM2=ALAM2*(Q2THR/ALAM2)**(2D0/(33D0-2D0*NF))
c          write(6,*) nf
          if(nf.eq.3) then
             ALAM2=ALAM2*(Q2THR/ALAM2)**(2./27.)
           if(Ialps.eq.2) ALAM2=ALAM2*(log(Q2THR/ALAM2))**(107./2025)
          endif
          GOTO 100
        ENDIF
      ENDIF
  110 IF(NF.LT.MIN(8,nfmax)) THEN
        Q2THR=mass(NF+1)**2
        IF(Q2EFF.GT.Q2THR) THEN
          NF=NF+1
          if(nf.eq.5) then 
             ALAM2=ALAM2*(ALAM2/Q2THR)**(2./23.)
             if(Ialps.eq.2) ALAM2=ALAM2*(log(Q2THR/ALAM2)**(-2.*963./13225.))
          endif
          if(nf.eq.6) then 
             ALAM2=ALAM2*(ALAM2/Q2THR)**(2./21.)
             if(Ialps.eq.2) ALAM2=ALAM2*(log(Q2THR/ALAM2)**(-2.*321./3381.)) 
          endif
         GOTO 110
        ENDIF
      ENDIF
	if(nf.gt.6.or.nf.le.0) write(6,*) nf,q2

      B0=(33D0-2D0*NF)/6D0
      ALGQ=LOG(MAX(1.0001D0,Q2EFF/ALAM2))
      IF(Ialps.EQ.1) THEN
        alphas=MIN(10.0,2.*pi/(B0*ALGQ))
      ELSE
        B1=(153D0-19D0*NF)/6D0
        alphas=MIN(10.0,2.*pi/(B0*ALGQ)*(1.-B1*LOG(ALGQ)/
     &  (B0**2*ALGQ)))
      ENDIF

      asbmy = 2.*Ca*alphas/2./pi ! hannes new
      if(asbfirst) then
         asbfirst=.false.
         write(6,*) ' asbmy: using ',ialps,' order of alphas'
c        asbmy = alphas
c        write(6,*) ' new alphas ',alphas,qcdlam,b0,algq,q2eff
      endif
	return
	end
	
      function g1p(y)
	
	Implicit None
	Real g1p,y
	Double Precision g1t,z1,z1min,z1max
	Double Precision zmax,asb
	Double Precision Qtem,Capb,QCDlam
	Common/g1int/Qtem,Capb,QCDlam,zmax
	Integer Ipgg,ns_sel
	Double Precision scal 
	common/Pggsel/Ipgg,ns_sel,scal

	Double Precision splitt,asbmy
	External splitt,asbmy
	Double Precision Cf,Ca
	Ca = 3d0
	Cf = 4d0/3d0
c z1= 1-z	
c z1min= = Q0/q, z1max=1
       z1min=1.-zmax
	 z1max=1.
	 z1= z1min*(z1max/z1min)**y
       g1t = asbmy(scal*qtem*z1) *Splitt(z1)	 
	 g1t=g1t*z1*log(z1max/z1min)
	 g1p=g1t
	 return
	 end
	 
      function g0p(y)
	
	Implicit None
	Real g0p,y
	Double Precision g0t,z1,z1min,z1max,asbold
	Double Precision Qtem,QQt,xmin,zmax,Qg
	Common/g0int/Qtem,QQt,xmin,asbold,zmax,Qg
	Integer Ipgg,ns_sel
	Double Precision scal 
	common/Pggsel/Ipgg,ns_sel,scal

	Double Precision Splitt0,Nonsud,asbmy
	External Splitt0,Nonsud,asbmy
	Double Precision as
	Double Precision Cf,Ca
	Ca = 3d0
	Cf = 4d0/3d0
      z1min=xmin
	
	
	
	z1max=zmax
	z1= z1min*(z1max/z1min)**y
	as=asbmy(scal*QQt)
      if(Ipgg.ge.2) then
	as= asbmy(scal*Qtem)
	else
	as=asbmy(scal*QQt)
	endif
      g0t = as *Splitt0(z1)*Nonsud(z1)	
	g0t=g0t*z1*log(z1max/z1min)
	g0p=g0t
	
	return
	end
	Function Splitt(z1)
	Implicit None
	Integer Ipgg,ns_sel
	Double Precision scal
	common/Pggsel/Ipgg,ns_sel,scal

	Double Precision Splitt,z1
	Double Precision BB
	Common /splitting/ BB
	If(Ipgg.eq.0.or.Ipgg.eq.2) then
c first the simple 1/(1-z)=1/z1	 
	  Splitt=1./z1
	  else
c now gavins version
       Splitt = ((1.-z1)/z1 +BB*(1.-z1)*z1)
	endif
	Return
	End 
	 
	
      Function Splitt0(z)
	Implicit None
	Integer Ipgg,ns_sel
	Double Precision scal
	common/Pggsel/Ipgg,ns_sel,scal
	Double Precision Splitt0,z,z1
	Double Precision BB
	Common /splitting/ BB
	z1 = z
c first the simple 1/z	 
	If(Ipgg.eq.0.or.Ipgg.eq.2) then
         Splitt0=1./z1
	   else
c now gavins version
         Splitt0 = (1.-z1)/z1 + (1.-BB)*z1*(1.-z1)
	 endif
	 Return
	 End 
	 
      Function Nonsud(zz)
	Implicit None
      Logical       onel,tevol,hspot,powin,afix
      Common/SMlcom/onel,tevol,hspot,powin,afix
	
	Double Precision Nonsud,z,zz
	Double Precision Nonsud1,z1,ktqt,asbold,z0,z0c,z01
	Double PRecision Nonsud2
	Double Precision BB
	Common /splitting/ BB
	Double Precision Qtem,QQt,xmin,zmax,Qg
	Common/g0int/Qtem,QQt,xmin,asbold,zmax,Qg
	Double Precision Qtemp,Capbp,QCDlamp,zmax1
	Common/g1int/Qtemp,Capbp,QCDlamp,zmax1
	Double Precision Q0	
	Double Precision term1,term2,term3,term4,term5,term6,term7
	Double Precision t0,B,k,q
	Double Precision tq,tk,tc,zc,zt
	Integer ncall
	Integer Ipgg,ns_sel
	Double Precision scal 
	common/Pggsel/Ipgg,ns_sel,scal
	Double Precision t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
	Double Precision t11,t12,t13,t14,t15,t16,t17,t18,t19
	Double Precision t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
	Double Precision t30,t31,t32,t33,t34,t35,t36,t37,t38,t39
	Double Precision t40,t42,t41,t43,t45,t46,t47,t48
	Double Precision t50,t51,t52,t53,t55,t56,t58,t59
	Double Precision t62,t63,t65,t66,t68
	Double Precision t70,t71,t74,t75,t76
	Double Precision t80,t82,t85,t92,t94,t95
	Double Precision t102,t107,t108
	Double Precision a0,b0,c0,d0
	Double Precision one
	Double Precision Ei
	External Ei

c ns_sel = -1: standard non sudakov old
c ns_sel = 0: standard non sudakov with cut
c ns_sel = 1: full splitting
c ns_sel = 2: alphas(q) sing. terms
c ns_sel = 3: alphas(q) full splitting
      Integer Insq
c     Insq = 0 non-sudakov with alphasq with hannes treatment
c	Insq = 1 non-sudakov with alphasq with gavin treatment
      Data insq/0/
	data ncall/0/
	Double Precision asbmy
	External asbmy
	Logical first
	Data first/.true./
	if(first) then 
	   write(6,*) ' non-sudakov selection : ns_sel = ',ns_sel
	   write(6,*) ' Insq = ',Insq
	   first=.false.
	endif
	if(onel) then
	  Nonsud = 1.
	  else 
	
	z1=zz
	ktqt=QQt/Qtem
	if(ktqt.ge.1.) then
	  z0= 1.
	Elseif(z1.lt.ktqt.and.ktqt.le.1.) then
	  z0=ktqt
	elseif(ktqt.le.z1) then
	  z0=z1
	else
	  write(6,*) ' Nonsud : this should never happen QQt,Qtem:',QQt,Qtem
	Endif
	if(z1.gt.z0) then
	   write(6,*) ' Nonsud z1,z0',z1,z0,qtem,qqt
	endif
      k=QQt
      q=Qtem
      Q0=Qg  
      Q0=0.2
      z=z1
      tk = log(k**2/QCDlamp**2)
	tq = log(q**2/QCDlamp**2)
	tc = log(Q0**2/QCDlamp**2)
	zc = Q0/q
	zt = max(z1,zc)
c now gavins version
	B=BB
	if(ns_sel.eq.-1) then
         Nonsud1 =  log(z0/z1)*log(k**2/z0/z1/q**2)
         Nonsud1 =  -asbmy(scal*k)*Nonsud1
	elseif(ns_sel.eq.0) then
c the simple non - sudakov
	   if(tk.lt.tc) then
	     Nonsud1 = 0.
	   else
           Nonsud1 =  log(z0/zt)*log(k**2/z0/zt/q**2)
	   if(zc.gt.z1) then      
              Nonsud1 = Nonsud1+(log(zc/z1))*log(k**2/Q0**2)
	   endif
         Nonsud1 =  -asbmy(scal*k)*Nonsud1
	   endif  
	
	elseif(ns_sel.eq.1) then
	   if(tk.lt.tc) then
	     t0 = 0.
	   else           
c
c
c here comes the new part
           t0 = -log(k/z0/q)**2+log(k/z0/q)*z0**2+z0**2/2-
     #          B*log(k/z0/q)*z0**2-
     #          z0**2*B/2-2*z0*log(k/z0/q)-
     #          2*z0-2.D0/3.D0*z0**3*log(k/z0/q)-2.D0/9.D0*z0**3+
     #          2.D0/3.D0*B*z0**3*log(k/z0/q)+2.D0/9.D0*z0**3*B+
     #          log(k/zt/q)**2-log(k/zt/q)*zt**2-
     #          zt**2/2+B*log(k/zt/q)*zt**2+zt**2*B/2+2*zt*log(k/zt/q)+
     #          2*zt+2.D0/3.D0*zt**3*log(k/zt/q)+2.D0/9.D0*zt**3-
     #          2.D0/3.D0*B*zt**3*log(k/zt/q)-2.D0/9.D0*zt**3*B
           
            If(Insq.eq.0) then
	        if(zc.gt.z1) then      
                t1 = -log(k/zc/q)**2+log(k/zc/q)*zc**2+zc**2/2-
     #          B*log(k/zc/q)*zc**2-
     #          zc**2*B/2-2*zc*log(k/zc/q)-
     #          2*zc-2.D0/3.D0*zc**3*log(k/zc/q)-2.D0/9.D0*zc**3+
     #          2.D0/3.D0*B*zc**3*log(k/zc/q)+2.D0/9.D0*zc**3*B+
     #          log(k/z1/q)**2-log(k/z1/q)*z1**2-
     #          z1**2/2+B*log(k/z1/q)*z1**2+z1**2*B/2+2*z1*log(k/z1/q)+
     #          2*z1+2.D0/3.D0*z1**3*log(k/z1/q)+2.D0/9.D0*z1**3-
     #          2.D0/3.D0*B*z1**3*log(k/z1/q)-2.D0/9.D0*z1**3*B
                t0 = t0 + t1
	        endif
     
 	      Elseif(Insq.eq.1) then
c now gavins version ( for full splitting fct )
	          if(zc.gt.z1) then      
                   t1 = 2*log(k/q0)*log(zc)-2*log(k/q0)*zc-
     #	       2.D0/3.D0*log(k/q0)*zc**3+
     #             log(k/q0)*zc**2+2.D0/3.D0*log(k/q0)*B*zc**3-
     #             log(k/q0)*B*zc**2-2*log(k/q0)*log(z1)+
     #             2*log(k/q0)*z1+2.D0/3.D0*log(k/q0)*z1**3-
     #             log(k/q0)*z1**2-2.D0/3.D0*log(k/q0)*B*z1**3+
     #             log(k/q0)*B*z1**2
                   t0 = t0 + t1
	          endif
	      endif
	   endif
         Nonsud1 =  -asbmy(scal*k)*t0
	elseif(ns_sel.eq.2) then
c the simple non - sudakov

c here non sudakov with alphas(q) 
	   if(tk.lt.tc) then
	      Nonsud1 = 0.
	   else
            Nonsud1 = 0.5*log(tk/(tq+2.*log(z0)))*tq+
     &          log(tk/(tq+2*log(z0)))*log(z0)+0.5*tq+log(z0) -
     &         (0.5*log(tk/(tq+2.*log(zt)))*tq+
     &          log(tk/(tq+2*log(zt)))*log(zt)+0.5*tq+log(zt))
            If(Insq.eq.0) then
               Nonsud1 =  Nonsud1 * Capbp 
	         if(zc.gt.z1) then      
   	            z01=zc
	            ktqt=QQt/Qtem
                  if(ktqt.ge.1.) then
	              z0c= 1.
	              z0c=z01
	            Elseif(z01.lt.ktqt.and.ktqt.le.1.) then
	              z0c=min(ktqt,z01)
	            elseif(ktqt.le.z01) then
	              z0c=z01
	            else
	              write(6,*) ' Nonsud : this should never happen '
	            Endif 
	            Nonsud1 = Nonsud1 + asbmy(q0)*
     &	       log(z0c/z1)*log(k**2/z0c/z1/q**2)
	         endif
 	      Elseif(Insq.eq.1) then
	         if(zc.gt.z1) then      
                 Nonsud1 = Nonsud1+(log(zc)-log(z1))*log(tk/tc)
	         endif
               Nonsud1 = Nonsud1 * Capbp  
	      endif
            if(Nonsud1.ne.Nonsud1) then
	        write(6,*) ' problem with new Non Sudakov '
	        write(6,*) ' z0,z1 ',z0,z1,k,q
	        write(6,*) ' log(tk/(tq+2*log(z0)))',
     #	     log(tk/(tq+2.*log(z0)))
	        write(6,*) ' log(tk/(tq+2*log(z1)))',
     #           log(tk/(tq+2.*log(z1)))
	        write(6,*) ' tq ',tq,' tk ',tk,log(z0),log(z1)
	        write(6,*) '  Non - Suadkov: ',Nonsud1
            endif
	   endif
         Nonsud1 = -Nonsud1  
	elseif(ns_sel.eq.3) then
c here non sudakov with fullspiltting and alphas(q) 
	   if(tk.lt.tc) then
	     t0 = 0.
	     Nonsud1 = 0.
	   else
           one = 1.d0
           t1 = log(tk)
           t3 = tq/2
           t4 = log(zt)
           t5 = t3+t4
           t6 = log(t5)
           t7 = zt**2
           t8 = t7*zt
           t9 = t6*t8
           t12 = t6*t7
           t14 = log(z0)
           t16 = 3./2.*tq
           t17 = exp(-t16)
           a0=-3*t4-t16
           t20 = Ei(one,a0)
           t24 = z0**2
           t25 = t24*z0
           t28 = 2*tq
           t31 = log(t28+4*t14)
	     b0=-3*t14-t16
           t36 = Ei(one,b0)
           t39 = exp(-t3)
           t40 = -t14-t3
           t41 = Ei(one,t40)
           t45 = log(-t40)
           t47 = exp(-tq)
	     c0=-2*t14-tq
           t50 = Ei(one,c0)
           t56 = t45*t25
           t58 = -z0*t1-t9/3-t6*zt-t4+t12/2+t1*t14+t14-t17*t20/3+
     #          t6*t4-t25*t1/3
     #          -tq*t31/2+t17*t36/3+t39*t41+t24*t1/2+t45*z0-t47*t50/2+
     #          zt*t1-t7*t1/2+t56/3
           t62 = log(t28+4*t4)
           t65 = t45*t24
           t70 = B*t17
           t80 = B*t47
	     d0=-2*t4-tq
           t92 = Ei(one,d0)
           t95 = Ei(one,-t5)
           t107 = -t1*t4+tq*t62/2-t65/2+t8*t1/3-t45*t14+t70*t20/3-
     #           t24*B*t1/2-
     #           t56*B/3-t12*B/2+t80*t50/2+t25*B*t1/3-t70*t36/3+t9*B/3+
     #           t47*t92/2-t39
     #           *t95+t7*B*t1/2+t65*B/2-t8*B*t1/3-t80*t92/2
           t108 = t58+t107

           t0 = t108
           Nonsud1 = -t0 * Capbp 
	     t0 = 0
           If(Insq.eq.0) then
	        if(zc.gt.z1) then      
   	           z01=zc
	           ktqt=QQt/Qtem
                 if(ktqt.ge.1.) then
	              z0c= 1.
	              z0c=z01
	           Elseif(z01.lt.ktqt.and.ktqt.le.1.) then
	              z0c=min(ktqt,z01)
	           elseif(ktqt.le.z01) then
	              z0c=z01
	           else
	             write(6,*) ' Nonsud : this should never happen '
	        Endif 
           
              t0 = -log(k/z0c/q)**2+log(k/z0c/q)*z0c**2+z0c**2/2-B*
     #	   log(k/z0c/q)*
     #         z0c**2-z0c**2*B/2-2*log(k/z0c/q)*z0c-2*z0c-2.D0/3.D0*
     #         z0c**3*log(k/z0c/q)
     #         -2.D0/9.D0*z0c**3+2.D0/3.D0*B*z0c**3*log(k/z0c/q)+
     #         2.D0/9.D0*z0c**3*B+log(k/z1/q)**2-log(k/z1/q)*z1**2-
     #         z1**2/2+B*log(k/z1/q)*z1**2+z1**2*B/2
     #         +2*log(k/z1/q)*z1+2*z1+2.D0/3.D0*z1**3*log(k/z1/q)+
     #         2.D0/9.D0*z1**3
     #         -2.D0/3.D0*B*z1**3*log(k/z1/q)-2.D0/9.D0*z1**3*B
              t0 = - asbmy(q0)* t0 
              if(t0.gt.0) write(6,*) 'exp(t0) > 0 ',t0
	        Nonsud1 = Nonsud1+ t0
              if(Nonsud1.gt.0) write(6,*) 'exp(Nonsud1) > 0 ',Nonsud1
		  endif
 	     Elseif(Insq.eq.1) then
	        if(zc.gt.z1) then      
                 t0 = 2*log(k/q0)*log(zc)-2*log(k/q0)*zc-
     #	     2.D0/3.D0*log(k/q0)*zc**3+
     #           log(k/q0)*zc**2+2.D0/3.D0*log(k/q0)*B*zc**3-
     #           log(k/q0)*B*zc**2-2*log(k/q0)*log(z1)+
     #           2*log(k/q0)*z1+2.D0/3.D0*log(k/q0)*z1**3-
     #           log(k/q0)*z1**2-2.D0/3.D0*log(k/q0)*B*z1**3+
     #           log(k/q0)*B*z1**2
                 t0=-asbmy(scal*k)*t0
                 Nonsud1 = Nonsud1 + t0 
		     
	        endif
	     Endif

           if(Nonsud1.ne.Nonsud1) then
	        write(6,*) ' problem with new Non Sudakov '
	        write(6,*) ' z0,z1 ',z0,z1,k,q
	        write(6,*) '  Non - Suadkov: ',Nonsud1,t0,t108
           endif
	   endif
	else 
	  write(6,*) '  Non - Suadkov: ',ns_sel,' not yet implemented '
        stop
      endif


      Nonsud1 = exp(Nonsud1)
     
     
	If(Nonsud1.lt.0.0.or.Nonsud1.gt.1.0001) Then
	  ncall=ncall+1
	  If(ncall.le.1000) then
	  write(6,*) ' Problem with new Non - Suadkov: ',Nonsud1,ns_sel
	  endif 
	  Nonsud1 = 1.
	Endif
      Nonsud = Nonsud1
	endif
	Return
	End 
      Function Ei(x1,x2)
	Implicit None
	Double Precision x1,x2,Dexpin,Ei
	Integer i
	External Dexpin
	if(x1.eq.1.) then
	  Ei=dexpin(x2)
	  else
	  write(6,*) ' Ei: i= ',x1,' not implemented '
	  Endif
	Return
	End 
      Function Nonsud0(z)
	Implicit None
      
      Logical       onel,tevol,hspot,powin,afix
      Common/SMlcom/onel,tevol,hspot,powin,afix
	Double Precision Nonsud0,z
	Double Precision Nonsud1,z1,ktqt,asbold,z0
	Double Precision BB
	Common /splitting/ BB
	Double Precision Qtem,QQt,xmin,zmax,Qg
	Common/g0int/Qtem,QQt,xmin,asbold,zmax,Qg
	if(onel) then 
	 Nonsud0 = 1.
	 else
	z1=z
	ktqt=QQt/Qtem

	if(ktqt.ge.1.) then
	  z0= 1.
	  Elseif(z1.lt.ktqt.and.ktqt.le.1.) then
	  z0=ktqt
	  elseif(ktqt.le.z1) then
	  z0=z1
	  else
	  write(6,*) ' Nonsud : this should never happen '
	  Endif
c the simple non - sudakov	 
	Nonsud1 =  -asbold*log(z0/z1)*log(QQt**2/z0/z1/Qtem**2) !hannes new 
      Nonsud1 = exp(Nonsud1)
	If(Nonsud1.lt.0.0.or.Nonsud1.gt.1.0) Then
	write(6,*) ' Problem with old Non - Suadkov: ',Nonsud1
	Endif
      Nonsud0 = Nonsud1
	endif
	 Return
	 End 
	
C------------------------------------------------------------------------
      Subroutine SMbran
C--Small x branching algorithm
      Include 'SMallx.inc'
	Integer Kid
	Common /SMbran1/Kid(0:1000)
	Double Precision Qtemgg,Ttemgg,Qtemqq,Ttemqq,Qtemqg,Ttemqg,Qtemgq,Ttemgq 
	Common/SMSUD2/ Qtemgg,Ttemgg,Qtemqq,Ttemqq,Qtemqg,Ttemqg,Qtemgq,Ttemgq 
	
      Double Precision sL,qtx,qty,QQtx,QQty,TQt,ynow,qt,ybar,ximax
	Double Precision QQt0
      Double Precision ktold,ktnew,k2new,k2old,qt2_1,qt2_2
      Double Precision tmax,xsitest,p_max,xsiold,xsi
      Double Precision wg_rej,wg_new,wgnew,wgn,qqtn
      Common/mygrid/tmax,p_max
      Integer I,j,Nbad
      Integer nobran,ikincut
      common/myvar/nobran,ikincut
	Double Precision Qbarmy,Qbar_min,Qbar_max
      Common/mglubran/Qbarmy,Qbar_min,Qbar_max
	Integer nrglu
	Common/mynrglu/nrglu
      Integer Ineg
      common/mycheck/Ineg	
	Integer nsidx
	common /myns/nsidx
	Integer ngl_tot,nsreg
	Real zgen,zminus,qp,qm,ktqt,wtg,ktg,qtg
	common /myzval/zgen(500),zminus(500),qp(500),qm(500),
     &	nsreg(500),ngl_tot,ktqt(500),wtg(500),ktg(500),qtg(500)
	Double Precision xgb_old,z_minus,z_plus
	Integer maxglu,nglu
	parameter (maxglu=100)
      Double Precision zp,zm,xkt,xqt,xglu
	Integer II,I1
	Logical QQtcut
ccc 	Data QQtcut/.false./
	Data QQtcut/.true./  ! new standard
	Double Precision QQtc
      Integer isat
	Data isat/0/ ! no saturation ktcut
c      Data isat/1/ ! sat. ktcut
	Double Precision ktcut
        Common/kcuts/ktcut 
        Data Nbad/0/
	Double Precision alp_dns,as
	Integer Ismallz
	Common /smallz/ Ismallz,alp_dns,as
      Double Precision weight,wtot1
	Double Precision SMrgen
	
	double precision xsim,ykt,pzt,pet,pem
	common/hzbval/xsim,ykt,pzt,pet,pem
	Double Precision Qtem_old
	Integer inew
      common/hannestest/Qtem_old,inew
	logical my_glu
	common/glumy/my_glu
	Double Precision xgbarp
	
	
	Double Precision kt2
	
      Double Precision zz
	double Precision xsat,lamb
	data xsat/0.41d-4/
cccc        data xsat/1./

	data lamb/0.277/
      Double Precision Qscut
	Integer Iqbar
	Data iqbar/0/
cc	Data iqbar/1/
      Logical First
	Data First/.true./
	Double Precision scal_qbar
	Data scal_qbar/1./
	
	Real Qmax_check,xg1,xg1b,xg2,xg2b
	Common/parton_cas/Qmax_check,xg1,xg1b,xg2,xg2b
	Integer Nbp
      Parameter (Nbp=50)
      Double Precision X3M(0:Nbp)
      common/gridtt/x3m
      Integer ix3,ixold
      
      
C Parameters for Sudakov form factors
C (NMXSUD= max no of entries in lookup table)
	Real ACCUR,QEV,SUD
	Integer INTER,NQEV,NSUD
	Integer NMXSUD
      PARAMETER (NMXSUD=1024)
      COMMON/SMUSUD/ACCUR,QEV(NMXSUD,6),SUD(NMXSUD,6),INTER,NQEV,NSUD
	Real Fdiv(NMXSUD),Adiv(NMXSUD)
	Real DIVDIF
	External DIVDIF
      Double Precision Sudakov(4)
      Real Qtest,Sudtest
      
      Integer Igq,Iqg
      Data Iqg/0/ ! Iqg=1 for g->qq splitting
      Data Igq/0/ ! Igq=1 for q->gq splitting
*	
      Integer Ipgg,ns_sel
      Double Precision scal 
      common/Pggsel/Ipgg,ns_sel,scal
      Double Precision Splittgq,Splittqg,Splittqq,Splitt,asbmy
      External Splittgq,Splittqg,Splittqq,Splitt, asbmy
      Double Precision wxgg,wxqg,wxqq,wxgq
      Double Precision wz
      Common/smweights/wz
      Double Precision Cf,Ca
      Double Precision nf_light
      
      nf_light = 3.
      Cf = 4.d0/3.d0
      Ca = 3.d0
      ktcut = 0.4
	QQtc = Qg   
      Qscut = 0.5
	if(first) then
	   first=.false.
	   write(6,*) ' smbran_new-all: QQtcut,ktcut,Isat',QQtcut,ktcut,Isat
	   write(6,*) ' smbran_new-all: iqbar',iqbar
	   write(6,*) ' smbran QQtc = ',QQtc,'Qg = ',Qg,' Qs = ',Qs
         write(6,*) ' smbran: g->qq splitting (Iqg) ',Iqg
         write(6,*) ' smbran: q->qq splitting (Igq) ',Igq
	endif   
	Ineg = 0
      ktold=0.
      ktnew =0.
      k2new = 0.
      k2old = 0.
	xsiold = 0.
      wg_rej = 1.
      wg_new = 1.
      Excode=0.
      I=0
      Told=Ts
      Nsoft=0
      nglu=I+1
      if(my_glu) then
	  yphton = 1.
	  yphbar = 0.
      endif

      ximax=(xgbar+yphton)/(xg+yphbar)
		
	
      xsim = 999999.
      QQtx=P(1,0) 
      QQty=P(2,0) 
      ktold = P(1,0)**2 + P(2,0)**2	
c new acc to cascade
      xgbarp =-ktold/(1.d0-xg)/Sbig
	
      ximax=(xgbarp+yphton)/(xg+yphbar)  

		
      tmax = 1.E-5
      tmax = max(1.E-5,ximax)
		
	Qtem_old = Qg

c--loop over branches
      ixold = 0
  100 Ifin=I
      Jfin=Ifin
C--Limit on Qtem
      Qbar = Qbarmy
      Qhard=Qbar 
      ktold = P(1,I)**2 + P(2,I)**2
      QQt0=Sqrt(QQtx**2+QQty**2)
      QQt0=Sqrt(ktold)
c      if(nobran.ge.1) then
cc no branching at all:
c         Qbar=-1.
c         tmax = 1.E12
cc         write(6,*) ' no-branching selected '	  
c      endif
	
C--Find Qtem
      Call SMqtem(.false.)

      If(Kid(Ifin).ne.21) then
c treat quark splittings
c only q -> g up to now
         Qtem = Qtemqg
 	   Ttem = Ttemqg 
      Else If(Kid(Ifin).eq.21) then
c treat gluon splittings
c only g -> g g 
	   Qtem = Qtemgg
	   Ttem = Ttemgg
      Endif   
      if(nobran.ge.1) then
c no branching at all:
         Qbar=-1.
         tmax = 1.E12
	 Qtem = Qmax
c         write(6,*) ' no-branching selected '	  
      endif
      if(Qtem.gt.Qmax) then
c	   write(6,*) ' Qtem = ',Qtem,' > ',Qmax 
	   go to 600
      endif
      ix3=ixold-1
345   ix3=ix3+1
      if(ix3.gt.50) then 
         write(6,*) ' FATAL ix3>50',qtem,ix3,x3m(ix3)
         stop
      endif
      if(qtem.lt.x3m(ix3)) then
         ixold = ix3
      else 
        p_max = x3m(ix3) 
        kt2 = P(1,Ifin)**2 + P(2,Ifin)**2
        call updfgrid(2)
        goto 345
      endif

      If (Qtem.gt.Qbar) go to 600
	if(Qtem.gt.Qmax) then
	   go to 600
	endif
c--Compute provisional qt and Qt
      Call SMrazm(Qtem,qtx,qty)
      QQtx=P(1,I)+qtx
      QQty=P(2,I)+qty
      ktold = P(1,I)**2 + P(2,I)**2
      QQt=Sqrt(QQtx**2+QQty**2)
      if(isat.eq.1) ktcut = max(0.3,Qscut*(xg/xsat)**(-lamb/2))

      if(QQtcut) then
         If (QQt.le.QQtc) then
            Told=Ttem
            go to 100
         Endif
	endif
      If (QQt.le.ktcut) then
         Told=Ttem
         go to 100
         else
      Endif
      TQt=Log(QQt/QCDLam)
c check on Sudakov factors
      Qtest = Qtem**2
      Do IS=1,4
c IS=1 -> gg; =2 -> gq; =3 -> qq; =4 -> qg   
      sudtest = DIVDIF(SUD(1,IS),QEV(1,IS),NQEV,Qtest,INTER)
      sudakov(IS) = exp(-sudtest)
c      write(6,*) ' smbran suda: ',IS,Qtem,sudtest,sudakov(IS)
      End Do
      
      
      weight = 0.
      wxqq = 0.
      wxqg = 0.
      wxgq = 0.
      wxgg = 0.
      Call Szvalnew
      ynow=1.-znow
      wxgg = wg  
      if(wg.lt.0 ) write(6,*) ' smbran_new-all-std: wg < 0 ',wg

      if(Iqg.eq.1) then 
         wxqg = asbmy(scal*qtem*ynow)* 2*nf_light*Splittqg(znow) *0.5/2./Ca *wz
         wxqg = wxqg * Sudakov(4)/Sudakov(1)
      endif
      wxqq = asbmy(scal*qtem*ynow)* Splittqq(znow)*Cf /2./Ca *wz
ccccccccccccccccccccccc      wxqq = wxqq * Sudakov(3)/Sudakov(1)
      if(Igq.eq.1) then 
         wxgq = asbmy(scal*qtem*ynow)* Splittgq(znow)*Cf /2./Ca *wz
         wxgq = wxgq * Sudakov(2)/Sudakov(1)
      Endif 
c call my modified smzval
      If(Kid(Ifin).ne.21) Then
         if(smrgen(3).le.wxqq/(wxqq+wxgq)) then
	      Kid(I+1) = Kid(Ifin)
 	      wg = wxqq
            else
cc            write(6,*) ' SBRAN FATAL: select g - splitting for quark'
            
	      Kid(I+1) = 21 
         Endif
	Elseif(Kid(Ifin).eq.21) Then 
         if(smrgen(4).le.wxgg/(wxgg+wxqg)) then
	      Kid(I+1) = 21
	      wg = wxgg
            else
cc            write(6,*) ' SBRAN FATAL: select q - splitting for gluon'
c	      Kid(I+1) = Kid(Ifin)
	      Kid(I+1) = 3 ! test
            wg=wxqg
         Endif
	Endif
c	wg = wxgg
      if(wg.eq.0) then
         write(6,*) 'smbran weight = 0',wg,wxqq, Qtemgg,Qtemqg,Qbar
      endif 
      if(wg.lt.0) write(6,*) ' smbran_new-all-std test ',wxgg,wxqg,weight
c--Compute true qt and Qt
      ynow=1.-znow
      qt=Qtem*ynow
      ybar=qt**2/(xg*ynow*Sbig)

      xsi  = ybar/(xg*ynow)
	if(xsi.lt.xsiold) then
	  write(6,*) ' non ordred xsi_new ',xsi ,' xsi_old ',xsiold
	  endif
	  xsiold = xsi

      qtx=qtx*ynow
      qty=qty*ynow
      QQtx=P(1,I)+qtx
      QQty=P(2,I)+qty
      QQt=Sqrt(QQtx**2+QQty**2)
	QQtn = QQt

      ktnew = QQt**2
	
c my check
      qt2_1 = qt**2
      qt2_2 = qtx**2 + qty**2
      If(ikincut.eq.1) then
c LDC/ Catani kin constraint
         if(ktnew.le.znow*qt2_1) then
            Told=Ttem
            go to 100
         endif
      Elseif(ikincut.eq.2) then
c martin kwiecinski
         if(ktnew.le.znow*qt2_1/(1-znow)) then
            Told=Ttem
            go to 100
         endif
      Elseif(ikincut.eq.3) then
c BFKL kin constraint
         if(ktnew/znow.le.ktold) then
            Told=Ttem
            go to 100
         endif
      Endif
      If (onel.and.QQt.gt.Qhard) then
	  go to 600
	endif
      if(QQtcut) then
         If (QQt.le.QQtc) then
            Told=Ttem
            go to 100
         Endif
      Endif
      If (QQt.le.ktcut) then
         Told=Ttem
         if(isat.eq.0) then
         	go to 100
         else
         	go to 100
	   endif
      Endif
c---Exit conditions for outside kinematic limits
      If (QQt.gt.Qmax) then
         if(my_glu) then
	      else
            Excode=10.
	      goto 600
	   endif
      Endif
      xg=znow*xg
	xgb_old=xgbar
      xgbar=xgbar-ybar
	z_minus=abs(xgb_old/xgbar)
      if(xgbar.gt.0) then
         write(6,*) ' xg = ',xg,' xgbar = ',xgbar
      endif
      ktold = QQt**2
c new acc to cascade
      xgbarp =-ktold/(1.d0-xg)/Sbig
	
      ximax=(xgbarp+yphton)/(xg+yphbar) ! test 1


c---Exit conditions for xg < xmin
      If (xg.le.xmin) then
         Excode=15.
	   goto 600
      Endif

c--add a branch
      I=I+1
      P(1,I)=QQtx
      P(2,I)=QQty
      P(3,I)=xg*p2(3)+xgbar*p1(3)
      P(4,I)=xg*p2(4)+xgbar*p1(4)
      Call SMumas(P(1,I))
      k2new = P(5,I)**2
      xsitest = xgbar/xg
	if(abs(xsitest).ge.abs(p2(3)/p1(3))) then
         Ineg = 1
	endif
      k2old = P(5,I)**2
C--Update weight
      Wt=P(6,I-1)*wg
	wtg(i)=wt
      P(6,I)=Wt
      P(7,I)=0.
      P(8,I)=qt
      
      
      

      If (onel) then
c--1-loop
         Told=Ttem
      Else
c--all loops
         Told=Ttem+Log(znow)
         If (Told.lt.Tg) then
            Told=Tg
            go to 100
         Endif
      Endif
        Inew = 0
      p_max = Qtem_old	
      If (QQt.le.QQtc) then
            write(6,*) ' QQtcut ',QQt,QQtc,Qtem_old
      endif 

      go to 100
C--End spacelike branching
  600 Ifin=I
  
      xsim = ximax
	ykt = sqrt(ktold)
      Qmax_check = Qbar
	xg1 = yphton
	xg1b = yphbar
	xg2 = xg
	xg2b =xgbarp
c new acc to cascade
      kt2 = P(1,Ifin)**2+P(2,Ifin)**2
cc      xgbarp =-kt2/(1.d0-xg)/Sbig

      P(3,Ifin)=xg*p2(3)+xgbar*p1(3)
      P(4,Ifin)=xg*p2(4)+xgbar*p1(4)


      Call SMumas(P(1,Ifin))
      
      pzt = P(3,Ifin)
	pet = P(4,Ifin)
	pem = P(5,Ifin)
      if(nglu.lt.0) write(6,*) ' nglu = ',nglu,Ifin
	

      ngl_tot = Ifin-1
      ngl_tot = Ifin

      End
	Function wgnew(qqtnew)
	
      Include 'SMallx.inc'
	
	Double Precision wgnew,qqtnew
	Double Precision as,D_ns
	Double Precision alp_dns,asold
	Integer Ismallz
	Common /smallz/ Ismallz,alp_dns,asold
	Integer Ipgg,ns_sel
	Double Precision scal 
	common/Pggsel/Ipgg,ns_sel,scal
	Double Precision Nonsud
	External Nonsud
	Double Precision asbmy
	External asbmy
	Double Precision Qtemn,QQtn,xminn,asbn,zmaxn,Qgn
	Common/g0int/Qtemn,QQtn,xminn,asbn,zmaxn,Qgn
	QQtn = QQtnew
	
	If(Ismallz.eq.0) then
	   wgnew=1.
	else
         if(Ipgg.ge.2) then
c	      write(6,*) ' asmy (q) Ipgg: ',Ipgg
	      as= asbmy(scal*Qtem)
	   else
c	      write(6,*) ' asmy (k) Ipgg: ',Ipgg
	      as=asbmy(scal*QQt)
	   endif
         D_ns = Nonsud(znow)
	   wgnew=as*D_ns/alp_dns
	endif
	Return
	End
