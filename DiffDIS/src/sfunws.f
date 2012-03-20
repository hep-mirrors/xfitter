C !!! suspicious lines 294-5: xcmax instead of xbmax?
C
C==========================================================
cws      subroutine sfun(x,q2,f2p,flp,f1p,rp,f2n,fln,f1n,rn,
cws     +f2c,flc,f1c,f2b,flb,f1b)
      subroutine sfunws(x,q2,f2p,flp,f2c,flc,f2b,flb)
cws      implicit real*8(a-h,o-z)
      implicit none
      DOUBLE PRECISION x,q2,f2p,flp,f2c,flc,f2b,flb,fflc,fflb
      DOUBLE PRECISION XI,WI,XX
      INTEGER NTERMS
      DOUBLE PRECISION alambda,flavor,qsct,qsdt
      INTEGER iord,IFL,iorder,IF3,I
      DOUBLE PRECISION ALPHAzo,pi,pi2,xlam2,al,cf,ca,enf,dpsi2
      DOUBLE PRECISION ZETA2,ZETA3
      DOUBLE PRECISION Y,XY,fpxy,fcxy,fbxy,gluxy
      DOUBLE PRECISION ww2,thcq,thcw,epsc,epsc4,fpsc4,FAC
      DOUBLE PRECISION thbq,thbw,epsb,epsb4,fpsb4,AL1
      DOUBLE PRECISION upv,dnv,usea,dsea,str,chm,bot,glu
      DOUBLE PRECISION fp,fc,fb,ffp,ffc,ffb,fflp,facc,facb
      DOUBLE PRECISION C22,C23,CG2,f1lq,f1lg,FSXY
      DOUBLE PRECISION Y1,Y2,Y3,Y4,Y5,YP1
      DOUBLE PRECISION DL,DL2,DL3,DL4,DLM1,DLP1,DLM2
      DOUBLE PRECISION ALI2,ALI2M,ALI21,ALI3,ALI3M,S12M
      DOUBLE PRECISION FASTLI2,FASTLI3,FASTS12,cheavyzo,coeff2
      DOUBLE PRECISION FACT1,FACT2,FACT3,FACT4,FACT5,FACT6,FACT7
      DOUBLE PRECISION FNS2LQ,FS2LQ,F2LG
      DOUBLE PRECISION del,xyu,dfcxy,fcxyu,fcxyl,fbxyu,fbxyl,xyl
      DOUBLE PRECISION c0c,cg21c,cg22c,clg2c,cl1c,xcmax
      DOUBLE PRECISION c0b,cg21b,cg22b,clg2b,cl1b,xbmax
      DOUBLE PRECISION dfbxy
      COMMON/GAUS96/XI(96),WI(96),XX(97),NTERMS
      COMMON/TRSFPAR/alambda,flavor,qsct,qsdt,iord
      DOUBLE PRECISION F2pg,FLpg,F2cg,FLcg,F2bg,FLbg
      DOUBLE PRECISION www,Fg2
      COMMON/GCONTR/F2pg,FLpg,F2cg,FLcg,F2bg,FLbg
c-- alambda = 4-flavour Lambda_QCD
c-- qsdt = 4*mc^2, qsct = 4*mb^2 -- really!
      parameter (pi=3.1415926535897932384626433832795D0, pi2=pi**2)
cws      data pi,pi2/3.14159,9.8696/
      IFL=4
      xlam2=alambda*alambda
c      t=dlog(q2/xlam2)
c      write(*,*)'t',t
      al=ALPHAzo(dlog(q2/xlam2))/(4.*pi)
c      argmin=qsdt/4./xlam2
c      scale=dsqrt(q2)

      cf=4./3.
      ca=3.
      enf=flavor
      iorder=iord
      dpsi2=2./9.
      ZETA2=1.64493406684823
      ZETA3=1.20205690315959

      ww2=(1.-x)*q2/x
      epsc4=qsdt/q2
c      write(*,*)'ww2,epsc4',ww2,epsc4
      fpsc4=qsdt/ww2
      epsc=epsc4/4.
      thcq=1.
      if(epsc.gt.1.) thcq=0.
      thcw=1.
      if(fpsc4.gt.1.) thcw=0.
      epsb4=qsct/q2
      fpsb4=qsct/ww2
      epsb=epsb4/4.
      thbq=1.
      if(epsb.gt.1.) thbq=0.
      thbw=1.
      if(fpsb4.gt.1.) thbw=0.
      call xpdfvs(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
c      write(*,*)'did we get here',x,q2,upv,dnv,glu
      if(epsc.gt.1.) chm=0.
      if(epsb.gt.1.) bot=0.
      fp=(4.*upv+dnv+8.*usea+2.*dsea+2.*str)/9.
      fc=8.*chm/9.
      fb=2.*bot/9.

      ffp=0.
      ffc=0.
      ffb=0.
      fflp=0.
      fflc=0.
      fflb=0.

      IF(IORD.LE.0.) THEN 
        GO TO 27
      ELSE
        GO TO 22
      ENDIF
   22 CONTINUE
      IF3=0

      FAC=FLAVOR
      facc=1.
      facb=1.

      IF(ifl.EQ.3.OR.ifl.EQ.4) then
      FAC=6./9.
      facc=4./9.
      facb=1./9.
      endif

      IF(IFL.EQ.1) IF3=1
      CF=4./3.
      AL1=dLOG(1.-X)

      ffp=Fp+Fp*AL*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))
c      ffc=ffc+fc*AL*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))
c      ffb=ffb+fb*AL*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))

      DO 23 I=1,NTERMS
      Y=0.5*(1.-X)*XI(I)+0.5*(1.+X)
      XY=X/Y
      AL1=dLOG(1.-Y)
      call xpdfvs(xy,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
c      write(*,*)'did we make it',xy,upv,glu
      if(epsc.gt.1.) chm=0.
      if(epsb.gt.1.) bot=0.
      fpxy=(4.*upv+dnv+8.*usea+2.*dsea+2.*str)/9.
      fcxy=8.*chm/9.
      fbxy=2.*bot/9.
      gluxy=glu
      C22=CF*(6.+4.*Y-2.*(1.+Y*Y)/(1.-Y)*dLOG(Y)-2.*(1.+Y)*AL1
     2-IF3*2.*(1.+Y))
      C23=CF*(-3.+4.*AL1)/(1.-Y)
      CG2=2.*FAC*(-1.+8.*Y*(1.-Y)+(1.-2.*Y+2.*Y*Y)*dLOG(1./Y-1.))
      f1lq=4.*cf*y
      f1lg=8.*enf*y*(1.-y)

      ffp=ffp+.5*(1.-X)*WI(I)*AL*(C22*fpxy+C23*(fpxy-fp))
c      ffc=ffc+.5*(1.-X)*WI(I)*AL*(C22*fcxy+C23*(fcxy-fc))
c      ffb=ffb+.5*(1.-X)*WI(I)*AL*(C22*fbxy+C23*(fbxy-fb))

      fflp=fflp+.5*(1.-x)*wi(i)*al*f1lq*fpxy

      IF(IFL-1) 23,23,24
   24 CONTINUE

      www = .5*(1.-X)*WI(I)*AL
      F2pg = www*CG2*gluxy
      ffp=ffp + F2pg
      FLpg = www*dpsi2*f1lg*gluxy
      fflp=fflp+ FLpg

      FSXY=DPSI2*(UPV+DNV+2.*USEA+2.*DSEA+2.*STR)
      Y1=1.-Y
      Y2=Y*Y
      Y3=Y*Y2
      Y4=Y3*Y
      Y5=Y4*Y
      DL=DLOG(Y)
      DL2=DL*DL
      DLM1=DLOG(Y1)
      YP1=Y+1.
      DLP1=DLOG(YP1)
      DLM2=DLOG(1.-Y2)
      DL3=DLOG(Y/Y1)
      DL4=DLOG(Y2/Y1)
      ALI2 =FASTLI2(Y)
      ALI2M=FASTLI2(-Y)
      ALI21=FASTLI2(Y1)
      ALI3 =FASTLI3(Y)
      ALI3M=FASTLI3(-Y)
      S12M =FASTS12(-Y)
      FACT1=4.*DL*(6.-3.*Y+47.*Y2-9.*Y3)/(15.*Y2)
     X-4.*ALI2M*(DL-2.*DLP1)-8.*ZETA3-2.*DL2*DLM2
     X+4.*DL*DLP1*DLP1-4.*DL*ALI2+0.4*(5.-3.*Y2)*DL2
     X-4.*(2.+10.*Y2+5.*Y3-3.*Y5)*(ALI2M+DL*DLP1)/(5.*Y3)
     X+4.*ZETA2*(DLM2-0.2*(5.-3.*Y2))+8.*S12M+4.*ALI3+4.*ALI3M
     X-23.*DLM1/3.-(144.+294.*Y-1729.*Y2+216.*Y3)/(90.*Y2)
      FACT2=ALI2+DL3*DL3-3.*ZETA2-(3.-22.*Y)*DL/(3.*Y)
     X+(6.-25.*Y)*DLM1/(6.*Y)-(78.-355.*Y)/(36.*Y)
      FACT3=DL4-(6.-25.*Y)/(6.*Y)
      FNS2LQ=4.*CF*(CA-2.*CF)*Y*FACT1
     X+8.*CF*CF*Y*FACT2-(8./3.)*CF*ENF*Y*FACT3
      FS2LQ=(16./9.)*CF*ENF*(3.*(1.-2.*Y-2.*Y2)*Y1*DLM1/Y
     X+9.*Y*(ALI2+DL2-ZETA2)+9.*(1.-Y-2.*Y2)*DL
     X-9.*Y*Y1-Y1*Y1*Y1/Y)
      FACT4=(1.-3.*Y-27.*Y2+29.*Y3)*DLM1/(3.*Y2)
     X-2.*Y1*DL*DLM1+2.*YP1*ALI2M+4.*ALI2+3.*DL2
     X+2.*(Y-2.)*ZETA2+Y1*DLM1*DLM1+2.*YP1*DL*DLP1
     X+(24.+192.*Y-317.*Y2)*DL/(24.*Y)+(-8.+24.*Y+501.*Y2-517.*Y3)
     X/(72.*Y2)
      FACT5=ALI2+2.*(5.+3.*Y2)*DL2/15.-(1.+3.*Y-4.*Y2)*DLM1/
     X(2.*Y)+(-2.+10.*Y3-12.*Y5)*(ALI2M+DL*DLP1)/(15.*Y3)
     X-(5.+12.*Y2)*ZETA2/15.+(4.+13.*Y+48.*Y2-36.*Y3)*DL/(30.*Y2)
     X-(4.-Y-198.*Y2+195.*Y3)/(30.*Y2)
      FACT6=-4.*ALI21+2.*YP1*(ALI2M+DL*DLP1)+Y1*DLM1*DLM1
     X+(-6.+2.*Y)*DL*DLM1+(-1./Y-9.+29./3.*Y+1./3./Y2)*DLM1
     X+3.*DL2+(1./Y+8.-13.*Y)*DL+2.*Y*ZETA2+1./3./Y+17./3.-53./9.*Y
     X-1./9./Y2
      FACT7=ALI21+DL*DLM1+(-2./3.+4./5.*Y2+2./15./Y3)*(ALI2M+DL*DLP1)
     X+(1./2./Y+3./2.-2.*Y)*DLM1-(2./3.+2./5.*Y2)*DL2
     X+(-13./2./Y-39.+18.*Y-2./Y2)/15.*DL
     X+(-2./3.+4./5.*Y2)*ZETA2-8./15./Y-19./5.+21./5.*Y+2./15./Y2
      F2LG=16.*CA*ENF*Y*FACT6+16.*CF*ENF*Y*FACT7
      www = 0.5*(1.-x)*WI(I)*AL*AL
      Fg2 = www*DPSI2*F2LG*gluxy
      FLpg = FLpg + Fg2
      fflp=fflp + www*(FNS2LQ*fpxy + FS2LQ*FSXY) + Fg2

   23 CONTINUE
   21 CONTINUE

c      r7=dsqrt(7.d0)
      xcmax=1./(1.+epsc4)
      if(xcmax.le.x) go to 321
      DO 323 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      call xpdfvs(xy,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsc.gt.1.) chm=0.
      gluxy=glu
      fcxy=8.*chm/9.
      del=0.01
      xyu=xy*(1.+0.5*del)
      if(xyu.gt.1.d0) then
      dfcxy=0.d0
      go to 1324
      endif
      call xpdfvs(xyu,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsc.gt.1.) chm=0.
      fcxyu=8.*chm/9.
      xyl=xy*(1.-0.5*del)
      call xpdfvs(xyl,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsc.gt.1.) chm=0.
      fcxyl=8.*chm/9.
      dfcxy=(fcxyu-fcxyl)/del
 1324 c0c=cheavyzo(2,y,epsc)
      if(epsc.gt.1.d0) c0c=0.d0
      cg21c=2.*facc*cheavyzo(1,y,epsc)
      cg22c=2.*facc*c0c*dlog(1./epsc)
      clg2c=2.*facc*cheavyzo(3,y,epsc)
      f1lq=cheavyzo(4,y,epsc)      

      www = 0.5*(xcmax-x)*wi(i)
      ffc=ffc + www*c0c*(-dfcxy+3.*fcxy)
      F2cg = www*al*(cg21c-cg22c)*gluxy
      ffc=ffc + F2cg
      fflc=fflc+www*al*f1lq*fcxy
      FLcg = www*al*clg2c*gluxy
      fflc=fflc + FLcg

c      fcextra=coeff2(y,epsc)
      ffc=ffc+www*fcxy*coeff2(y,epsc)
c      c1c=cheavyzo(5,y,epsc)
      fflc=fflc+www*al*al*cheavyzo(5,y,epsc)*fcxy
      cl1c=4.*cf*(1.+y-2.*y*y+2.*y*dlog(y))
c      clg22c=2.*facc*cl1c*dlog(1./epsc)
      Fg2 = -www*al*al *2.*facc*cl1c*dlog(1./epsc)*gluxy
      FLcg = FLcg + Fg2
      fflc=fflc + Fg2 

  323 CONTINUE
  321 CONTINUE

      xbmax=1./(1.+epsb4)
      if(xbmax.le.x) go to 421
      DO 423 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      call xpdfvs(xy,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsb.gt.1.) bot=0.
      gluxy=glu
      fbxy=2.*bot/9.
      del=0.01
      xyu=xy*(1.+0.5*del)
      if(xyu.gt.1.d0) then
      dfbxy=0.d0
      go to 1424
      endif
      call xpdfvs(xyu,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsb.gt.1.) bot=0.
      fbxyu=2.*bot/9.
      xyl=xy*(1.-0.5*del)
      call xpdfvs(xyl,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsb.gt.1.) bot=0.
      fbxyl=2.*bot/9.
      dfbxy=(fbxyu-fbxyl)/del

 1424 c0b=cheavyzo(2,y,epsb)
      if(epsb.gt.1.d0) c0b=0.d0
      cg21b=2.*facb*cheavyzo(1,y,epsb)
      cg22b=2.*facb*c0b*dlog(1./epsb)
      clg2b=2.*facb*cheavyzo(3,y,epsb)
      f1lq=cheavyzo(4,y,epsb)      

      www = 0.5*(xbmax-x)*wi(i)
      ffb=ffb+www*c0b*(-dfbxy+3.*fbxy)
      F2bg = www*al*(cg21b-cg22b)*gluxy
      ffb=ffb+ F2bg

C !!! suspicious: xcmax instead of xbmax?
cws      fflb=fflb+0.5*(xcmax-x)*wi(i)*al*f1lq*fbxy
      fflb=fflb+www*al*f1lq*fbxy
cws      FLbg = 0.5*(xcmax-x)*wi(i)*al*clg2b*gluxy
      FLbg = www*al*clg2b*gluxy

      fflb=fflb+ FLbg

c      fbextra=coeff2(y,epsb)
      ffb=ffb+www*fbxy*coeff2(y,epsb)
c      c1b=cheavyzo(5,y,epsb)
      fflb=fflb+www*al*al*cheavyzo(5,y,epsb)*fbxy
      cl1b=4.*cf*(1.+y-2.*y*y+2.*y*dlog(y))
c      clg22b=2.*facb*cl1b*dlog(1./epsb)
      Fg2 = -www*al*al*2.*facb*cl1b*dlog(1./epsb)*gluxy
      FLbg = FLbg + Fg2
      fflb=fflb + Fg2

  423 CONTINUE
  421 CONTINUE

      if(ffc.lt.0.) ffc=0.
      if(ffb.lt.0.) ffb=0.

      f2p=ffp+ffc+ffb
      F2pg = F2pg + F2cg + F2bg
      f2c=ffc
      f2b=ffb

      flp=fflp+fflc+fflb
      FLpg = FLpg + FLcg + FLbg
      flc=fflc
      flb=fflb

   27 RETURN
      END

