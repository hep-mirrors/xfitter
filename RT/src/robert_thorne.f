**********************************************************************
*                                                                     *
*    Program for generating Electromagnetic Structure Functions using *
*    consistent treatment of charm and bottom structure functions     *
*    Not included are the effects due to NLO corrections to photon-   *
*    gluon fusion.   Charm mass = 1.4 GeV   Bottom mass = 4.75 GeV    *
*                                                                     *
*    The program should be run only with iord set to 2                *
*    The calculation of F_L includes the full order(alpha_s^2)        *
*    contribution                                                     *
*    The program is self contained, only requiring the subroutine     *
*    mrst2002.f and the grid file mrst2002nlo.dat.dat to be accessible*
*                                                                     *
***********************************************************************

      subroutine sfun(x,q2,mode,f2p,flp,f1p,rp,f2n,fln,f1n,rn,
     xf2c,flc,f1c,f2b,flb,f1b)
      implicit real*8(a-h,o-z)
      COMMON/GAUS96/XI(96),WI(96),XX(97),NTERMS
      COMMON/INPUT/alphaS0,alambda,flavor,qsct,qsdt,iord,inull
      common/alfscl/sclfac
      common/comptc/chmq,chmg,chmgsw,chmdcg
      common/comptcl/chmlq,chmlg
      common/comptb/botq,botg,botgsw,botdcg
      dimension pdfsf(-6:6)
      data pi,pi2,z2,z3/3.14159,9.8696,1.6449,1.20205/
      CHARACTER prefix*50
      XLAM=alambda
      IFL=4
      xlam2=xlam*xlam
      t=dlog(q2/xlam2)
      tchm=dlog(qsdt/(4.*xlam2))
      tbot=dlog(qsct/(4.*xlam2)) 
      al=alpha(t)/(4.*pi)
      alchm=alpha(tchm)/(4.*pi)
      albot=alpha(tbot)/(4.*pi)

      argmin=qsdt/4./xlam2
      q=dsqrt(q2)

      cf=4./3.
      ca=3.
      enf=flavor
      nf=enf
      iorder=iord
      dpsi2=2./9.
      z2=1.6449
      z3=1.20205

      ww2=(1.-x)*q2/x
      epsc4=qsdt/q2
      fpsc4=qsdt/ww2
      epsc=epsc4/4.
      epsc4=1.00*qsdt/q2      
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
c      call mrst2004(x,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
c      CALL GetAllPDFs(prefix,iset,x,q,upv,dnv,usea,dsea,str,sbar,
c     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,x,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)
          write(42,*) upv,dnv,usea,dsea,str,chm,bot,glu
      if(epsc.gt.1.) chm=0.
      if(epsb.gt.1.) bot=0.
      fp=(4.*upv+dnv+8.*usea+2.*dsea+2.*str)/9.
      fn=(4.*dnv+upv+2.*usea+8.*dsea+2.*str)/9.

      ffp=0.
      ffn=0.
      ffc=0.
      ffb=0.
      chmq=0.
      chmg=0.
      chmgsw=0.
      chmdcg=0.
      chmlq=0.
      chmlg=0.
      botq=0.
      botg=0.
      botgsw=0.
      botdcg=0.
      fflp=0.
      ffln=0.
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
      sfac=2./9.
      endif

      IF(IFL.EQ.1) IF3=1
      CF=4./3.
      AL1=dLOG(1.-X)

      ffp=Fp+Fp*AL*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))
      ffn=Fn+Fn*AL*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))
      fflp=fflp+fp*al**2*(-0.012)
      ffln=ffln+fn*al**2*(-0.012)




      DO 23 I=1,NTERMS
      Y=0.5*(1.-X)*XI(I)+0.5*(1.+X)
      XY=X/Y
      AL1=dLOG(1.-Y)
c      call mrst2004(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
c      CALL GetAllPDFs(prefix,iset,xy,q,upv,dnv,usea,dsea,str,sbar,
c    &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xy,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)
      if(epsc.gt.1.) chm=0.
      if(epsb.gt.1.) bot=0.
      fpxy=(4.*upv+dnv+8.*usea+2.*dsea+2.*str)/9.
      fnxy=(4.*dnv+upv+2.*usea+8.*dsea+2.*str)/9.
      fcxy=8.*chm/9.
      fbxy=2.*bot/9.
      gluxy=glu
      singxy=(upv+dnv+2.*usea+2.*dsea+2.*str)
      C22=CF*(6.+4.*Y-2.*(1.+Y*Y)/(1.-Y)*dLOG(Y)-2.*(1.+Y)*AL1
     2-IF3*2.*(1.+Y))
      C23=CF*(-3.+4.*AL1)/(1.-Y)
      CG2=2.*FAC*(-1.+8.*Y*(1.-Y)+(1.-2.*Y+2.*Y*Y)*dLOG(1./Y-1.))
      f1lq=4.*cf*y
      f1lg=8.*enf*y*(1.-y)

      ffp=ffp+.5*(1.-X)*WI(I)*AL*(C22*fpxy+C23*(fpxy-fp))
      ffn=ffn+.5*(1.-X)*WI(I)*AL*(C22*fnxy+C23*(fnxy-fn))


      fflp=fflp+.5*(1.-x)*wi(i)*(al*f1lq-2*0.0*y**2/q2)*fpxy
      ffln=ffln+.5*(1.-x)*wi(i)*al*f1lq*fnxy

      IF(IFL-1) 23,23,24
   24 CONTINUE

      ffp=ffp+.5*(1.-X)*WI(I)*AL*CG2*gluxy
      ffn=ffn+.5*(1.-X)*WI(I)*AL*CG2*gluxy

      fflp=fflp+.5*(1.-x)*wi(i)*al*dpsi2*f1lg*gluxy
      ffln=ffln+.5*(1.-x)*wi(i)*al*dpsi2*f1lg*gluxy

 
      FSXY=DPSI2*(UPV+DNV+2.*USEA+2.*DSEA+2.*STR)
      Y1=1.-Y
      DL=DLOG(Y)
      DL2=DL*DL
      DLM1=DLOG(Y1)
      DLM2=DLM1*DLM1
      DLM3=DLM2*DLM1
      DLM4=DLM3*DLM1
      FNS2LQ=128./9.*y*DLM2-46.50*y*DLM1-84.094*DL*DLM1-37.338
     X+89.53*y
     X+33.82*y**2+x*DL*(32.90+18.41*DL)-128./9.*DL
     X+16./27.*ENF*(6.*y*DLM1-12*y*DL-25.*y+6.)
      FS2LQ=ENF*((15.94-5.212*y)*Y1*Y1*DLM1+(0.421+1.520*y)*DL*DL
     X+28.09*Y1*DL-(2.370/Y-19.27)*Y1**3)
      F2LG=ENF*((94.74-49.20*y)*y1*DLM2+864.8*Y1*DLM1
     X+1161.*y*DLM1*DL
     X+60.06*y*DL*DL+39.66*Y1*DL-5.333*(1./Y-1.))
      fflp=fflp+0.5*(1.-x)*WI(I)*AL*AL*
     X(FNS2LQ*fpxy+FS2LQ*FSXY+DPSI2*F2LG*gluxy)
      ffln=ffln+0.5*(1.-x)*WI(I)*AL*AL*
     X(FNS2LQ*fnxy+FS2LQ*FSXY+DPSI2*F2LG*gluxy)



   23 CONTINUE
   
      if(epsc.gt.1.) go to 81 
      if(epsc.lt.1.) go to 82

  81  xcmax=1./(1.+epsc4)
      eps=epsc
      if(xcmax.le.x) go to 821
      DO 84 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
c      call mrst2004(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
c      CALL GetAllPDFs(prefix,iset,xy,q,upv,dnv,usea,dsea,str,sbar,
c     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xy,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)
      fpxy=(4.*upv+dnv+8.*usea+2.*dsea+2.*str)/9.
      cllff2l=cllffnsl(y,eps)
      fflp=fflp+0.5*(xcmax-x)*wi(i)*al**2.*(cllff2l*fpxy)
   84 CONTINUE
      go to 821
      
  82  xcmax=(1./(1.+epsc4)) 
      eps=epsc
      if(xcmax.le.x) go to 821
      DO 85 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
c      call mrst2004(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
c      CALL GetAllPDFs(prefix,iset,xy,q,upv,dnv,usea,dsea,str,sbar,
c     &        chm,cbar,bot,bbar,glu,phot)
         CALL FPDFXQ(1,xy,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)


      fpxy=(4.*upv+dnv+8.*usea+2.*dsea+2.*str)/9.
      cllff2h=cllffnsh(y,eps)
      fflp=fflp+0.5*(xcmax-x)*wi(i)*al**2.*(cllff2h*fpxy)
  85  CONTINUE
      go to 821 

  821 CONTINUE

      if(epsb.gt.1.) go to 91 
      if(epsb.lt.1.) go to 92

  91  xbmax=1./(1.+epsb4)
      eps=epsb
      if(xbmax.le.x) go to 921
      DO 94 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
c      call mrst2004(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
c      CALL GetAllPDFs(prefix,iset,xy,q,upv,dnv,usea,dsea,str,sbar,
c     &        chm,cbar,bot,bbar,glu,phot)
         CALL FPDFXQ(1,xy,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)

      fpxy=(4.*upv+dnv+8.*usea+2.*dsea+2.*str+8.*chm)/9.
      cllff2l=cllffnsl(y,eps)
      fflp=fflp+0.5*(xbmax-x)*wi(i)*al**2.*(cllff2l*fpxy)
   94 CONTINUE
      go to 921
      
  92  xbmax=(1./(1.+epsb4)) 
      eps=epsb
      if(xbmax.le.x) go to 921
      DO 95 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
c      call mrst2004(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
c      CALL GetAllPDFs(prefix,iset,xy,q,upv,dnv,usea,dsea,str,sbar,
c     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xy,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)
      fpxy=(4.*upv+dnv+8.*usea+2.*dsea+2.*str+8.*chm)/9.
      cllff2h=cllffnsh(y,eps)
      fflp=fflp+0.5*(xbmax-x)*wi(i)*al**2.*(cllff2h*fpxy)
  95  CONTINUE
      go to 921 

  921 CONTINUE


   21 CONTINUE

      xcmax=1./(1.+epsc4)
      if(xcmax.le.x) go to 321
      xcmup=x/xcmax
c      call mrst2004(xcmup,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
c      CALL GetAllPDFs(prefix,iset,xcmup,q,upv,dnv,usea,dsea,str,sbar,
c     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xcmup,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)

      if(epsc.gt.1.) chm=0.
      ffc=8.*chm/9. 
      continue

      xcmax=1./(1.+epsc4)
      if(xcmax.le.x) go to 321
      xcmup=x/xcmax
c      call mrst2004(xcmup,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
C      CALL GetAllPDFs(prefix,iset,xcmup,q,upv,dnv,usea,dsea,str,sbar,
C     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xcmup,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)

      if(epsc.gt.1.) chm=0.
      fc=8.*chm/9.
      AL1c=dLOG(1.-xcmup)
      ffc=Ffc+Fc*AL*CF*(-9.-2.*PI2/3.+AL1c*(-3.+2.*AL1c))

      DO 323 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
c      call mrst2004(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
C      CALL GetAllPDFs(prefix,iset,xy,q,upv,dnv,usea,dsea,str,sbar,
C     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xy,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)

      if(epsc.gt.1.) chm=0.
      gluxy=glu
      fcxy=8.*chm/9.
      ycmup=y/xcmax
      if(ycmup.gt.0.9999) ycmup=0.9999
      c0c=1.
      p0qg=ycmup**2+(1.-ycmup)**2
      C22c=CF*(6.+4.*Ycmup-2.*(1.+Ycmup*Ycmup)/(1.-Ycmup)*
     .dLOG(Ycmup)-2.*(1.+Ycmup)*dlog(1.-ycmup)
     2-IF3*2.*(1.+Ycmup))
      C23c=CF*(-3.+4.*dlog(1.-ycmup))/(1.-Ycmup)
      if(epsc.gt.1.d0) c0c=0.d0
      cg21c=2.*facc*cheavy(1,y,epsc)
      cg22c=2.*facc*c0c*p0qg*dlog(1./epsc)
      clg2c=2.*facc*cheavy(3,y,epsc)
      f1lq=cheavy(4,ycmup,epsc)
      ffc=ffc+0.5/xcmax*(xcmax-X)*WI(I)*AL*(C22c*fcxy+C23c*(fcxy-fc))      
      ffc=ffc+0.5*(xcmax-x)*wi(i)*al*(cg21c-cg22c/xcmax)*gluxy
      fflc=fflc+0.5/xcmax*(xcmax-x)*wi(i)*al*f1lq*fcxy
      fflc=fflc+0.5*(xcmax-x)*wi(i)*al*clg2c*gluxy



  323 CONTINUE

      if(epsc.gt.1.) then 
      xcmax=1./(1.+epsc4)
      eps=epsc
      if(xcmax.le.x) go to 321
      DO 324 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
c      call mrst2004(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
C      CALL GetAllPDFs(prefix,iset,xy,q,upv,dnv,usea,dsea,str,sbar,
C     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xy,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)
      gluxy=glu
      singxy=upv+2.*usea+dnv+2.*dsea+2.*str
      cgff2=facc*(c2gffnsl(y,eps)*(1-0.5*dexp(1-eps**2))+
     .c2gffnsh(y,eps)*0.5*dexp(1-eps**2))
      cqff2=facc*(c2qffnsl(y,eps)*(1-0.5*dexp(1-eps**2))+
     .c2qffnsh(y,eps)*0.5*dexp(1-eps**2))
      ffc=ffc+0.5*(xcmax-x)*wi(i)*al**2.*(cgff2*gluxy+cqff2*singxy)

  324 CONTINUE

      else 
      xcmax=1./(1.+ 4.)
      eps=1.
      if(xcmax.le.x+0.00001) go to 321
      DO 325 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      scalef=dsqrt(qsdt/4.)
      q2f=qsdt/4.
c      call mrst2004(xy,scalef,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C      CALL GetAllPDFs(prefix,iset,xy,scalef,upv,dnv,usea,dsea,str,sbar,
C     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xy,q2f,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)

      gluxy=glu
      singxy=upv+2.*usea+dnv+2.*dsea+2.*str
      cgff2=facc*c2gffnsl(y,eps)
      cqff2=facc*c2qffnsl(y,eps)
      ffc=ffc+0.5*(xcmax-x)*wi(i)*alchm**2*(cgff2*gluxy+cqff2*singxy)



  325 CONTINUE
      endif


  321 CONTINUE

      xbmax=1./(1.+epsb4)
      if(xbmax.le.x) go to 421
      xbmup=x/xbmax
c      call mrst2004(xbmup,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
C      CALL GetAllPDFs(prefix,iset,xbmup,q,upv,dnv,usea,dsea,str,sbar,
C     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xbmup,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)

      if(epsb.gt.1.) bot=0.
      fb=2.*bot/9.
      AL1b=dLOG(1.-xbmup)
      ffb=Fb+Fb*AL*CF*(-9.-2.*PI2/3.+AL1b*(-3.+2.*AL1b))

      DO 423 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
c      call mrst2004(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
C      CALL GetAllPDFs(prefix,iset,xy,q,upv,dnv,usea,dsea,str,sbar,
C     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xy,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)

      if(epsb.gt.1.) bot=0.
      gluxy=glu
      fbxy=2.*bot/9.
      ybmup=y/xbmax
      if(ybmup.gt.0.99999) ybmup=0.99999
      c0b=1.
      p0qg=ybmup**2+(1.-ybmup)**2
      C22b=CF*(6.+4.*Ybmup-2.*(1.+Ybmup*Ybmup)/(1.-Ybmup)*
     .dLOG(Ybmup)-2.*(1.+Ybmup)*dlog(1.-ybmup)
     2-IF3*2.*(1.+Ybmup))
      C23b=CF*(-3.+4.*dlog(1.-ybmup))/(1.-Ybmup)
      if(epsb.gt.1.d0) c0b=0.d0
      cg21b=2.*facb*cheavy(1,y,epsb)
      cg22b=2.*facb*c0b*p0qg*dlog(1./epsb)
      clg2b=2.*facb*cheavy(3,y,epsb)
      f1lq=cheavy(4,ybmup,epsb)
      ffb=ffb+0.5/xbmax*(xbmax-X)*WI(I)*AL*(C22b*fbxy+C23b*(fbxy-fb))      
      ffb=ffb+0.5*(xbmax-x)*wi(i)*al*(cg21b-cg22b/xbmax)*gluxy
      fflb=fflb+0.5*(xbmax-x)*wi(i)*al*f1lq*fbxy
      fflb=fflb+0.5*(xbmax-x)*wi(i)*al*clg2b*gluxy



  423 CONTINUE

      if(epsb.gt.1.) then 
      xbmax=1./(1.+epsb4)
      eps=epsb
      if(xbmax.le.x) go to 421
      DO 424 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
c      call mrst2004(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
C      CALL GetAllPDFs(prefix,iset,xy,q,upv,dnv,usea,dsea,str,sbar,
C     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xy,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)


      gluxy=glu
      singxy=upv+2.*usea+dnv+2.*dsea+2.*str+2.*chm
      cgff2=facb*(c2gffnsl(y,eps)*(1-0.5*dexp(1-eps**2))+
     .c2gffnsh(y,eps)*0.5*dexp(1-eps**2))
      cqff2=facb*(c2qffnsl(y,eps)*(1-0.5*dexp(1-eps**2))+
     .c2qffnsh(y,eps)*0.5*dexp(1-eps**2))
      ffb=ffb+0.5*(xbmax-x)*wi(i)*al**2.*(cgff2*gluxy+cqff2*singxy)

  424 CONTINUE

      else 
      xbmax=1./(1.+ 4.)
      eps=1.
      if(xbmax.le.x+0.00001) go to 421
      DO 425 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      scalef=dsqrt(qsct/4.)
      q2f=qsct/4.
c      call mrst2004(xy,scalef,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C      CALL GetAllPDFs(prefix,iset,xy,scalef,upv,dnv,usea,dsea,str,sbar,
C     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xy,q2f,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)

      gluxy=glu
      singxy=upv+2.*usea+dnv+2.*dsea+2.*str+2.*chm
      cgff2=facb*c2gffnsl(y,eps)
      cqff2=facb*c2qffnsl(y,eps)
      ffb=ffb+0.5*(xbmax-x)*wi(i)*albot**2*(cgff2*gluxy+cqff2*singxy)


  425 CONTINUE
      endif



  421 CONTINUE

      xcmax=1./(1.+epsc4)
      if(xcmax.le.x) go to 521
      xcmup=x/xcmax
c      call mrst2004(xcmup,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
C      CALL GetAllPDFs(prefix,iset,xcmup,q,upv,dnv,usea,dsea,str,sbar,
C     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xcmup,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)

      if(epsc.gt.1.) chm=0.
      fc=8.*chm/9.
      fflc=fflc+Fc*(AL**2*(-0.0012))
     .*1.25*(1/(1+4.*epsc)-0.2)


  523 continue
      if(epsc.gt.1.) then 
      xcmax=1./(1.+epsc4)
      eps=epsc
      if(xcmax.le.x) go to 521
      DO 524 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
c      call mrst2004(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
C      CALL GetAllPDFs(prefix,iset,xy,q,upv,dnv,usea,dsea,str,sbar,
C     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xy,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)

      gluxy=glu
      singxy=upv+2.*usea+dnv+2.*dsea+2.*str
      cgffl=facc*(clgffnsl(y,eps)*(1-0.5*dexp(1-eps**2))+
     .clgffnsh(y,eps)*0.5*dexp(1-eps**2)) 
      cqffl=facc*(clqffnsl(y,eps)*(1-0.5*dexp(1-eps**2))+
     .clqffnsh(y,eps)*0.5*dexp(1-eps**2))
      fflc=fflc+0.5*(xcmax-x)*wi(i)*al**2.*(cgffl*gluxy+cqffl*singxy)

  524 CONTINUE
      else 
      xcmax=(1./(1.+epsc4)) 
      eps=epsc
      if(xcmax.le.x) go to 521
      DO 525 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
c      call mrst2004(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
C      CALL GetAllPDFs(prefix,iset,xy,q,upv,dnv,usea,dsea,str,sbar,
C     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xy,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)

      if(epsc.gt.1.) chm=0.
      gluxy=glu
      fcxy=8./9.*chm
      singxy=upv+2.*usea+dnv+2.*dsea+2.*str
      ymul=y*(1+epsc4)
      Y1mul=1.-Ymul
      DL=DLOG(Ymul)
      DL2=DL*DL
      DLM1=DLOG(Y1mul)
      DLM2=DLM1*DLM1
      DLM3=DLM2*DLM1
      DLM4=DLM3*DLM1
      FNS2LQmul=128./9.d0*ymul*DLM1**2-46.50*ymul*DLM1-84.094*DL*DLM1
     x-37.338 +89.53*ymul
     X+33.82*ymul**2+ymul*DL*(32.90+18.41*DL)-128./9.d0*DL
     X+16./27.d0*enf*(6.*ymul*DLM1-12*ymul*DL-25.*ymul+6.)
      FS2LQmul=((15.94-5.212*ymul)*Y1mul*Y1mul*DLM1+(0.421+1.520*ymul)
     x*DL*DL+28.09*Y1mul*DL-(2.370/Ymul-19.27)*Y1mul**3)
      cgvfl=facc*((clgffnsh(y,eps)*(1-0.5*dexp(1-1/eps**2))+
     .clgffnsl(y,eps)*0.5*dexp(1-1/eps**2))-clgvfsub(ymul,eps)/xcmax)
      cqvfl=facc*((clqffnsh(y,eps)*(1-0.5*dexp(1-1/eps**2))+
     .clqffnsl(y,eps)*0.5*dexp(1-1/eps**2)))
      fflc=fflc+0.5/xcmax*(xcmax-x)*wi(i)*al**2.*(fcxy*FNS2LQmul+
     xfcxy*FS2LQmul)*1.25*(1/(1+4.*eps)-0.2)
      fflc=fflc+0.5*(xcmax-x)*wi(i)*al**2.*(cgvfl*gluxy+cqvfl*singxy)



  525 CONTINUE
      endif

  526 CONTINUE


  521 continue

      xbmax=1./(1.+epsb4)
      if(xbmax.le.x) go to 621
      xbmup=x/xbmax
c      call mrst2004(xbmup,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
C      CALL GetAllPDFs(prefix,iset,xbmup,q,upv,dnv,usea,dsea,str,sbar,
C     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xbmup,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)

      if(epsb.gt.1.) bot=0.
      fb=2.*bot/9.
      fflb=fflb+Fb*(AL**2*(-0.0012))
     .*1.25*(1/(1+4.*epsb)-0.2)


  623 continue
      if(epsb.gt.1.) then 
      xbmax=1./(1.+epsb4)
      eps=epsb
      if(xbmax.le.x) go to 621
      DO 624 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
c      call mrst2004(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
c      CALL GetAllPDFs(prefix,iset,xy,q,upv,dnv,usea,dsea,str,sbar,
c     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xy,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)

      gluxy=glu
      singxy=upv+2.*usea+dnv+2.*dsea+2.*str+2*chm
      cgffl=facb*(clgffnsl(y,eps)*(1-0.5*dexp(1-eps**2))+
     .clgffnsh(y,eps)*0.5*dexp(1-eps**2)) 
      cqffl=facb*(clqffnsl(y,eps)*(1-0.5*dexp(1-eps**2))+
     .clqffnsh(y,eps)*0.5*dexp(1-eps**2))
      fflb=fflb+0.5*(xbmax-x)*wi(i)*al**2.*(cgffl*gluxy+cqffl*singxy)

  624 CONTINUE
      else 
      xbmax=(1./(1.+epsb4)) 
      eps=epsb
      if(xbmax.le.x) go to 621
      DO 625 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
c      call mrst2004(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
C      CALL GetAllPDFs(prefix,iset,xy,q,upv,dnv,usea,dsea,str,sbar,
C     &        chm,cbar,bot,bbar,glu,phot)
        CALL FPDFXQ(1,xy,q2,PDFSF,inull)
          glu=pdfSF(0)
          upv=pdfSF(2)-pdfSF(-2)
          dnv=pdfSF(1)-pdfSF(-1)
          usea=pdfSF(-2)
          str=pdfSF(-3)
          dsea=pdfSF(-1)
          chm=pdfSF(-4)
          bot=pdfSF(-5)

      if(epsb.gt.1.) bot=0.
      gluxy=glu
      fbxy=2./9.*bot
      singxy=upv+2.*usea+dnv+2.*dsea+2.*str+2.*chm
      ymul=y*(1+epsb4)
      Y1mul=1.-Ymul
      DL=DLOG(Ymul)
      DL2=DL*DL
      DLM1=DLOG(Y1mul)
      DLM2=DLM1*DLM1
      DLM3=DLM2*DLM1
      DLM4=DLM3*DLM1
      FNS2LQmul=128./9.d0*ymul*DLM1**2-46.50*ymul*DLM1-84.094*DL*DLM1
     x-37.338 +89.53*ymul
     X+33.82*ymul**2+ymul*DL*(32.90+18.41*DL)-128./9.d0*DL
     X+16./27.d0*enf*(6.*ymul*DLM1-12*ymul*DL-25.*ymul+6.)
      FS2LQmul=((15.94-5.212*ymul)*Y1mul*Y1mul*DLM1+(0.421+1.520*ymul)
     x*DL*DL+28.09*Y1mul*DL-(2.370/Ymul-19.27)*Y1mul**3)

  625 CONTINUE
      endif


  621 continue




      if(ffc.lt.0.) ffc=0.
c      if(ffb.lt.0.) ffb=0.
      f2p=ffp+ffc+ffb
      f2n=ffn+ffc+ffb
      f2c=ffc
      f2b=ffb
      chmq=ffc-chmg
      chmlq=fflc-chmlg
      botq=ffb-botg

      flp=fflp+fflc+fflb+0.0/q2*f2p
      fln=ffln+fflc+fflb
      flc=fflc
      flb=fflb
      f1p=f2p-flp
      f1n=f2n-fln
      f1c=f2c-flc
      f1b=f2b-flb
      rp=flp/f1p
      rn=fln/f1n

   27 RETURN
      END

      FUNCTION ALPHA(T)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) ! G.W. 15/02/2007


      COMMON/INPUT/alphaS0,alambda,flavor,qsct,qsdt,iord,inull
C--   G.W. 27/06/2007 Initialise new alpha_S routine from PEGASUS.
C--   Input parameter 62 is alpha_S(Q02).
C--   Set top quark mass to be twice maximum Q2 grid point.
         CALL INITALPHAS( IORD, 1.D0, sqrt(1.d0), alphaS0,
     &        sqrt(QSDT/4.D0), sqrt(QSCT/4.D0),
     &        sqrt(2.D0*1.D9) )


C--   G.W. 15/06/2007 Use new routine from PEGASUS.
      QS = alambda**2*EXP(T)
      ALPHA = ALPHASRT(sqrt(QS))

      RETURN
      END






      SUBROUTINE WATE96
C*******************************************************************
C*****							       *****
C***** THE X(I)	AND W(I) ARE THE DIRECT	OUTPUT FROM A PROGRAM  *****
C***** USING NAG ROUTINE D01BCF	TO CALCULATE THE	       *****
C***** GAUSS-LEGENDRE WEIGHTS FOR 96 POINT INTEGRATION.	       *****
C***** THEY AGREE TO TYPICALLY 14 DECIMAL PLACES WITH THE      *****
C***** TABLE IN	ABRAMOWITZ & STEGUN, PAGE 919.		       *****
C*****							       *****
C***** ---->   PETER HARRIMAN, APRIL 3RD 1990.		       *****
C*****							       *****
C*******************************************************************
      IMPLICIT REAL*8 (A-H,o-Z)
      DIMENSION	X(48),W(48)
      COMMON/GAUS96/XI(96),WI(96),XX(97),NTERMS
      NTERMS=96

      X( 1)=   0.01627674484960183561
      X( 2)=   0.04881298513604856015
      X( 3)=   0.08129749546442434360
      X( 4)=   0.11369585011066471632
      X( 5)=   0.14597371465489567682
      X( 6)=   0.17809688236761733390
      X( 7)=   0.21003131046056591064
      X( 8)=   0.24174315616383866556
      X( 9)=   0.27319881259104774468
      X(10)=   0.30436494435449495954
      X(11)=   0.33520852289262397655
      X(12)=   0.36569686147231213885
      X(13)=   0.39579764982890709712
      X(14)=   0.42547898840729897474
      X(15)=   0.45470942216774136446
      X(16)=   0.48345797392059470382
      X(17)=   0.51169417715466604391
      X(18)=   0.53938810832435567233
      X(19)=   0.56651041856139533470
      X(20)=   0.59303236477757022282
      X(21)=   0.61892584012546672523
      X(22)=   0.64416340378496526886
      X(23)=   0.66871831004391424358
      X(24)=   0.69256453664216964528
      X(25)=   0.71567681234896561582
      X(26)=   0.73803064374439816819
      X(27)=   0.75960234117664555964
      X(28)=   0.78036904386743123629
      X(29)=   0.80030874413913884180
      X(30)=   0.81940031073792957139
      X(31)=   0.83762351122818502758
      X(32)=   0.85495903343459936363
      X(33)=   0.87138850590929436968
      X(34)=   0.88689451740241818933
      X(35)=   0.90146063531585023110
      X(36)=   0.91507142312089592706
      X(37)=   0.92771245672230655266
      X(38)=   0.93937033975275308073
      X(39)=   0.95003271778443564022
      X(40)=   0.95968829144874048809
      X(41)=   0.96832682846326217918
      X(42)=   0.97593917458513455843
      X(43)=   0.98251726356301274934
      X(44)=   0.98805412632962202890
      X(45)=   0.99254390032376081654
      X(46)=   0.99598184298720747465
      X(47)=   0.99836437586317963722
      X(48)=   0.99968950388322870559
      W( 1)=   0.03255061449236316962
      W( 2)=   0.03251611871386883307
      W( 3)=   0.03244716371406427668
      W( 4)=   0.03234382256857594104
      W( 5)=   0.03220620479403026124
      W( 6)=   0.03203445623199267876
      W( 7)=   0.03182875889441101874
      W( 8)=   0.03158933077072719007
      W( 9)=   0.03131642559686137819
      W(10)=   0.03101033258631386231
      W(11)=   0.03067137612366917839
      W(12)=   0.03029991542082762553
      W(13)=   0.02989634413632842385
      W(14)=   0.02946108995816795100
      W(15)=   0.02899461415055528410
      W(16)=   0.02849741106508543861
      W(17)=   0.02797000761684838950
      W(18)=   0.02741296272602931385
      W(19)=   0.02682686672559184485
      W(20)=   0.02621234073567250055
      W(21)=   0.02557003600534944960
      W(22)=   0.02490063322248370695
      W(23)=   0.02420484179236479915
      W(24)=   0.02348339908592633665
      W(25)=   0.02273706965832950717
      W(26)=   0.02196664443874448477
      W(27)=   0.02117293989219144572
      W(28)=   0.02035679715433347898
      W(29)=   0.01951908114014518992
      W(30)=   0.01866067962741165898
      W(31)=   0.01778250231604547316
      W(32)=   0.01688547986424539715
      W(33)=   0.01597056290256253144
      W(34)=   0.01503872102699521608
      W(35)=   0.01409094177231515264
      W(36)=   0.01312822956696188190
      W(37)=   0.01215160467108866759
      W(38)=   0.01116210209983888144
      W(39)=   0.01016077053500880978
      W(40)=   0.00914867123078384552
      W(41)=   0.00812687692569928101
      W(42)=   0.00709647079115442616
      W(43)=   0.00605854550423662775
      W(44)=   0.00501420274292825661
      W(45)=   0.00396455433844564804
      W(46)=   0.00291073181793626202
      W(47)=   0.00185396078894924657
      W(48)=   0.00079679206555731759
      DO 1 I=1,48
      XI(I)=-X(49-I)
      WI(I)=W(49-I)
      XI(I+48)=X(I)
      WI(I+48)=W(I)
    1 CONTINUE
      DO 2 I=1,96
    2 XX(I)=0.5*(XI(I)+1.)
      XX(97)=1.0
      EXPON=1.5
      DO 3 I=1,96
      YI=2.*(0.5*(1.+XI(I)))**EXPON-1.
      WI(I)=WI(I)/(1.+XI(I))*(1.+YI)*EXPON
      XI(I)=YI
      XX(I)=0.5*(1.+YI)
    3 CONTINUE
      RETURN
      END



C----------------------------------------------------------------------

      integer function locx(xx,nx,x)
C--   returns an integer j such that x lies inbetween xx(j) and xx(j+1).
C--   nx is the length of the array with xx(nx) the highest element.
      implicit none
      integer nx,jl,ju,jm
      double precision x,xx(nx)
      if(x.eq.xx(1)) then
         locx=1
         return
      endif
      if(x.eq.xx(nx)) then
         locx=nx-1  
         return
      endif
      ju=nx+1
      jl=0
    1 if((ju-jl).le.1) go to 2
      jm=(ju+jl)/2
      if(x.ge.xx(jm)) then
         jl=jm
      else
         ju=jm
      endif
      go to 1
    2 locx=jl
      return
      end

C----------------------------------------------------------------------

      double precision function polderiv1(x1,x2,x3,y1,y2,y3)
C--   returns the estimate of the derivative at x1 obtained by a
C--   polynomial interpolation using the three points (x_i,y_i).
      implicit none
      double precision x1,x2,x3,y1,y2,y3
      polderiv1=(x3*x3*(y1-y2)+2.d0*x1*(x3*(-y1+y2)+x2*(y1-y3))
     &     +x2*x2*(-y1+y3)+x1*x1*(-y2+y3))/((x1-x2)*(x1-x3)*(x2-x3))
      return
      end

      double precision function polderiv2(x1,x2,x3,y1,y2,y3)
C--   returns the estimate of the derivative at x2 obtained by a
C--   polynomial interpolation using the three points (x_i,y_i).
      implicit none
      double precision x1,x2,x3,y1,y2,y3
      polderiv2=(x3*x3*(y1-y2)-2.d0*x2*(x3*(y1-y2)+x1*(y2-y3))
     &     +x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))
      return
      end

      double precision function polderiv3(x1,x2,x3,y1,y2,y3)
C--   returns the estimate of the derivative at x3 obtained by a
C--   polynomial interpolation using the three points (x_i,y_i).
      implicit none
      double precision x1,x2,x3,y1,y2,y3
      polderiv3=(x3*x3*(-y1+y2)+2.d0*x2*x3*(y1-y3)+x1*x1*(y2-y3)
     &     +x2*x2*(-y1+y3)+2.d0*x1*x3*(-y2+y3))/
     &     ((x1-x2)*(x1-x3)*(x2-x3))
      return
      end

C----------------------------------------------------------------------


      real*8 function cheavy(i,z,eps)
      implicit real*8(a-h,o-z)

c     this function returns the values of C_g(z,Q^2) 
c     and the deriv. wrt log Q^2. Here eps=m^2/Q^2.
c     If i=1  C_g for F2.  If i=2 deriv of C_g for F2 
c     If i=3  C_g for FL.  If i=4 (1-m2/Q2)*beta*Clq(massless)

      if(i.gt.4) stop
      z1=1.-z
      z2=z*z
      eps2=eps*eps
      beta2=1.-4.*eps*z/z1
      beta2alt=1.-4.*eps*z/z1
C      if(beta2.lt.0.) go to 10

C
C Jan 24 2010, SG: avoid negative beta2:
C
      if (beta2.lt.0) then
         beta2 = 0.
         beta2alt = 0.
      endif
      beta=dsqrt(beta2)
      betaalt=dsqrt(beta2alt)
      a=z2+z1*z1
      b=4.*z*(1.-3.*z)
      c=-8.*z2
      aa=8.*z*z1-1.
      bb=-4.*z*z1
      arg=(1.+beta)/(1.-beta)
      fac=dlog(arg)
      cf=4./3.
      go to (1,2,3,4) i
    1 cheavy=betaalt/beta*((a+b*eps+c*eps2)*fac+(aa+bb*eps)*beta)
      return
    2 cheavy=(-b*eps-2.*c*eps2)*fac+(a+b*eps+c*eps2)/beta
     .      +(-bb*eps)*beta +(aa+bb*eps)*2.*z*eps/z1/beta
      return
    3 cheavy=-bb*beta+c*eps*fac
      return
    4 cheavy=4.*cf*z*1.25*(1/(1+4.*eps)-0.2)
      return
   10 print 99
   99 format(1x,'x > x0')
      stop
      end

c Simplified version of NLO c2g in FFNS Q^2<M^2
      real*8 function c2gffnsl(z,eps)
      implicit real*8(a-h,o-z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)
      xi=1./eps
      pi = 3.14159265359d0
      rho=1.-beta2 
 
      term1=rho*(1.5*0.25d0/pi*1.d0/(1.d0 
     .+ 0.25d0*(xi))*
     .(beta*(dlog(8.d0*beta*beta))**2- 5.d0*beta*dlog(8.d0*beta*beta) 
     .- 0.25d0*pi*pi)+2./3.*0.25d0/pi*1.d0/(1.d0 
     .+ 0.25d0*(xi))*pi*pi/2.d0)*xi*16.*pi/z

      term2 = beta*(-0.747*(xi)*(z*(dlog(z))**4)*(1+
     .1.626*dlog(xi)+1.008*(dlog(xi))**2-0.205*(dlog(xi))**3)
     .-39.34*(xi)*((dlog(1.-z)))
     .*(1-5.521*dlog(xi)+1.934*(dlog(xi))**2)
     .+867.6*(xi)*((dlog(1.-z))**2)
     .*(1-0.6569*dlog(xi)+0.1751*(dlog(xi))**2)
     .+1.978*(xi)*(z*(dlog(z))**3)
     .*(1-5.107*dlog(xi)-1.874*(dlog(xi))**2) 
     .+8.008*(xi)*(1+0.0333*dlog(1./z))*beta**14.238*(1-
     .0.5042*dlog(xi)-0.1053*(dlog(xi))**2+0.01987*(dlog(xi))**3)
     .+6.541*(xi)**0.4252*beta)/z
 
      c2gffnsl=term1+term2 
  
      return
      end

c Simplified version of NLO c2q in FFNS Q^2<M^2
      real*8 function c2qffnsl(z,eps)
      implicit real*8(a-h,o-z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)
      xi=1./eps

      term1=beta*(320.54*(xi)*z**2.*(1+
     .2.257*dlog(xi)+3.104*(dlog(xi))**2-0.5681*(dlog(xi))**3)
     .+11.37*(xi)*(z*(dlog(z)))
     .*(1+4.427*dlog(xi)-18.199*(dlog(xi))**2)
     .-9.518*(xi)*(z*(dlog(z))**2)
     .*(1-2.815*dlog(xi)+2.935*(dlog(xi))**2)
     .-12.684*(xi)*(z)
     .*(1+5.0125*dlog(xi)+34.086*(dlog(xi))**2) 
     .+4.654*(xi)*(1+0.05894*dlog(1./z))*(1-
     .0.6226*dlog(xi)-0.0299*(dlog(xi))**2
     .-0.00108*(dlog(xi))**3))/z
 
      c2qffnsl=term1 
  
      return
      end

c Simplified version of NNLO l2q in FFNS Q^2<M^2
      real*8 function cl2ffnsl(z,eps)
      implicit real*8(a-h,o-z)
      xi=1./eps
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)

      term1=beta**2*xi*(-78.24*z**2*(1+
     .0.388*dlog(xi)-5.232*(dlog(xi))**2-0.197*(dlog(xi))**3)
     .-1.458*(z*dlog(z))
     .*(1-6.136*dlog(xi)+8.047*(dlog(xi))**2)
     .+153.1*(z**3)
     .*(1-1.203*dlog(xi)-7.212*(dlog(xi))**2)
     .+2.8375*(z)
     .*(1+9.749*dlog(xi)-19.823*(dlog(xi))**2) 
     .+0.262*(1-0.1015*dlog(1./z))*(1-
     .0.0239*dlog(xi)+0.07244*(dlog(xi))**2+0.0435*(dlog(xi))**3))/z
 
      cl2ffnsl=term1 
  
      return
      end



c Simplified version of NLO c2g in FFNS Q^2>M^2
      real*8 function c2gffnsh(z,eps)
      implicit real*8(a-h,o-z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)
      xi=1./eps

      term1 = (2.d0/3.d0*((8.d0-16.*z+16.*z**2)*dlog(1.-z)
     x-(4.d0-8.*z+16.*z**2)*dlog(z)-2.d0+8.*z)
     x-3.d0/2.d0*((8.d0-16.*z+16.*z**2)*dlog(1.-z)
     x+(8.d0+32.*z)*dlog(z)+16./3./z+4.d0+32.*z-124./3.*z**2))
     x*beta*(dlog(xi))**2.

      term2 = (80./3./z +304.4*z-378.3*z**2+5.96*dlog(1.-z)
     x+33.9*dlog(z)
     x-0.277*(dlog(1.-z))**2-4.47*(dlog(z))**2+80.15*z**3
     x-153.2*dlog(z)*dlog(1.-z))*beta*(dlog(xi))

       DL  = dlog(z)
       DL1 = dlog(1.-z)

      term3 =   ( 1./z * (11.90 - 42.77* DL1) + 6.10 * DL**3  
     1         - 49.92 * DL**2 - 251.02 * DL - 1123.1 - 665.7* DL1
     2         + (214.4 - 215.1*z) * DL1**3 - 145.75 * DL1**2
     3         - 770.5 * DL**2 * DL1 - 935.04 * DL * DL1**2 )*beta

      term4 = (-224./9./z -10./9.*(dlog(1.-z))**3-316.15*z
     x+200.0*z**2-27.24*dlog(1.-z)-14.52*dlog(z)
     x-2.28*(dlog(1.-z))**2+13.21*(dlog(z))**2+96.77*z**3
     x+217.06*dlog(z)*dlog(1.-z))*beta

      term5 =(9.540*((dlog(1.-z))**4)*(1+
     .0.8336*dlog(xi)+1.072*(dlog(xi))**2
     .-0.06285*(dlog(xi))**3)
     .+3.853*(-(dlog(1.-z)))
     .*(1+9.5591*dlog(xi)+14.252*(dlog(xi))**2)
     .+365.18*((dlog(1.-z))**2)
     .*(1-0.202*dlog(xi)-0.1031*(dlog(xi))**2)
     .+3.125*((dlog(1.-z))**6)
     .*(1-1.578*dlog(xi)+0.2167*(dlog(xi))**2) 
     .+31.652*(1-0.002965*dlog(1./z))*beta**1.2823*(1+
     .0.31946*dlog(xi)+0.13524*(dlog(xi))**2
     .+0.06488*(dlog(xi))**3)-0.937*beta*xi**(-0.9556))/(z*xi)

      c2gffnsh=(term1+term2+term3+term4+term5)

      return
      end

c Simplified version of NLO c2q in FFNS Q^2>M^2
      real*8 function c2qffnsh(z,eps)
      implicit real*8(a-h,o-z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)
      xi=1/eps

      term1 = -(2.d0/3.d0*(8.*(1.+z)*dlog(z)+16./3./z+4.-4.*z
     .-16./3.*z**2))*beta*(dlog(xi))**2.

      term2 = -2./3.*(8.*(1.+z)*(dlog(z))**2
     .-(8.+40.*z+64./3.*z**2)*dlog(z)-160./9./z+16.-48.*z
     .+448./9.*z**2)*beta*(dlog(xi))

      term3 = (-896./81./z -88.76*z
     x+59.07*z**2+3.41*dlog(z)+8.078*(dlog(z))**2
     x+16.42*dlog(z)*(dlog(1.-z))**2+40.98*z**3 
     x+97.04*dlog(z)*dlog(1.-z))*beta

       DL  = LOG (z)
       DL1 = LOG (1.-z)

       term4 =   ( 5.290 * (1./z-1.) + 4.310 * DL**3   
     1         - 2.086 * DL**2 + 39.78 * DL - 0.101 * (1.-z) * DL1**3 
     2- (24.75 - 13.80 * z) * DL**2 * DL1 + 30.23 * DL * DL1)*beta

      term5 =(0.219*(1+
     .15.57*dlog(xi)+13.89*(dlog(xi))**2
     .+18.17*(dlog(xi))**3)
     .+0.0522*(-(dlog(z)))
     .*(1-12.47*dlog(xi)-30.14*(dlog(xi))**2)
     .-0.0075*((dlog(z))**2)
     .*(1-4.404*dlog(xi)-14.44*(dlog(xi))**2)
     .+1.86*(z)
     .*(1-5.697*dlog(xi)-8.0884*(dlog(xi))**2) 
     .+13.738*(1-0.00555*dlog(1./z))*beta*(1+
     .0.3379*dlog(xi)+0.3241*(dlog(xi))**2
     .-0.2286*(dlog(xi))**3))/(z*xi)

      c2qffnsh=(term1+term2+term3+term4+term5)

      return
      end

c Simplified version of NLO l2 in FFNS Q^2>M^2
      real*8 function cl2ffnsh(z,eps)
      implicit real*8(a-h,o-z)
      pi = 3.14159265359d0
      xi=1/eps

      term1 = 2.d0/3.d0*(4./3.*(1.+z**2)/(1.-z))*(dlog(xi))**2.

      term2 =2.d0/3.d0*((1.+z**2)/(1.-z)*(8./3.*dlog(1.-z)
     .-16./3.*dlog(z)-58./9.)+2./3.+26./3.*z)*(dlog(xi))

      term3 =2.d0/3.d0*((1.+z**2)/(1.-z)*(-8./3.*(1.0379*(1.-z)
     .+0.509*(1.-z)**3.-0.2713*(1.-z)**7.+0.3676*(1.-z)**10.)
     .-8./3.*pi**2./6.-16./3.*dlog(z)*dlog(1.-z)
     .+4./3.*(dlog(1.-z))**2.+4.*(dlog(z))**2.-58./9.*dlog(1.-z)
     .+134./9.*dlog(z)+359./27.)+(2./3.+26./3.*z)*dlog(1.-z)
     .-(2.+46./3.*z)*dlog(z)+29./9.-295./9.*z) 

      term4 =1/(z*xi)*(7.455*(dlog(1.-z)/(1.-z))*(1-
     .18.32*dlog(xi)+1.741*(dlog(xi))**2-0.235*(dlog(xi))**3)
     .+0.4824*(z/(1.-z))
     .*(1+49.627*dlog(xi)-11.47*(dlog(xi))**2)
     .-15.76*((dlog(1.-z)))
     .*(1-10.16*dlog(xi)+0.843*(dlog(xi))**2)
     .-9.023*(((dlog(1.-z))**2)/(1.-z))
     .*(1+5.997*dlog(xi)-0.84*(dlog(xi))**2+0.0688*(dlog(xi))**3) 
     .-17.11*((dlog(1.-z))**2.)*(1+
     .2.768*dlog(xi)+0.142*(dlog(xi))**2+0.104*(dlog(xi))**3))

      cl2ffnsh=(term1+term2+term3+term4)

      return
      end

c Subtraction term for NNLO c2g in VFNS
      real*8 function c2gvfsub(z,eps)
      implicit real*8(a-h,o-z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)
      xi=1./eps

      term1 = (2.d0/3.d0*((8.d0-16.*z+16.*z**2)*dlog(1.-z)
     x-(4.d0-8.*z+16.*z**2)*dlog(z)-2.d0+8.*z)
     x-3.d0/2.d0*((8.d0-16.*z+16.*z**2)*dlog(1.-z)
     x+(8.d0+32.*z)*dlog(z)+16./3./z+4.d0+32.*z-124./3.*z**2))
     x*(dlog(xi))**2.

      term2 = (80./3./z +304.4*z-378.3*z**2+5.96*dlog(1.-z)
     x+33.9*dlog(z)
     x-0.277*(dlog(1.-z))**2-4.47*(dlog(z))**2+80.15*z**3
     x-153.2*dlog(z)*dlog(1.-z))*(dlog(xi))

      term3 = (-224./9./z -10./9.*(dlog(1.-z))**3-316.15*z
     x+200.0*z**2-27.24*dlog(1.-z)-14.52*dlog(z)
     x-2.28*(dlog(1.-z))**2+13.21*(dlog(z))**2+96.77*z**3
     x+217.06*dlog(z)*dlog(1.-z))

      c2gvfsub=(term1+term2+term3)

      return
      end

c Subraction term for NNLO c2q in VFNS
      real*8 function c2qvfsub(z,eps)
      implicit real*8(a-h,o-z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)
      xi=1/eps

      term1 = -(2.d0/3.d0*(8.*(1.+z)*dlog(z)+16./3./z+4.-4.*z
     .-16./3.*z**2))*(dlog(xi))**2.

      term2 = -2./3.*(8.*(1.+z)*(dlog(z))**2
     .-(8.+40.*z+64./3.*z**2)*dlog(z)-160./9./z+16.-48.*z
     .+448./9.*z**2)*(dlog(xi))

      term3 = (-896./81./z -88.76*z
     x+59.07*z**2+3.41*dlog(z)+8.078*(dlog(z))**2
     x+16.42*dlog(z)*(dlog(1.-z))**2+40.98*z**3 
     x+97.04*dlog(z)*dlog(1.-z))


      c2qvfsub=(term1+term2+term3)

      return
      end

c Simplified version of NLO clg in FFNS Q^2>M^2
      real*8 function clgffnsh(z,eps)
      implicit real*8(a-h,o-z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)
      xi=1./eps


      term1 = 2./3.*(32.*z*dlog(z)+16.+16.*z-32.*z**2)*beta*(dlog(xi))

       DL  = dlog(z)
       DL1 = dlog(1.-z)

      term2 = ((94.74-49.2*z)*(1.-z)*DL1**2. + 864.8*(1.-z)*DL1 
     1 +1161.*z*DL1*DL+60.06*z*DL**2.+39.66*(1.-z)*DL
     2 -5.333*(1./z-1))*beta**3.

      term3 =beta**2.3/(z*xi)*(-155.41*z**2*(1+
     .2.057*dlog(xi)-0.017*(dlog(xi))**2-0.0639*(dlog(xi))**3)
     .-14.67*(z*dlog(z))
     .*(1-1.378*dlog(xi)-2.755*(dlog(xi))**2)
     .+948.72*(z**3)
     .*(1+1.144*dlog(xi)-0.400*(dlog(xi))**2)
     .+35.24*(z)
     .*(1+3.634*dlog(xi)+2.395*(dlog(xi))**2) 
     .+7.134*(1-0.0199*dlog(1./z))*(1+
     .1.144*dlog(xi)+0.8324*(dlog(xi))**2+0.4316*(dlog(xi))**3))

      clgffnsh=(term1+term2+term3)

      return
      end

c Simplified version of NLO clq in FFNS Q^2>M^2
      real*8 function clqffnsh(z,eps)
      implicit real*8(a-h,o-z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)
      xi=1/eps

       DL  = LOG (z)
       DL1 = LOG (1.-z)

      term1 = ((15.94-5.212*z)*(1.-z)**2.*DL1 + 
     1(0.421+1.520*z)*DL**2 +28.09*(1.-z)*DL
     2-(2.370/z-19.27)*(1.-z)**3. )*beta**3.

      term2 =beta*(-1.198*z**2*(1+
     .4.146*dlog(xi)-8.22*(dlog(xi))**2+1.236*(dlog(xi))**3)
     .+0.00335*(z*(dlog(z)))
     .*(1-7.135*dlog(xi)+12.282*(dlog(xi))**2)
     .-0.188*(z*(dlog(z))**2)
     .*(1+6.095*dlog(xi)-17.989*(dlog(xi))**2)
     .-4.8214*(z)
     .*(1-1.686*dlog(xi)+3.069*(dlog(xi))**2) 
     .+4.214*(1-0.03673*dlog(1./z))*(1+
     .1.1373*dlog(xi)+0.53361*(dlog(xi))**2
     .+0.37647*(dlog(xi))**3))/(z*xi)

      clqffnsh=(term1+term2)

      return
      end

c Simplified version of NLO clg in FFNS Q^2<M^2
      real*8 function clgffnsl(z,eps)
      implicit real*8(a-h,o-z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)
      xi=1./eps
      pi = 3.14159265359d0
 
      term1 =beta**2*xi*(-22.06*z**2*(1+
     .2.24*dlog(xi)+1.298*(dlog(xi))**2+6.15*(dlog(xi))**3)
     .+2.091*(z*dlog(z))
     .*(1+3.175*dlog(xi)+0.7*(dlog(xi))**2)
     .+1566.8*(z**3)
     .*(1-1.923*dlog(xi)+1.931*(dlog(xi))**2)
     .-10.41*(z)
     .*(1-2.791*dlog(xi)-0.181*(dlog(xi))**2) 
     .+0.679*(1-0.00156*dlog(1./z))*(1+
     .0.2368*dlog(xi)-0.1357*(dlog(xi))**2-0.0484*(dlog(xi))**3))/z
 
      clgffnsl=term1 
  
      return
      end

c Simplified version of NLO clq in FFNS Q^2<M^2
      real*8 function clqffnsl(z,eps)
      implicit real*8(a-h,o-z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)
      xi=1./eps

      term1=beta**2*xi*(24.83*z**2*(1+
     .0.324*dlog(xi)+0.429*(dlog(xi))**2+0.0923*(dlog(xi))**3)
     .+1.091*(z*dlog(z))
     .*(1-1.251*dlog(xi)-0.3127*(dlog(xi))**2)
     .-5.358*(z**3)
     .*(1-1.684*dlog(xi)-5.652*(dlog(xi))**2)
     .-3.111*(z)
     .*(1+1.553*dlog(xi)+0.233*(dlog(xi))**2) 
     .+0.2443*(1+0.02402*dlog(1./z))*(1+
     .0.2075*dlog(xi)-0.1234*(dlog(xi))**2-0.0314*(dlog(xi))**3))/z
 
      clqffnsl=term1 
  
      return
      end

c Subtraction term for NLO clg in VFNS
      real*8 function clgvfsub(z,eps)
      implicit real*8(a-h,o-z)
      beta2=1.-4.*eps*z/(1.-z)

C
C Jan 24, 2010 SG: protect against negative beta2:
C
      beta=dsqrt(max(beta2,0))
      xi=1./eps

      term1 = 2./3.*(32.*z*dlog(z)+16.+16.*z-32.*z**2)*(dlog(xi))

      clgvfsub=(term1)*1.25*(1/(1+4.*eps)-0.2)

      return
      end

c Simplified version of NLO llq in FFNS Q^2<M^2
      real*8 function cllffnsl(z,eps)
      implicit real*8(a-h,o-z)
      xi=1./eps
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)

      term1=beta**5*xi*(3.237*z**2*(1-
     .3.253*dlog(xi)-1.929*(dlog(xi))**2+4.693*(dlog(xi))**3)
     .-0.234*(z*dlog(z))
     .*(1-1.237*dlog(xi)+3.638*(dlog(xi))**2)
     .-150.2*(z**3)
     .*(1-2.321*dlog(xi)+0.138*(dlog(xi))**2)
     .+1.723*(z)
     .*(1-0.786*dlog(xi)-0.1119*(dlog(xi))**2) 
     .+1.473*z*(1-0.09154*dlog(1./z))*(1-
     .0.8051*dlog(xi)-1.567*(dlog(xi))**2-0.7084*(dlog(xi))**3))/z
 
      cllffnsl=term1 
  
      return
      end

c Simplified version of NLO llq in FFNS Q^2>M^2
      real*8 function cllffnsh(z,eps)
      implicit real*8(a-h,o-z)
      xi=1./eps
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)

      term1=2.d0/3.d0*(16./3.*z*(dlog(1.-z)-2*dlog(z))
     .+16./3.-400./18.*z)
 
      term2=1/xi*(2.473*((dlog(1-z))**4)*(1-
     .2.98*dlog(xi)+0.694*(dlog(xi))**2-0.047*(dlog(xi))**3)
     .+6.924*(z**2*(dlog(z)))
     .*(1+0.628*dlog(xi)+0.496*(dlog(xi))**2)
     .+10.24*(z*(dlog(1-z))**2)
     .*(1+2.51*dlog(xi)-0.658*(dlog(xi))**2)
     .-6.869*(z*(dlog(1-z))**3)
     .*(1+0.413*dlog(xi)+0.166*(dlog(xi))**2) 
     .-2.46*(1-0.932*dlog(1./z))*beta**(-0.671)*z**2*(1+
     .0.779*dlog(xi)+1.966*(dlog(xi))**2-0.439*(dlog(xi))**3)
     .-26.84*z**2*(xi)**0.14*beta)/z     
 
      cllffnsh=term1+term2 
  
      return
      end

c Simplified version of LO clg convoluted with a2gg
      real*8 function clgconagg(z,eps)
      implicit real*8(a-h,o-z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)
      xi=1./eps

      term1 = 3.85*(0.2-z)+16525*z*(z-0.2)*dlog(5*z)
     .+457.82*dlog(5*z)*dlog(1-5*z)-54222*(1-5.*z)*z**3
     .+17745*z**2*(dlog(5*z))**2*dlog(1-5*z)+5619*z*dlog(5.*z)
     .+282565*z**2*(0.2-z)+138.7*z*(1-5*z)-18965*z**2*(1-5*z)
     .-11346*z**2*dlog(5*z)*dlog(1-5*z)

      clgconagg=(term1)/z

      return
      end





*
* =======================================================================
*
* ..This is the regular non-singlet piece for the electromagnetic F2, 
*    corresponding to C2NSP+C2NSN in W. van Neerven's program. The 
*    (10+8) numerical coefficients are fitted to his results, using x 
*    values between 10^-6 and 1-10^-6. 
*
       DOUBLE PRECISION FUNCTION C2NN2ART (Y, NF)
      IMPLICIT REAL*8(A-H,O-Z)
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2NN2ART = 
     1          - 69.59 - 1008.* Y
     2          - 2.835 * DL**3 - 17.08 * DL**2 + 5.986 * DL 
     3          - 17.19 * DL1**3 + 71.08 * DL1**2 - 660.7 * DL1
     4          - 174.8 * DL * DL1**2 + 95.09 * DL**2 * DL1
     5        + NF * ( - 5.691 - 37.91 * Y 
     6          + 2.244 * DL**2 + 5.770 * DL 
     7          - 1.707 * DL1**2  + 22.95 * DL1
     8          + 3.036 * DL**2 * DL1 + 17.97 * DL * DL1 )     
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the regular non-singlet piece for the charged-current F2, 
*    corresponding to C2NSP-C2NSN in WvN's program. For the NF^0 piece
*    8 numerical coefficients are fitted to his results, the ones of
*    ln^3(1-y) and ln^2(1-y) are taken over from C3NN2A. The NF-piece
*    is also the same as in C3NN2A. (I think he means C2NN2A (RGR))
*
       DOUBLE PRECISION FUNCTION C2NC2ART (Y, NF)
      IMPLICIT REAL*8(A-H,O-Z)
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2NC2ART = 
     1          - 84.18 - 1010.* Y
     2          - 3.748 * DL**3 - 19.56 * DL**2 - 1.235 * DL 
     3          - 17.19 * DL1**3 + 71.08 * DL1**2 - 663.0 * DL1
     4          - 192.4 * DL * DL1**2 + 80.41 * DL**2 * DL1
     5        + NF * ( - 5.691 - 37.91 * Y 
     6          + 2.244 * DL**2 + 5.770 * DL 
     7          - 1.707 * DL1**2  + 22.95 * DL1
     8          + 3.036 * DL**2 * DL1 + 17.97 * DL * DL1 )     
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular NS piece, denoted by SOFT2 in WvN's program. 
*    It is the same for all F2 and F3 cases. The numerical coefficients 
*    are exact, but truncated.
*
       DOUBLE PRECISION FUNCTION C2NS2BRT (Y, NF)
      IMPLICIT REAL*8(A-H,O-Z)
*
       DL1 = LOG (1.-Y)
       DM  = 1./(1.-Y)
*
       C2NS2BRT = 
     1          + 14.2222 * DL1**3 - 61.3333 * DL1**2 - 31.105 * DL1 
     2          + 188.64 
     3        + NF * ( 1.77778 * DL1**2 - 8.5926 * DL1 + 6.3489 ) 
       C2NS2BRT = DM * C2NS2BRT
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' NS piece for the e.m. F2, denoted by COR2 in 
*    WvN's program. The numerical coefficients of the logs are exact,
*    but truncated, the constant one (from the delta-function) is 
*    slightly adjusted (+ 0.485 - 0.0035 NF) using the lowest moments.
*
       DOUBLE PRECISION FUNCTION C2NN2CRT (Y, NF)
      IMPLICIT REAL*8(A-H,O-Z)

       DL1 = LOG (1.-Y)
*
       C2NN2CRT = 
     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
     2          + 188.64 * DL1 - 338.531 + 0.485 
     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
     4          + 6.3489 * DL1 + 46.844 - 0.0035)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' NS piece for the CC F2, also given by COR2 in 
*    WvN's program. The numerical coefficients of the logs are exact,
*    but truncated, the constant one is adjusted (- 0.2652 - 0.0035 NF) 
*    using the first moment.
*
       DOUBLE PRECISION FUNCTION C2NC2CRT (Y, NF)
      IMPLICIT REAL*8(A-H,O-Z)

       DL1 = LOG (1.-Y)
*
       C2NC2CRT = 
     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
     2          + 188.64 * DL1 - 338.531 + 0.537
     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
     4          + 6.3489 * DL1 + 46.844 - 0.0035)
*
       RETURN
       END
*
* ----------------------------------------------------------------------
*
* ..This is the regular non-singlet piece for the charged current F3,
* .. typed in by RGR
* 
       DOUBLE PRECISION FUNCTION C3NC2ART (Y, NF)
      IMPLICIT REAL*8(A-H,O-Z)
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C3NC2ART = 
     1          - 206.1 - 576.8 * Y
     2          - 3.922 * DL**3 - 33.31 * DL**2 -67.60 * DL 
     3          - 15.20 * DL1**3 + 94.61 * DL1**2 - 409.6 * DL1
     4          - 147.9 * DL * DL1**2 
     5        + NF * ( - 6.337 - 14.97 * Y 
     6          + 2.207 * DL**2 + 8.683 * DL 
     7          + 0.042 * DL1**3 - 0.808 * DL1**2  + 25.00 * DL1
     8          + 9.684 * DL * DL1 )     
*
       RETURN
       END
*
* -------------------------------------------------------------------------
*
*
* ..This is the singular non-singlet piece for the charged current F3,
* .. precisely the same as C2NS2B from vN+V paper
*
       DOUBLE PRECISION FUNCTION C3NC2B (Y, NF)
      IMPLICIT REAL*8(A-H,O-Z)
*
       DL1 = LOG (1.-Y)
       DM  = 1./(1.-Y)
*
       C3NC2B = 
     1          + 14.2222 * DL1**3 - 61.3333 * DL1**2 - 31.105 * DL1 
     2          + 188.64 
     3        + NF * ( 1.77778 * DL1**2 - 8.5926 * DL1 + 6.3489 ) 
       C3NC2B = DM * C3NC2B
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
* ..This is the 'local' NS piece for the charged current F3,
* .. smae as C2NN2C from vN+V paper
*
       DOUBLE PRECISION FUNCTION C3NC2CRT(Y, NF)
      IMPLICIT REAL*8(A-H,O-Z)

       DL1 = LOG (1.-Y)
*
       C3NC2CRT = 
     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
     2          + 188.64 * DL1 - 338.531 + 0.485 
     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
     4          + 6.3489 * DL1 + 46.844 - 0.0035)
*
       RETURN
       END
*
* -----------------------------------------------------------------------
*
*
* ..This is the pure singlet piece, denoted by C2S in WvN's program. 
*    Seven numerical coefficients (all but the one of 1/y, which is 
*    exact up to truncation) are fitted to his results, using x values
*    between 10^-6 and 1-10^-6.
*
       DOUBLE PRECISION FUNCTION C2S2ART (Y, NF)
      IMPLICIT REAL*8(A-H,O-Z)
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2S2ART =   NF * ( 5.290 * (1./Y-1.) + 4.310 * DL**3   
     1         - 2.086 * DL**2 + 39.78 * DL - 0.101 * (1.-Y) * DL1**3 
     2         - (24.75 - 13.80 * Y) * DL**2 * DL1 + 30.23 * DL * DL1 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the regular gluon piece, denoted by C2G2 in WvN's program. 
*    Nine numerical coefficients are fitted as above, the ones of 1/y, 
*    ln^3(1-y), and ln^2(1-y) are exact up to truncation.
*
       DOUBLE PRECISION FUNCTION C2G2ART (Y, NF)
      IMPLICIT REAL*8(A-H,O-Z)
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2G2ART =   NF * ( 1./Y * (11.90 + 1494.* DL1) + 5.319 * DL**3  
     1         - 59.48 * DL**2 - 284.8 * DL + 392.4 - 1483.* DL1
     2         + (6.445 + 209.4 * (1.-Y)) * DL1**3 - 24.00 * DL1**2
     3         - 724.1 * DL**2 * DL1 - 871.8 * DL * DL1**2 )
*
       RETURN
       END
* 
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' gluon piece, which has no counterpart in WvN's
*    program, as it does not exist in the exact expressions. Here it 
*    is, however, relevant for achieving a high accuracy of the convo-
*    lution, as are the adjustments of the constant in the non-singlet
*    quark coefficient functions. The value is fixed from the lowest 
*    even-integer moments. 
*
       DOUBLE PRECISION FUNCTION C2G2CRT (Y, NF)
      IMPLICIT REAL*8(A-H,O-Z)
*
       C2G2CRT = - NF * 0.28  
*
       RETURN
       END

c Approx version of NNLO c2g in FFNS Q^2<M^2
      real*8 function c2gffns3(z,eps)
      implicit real*8(a-h,o-z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)
      xi=1./eps
      eta=xi*(1.-z)/(4.*z)-1
      xmax=(1./(1.+eps*4))
      pi=3.14159265359d0

      term1 = 96*beta*(13.073*xi-23.827*xi**2+24.107*xi**3
     .-9.173*xi**4)/z*(dlog(xmax/z)-4.)*(1.-z/xmax)**20

      term2=256./z*pi**3*xi*(1./(1.+xi/4))*2.1/(1.+eta)*
     .((0.0144-0.0273*dlog(eta)+0.00235*(dlog(eta))**2
     .-0.0001033*(dlog(eta))**3-0.0000478*(dlog(eta))**4)
     .+dlog(xi)*(0.0205-0.00373*dlog(eta)+0.00339*(dlog(eta))**2
     .+0.000128*(dlog(eta))**3-0.000044*(dlog(eta))**4)  
     .+(dlog(xi))**2*(-0.00065-0.0003*dlog(eta)
     .+0.000178*(dlog(eta))**2+0.0000206*(dlog(eta))**3))

      c2gffns3=(term1+term2)

      return
      end

c Approx version of NNLO c2q in FFNS Q^2<M^2
      real*8 function c2qffns3(z,eps)
      implicit real*8(a-h,o-z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)
      xi=1./eps
      eta=1.-beta2
      xmax=(1./(1.+eps*4.))
      pi=3.14159265359d0

      term1 = 4./9.*96*beta*(13.073*xi-23.827*xi**2+24.107*xi**3
     .-9.173*xi**4)/z*(dlog(xmax/z)-4.)*(1.-z/xmax)**20

      c2qffns3=term1

      return
      end

c Approx version of NNLO clg in FFNS Q^2<M^2
      real*8 function clgffns3(z,eps)
      implicit real*8(a-h,o-z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)
      xi=1./eps
      eta=xi*(1.-z)/(4.*z)-1
      xmax=(1./(1.+eps*4))
      pi=3.14159265359d0

      term1 = 96*beta**3*(0.484*xi**2-0.567*xi**3
     .+0.239*xi**4)/z*(dlog(xmax/z)-4.)*(1.-z/xmax)**20

      clgffns3=term1

      return
      end

c Approx version of NNLO clq in FFNS Q^2<M^2
      real*8 function clqffns3(z,eps)
      implicit real*8(a-h,o-z)
      beta2=1.-4.*eps*z/(1.-z)
      beta=dsqrt(beta2)
      xi=1./eps
      eta=1.-beta2
      xmax=(1./(1.+eps*4.))
      pi=3.14159265359d0

      term1 = 4./9.*96*beta**3*(0.484*xi**2-0.567*xi**3
     .+0.239*xi**4)/z*(dlog(xmax/z)-4.)*(1.-z/xmax)**20

      clqffns3=term1

      return
      end



* =======================================================================
*
* ..This is the regular non-singlet piece for the electromagnetic F2, 
*    corresponding to C2NSP+C2NSN in W. van Neerven's program. The 
*    (10+8) numerical coefficients are fitted to his results, using x 
*    values between 10^-6 and 1-10^-6. Heavy flvour part. 
*
       DOUBLE PRECISION FUNCTION C2NN2AH (Y)
      IMPLICIT REAL*8(A-H,O-Z)
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2NN2AH = ( - 5.691 - 37.91 * Y 
     .          + 2.244 * DL**2 + 5.770 * DL 
     .          - 1.707 * DL1**2  + 22.95 * DL1
     .          + 3.036 * DL**2 * DL1 + 17.97 * DL * DL1 )     
*
       RETURN
       END

* ---------------------------------------------------------------------
*
*
* ..This is the 'local' NS piece for the e.m. F2, denoted by COR2 in 
*    WvN's program. The numerical coefficients of the logs are exact,
*    but truncated, the constant one (from the delta-function) is 
*    slightly adjusted (+ 0.485 - 0.0035 NF) using the lowest moments.
*  Heavy flavour piece.
       DOUBLE PRECISION FUNCTION C2NN2CH (Y)
      IMPLICIT REAL*8(A-H,O-Z)

       DL1 = LOG (1.-Y)
*
       C2NN2CH = (0.592593 * DL1**3 - 4.2963 * DL1**2 
     .          + 6.3489 * DL1 + 46.844 - 0.0035)
*
       RETURN
       END
*

* ---------------------------------------------------------------------
*
*
* ..This is the singular NS piece, denoted by SOFT2 in WvN's program. 
*    It is the same for all F2 and F3 cases. The numerical coefficients 
*    are exact, but truncated. Heavy flavour piece
*
       DOUBLE PRECISION FUNCTION C2NS2BH (Y)
      IMPLICIT REAL*8(A-H,O-Z)
*
       DL1 = LOG (1.-Y)
       DM  = 1./(1.-Y)
*
       C2NS2BH = ( 1.77778 * DL1**2 - 8.5926 * DL1 + 6.3489 ) 
       C2NS2BH = DM * C2NS2BH
*
       RETURN
       END
*

C----------------------------------------------------------------------
C--   G.Watt 11/06/2007.
C----------------------------------------------------------------------
C--   Stand-alone code for alpha_s cannibalised from A.Vogt's
C--   QCD-PEGASUS package (hep-ph/0408244).  The running coupling
C--   alpha_s is obtained at N^mLO (m = 0,1,2,3) by solving the
C--   renormalisation group equation in the MSbar scheme by a
C--   fourth-order Runge-Kutta integration.  Transitions from
C--   n_f to n_f+1 flavours are made when the factorisation scale
C--   mu_f equals the pole masses m_h (h = c,b,t).  At exactly
C--   the thresholds m_{c,b,t}, the number of flavours n_f = {3,4,5}.
C--   The top quark mass should be set to be very large to evolve with
C--   a maximum of five flavours.  The factorisation scale mu_f may be
C--   a constant multiple of the renormalisation scale mu_r.  However,
C--   the input factorisation scale mu_(f,0) must be less than or equal
C--   to the charm quark mass.
C--
C--   Example of usage.
C--   First call the initialisation routine (only needed once):
C--
C--    IORD = 2                  ! perturbative order (N^mLO,m=0,1,2,3)
C--    FR2 = 1.D0                ! ratio of mu_f^2 to mu_r^2
C--    R0 = 1.D0                 ! input mu_r in GeV
C--    ASI = 0.5D0               ! input value of alpha_s at mu_r = R0
C--    MC = 1.43D0               ! charm quark mass
C--    MB = 4.3D0                ! bottom quark mass
C--    MT = 1.D10                ! top quark mass
C--    CALL INITALPHAS(IORD, FR2, R0, ASI, MC, MB, MT)
C--
C--   Then get alpha_s at a renormalisation scale mu_r with:
C--
C--    MUR = 100.D0              ! renormalisation scale in GeV
C--    ALFAS = ALPHAS(MUR)
C--
C----------------------------------------------------------------------

      SUBROUTINE INITALPHAS(IORD, FR2, R0, ASI, MC, MB, MT)
C--   IORD = 0 (LO), 1 (NLO), 2 (NNLO), 3 (NNNLO).
C--   FR2 = ratio of mu_f^2 to mu_r^2 (must be a fixed value).
C--   R0 = input renormalisation scale (in GeV) for alpha_s.
C--   ASI = input value of alpha_s at the renormalisation scale R0.
C--   MC,MB,MT = heavy quark masses in GeV.
      IMPLICIT NONE

      INTEGER IORD,NAORD,NASTPS,IVFNS,NFF
      DOUBLE PRECISION FR2,R0,ASI,MC,MB,MT,LOGFR,R20,
     &     PI,ZETA,CF,CA,TR,AS0,M20,MC2,MB2,MT2
      PARAMETER(PI = 3.1415 92653 58979 D0)

      COMMON / RZETA  / ZETA(6)
      COMMON / COLOUR / CF, CA, TR
      COMMON / ASINP  / AS0, M20
      COMMON / ASPAR  / NAORD, NASTPS
      COMMON / VARFLV / IVFNS
      COMMON / NFFIX  / NFF
      COMMON / FRRAT  / LOGFR

*
* ..QCD colour factors
*
      CA = 3.D0
      CF = 4./3.D0
      TR = 0.5 D0
*
* ..The lowest integer values of the Zeta function
*
      ZETA(1) = 0.5772 15664 90153 D0
      ZETA(2) = 1.64493 40668 48226 D0
      ZETA(3) = 1.20205 69031 59594 D0
      ZETA(4) = 1.08232 32337 11138 D0
      ZETA(5) = 1.03692 77551 43370 D0
      ZETA(6) = 1.01734 30619 84449 D0

      IVFNS = 1                 ! variable flavour-number scheme (VFNS)
C      IVFNS = 0                 ! fixed flavour-number scheme (FFNS)
      NFF = 4                   ! number of flavours for FFNS
      NAORD = IORD              ! perturbative order of alpha_s
      NASTPS = 20               ! num. steps in Runge-Kutta integration
      R20 = R0**2               ! input renormalisation scale
      MC2 = MC**2               ! mu_f^2 for charm threshold
      MB2 = MB**2               ! mu_f^2 for bottom threshold
      MT2 = MT**2               ! mu_f^2 for top threshold
      LOGFR = LOG(FR2)          ! log of ratio of mu_f^2 to mu_r^2
      M20 = R20 * FR2           ! input factorisation scale

*
* ..Stop some nonsense
*
      IF ( (IVFNS .EQ. 0) .AND. (NFF .LT. 3) ) THEN
         WRITE (6,*) 'Wrong flavour number for FFNS evolution. STOP'
         STOP
      END IF
      IF ( (IVFNS .EQ. 0) .AND. (NFF .GT. 5) ) THEN
         WRITE (6,*) 'Wrong flavour number for FFNS evolution. STOP'
         STOP
      END IF
*     
      IF ( NAORD .GT. 3 ) THEN
         WRITE (6,*) 'Specified order in a_s too high. STOP' 
         STOP
      END IF
*
      IF ( (IVFNS .NE. 0) .AND. (FR2 .GT. 4.001D0) ) THEN
         WRITE (6,*) 'Too low mu_r for VFNS evolution. STOP'
         STOP
      END IF
*
      IF ( (IVFNS .EQ. 1) .AND. (M20 .GT. MC2) ) THEN
         WRITE (6,*) 'Too high mu_0 for VFNS evolution. STOP'
         STOP
      END IF
*     
      IF ( (ASI .GT. 2.D0) .OR. (ASI .LT. 2.D-2) ) THEN
         WRITE (6,*) 'alpha_s out of range. STOP'
         STOP
      END IF
*     
      IF ( (IVFNS .EQ. 1) .AND. (MC2 .GT. MB2) ) THEN
         WRITE (6,*) 'Wrong charm-bottom mass hierarchy. STOP'
         STOP
      END IF
      IF ( (IVFNS .EQ. 1) .AND. (MB2 .GT. MT2) ) THEN
         WRITE (6,*) 'Wrong bottom-top mass hierarchy. STOP'
         STOP
      END IF
*

C--   Store the beta function coefficients in a COMMON block.
      CALL BETAFCT

C--   Store a_s = alpha_s(mu_r^2)/(4 pi) at the input scale R0.
      AS0 = ASI / (4.D0* PI)

C--   Store alpha_s at the heavy flavour thresholds in a COMMON block.
       IF (IVFNS .NE. 0) THEN
          CALL EVNFTHR (MC2, MB2, MT2)
       END IF

      RETURN
      END

C----------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION ALPHASRT(MUR)
      IMPLICIT NONE
      INTEGER NFF,IVFNS,NF
      DOUBLE PRECISION PI,LOGFR,AS0,M20,ASC,M2C,ASB,M2B,AST,M2T,M2,MUR,
     &     R2SAVE,ASSAVE,R2,ASI,ASF,R20,R2T,R2B,R2C,AS
      PARAMETER ( PI = 3.1415 92653 58979 D0 )
*
* ..Input common blocks 
* 
       COMMON / NFFIX  / NFF
       COMMON / VARFLV / IVFNS 
       COMMON / FRRAT  / LOGFR
       COMMON / ASINP  / AS0, M20
       COMMON / ASFTHR / ASC, M2C, ASB, M2B, AST, M2T
*
* ..Variables to be saved for the next call
*
       SAVE R2SAVE, ASSAVE
*
* ..Check for new scale mu_r^2
*
       R2 = MUR**2
       IF ( R2.EQ.R2SAVE ) THEN
          ALPHASRT = ASSAVE
          RETURN
       END IF

       M2 = R2 * EXP(+LOGFR)
       IF (IVFNS .EQ. 0) THEN
*
*   Fixed number of flavours
*
          NF  = NFF
          R20 = M20 * R2/M2
          ASI = AS0
          ASF = AS (R2, R20, AS0, NF)
*
       ELSE
*
* ..Variable number of flavours
*
          IF (M2 .GT. M2T) THEN
             NF = 6
             R2T = M2T * R2/M2
             ASI = AST
             ASF = AS (R2, R2T, AST, NF)
*
          ELSE IF (M2 .GT. M2B) THEN
             NF = 5
             R2B = M2B * R2/M2
             ASI = ASB
             ASF = AS (R2, R2B, ASB, NF)
*     
          ELSE IF (M2 .GT. M2C) THEN
             NF = 4
             R2C = M2C * R2/M2
             ASI = ASC
             ASF = AS (R2, R2C, ASC, NF)
*     
          ELSE
             NF = 3
             R20 = M20 * R2/M2
             ASI = AS0
             ASF = AS (R2, R20, AS0, NF)
*       
          END IF
*
       END IF
*
* ..Final value of alpha_s
*
       ALPHASRT = 4.D0*PI*ASF
*
* ..Save variables for the next call
*
       R2SAVE = R2
       ASSAVE = ALPHASRT
*
       RETURN
       END
*
* =================================================================av==


* =====================================================================
*
* ..The threshold matching of the QCD coupling in the MS(bar) scheme,  
*    a_s = alpha_s(mu_r^2)/(4 pi),  for  NF -> NF + 1  active flavours 
*    up to order a_s^4 (NNNLO).
*
* ..The value  ASNF  of a_s for NF flavours at the matching scale, the 
*    logarithm  LOGRH = ln (mu_r^2/m_H^2) -- where m_H is the pole mass
*    of the heavy quark -- and  NF  are passed as arguments to the 
*    function  ASNF1.  The order of the expansion  NAORD  (defined as 
*    the 'n' in N^nLO) is provided by the common-block  ASPAR.
*
* ..The matching coefficients are inverted from Chetyrkin, Kniehl and
*    Steinhauser, Phys. Rev. Lett. 79 (1997) 2184. The QCD colour
*    factors have been hard-wired in these results. The lowest integer 
*    values of the Zeta function are given by the common-block  RZETA.
*
* =====================================================================
*
*
      DOUBLE PRECISION FUNCTION ASNF1 (ASNF, LOGRH, NF)
*
      IMPLICIT NONE
      INTEGER NF, NAORD, NASTPS, PRVCLL, K1, K2
      DOUBLE PRECISION ASNF,LOGRH,ZETA,CMC,CMCI30,CMCF30,CMCF31,
     &     CMCI31,ASP,LRHP

      DIMENSION CMC(3,0:3)
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks 
*
      COMMON / ASPAR  / NAORD, NASTPS
      COMMON / RZETA  / ZETA(6)
*
* ..Variables to be saved for the next call
*
      SAVE CMC, CMCI30, CMCF30, CMCF31, CMCI31, PRVCLL
*
* ---------------------------------------------------------------------
*
* ..The coupling-constant matching coefficients (CMC's) up to NNNLO 
*   (calculated and saved in the first call of this routine)
*
       IF (PRVCLL .NE. 1) THEN
*
         CMC(1,0) =  0.D0
         CMC(1,1) =  2./3.D0
*
         CMC(2,0) = 14./3.D0
         CMC(2,1) = 38./3.D0
         CMC(2,2) =  4./9.D0  
*
         CMCI30 = + 80507./432.D0 * ZETA(3) + 58933./1944.D0 
     1            + 128./3.D0 * ZETA(2) * (1.+ DLOG(2.D0)/3.D0)
         CMCF30 = - 64./9.D0 * (ZETA(2) + 2479./3456.D0)
         CMCI31 =   8941./27.D0
         CMCF31 = - 409./27.D0
         CMC(3,2) = 511./9.D0
         CMC(3,3) = 8./27.D0
*
         PRVCLL = 1
*
       END IF
*
* ---------------------------------------------------------------------
*
* ..The N_f dependent CMC's, and the alpha_s matching at order NAORD 
*
       CMC(3,0) = CMCI30 + NF * CMCF30
       CMC(3,1) = CMCI31 + NF * CMCF31
*
       ASNF1 = ASNF
       IF (NAORD .EQ. 0) GO TO 1
       ASP   = ASNF
*
       DO 11 K1 = 1, NAORD 
         ASP = ASP * ASNF
         LRHP = 1.D0
*
       DO 12 K2 = 0, K1
         ASNF1 = ASNF1 + ASP * CMC(K1,K2) * LRHP
         LRHP = LRHP * LOGRH
*
  12   CONTINUE
  11   CONTINUE
*
* ---------------------------------------------------------------------
*
   1   RETURN
       END

*
* =================================================================av==
*
* ..The subroutine  EVNFTHR  for the evolution of  a_s = alpha_s/(4 pi)
*    from a three-flavour initial scale to the four- to six-flavour
*    thresholds (identified with the squares of the corresponding quark
*    masses).  The results are written to the common-block  ASFTHR.
*
* ..The input scale  M20 = mu_(f,0)^2  and the corresponding value 
*    AS0  of a_s  are provided by  ASINP.  The fixed scale logarithm
*    LOGFR = ln (mu_f^2/mu_r^2) is specified in  FRRAT.  The alpha_s
*    matching is done by the function ASNF1.
*
* =====================================================================
*
*
       SUBROUTINE EVNFTHR (MC2, MB2, MT2)
*
       IMPLICIT NONE
       DOUBLE PRECISION MC2, MB2, MT2, M20, M2C, M2B, M2T, R20, R2C, 
     1                  R2B, R2T, AS, ASNF1, AS0, ASC, ASB, AST,
     2                  ASC3, ASB4, AST5, LOGFR, SC, SB, ST
*
* ---------------------------------------------------------------------
* 
* ..Input common blocks
*  
       COMMON / ASINP  / AS0, M20
       COMMON / FRRAT  / LOGFR
*
* ..Output common blocks
*
       COMMON / ASFTHR / ASC, M2C, ASB, M2B, AST, M2T

* ---------------------------------------------------------------------
*
* ..Coupling constants at and evolution distances to/between thresholds
* 
       R20 = M20 * EXP(-LOGFR)
*
* ..Charm
*
       M2C  = MC2
       R2C  = M2C * R20/M20
       ASC3 = AS (R2C, R20, AS0, 3)
       SC   = LOG (AS0 / ASC3)
       ASC  = ASNF1 (ASC3, -LOGFR, 3)
*
* ..Bottom 
*
       M2B  = MB2
       R2B  = M2B * R20/M20
       ASB4 = AS (R2B, R2C, ASC, 4)
       SB   = LOG (ASC / ASB4)
       ASB  = ASNF1 (ASB4, -LOGFR, 4)
*
* ..Top
*
       M2T  = MT2
       R2T  = M2T * R20/M20
       AST5 = AS (R2T, R2B, ASB, 5)
       ST   = LOG (ASB / AST5)
       AST  = ASNF1 (AST5, -LOGFR, 5)

       RETURN
       END

*
* =================================================================av==
*
* ..The running coupling of QCD,  
*
*         AS  =  a_s  =  alpha_s(mu_r^2)/(4 pi),
*
*    obtained by integrating the evolution equation for a fixed number
*    of massless flavours  NF.  Except at leading order (LO),  AS  is 
*    obtained using a fourth-order Runge-Kutta integration.
*
* ..The initial and final scales  R20  and  R2,  the value  AS0  at
*    R20, and  NF  are passed as function arguments.  The coefficients 
*    of the beta function up to  a_s^5 (N^3LO)  are provided by the 
*    common-block  BETACOM.  The order of the expansion  NAORD (defined
*    as the 'n' in N^nLO) and the number of steps  NASTPS  for the 
*    integration beyond LO are given by the common-block  ASPAR.
*
* =====================================================================
*
*
      DOUBLE PRECISION FUNCTION AS (R2, R20, AS0, NF)
*
      IMPLICIT NONE
      INTEGER NFMIN, NFMAX, NF, NAORD, NASTPS, K1
      DOUBLE PRECISION R2, R20, AS0, SXTH, BETA0, BETA1, BETA2, BETA3,
     &     FBETA1,FBETA2,FBETA3,A,LRRAT,DLR,XK0,XK1,XK2,XK3
      PARAMETER (NFMIN = 3, NFMAX = 6)
      PARAMETER ( SXTH = 0.16666 66666 66666 D0 )
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks 
*
       COMMON / ASPAR  / NAORD, NASTPS
       COMMON / BETACOM   / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     ,                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
*
* ..The beta functions FBETAn at N^nLO for n = 1, 2, and 3
*
       FBETA1(A) = - A**2 * ( BETA0(NF) + A *   BETA1(NF) )
       FBETA2(A) = - A**2 * ( BETA0(NF) + A * ( BETA1(NF)
     ,                        + A * BETA2(NF) ) )
       FBETA3(A) = - A**2 * ( BETA0(NF) + A * ( BETA1(NF)
     ,                        + A * (BETA2(NF) + A * BETA3(NF)) ) )
*
* ---------------------------------------------------------------------
*
* ..Initial value, evolution distance and step size
*
       AS = AS0
       LRRAT = LOG (R2/R20)
       DLR = LRRAT / NASTPS
*
* ..Solution of the evolution equation depending on  NAORD
*   (fourth-order Runge-Kutta beyond the leading order)
*
       IF (NAORD .EQ. 0) THEN
*
         AS = AS0 / (1.+ BETA0(NF) * AS0 * LRRAT)
*
       ELSE IF (NAORD .EQ. 1) THEN
*
       DO 2 K1 = 1, NASTPS
         XK0 = DLR * FBETA1 (AS)
         XK1 = DLR * FBETA1 (AS + 0.5 * XK0)
         XK2 = DLR * FBETA1 (AS + 0.5 * XK1)
         XK3 = DLR * FBETA1 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  2    CONTINUE
*
       ELSE IF (NAORD .EQ. 2) THEN
*
       DO 3 K1 = 1, NASTPS
         XK0 = DLR * FBETA2 (AS)
         XK1 = DLR * FBETA2 (AS + 0.5 * XK0)
         XK2 = DLR * FBETA2 (AS + 0.5 * XK1)
         XK3 = DLR * FBETA2 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  3    CONTINUE
*  
       ELSE IF (NAORD .EQ. 3) THEN
*
       DO 4 K1 = 1, NASTPS
         XK0 = DLR * FBETA3 (AS)
         XK1 = DLR * FBETA3 (AS + 0.5 * XK0)
         XK2 = DLR * FBETA3 (AS + 0.5 * XK1)
         XK3 = DLR * FBETA3 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  4    CONTINUE
       END IF
*
* ---------------------------------------------------------------------
*
       RETURN
       END

*
* =================================================================av==
*
* ..The subroutine BETAFCT for the coefficients  BETA0...BETA3  of the 
*    beta function of QCD up to order alpha_s^5 (N^3LO), normalized by 
*
*        d a_s / d ln mu_r^2  =  - BETA0 a_s^2 - BETA1 a_s^3 - ... 
*
*    with  a_s = alpha_s/(4*pi). 
*
* ..The MSbar coefficients are written to the common-block BETACOM for 
*   NF = 3...6  (parameters NFMIN, NFMAX) quark flavours.
*
* ..The factors CF, CA and TF  are taken from the common-block  COLOUR.
*    Beyond NLO the QCD colour factors are hard-wired in this routine,
*    and the numerical coefficients are truncated to six digits.
*
* =====================================================================
*
*
       SUBROUTINE BETAFCT
*
       IMPLICIT DOUBLE PRECISION (A - Z)
       INTEGER NFMIN, NFMAX, NF
       PARAMETER (NFMIN = 3, NFMAX = 6)
*
* ---------------------------------------------------------------------
*
* ..Input common-block
*
       COMMON / COLOUR / CF, CA, TR
*
* ..Output common-block
*
       COMMON / BETACOM   / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     1                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)

*
* ---------------------------------------------------------------------
*
* ..The full LO and NLO coefficients 
*
       B00 =  11./3.D0 * CA
       B01 =  -4./3.D0 * TR
       B10 =  34./3.D0 * CA**2
       B11 = -20./3.D0 * CA*TR - 4.* CF*TR
*
* ..Flavour-number loop and output to the array
*
       DO 1 NF = NFMIN, NFMAX
*
       BETA0(NF) = B00 + B01 * NF
       BETA1(NF) = B10 + B11 * NF
*
       BETA2(NF) = 1428.50 - 279.611 * NF + 6.01852 * NF**2
       BETA3(NF) = 29243.0 - 6946.30 * NF + 405.089 * NF**2 
     1             + 1.49931 * NF**3
*
* ---------------------------------------------------------------------
*
  1    CONTINUE
*
       RETURN
       END
*
* =================================================================av==
