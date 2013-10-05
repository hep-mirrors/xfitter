
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
