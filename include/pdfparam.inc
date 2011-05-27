* modified by JF to fit a la ZEUS
* added fg for SPECIAL TEST

* added for negative gluon
        double precision
     +  fcharm,fstrange,
     +  ag,bg,cg,dg,fg, apg,bpg,cpg,
     +  aU,bU,cU,dU,eU,fU,
     +  aD,bD,cD,dD,fD,
     +  aUbar,bUbar,cUbar,dUbar,
     +  aDbar,bDbar,cDbar,dDbar,
     +  aUv,bUv,cUv,dUv,eUv,fUv,
     +  aDv,bDv,cDv,dDv,fDv,
     +  Asea,Bsea,Csea,Dsea,
     +  Adel,Bdel,Cdel

C VR:  Add HT + Fudge FL parameters
        double precision
     +  f2ht1,f2ht2,afudge,Rfudge
 	
C VR:  Add CI (hayk)
        double precision
     +  etaull,etadll,etaulr,etadlr,
     +  etaurr,etadrr,etaurl,etadrl

    
C
C SG:  Add Chebyshev polynoms for the gluon:
C
      integer nchebmax
      parameter (nchebmax=100)
C Gluon:
      double precision chebpars(nchebmax),polypars(nchebmax)
C Sea, for ZEUS param. only:
      double precision chebparsSea(nchebmax),polyparsSea(nchebmax)
      double precision chebxminlog

       common/pdfparam/
     +  fcharm,fstrange,
     +  ag,bg,cg,dg,fg,apg,bpg,cpg,
     +  aU,bU,cU,dU,eU,fU,
     +  aD,bD,cD,dD,fD,
     +  aUbar,bUbar,cUbar,dUbar,
     +  aDbar,bDbar,cDbar,dDbar,
     +  aUv,bUv,cUv,dUv,eUv,fUv,
     +  aDv,bDv,cDv,dDv,fDv,
     +  Asea,Bsea,Csea,Dsea,
     +  Adel,Bdel,Cdel,
     +  f2ht1,f2ht2,afudge, Rfudge


C VR:  Add CI (hayk)
	
       common/ciparam/
     +  etaull,etadll,etaulr,etadlr,
     +  etaurr,etadrr,etaurl,etadrl

	

       common/pdfcheb/chebpars,polypars,chebxminlog,
     $     chebparsSea,polyparsSea


C SG, 25 Jan 2011: add poly param. for valence quarks
       integer MaxPolyVal
       parameter (MaxPolyVal=50)

C SG: internal order of the polynom:
       integer NPOLYVALINT

       double precision PolyUval(MaxPolyVal),PolyDval(MaxPolyVal)
       common/pdfpolyval/PolyUval,PolyDval,NPOLYVALINT

C SG: Add cteq-like parameterisation
       double precision ctuval(0:5), ctdval(0:5), ctubar(0:5), ctdbar(0:5), 
     $      ctglue(0:5)
       common/ctpar/ctuval,ctdval,ctubar,ctdbar,ctglue
