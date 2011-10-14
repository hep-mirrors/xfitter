* modified by JF to fit a la ZEUS
* added fg for SPECIAL TEST

* added for negative gluon
        double precision
     +  fcharm,fstrange,
     +  ag,bg,cg,dg,eg,fg, apg,bpg,cpg,
     +  aU,bU,cU,dU,eU,fU,
     +  aD,bD,cD,dD,eD,fD,
     +  aUbar,bUbar,cUbar,dUbar,
     +  aDbar,bDbar,cDbar,dDbar,
     +  aUv,bUv,cUv,dUv,eUv,fUv,
     +  aDv,bDv,cDv,dDv,eDv,fDv,
     +  Asea,Bsea,Csea,Dsea,
     +  Adel,Bdel,Cdel,Ddel, Astr, Bstr,Cstr

    
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
     +  ag,bg,cg,dg,eg,fg,apg,bpg,cpg,
     +  aU,bU,cU,dU,eU,fU,
     +  aD,bD,cD,dD,eD,fD,
     +  aUbar,bUbar,cUbar,dUbar,
     +  aDbar,bDbar,cDbar,dDbar,
     +  aUv,bUv,cUv,dUv,eUv,fUv,
     +  aDv,bDv,cDv,dDv,eDv,fDv,
     +  Asea,Bsea,Csea,Dsea,
     +  Adel,Bdel,Cdel,Ddel,Astr,Bstr,Cstr



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
       double precision ctuval(1:6), ctdval(1:6), ctubar(1:6), ctdbar(1:6), 
     $      ctglue(1:6),ctother(1:6)
       common/ctpar/ctuval,ctdval,ctubar,ctdbar,ctglue,ctother

C VR: Add standard parameterisation
       double precision paruval(1:10), pardval(1:10), parubar(1:10), pardbar(1:10), 
     $      parglue(1:10),paru(1:10), pard(1:10),pardel(1:10),parsea(1:10),parother(1:10)
       common/stpar/paruval,pardval,parubar,pardbar,
     $      parglue,paru,pard,pardel,parsea,parother

