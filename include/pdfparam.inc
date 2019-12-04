 !> PDF parameterisations


C !> Standard parameterisation

       double precision paruval(1:10), pardval(1:10)
     $     ,parubar(1:10), pardbar(1:10)
     $     ,parglue(1:10),paru(1:10), pard(1:10)
     $     ,parstr(1:10),parsea(1:10),parother(1:10)
     $     ,parphoton(1:10)
       common/stpar/paruval,pardval,parubar,pardbar,
     $      parglue,paru,pard,parstr,parsea,parother,parphoton

C---
       double precision fcharm,fstrange !> charm, strange fractions
       double precision kappas		!> Marina Walt, nPDFs s+sbar=kappas*(ubar+dbar)
       common/pdfparam/fcharm,fstrange,kappas
    
C
C !>  Chebyshev polynoms for the gluon, sea
C
      integer nchebmax
      parameter (nchebmax=100)
C Gluon:
      double precision chebpars(nchebmax),polypars(nchebmax)
C Sea, for ZEUS param. only:
      double precision chebparsSea(nchebmax),polyparsSea(nchebmax)
      double precision chebxminlog


       common/pdfcheb/chebpars,polypars,chebxminlog,
     $     chebparsSea,polyparsSea


C !> 25 Jan 2011: add poly param. for valence quarks
       integer MaxPolyVal
       parameter (MaxPolyVal=50)

C !> internal order of the polynom:
       integer NPOLYVALINT

       double precision PolyUval(MaxPolyVal),PolyDval(MaxPolyVal)
       common/pdfpolyval/PolyUval,PolyDval,NPOLYVALINT

C !> CTEQ-like parameterisation
       
       integer idxADataSet, ratiostep          !> 15/03/2017 Marina Walt, new parameters
       double precision ratio
  
       !> 15/03/2017 dimension of ctparams changed to 27 - to reflect parameters for nuclear PDFs

       double precision ctuval(1:27), ctdval(1:27),ctubar(1:27),ctdbar(1:27), 
     $      ctglue(1:27),ctother(1:27),ctstr(1:27), ctphoton(1:27)
       common/ctpar/ctuval,ctdval,ctubar,ctdbar,ctglue,ctother,ctstr,ctphoton,
     $      idxADataSet, ratiostep, ratio

C !> Add AS parameterisation
       double precision asuval(1:5), asdval(1:5), asubar(1:5), asdbar(1:5),
     $      asglue(1:5),asother(1:5)
       common/aspar/asuval,asdval,asubar,asdbar,
     $      asglue,asother

C !> Results of sum-rule integrals
       double precision uv_sum, dv_sum, p_sum
       common/csumrules/uv_sum, dv_sum, p_sum

C !> Temperature paramter
       double precision Temperature
       Common/CTemperature/Temperature
