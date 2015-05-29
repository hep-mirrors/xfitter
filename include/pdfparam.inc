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
       common/pdfparam/fcharm,fstrange
    
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
       double precision ctuval(1:9), ctdval(1:9),ctubar(1:9),ctdbar(1:9), 
     $      ctglue(1:9),ctother(1:9),ctstr(1:9), ctphoton(1:9)
       common/ctpar/ctuval,ctdval,ctubar,ctdbar,ctglue,ctother,ctstr,ctphoton

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
