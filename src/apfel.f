************************************************************************
*
*     Initialization routine for APFEL
*
************************************************************************
      subroutine apfel_ini
*
      implicit none
*
#include "steering.inc"
#include "couplings.inc"
#include "thresholds.inc"
#include "extrapars.inc"
*
      integer GetParameterIndex
      double precision alphas
*
*     Reference value of alphas taken from the extraparameters
*
      alphas = ExtraParamValue(GetParameterIndex('alphas'))
*
*     Define basic settings
*
      if(iTheory.eq.35)then
         call SetTheory("QUniD")                 ! Set QCD+QED evolution (default)
         call SetPDFEvolution("exactalpha")      ! Use DGLAP evolution in terms of muF
      else
         call SetTheory("QCD")                   ! Set QCD evolution (default)
         call SetPDFEvolution("exactalpha")      ! Use DGLAP evolution in terms of alphas (rather than muF => faster for short steps)
      endif
      call SetFastEvolution(.true.)              ! Use fast evolution (default)
      call SetAlphaEvolution("exact")            ! Use exact solution on the beta functions (default)
      call SetQLimits(0.5d0,20000d0)             ! Evolution limits
      call SetNumberOfGrids(3)                   ! x-space grid settings
      call SetGridParameters(1,40,3,9.8d-7)
      call SetGridParameters(2,30,3,1d-2)
      call SetGridParameters(3,20,3,7d-1)
      call LockGrids(.true.)                     ! Lock subgrids
      call SetPerturbativeOrder(I_FIT_ORDER-1)   ! Set perturbative order
      if(HF_SCHEME(1:2).eq."FF")then             ! Set mass scheme for the PDF evolution
         call SetFFNS(3)
      else
         call SetVFNS
      endif
      call SetAlphaqcdRef(alphas,Mz)
      if(HF_SCHEME(9:12).eq."RUNM")then
         call SetMSbarMasses(mch,mbt,mtp) ! Heavy-quark thresholds in the MSbar scheme
         call EnableMassRunning(.false.)
         if(HF_SCHEME(14:15).eq."ON") call EnableMassRunning(.true.)
      else
         call SetPoleMasses(mch,mbt,mtp)  ! Heavy-quark thresholds in the Pole scheme
      endif
      call SetMassMatchingScales(kmuc,kmub,kmut)
*
*     Initialize APFEL
*
      call InitializeAPFEL
*
      return
      end
