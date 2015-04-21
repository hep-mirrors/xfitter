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
#include "alphas.inc"
#include "couplings.inc"
*
*     Define basic settings
*
      call SetTheory("QCD")                      ! Set QCD evolution (default)
      call SetFastEvolution(.true.)              ! Use fast evolution (default)
      call SetAlphaEvolution("exact")            ! Use exact solution on the beta functions (default)
      call SetPDFEvolution("exactalpha")         ! Use DGLAP evolution in terms of alphas (rather than muF => faster for short steps)
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
      if(HF_SCHEME(9:12).eq."RUNM")then
         call SetMSbarMasses(dble(HF_MASS(1)),dble(HF_MASS(2)), ! Heavy-quark thresholds in the MSbar scheme
     1                       dble(HF_MASS(3)))
         call EnableMassRunning(.false.)
         if(HF_SCHEME(14:15).eq."ON") call EnableMassRunning(.true.)
      else
         call SetPoleMasses(dble(HF_MASS(1)),dble(HF_MASS(2)),  ! Heavy-quark thresholds in the Pole scheme
     1                      dble(HF_MASS(3)))
      endif
*
*     Initialize APFEL
*
      call InitializeAPFEL
*
      return
      end
