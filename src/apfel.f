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
      if(HFSCHEME.eq.3)then                      ! Set mass scheme for the PDF evolution
         call SetFFNS(3)
      else
         call SetVFNS
      endif
      call SetPoleMasses(dble(HF_MASS(1)),dble(HF_MASS(2)),  ! Heavy-quark thresholds in the Pole scheme
     1                   dble(HF_MASS(3)))
c      call SetMSbarMasses(dble(HF_MASS(1)),dble(HF_MASS(2)),  ! Heavy-quark thresholds in the MSbar scheme
c     1                    dble(HF_MASS(3)))
*
*     Initialize APFEL
*
      call InitializeAPFEL
*
      return
      end
