c=====================================================
      subroutine DDIS_FixModelParams(p)
C-
C (to be called fron fcn. parsminuit is a vector of minuit parameters.
C
      Implicit none 
#include "extrapars.inc"

      integer GetParameterIndex
      double precision p(*)
      double precision pars(3)
      integer idxPomeron_a0,idxReggeon_factor,idxReggeon_a0

      idxPomeron_a0 = GetParameterIndex('Pomeron_a0')
      if ( idxPomeron_a0.eq.0) then
         print *,'Did not find Pomeron_a0 parameter'
         print *,'Add to ExtraParamters with the name  Pomeron_a0'
         call HF_stop
      else
         idxPomeron_a0 = iExtraParamMinuit(idxPomeron_a0)
      endif
      pars(1)=p(idxPomeron_a0) 
      
      idxReggeon_factor = GetParameterIndex('Reggeon_factor')
      if ( idxReggeon_factor.eq.0) then
         print *,'Did not find Reggeon_factor parameter'
         print *,'Add to ExtraParamters with the name Reggeon_factor'
            call HF_stop
      else
         idxReggeon_factor = iExtraParamMinuit(idxReggeon_factor)
      endif
      pars(2)=p(idxReggeon_factor) 


      idxReggeon_a0 = GetParameterIndex('Reggeon_a0')
      if ( idxReggeon_a0.eq.0) then
         print *,'Did not find Reggeon_a0 parameter'
         print *,'Add to ExtraParamters with the name Reggeon_a0'
            call HF_stop
      else
         idxReggeon_a0 = iExtraParamMinuit(idxReggeon_a0)
      endif
      pars(3)=p(idxReggeon_a0) 

      call DDISsetParams(pars)
      
      end
