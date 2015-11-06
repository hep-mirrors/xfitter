!> @brief Extra minimization parameters to be passed to Minuit
#include "dimensions.h"
      integer nExtraParamMax   !> Maximum number of parameters
      parameter (nExtraParamMax=NEXTRAPARAMMAX_C)

      integer nExtraParam      !> Actual number of parameters
      character*32 ExtraParamNames(nExtraParamMax)     !> Names of extra pars
      double precision ExtraParamValue(nExtraParamMax) !> Initial values
      double precision ExtraParamStep (nExtraParamMax) !> Initial step size
      double precision ExtraParamMin  (nExtraParamMax) !> Min value 
      double precision ExtraParamMax  (nExtraParamMax) !> Max value
      double precision ExtraParamConstrVal(nExtraParamMax) !> External constraint (central value)
      double precision ExtraParamConstrUnc(nExtraParamMax) !> External constraint (uncertainty)
      integer iExtraParamMinuit       (nExtraParamMax) !> Minuit param. index

      common/ExtraPars/ExtraParamNames,ExtraParamValue,ExtraParamStep,
     $     ExtraParamMin,ExtraParamMax,ExtraParamConstrVal,
     $     ExtraParamConstrUnc,iExtraParamMinuit,nExtraParam
