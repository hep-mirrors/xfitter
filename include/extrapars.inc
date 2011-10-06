!> @brief Extra minimization parameters to be passed to Minuit
      integer nExtraParamMax   !> Maximum number of parameters
      parameter (nExtraParamMax=50)

      integer nExtraParam      !> Actual number of parameters
      character*32 ExtraParamNames(nExtraParamMax)     !> Names of extra pars
      double precision ExtraParamValue(nExtraParamMax) !> Initial values
      double precision ExtraParamStep (nExtraParamMax) !> Initial step size
      double precision ExtraParamMin  (nExtraParamMax) !> Min value 
      double precision ExtraParamMax  (nExtraParamMax) !> Max value
      integer iExtraParamMinuit       (nExtraParamMax) !> Minuit param. index

      common/ExtraPars/ExtraParamNames,ExtraParamValue,ExtraParamStep,
     $     ExtraParamMin,ExtraParamMax,iExtraParamMinuit,nExtraParam
