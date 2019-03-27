#include "dimensions.h"

      integer NTOT        !> Maximal number of data points
      parameter(NTOT = NTOT_C)

      integer NCovarMax   !> Maximal number of data points, treated by covariance matrix
      parameter (NCovarMax = 1000)

      integer nset            !> Maximal number of datasets
      parameter (NSET = NSET_C)

      integer NKFactMax          !> Max. number of k-factors
      parameter (NKFactMax = 10)

      integer NapplgridMax          !> Max. number of applgrids
      parameter (NapplgridMax = 10)

      integer NBinDimensionMax   !> MAx. number of abstract bins
      parameter (NBinDimensionMax=10)

      integer NPlotMax   !> Max. number of plots to be prepared
      parameter (NPlotMax = 500)
      integer NXGridMax
      parameter (NXGridMax = 500)
