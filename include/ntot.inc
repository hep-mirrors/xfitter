C
C Basic definitions for 
C

      integer NTOT        !> Total number of data points
      parameter(NTOT = 2500)


      integer NCovarMax   !> Total number of data points, treated by covariance matrix
      parameter (NCovarMax = 1000)

      integer nset            !> Total number of datasets   
      parameter (NSET = 150)

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