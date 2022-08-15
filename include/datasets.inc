!>
!>  Declarations for datasets
!>

C Namelist for data set declarations
      integer NInputFiles                 !> Number of input files
      character*256 InputFileNames(NSET)  !> Names of input files

C Namelist for statistical correlation declarations
      integer NCorrFiles                    !> Number of correlation files
      character*256 CorrFileNames(NSET**2)   !> Names of correlation files

      logical UseFixedTheory(NSET)          !> Used fixed theory predictions

C General dataset info
      integer 
     $      NDATASETS                 !> Number of datasets 
     $     ,NDATAPOINTS(NSET)         !> Number of data points per dataset
     $     ,DATASETIDX(NSET,NTOT)     !> Index from dataset point to global data array
      integer DATASETNUMBER(NSET)     !> xFitter dataset number

      character*80 DATASETLABEL(NSET) !> Dataset name
C SG: reaction type:
      character*80 DATASETReaction(NSET)    !> Dataset reaction type
C SG: extra theory file:
      integer nchar_theory
      parameter (nchar_theory=200)

      character*(nchar_theory) DATASETTheoryFile(NSet)  !> File containing extra info for theory predictions
      character*80 DATASETTheoryType(NSET)  !> Theory type. Currently: kfactor, applgrid, fastnlo    
      integer      DATASETTheoryIndex(NSET) !> Theory index: utility for FastNLO/applgrid interface

C Dimension of k-factors
      integer DATASETNKfactors(NSET)                      !> Number of k-factors
      character*80 DATASETKFactorNames(NKFactMax,NSET)    !> Names of k-factors

C  Try to read two or more applgrid
      integer DATASETNapplgrid(NSET)                      !> Number of applgrids
      character*80 DATASETapplgridNames(NapplgridMax,NSET)    !> Names of applgrids
      integer DATASETApplgridTheoryIndex(NapplgridMax,NSET)

C SG: relation to bins:
      integer DATASETBinningDimension(NSET)                        !> Number of abstract bins
      character*80 DATASETBinNames(NBinDimensionMax,NSET)          !> Names of abstract bins

C KN: For plotting
      character*64 GPlotVarCol(NSET)
      character*256 GPlotOptions(NPlotMax)
      integer GNPlots(NSET)

C SG: Extra info:
C
      integer    NInfoDimensionMax                                 !> Max. number of info entries per dataset
      parameter (NInfoDimensionMax=10)
C
      integer DATASETInfoDimension(NSET)                           !> Number of info entries
      character*80 DATASETInfoNames(NInfoDimensionMax,NSET)        !> Names of info entries
      double precision DATASETInfo(NInfoDimensionMax,NSET)         !> info values.

      double precision chi2_poi_tot
      double precision chi2_poi(NSET)
      double precision chi2_poi_data(NTOT)

C VB: For the hybrid VFNS
      integer DATASETMaxNF(NSET)                                   !> Maximum number of active flavours per dataset
      double precision DATASETSwitchScales(4:6,NSET)               !> Charm, bottom, and top switching scales per dataset
      logical UseHVFNS                                             !> Whether the H-VFNS has to be used

      common/datasets/
     $     DATASETInfo
     $ 	   ,DATASETSwitchScales
     $     ,DATASETKFactorNames
     $     ,DATASETTheoryFile
     $     ,DATASETTheoryType
     $     ,DATASETLABEL
     +     ,NDATAPOINTS,DATASETIDX
     +     ,DATASETNUMBER
     $     ,NDATASETS
     $     ,DATASETREACTION
     $     ,DATASETBinNames
     $     ,DATASETInfoNames
     $     ,DATASETInfoDimension
     $     ,DATASETBinningDimension
     $     ,NInputFiles
     $     ,InputFileNames
     $     ,DATASETNKfactors
     $     ,DATASETMaxNF
     $     ,DATASETTheoryIndex
     $     ,NCorrFiles
     $     ,CorrFileNames
     $     ,DATASETNapplgrid
     $     ,DATASETapplgridNames
     $     ,DATASETApplgridTheoryIndex
     $     ,GPlotVarCol
     $     ,GNPlots
     $     ,GPlotOptions
     $     ,UseFixedTheory
     $     ,UseHVFNS

      Common/cdatapoi/chi2_poi_tot,chi2_poi,chi2_poi_data

C
C Keep track of uncertainty types per dataset
C
      integer ibConst,ibPoisson,ibLinear
      parameter (ibConst   = 1)
      parameter (ibPoisson = 2)
      parameter (ibLinear  = 4)

      integer iStatTypesBitMask(NSET), iUncorTypesBitMask(NSET)
      Common/CdatasetErrorTypes/iStatTypesBitMask, iUncorTypesBitMask
