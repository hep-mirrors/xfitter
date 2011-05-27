      integer NDATASETS, NDATAPOINTS,DATASETIDX,NQ2BINS,NXBINS
      integer DATASETNUMBER
C SG: add DYZ and DYW
      logical ISNC,ISEPLUS,ISJET,ISCHARM,ISBOTTOM,ISDYZ,ISDYW
      logical ISCHARMANY,ISBOTTOMANY
      double precision POLAR,EPOL
      character*80 DATASETLABEL
C SG: reaction type:
      character*80 DATASETReaction
      
C SG: relation to bins:
      integer DATASETBinningDimension
      character*80 DATASETBinNames 
      character*80 DATASETInfoNames
      integer DATASETInfoDimension
      
      double precision DATASETInfo

      INTEGER MAXQ2BINS
      parameter(MAXQ2BINS=21)

      integer NInfoDimensionMax
      parameter (NInfoDimensionMax=10)

      integer NInputFiles
      character*80 InputFileNames

C Dimension of k-factors
      integer DATASETNKfactors
      character*80 DATASETKFactorNames

      common/datasets/
     $     DATASETInfo(NInfoDimensionMax,NSET)
     $     ,DATASETKFactorNames(NKFactMax,NSET)
     $     ,DATASETLABEL(NSET),POLAR(NSET),EPOL(NSET),
     +     NDATAPOINTS(NSET),NQ2BINS(NSET),DATASETIDX(NSET,NMAX),
     +     NXBINS(NSET,MAXQ2BINS),DATASETNUMBER(NSET),
     +     ISNC(NSET),ISEPLUS(NSET),ISJET(NSET),
     +     ISCHARM(NSET),ISBOTTOM(NSET),NDATASETS,
     +     ISCHARMANY,ISBOTTOMANY,ISDYZ(NSET),ISDYW(NSET),
     $     DATASETREACTION(NSET),DATASETBinNames(NBinDimensionMax,NSET)
     $     ,DATASETInfoNames(NInfoDimensionMax,NSET)
     $     ,DATASETInfoDimension(NSET)
     $     ,DATASETBinningDimension(NSET)
     $     ,InputFileNames(NSET),DATASETNKfactors(NSET),NInputFiles

