
      integer NBANDS
      parameter (NBANDS=40)

      integer I_FIT_ORDER        !>  Fit order (LO=0, NLO=1, NNLO=2)
      integer iDH_MOD            !>  Correction factor to the chi2, scaling
                                 !>  of statistical errors

      double precision Chi2MaxError !> Switch to Poisson stat. for larger errros
      
      integer STATYPE,SYSTYPE  !> For MC method, type of errors.
      logical lDEBUG           !> Debug flag
      logical lRAND            !> MC method activation flag 
      logical lrandData        !>  MC method using data (true) or prediction (false)
      logical WriteLHAPDF5     

      real*4 outxrange(2),     !> x-range for PDF text output 
     $     Q2VAL(NBANDS)       !> Q2 values for PDF printout
      integer  outnx           !> number of x points for PDF text output

c WS: optionally use previous fit results present in the output dir.
       !    0 - never use (default) 
       !    1 - use *.save*.txt as input to the current fit
       !    2 - do not run the fit again 
      integer UsePrevFit       !> optionally use previous fit results

C-------------------------------------------------

C Random number seed for lRand = true
      integer iSeeDmc

C parameters to set HQ scale      
      double precision hqscale2in,hqscale2inmass, hqscale1in, scalea1, scaleb1, scalec1, aq2, bq2
      integer MASSH

C  Offset method parameters
      logical CorrSystByOffset
      integer CorSysIndex ! 1-based index of corr. syst. error

C Applgrid flags
      logical useapplg
      integer napplgrids
      logical LUseAPPLgridCKM     !> Flag to define source of CKM vaules

C
C Fast/slow mode for APPLGRID PDF convolution. 
C
      integer IPDFSET          !>  QCDNUM PDF set 5-for LHAPDF 
      integer vIPDFSET         !>  value of IPDFSET for each dataset
      integer cIDataSet        !>  dataset index to be used to dermine the number of active flavour in the H-VFNS
      integer ICHECK_QCDNUM    !>  Flag to perform out-of-grid check
      logical FLAGRW
      logical DORWONLY
      logical pdfrotate   ! Special branch to perform PDF re-diagonalization
      
      logical ExtraPdfs   ! additional PDFs (e.g. photon)

      character*32 StatScale, UncorSysScale, CorSysScale,UncorChi2Type,CorChi2Type
      logical Chi2PoissonCorr, Chi2FirstIterationRescale, Chi2ExtraSystRescale

      integer AsymErrorsIterations !> Number of iterations for asymmetric errors. 0 -- off 

      logical UseDataSetIndex !> use index from the file

      integer steering_check
      common/STEERING/
     &     Q2VAL,
     $     Chi2MaxError,
     $     lDEBUG,
     $     I_FIT_ORDER,
     $     lRAND,
     $     STATYPE,
     $     SYSTYPE,
     $     outxrange,
     $     outnx,
     $     iSeedMc,
     $     useapplg,
     $     napplgrids,
     $     Lranddata,
     $     iDH_MOD,
     $     IPDFSET,
     $     vIPDFSET,
     $     cIDataSet,
     $     ICHECK_QCDNUM,
     $     UsePrevFit,
     $     CorrSystByOffset,CorSysIndex,
     $     StatScale, UncorSysScale, CorSysScale,UncorChi2Type,
     $     CorChi2Type, AsymErrorsIterations,
     $     pdfrotate, ExtraPdfs, WriteLHAPDF5, 
     $     steering_check ! Keep this at the end of the common block
      common/chi2Options/Chi2PoissonCorr, Chi2FirstIterationRescale
     $     , Chi2ExtraSystRescale

      common/STEERING2/FLAGRW,DORWONLY, UseDataSetIndex

c     parameters for QCDNUM x grid (unified with lower x values in applgrid)
      integer NMXGRID           !> Total number of X-grid sub-ranges
      parameter (NMXGRID=5)
      double precision xmin_grid(NMXGRID) !> X-grid sub ranges

      common/GRIDS/xmin_grid

      common/heavyflav/
     $     aq2,bq2,scalea1,scaleb1,scalec1,hqscale1in,hqscale2in,hqscale2inmass,MASSH



C     Value for selection of a dipole model fit

      integer DipoleModel

      common/DipoleModelType/DipoleModel

      integer DipCsModel       !> to use saturation or not use this is the answer
      common/DipCsModelType/DipCsModel

      character*128  LHAPDFSET        !> LHAPDF set name
      character*128  LHAPDFVARSET     !> LHAPDF VAR set name
      integer ILHAPDFSET              !> LHAPDF set number
      integer nLHAPDF_Sets 	      !> Number of members in LHAPDF set
      integer NPARVAR                 !> For LHAPDF VAR PDF, number of parameter variations
      logical LHAPDFErrors	      !> Flag for the LHAPDF analysis
      logical lhapdfprofile           !> Loop over eigensets and add PDF uncertainties to the chi2
      logical lhascaleprofile         !> Add QCD scale variations to the chi2
      logical Scale68                 !> Scale LHAPDF uncertainties from 90cl to 68cl
      integer nremovepriors           !> remove prior of last <n> PDF nuisance parameters
      common/CLHAPDF/LHAPDFSET,LHAPDFVARSET,ILHAPDFSET,nLHAPDF_Sets,NPARVAR,
     $	 LHAPDFErrors,Scale68,nremovepriors,lhapdfprofile,lhascaleprofile

      character*256 OutDirName   !> output dir name
      common/COutDirName/OutDirName

      character*32 RunningMode
      common/CRunningMode/RunningMode

C Define number of charged+gluon and addtional neutral PDFs:
      integer N_CHARGE_PDF
      integer N_NEUTRAL_PDF
      parameter ( N_CHARGE_PDF = 6)
      parameter ( N_NEUTRAL_PDF = 1)


C  Extra handling of LHAPDF:
      logical DataToTheo  ! if true, reset data to theory for sensitivity studies using hessian errors
      common/CDataToTheo/DataToTheo

C     30/08/2015 KK
C     WS 2015-10-02 parameters, errors and limits moved to ExtraParams
C     Namelist for Higher Twists

      logical doHiTwist
      character*32 HiTwistType
      character*64 HiTwistSubType
      common/Twistvars/doHiTwist,HiTwistType,HiTwistSubType

      logical ReadParsFromFile
      character*128 ParsFileName, CovFileName  ! Names of 
      common/CparsFromFile/ ParsFileName, CovFileName, ReadParsFromFile
