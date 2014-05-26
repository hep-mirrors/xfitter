
      integer NBANDS
      parameter (NBANDS=40)

      integer Itheory          !>  Theory type
     $     ,I_FIT_ORDER        !>  Fit order (LO=0, NLO=1, NNLO=2)
     $     ,IPARAM             !>  Internal flag for parameterisation type
     $     ,EWFIT              !>  EWfit (=0 NO EWFIT)
     $     ,iDH_MOD            !>  Correction factor to the chi2, scaling
                               !>  of statistical errors
      character*32 HF_SCHEME      
      integer HFSCHEME          !>  Heavy flavour scheme (0 - zero mass, 2 - RT, 3 - FFNS)

      double precision Chi2MaxError !> Switch to Poisson stat. for larger errros
      
      integer STATYPE,SYSTYPE  !> For MC method, type of errors.
      logical lDEBUG           !> Debug flag
      logical DOBANDS            !> Produce PDF error bands 
      logical lRAND            !> MC method activation flag 
      logical lrandData        !>  MC method using data (true) or prediction (false)
      logical H1QCDFUNC        !> flag to store functions for grids in h1pdfs  		
      logical UseGridLHAPDF5     !> store PDFs in LHAPDF5 format
      logical WriteLHAPDF6     !> store PDFs in LHAPDF6 format
      real*4 
     $     starting_scale,     !> Evolution starting scale
     $     strange_frac, 
     $     HF_MASS(3),         !> masses of charm, bottom and top
     $     charm_frac

      real*4 outxrange(2),     !> x-range for PDF text output 
     $     Q2VAL(NBANDS)       !> Q2 values for PDF printout
      integer  outnx           !> number of x points for PDF text output

      logical LFitDY           !> k-factor fit to DY data

c WS: optionally use previous fit results present in the output dir.
       !    0 - never use (default) 
       !    1 - use *.save*.txt as input to the current fit
       !    2 - do not run the fit again 
      integer UsePrevFit       !> optionally use previous fit results

C-------------------------------------------------

C Flag for PDF length, chebyshev param. for gluon      
      integer ILENPDF
      real WMNLEN,WMXLEN   ! W range used for pdf length calculation

C Number of checbyshev parameters for gluon and sea (0 - standard param.)
      integer nchebGLU,nchebSEA,IOFFSETCHEBSEA
C Type of cheb. parameterization:
C 0 - default
C 1 - (1-x) suppression.
C
      integer ichebtypeGlu, ichebtypeSea
 
      real chebxmin        !> low X-range to map with Cheb. polynoms

C X-dependent fstrange:
      integer ifsttype  !> 0 - default, 1 - hermes-like

C Random number seed for lRand = true
      integer iSeeDmc

C Polynomial param for valence quarks
      integer NPOLYVAL,IZPOPOLY
C Use forced positive parameterisation for IPOLYSQR=1
      integer IPOLYSQR

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
      logical LFastAPPLGRID

      integer IPDFSET          !>  QCDNUM PDF set 5-for LHAPDF 
      integer ICHECK_QCDNUM    !>  Flag to perform out-of-grid check
      logical lead             !> Flag to trigger lead PDF	
      logical deuteron             !> Flag to trigger deuteron PDF	
      logical CachePDFs        !> Cache PDF calls to speedup
      logical FLAGRW
      logical DORWONLY
      logical FlexibleGluon    !> Flexible gluon parameterisation
      logical Read_QCDNUM_Tables
      logical WriteAlphaSToMemberPDF    !> write alphas header to LHAPDF6 member files

      character*32 StatScale, UncorSysScale, CorSysScale,UncorChi2Type,CorChi2Type
      logical Chi2PoissonCorr, Chi2FirstIterationRescale, Chi2ExtraSystRescale

      integer AsymErrorsIterations !> Number of iterations for asymmetric errors. 0 -- off 
      common/STEERING/
     &     Q2VAL,starting_scale,strange_frac, Chi2MaxError,
     $     HF_MASS,charm_frac,
     &     lDEBUG,DOBANDS,UseGridLHAPDF5,WriteLHAPDF6,h1qcdfunc,WriteAlphaSToMemberPDF,
     &     Itheory,I_FIT_ORDER,IPARAM,HFSCHEME,
     $     lRAND,STATYPE,SYSTYPE,
     $     outxrange, outnx,ILENPDF,nchebGLU,chebxmin,
     $     nchebSEA,wmnlen,wmxlen
     $     ,ichebtypeGlu, ichebtypeSea
     $     ,ifsttype,iSeeDmc
     $     ,IOFFSETCHEBSEA,EWFIT
     $     ,npolyval, lead
     $     ,IZPOPOLY,IPOLYSQR, useapplg, napplgrids,lfitdy
     $     ,LFastAPPLGRID,Lranddata,iDH_MOD
     $     ,IPDFSET, ICHECK_QCDNUM
     $     ,UsePrevFit
     $     ,CorrSystByOffset,CorSysIndex
     $     ,StatScale, UncorSysScale, CorSysScale,UncorChi2Type,
     $     CorChi2Type, hf_scheme, AsymErrorsIterations,
     $     LUseAPPLgridCKM
      
      common/chi2Options/Chi2PoissonCorr, Chi2FirstIterationRescale
     $     , Chi2ExtraSystRescale

      common/STEERING2/FLAGRW, DORWONLY, CachePDFs, FlexibleGluon,
     $     read_qcdnum_tables

c     parameters for QCDNUM x grid (unified with lower x values in applgrid)
      integer NMXGRID           !> Total number of X-grid sub-ranges
      parameter (NMXGRID=5)
      double precision xmin_grid(NMXGRID) !> X-grid sub ranges

      common/GRIDS/xmin_grid

      character*64 PDF_DECOMPOSITION      !> Style of PDF decompisition, e.g U_D_Ubar_Dbar
      common/STEERING3/PDF_DECOMPOSITION     

      common/heavyflav/
     $     aq2,bq2,scalea1,scaleb1,scalec1,hqscale1in,hqscale2in,hqscale2inmass,MASSH


C
C     Value for selection of a dipole model fit
C
      integer DipoleModel

      common/DipoleModelType/DipoleModel

      character*32 PDFStyle           !> PDF style definition
      common/CPdfStyle/PDFStyle

      character*32 PDFType           !> PDF type definition
      common/CPdfStyle/PDFType

      character*128  LHAPDFSET        !> LHAPDF set name
      character*128  LHAPDFVARSET     !> LHAPDF VAR set name
      integer ILHAPDFSET              !> LHAPDF set number
      integer nLHAPDF_Sets 	      !> Number of members in LHAPDF set
      integer NPARVAR                 !> For LHAPDF VAR PDF, number of parameter variations
      logical LHAPDFErrors	!> For LHAPDF PDFtype, loop over eighensets
      logical Scale68           !> Scale LHAPDF uncertainties from 90cl to 68cl
      common/CLHAPDF/LHAPDFSET,LHAPDFVARSET,ILHAPDFSET,nLHAPDF_Sets,NPARVAR,LHAPDFErrors,Scale68

      character*128 OutDirName   !> outout dir name
      character*128 LHAPDF6OutDir
      common/COutDirName/OutDirName, LHAPDF6OutDir
