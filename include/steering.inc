
      integer NBANDS
      parameter (NBANDS=40)

      integer Itheory          !>  Theory type
     $     ,I_FIT_ORDER        !>  Fit order (LO=0, NLO=1, NNLO=2)
     $     ,IPARAM             !>  Internal flag for parameterisation type
     $     ,ICHI2              !>  Internal flag for chi2-type
     $     ,EWFIT              !>  EWfit (=0 NO EWFIT)
     $     ,iDH_MOD            !>  Correction factor to the chi2, scaling
                               !>  of statistical errors
      integer HFSCHEME         !>  Heavy flavour scheme (0 - zero mass, 2 - RT, 3 - FFNS)

      double precision Chi2MaxError !> Switch to Poisson stat. for larger errros
      
      integer STATYPE,SYSTYPE  !> For MC method, type of errors.
      logical lDEBUG           !> Debug flag
      logical DOBANDS            !> Produce PDF error bands 
      logical lRAND            !> MC method activation flag 
      logical lrandData        !>  MC method using data (true) or prediction (false)
      logical H1QCDFUNC        !> flag to store functions for grids in h1pdfs  		
      real asatur,             !> Saturation inspired cut for DIS data
     $     lsatur              !> Saturation inspired cut for DIS data

      real*4 
     $     starting_scale,     !> Evolution starting scale
     $     strange_frac, 
     $     HF_MASS(3),         !> masses of charm, bottom and top
     $     charm_frac

      real*4 outxrange(2),     !> x-range for PDF text output 
     $     Q2VAL(NBANDS)       !> Q2 values for PDF printout
      integer  outnx           !> number of x points for PDF text output

      logical LFitDY           !> k-factor fit to DY data

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
      double precision aq2,bq2
      integer MASSH

C Applgrid flags
      logical useapplg
      integer napplgrids

C
C Fast/slow mode for APPLGRID PDF convolution. 
C
      logical LFastAPPLGRID

      integer IPDFSET

      logical FLAGNNPDF
      logical DONNPDFONLY


      common/STEERING/
     &     Q2VAL,starting_scale,strange_frac, Chi2MaxError,
     $     HF_MASS,charm_frac,
     &     lDEBUG,DOBANDS,h1qcdfunc,
     &     Itheory,I_FIT_ORDER,IPARAM,ICHI2,HFSCHEME,
     $     lRAND,STATYPE,SYSTYPE,
     $     outxrange, outnx,ILENPDF,nchebGLU,chebxmin,
     $     nchebSEA,wmnlen,wmxlen
     $     ,ichebtypeGlu, ichebtypeSea
     $     ,ifsttype,iSeeDmc
     $     ,IOFFSETCHEBSEA,EWFIT
     $     ,asatur,lsatur,npolyval
     $     ,IZPOPOLY,IPOLYSQR, useapplg, napplgrids,lfitdy
     $     ,LFastAPPLGRID,Lranddata,iDH_MOD, IPDFSET

      common/STEERING2/FLAGNNPDF, DONNPDFONLY

     
      common/heavyflav/
     $     aq2,bq2,MASSH


C
C     Value for selection of a dipole model fit
C
      integer DipoleModel

      common/DipoleModelType/DipoleModel
