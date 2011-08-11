
      integer NBANDS
      parameter (NBANDS=40)

      integer LDOFIT1(30),Itheory,Iseed,I_FIT_ORDER,IPARAM,ICHI2
      integer isdrn
      integer LDOFIT2(30), LDOFIT3(30), LDOFIT4(30)    
      integer HFSCHEME,STATYPE,SYSTYPE, vfnsINDX
      logical lfirst,lcorr,lONLINE,lDEBUG,lNORMA,lTHEO,DOBANDS
      logical lRAND, LTOY
      real pxmax,pxmin,pq2min,pq2max
      real asatur,lsatur

      real*4 starting_scale,rtalphas,strange_frac,
     $     HF_MASS(2),Q2VAL(NBANDS),charm_frac

      real*4 outxrange(2)
      integer outform, outnx

      logical LFitDY

C-------------------------------------------------

C SG: add flag for PDF length, chebyshev param. for gluon      
      integer ILENPDF
C Number of checbyshev parameters for gluon and sea (0 - standard param.)
      integer nchebGLU,nchebSEA,IOFFSETCHEBSEA
C Type of cheb. parameterization:
C 0 - default
C 1 - (1-x) suppression.
C
      integer ichebtypeGlu, ichebtypeSea
 
      real chebxmin        ! low X-range to map with Cheb. polynoms
      real WMNLEN,WMXLEN   ! W range used for pdf length calculation

C SG, 31 Oct 2009: Add x-dependent fstrange:
      integer ifsttype  ! 0 - default, 1 - hermes-like

C SG, 23 Jan 2010: Select random number seed for lRand = true
      integer iSeeDmc

      logical IDDIPOLE,IDFLPRINT,H1QCDFUNC,idsigdq2

C SG, 25 Jan 2011, add polynomial param for valence quarks
      integer NPOLYVAL,IZPOPOLY
C SG, 1 Feb 2011, use forced positive parameterisation for IPOLYSQR=1
      integer IPOLYSQR

C 26/07/2010 added Q2max cut variable: pq2max


C AS, 07 APR 2011 applgrid flags
      logical useapplg
      integer napplgrids

      common/STEERING/
     &     Q2VAL,starting_scale,rtalphas,strange_frac,
     $     HF_MASS,pxmax,pxmin,pq2min,pq2max,charm_frac,
     &     lfirst,lcorr,lONLINE,lDEBUG,lNORMA,lTHEO,DOBANDS,
     &     LDOFIT1,Itheory,Iseed,I_FIT_ORDER,IPARAM,ICHI2,HFSCHEME,
     $     isdrn,lRAND,STATYPE,SYSTYPE,
     $     outxrange, outform, outnx,ILENPDF,nchebGLU,chebxmin,
     $     nchebSEA,wmnlen,wmxlen
     $     ,ichebtypeGlu, ichebtypeSea, ltoy
     $     ,ifsttype,iSeeDmc, LDOFIT2, LDOFIT3, LDOFIT4
     $     ,IOFFSETCHEBSEA,vfnsINDX
     $     ,IDFLprint,h1qcdfunc,IDSIGDQ2
     $     ,asatur,lsatur,npolyval
     $     ,IZPOPOLY,IPOLYSQR, useapplg, napplgrids,lfitdy

C
C Fast/slow mode for APPLGRID PDF convolution. Since the gain in speed is 
C small, set false by default
C
      logical LFastAPPLGRID
      data LFastAPPLGRID/.false./
