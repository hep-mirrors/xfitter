C> Common block for theory expression
      integer NTermsMax
      parameter (NTermsMax = 128)

      double precision dynscale
      integer NTerms
      character*32 TermName(NTermsMax)
      character*80 TermType(NTermsMax)
      character*4096 TermInfo(NTermsMax)
      character*256 TermSource(NTermsMax)
      character*10000 TheorExpr
      integer normalised
      double precision normalisation
      !TODO: remove the following 3 variables
      integer ppbar_collisions
      integer murdef
      integer mufdef

C Also store dataset info here
      integer ninfoMax
      parameter (ninfoMax=100)
      integer ninfo
      double precision datainfo(ninfoMax)
      character *80 CInfo(ninfoMax)
      character *80 dsname
      integer ds_index

C And some basic vars:

      common/theorexpr/ dynscale, NTerms, TermName, TermType, TermInfo,
     &     TermSource, TheorExpr, ppbar_collisions, normalised, normalisation,
     $     murdef, mufdef,
     $     ninfo, datainfo, CInfo, dsname, ds_index
