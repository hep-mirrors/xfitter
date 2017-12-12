C> Common block for theory expression
      integer NTermsMax
      parameter (NTermsMax = 16)

      double precision dynscale
      integer NTerms
      character*8 TermName(NTermsMax)
      character*80 TermType(NTermsMax)
      character*2048 TermInfo(NTermsMax)
      character*256 TermSource(NTermsMax)
      character*1000 TheorExpr
      integer ppbar_collisions
      integer normalised
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
     &     TermSource, TheorExpr, ppbar_collisions, normalised, 
     $     murdef, mufdef,
     $     ninfo, datainfo, CInfo, dsname, ds_index
