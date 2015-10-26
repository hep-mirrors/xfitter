C> Common block for theory expression
      integer NTermsMax
      parameter (NTermsMax = 16)

      double precision dynscale
      integer NTerms
      character*8 TermName(NTermsMax)
      character*80 TermType(NTermsMax)
      character*80 TermInfo(NTermsMax)
      character*1000 TermSource(NTermsMax)
      character*1000 TheorExpr
      integer ppbar_collisions
      integer normalised
      integer murdef
      integer mufdef

      common/theorexpr/ dynscale, NTerms, TermName, TermType, TermInfo, 
     & TermSource, TheorExpr, ppbar_collisions, normalised, murdef, mufdef
