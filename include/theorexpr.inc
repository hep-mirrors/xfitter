C> Common block for theory expression
      integer NTermsMax
      parameter (NTermsMax = 16)

      integer NTerms
      character*8 TermName(NTermsMax)
      character*80  TermType(NTermsMax)
      character*1000 TermSource(NTermsMax)
      character*1000 TheorExpr

      common/theorexpr/ NTerms, TermName, TermType, TermSource, TheorExpr
