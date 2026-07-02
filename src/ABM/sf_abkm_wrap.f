C
C OZ 10.10.17 updated to openqcdrad-2.1
C New user routines for PDFs and alpha_S have to be provided in this version

      FUNCTION useralphas(q2,kschemepdf,kordpdf,kpdfset)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'APSCOM6.'
      q = dsqrt(q2)
      useralphas = alphas_wrapper(q)
      RETURN
      END


      FUNCTION userpdfs(xb,q2,i,kschemepdf,kordpdf,kpdfset)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer i
      dimension pdfsff(-6:6)
      q = dsqrt(q2)
      CALL pdf_xfxq_wrapper(xb,q,PDFSFF)
      userpdfs = PDFSFF(i)
      RETURN
      END


      DOUBLE PRECISION FUNCTION DiLog(X)
      double precision x
      double precision ddilog
      dilog = ddilog(x)
      return
      end
