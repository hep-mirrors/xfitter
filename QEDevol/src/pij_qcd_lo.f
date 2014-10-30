C     =============================================
      double precision function dqcP111R(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP111R = dqcP0FFR(x,nf)

      return
      end

C     =============================================
      double precision function dqcP111S(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP111S = dqcP0FFS(x,nf)

      return
      end

C     =============================================
      double precision function dqcP111D(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      xx       = x       !avoid compiler warning
      qq       = qmu2    !avoid compiler warning
      jf       = nf      !avoid compiler warning
      dqcP111D = 2.D0

      return
      end

C     =============================================
      double precision function dqcP131A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP131A = (int(nf/2)-int((nf+1)/2))
     $*2.D0*dqcP0FGA(x,nf)

      return
      end

C     =============================================
      double precision function dqcP231A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP231A = 2.D0*nf*dqcP0FGA(x,nf)

      return
      end

C     =============================================
      double precision function dqcP321A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP321A = dqcP0GFA(x,nf)

      return
      end

C     =============================================
      double precision function dqcP331A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP331A = dqcP0GGA(x,nf)

      return
      end

C     =============================================
      double precision function dqcP331R(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP331R = dqcP0GGR(x,nf)

      return
      end

C     =============================================
      double precision function dqcP331S(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP331S = dqcP0GGS(x,nf)

      return
      end

C     =============================================
      double precision function dqcP331D(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      xx       = x       !avoid compiler warning
      qq       = qmu2    !avoid compiler warning
      jf       = nf      !avoid compiler warning

      dqcP331D = 6.D0*(11.D0/12.D0 - nf/18.D0)

      return
      end
