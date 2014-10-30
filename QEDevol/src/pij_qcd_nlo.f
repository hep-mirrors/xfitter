C     =============================================
      double precision function dqcP112A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP112A = pp1sfunc(x,nf)-pm1sfunc(x,nf)

      return
      end

C     =============================================
      double precision function dqcP112B(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP112B = pm1sfunc(x,nf)

      return
      end

C     =============================================
      double precision function dqcP122A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP122A = (ff1sfunc(x,nf)-xf1tfunc(x,nf)
     $-pp1sfunc(x,nf)+pm1sfunc(x,nf))
     $*(int(nf/2)-int((nf+1)/2))/real(nf)

      return
      end

C     =============================================
      double precision function dqcP122B(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP122B = (xf1tfunc(x,nf)-pm1sfunc(x,nf))
     $*(int(nf/2)-int((nf+1)/2))/real(nf)

      return
      end

C     =============================================
      double precision function dqcP132A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP132A = (int(nf/2)-int((nf+1)/2))/real(nf)
     $*gf1sfunc(x,nf)

      return
      end

C     =============================================
      double precision function dqcP222A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP222A = ff1sfunc(x,nf)-xf1tfunc(x,nf)

      return
      end

C     =============================================
      double precision function dqcP222B(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP222B = xf1tfunc(x,nf)

      return
      end

C     =============================================
      double precision function dqcP232A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP232A = gf1sfunc(x,nf)

      return
      end

C     =============================================
      double precision function dqcP322A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP322A = fg1sfunc(x,nf)

      return
      end

C     =============================================
      double precision function dqcP332A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP332A = gg1sfunc(x,nf)-xg1tfunc(x,nf)

      return
      end

C     =============================================
      double precision function dqcP332B(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP332B = xg1tfunc(x,nf)

      return
      end

C     =============================================
      double precision function dqcP552B(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP552B = pm1sfunc(x,nf)

      return
      end
