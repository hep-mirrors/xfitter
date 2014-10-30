C     =============================================
      double precision function dqcP113A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP113A = p2nspa(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP113B(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP113B = p2nsb(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP113D(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      xx       = x       !avoid compiler warning
      qq       = qmu2    !avoid compiler warning
      jf       = nf      !avoid compiler warning
      dqcP113D = p2nspc(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP123A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP123A = p2psa(x,nf)/8.D0
     $*(int(nf/2)-int((nf+1)/2))/real(nf)

      return
      end

C     =============================================
      double precision function dqcP133A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP133A = (int(nf/2)-int((nf+1)/2))/real(nf)
     $*p2qga(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP223A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP223A = p2nspa(x,nf)/8.D0+p2psa(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP223B(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP223B = p2nsb(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP223D(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      xx       = x       !avoid compiler warning
      qq       = qmu2    !avoid compiler warning
      jf       = nf      !avoid compiler warning
      dqcP223D = p2nspc(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP233A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP233A = p2qga(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP323A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP323A = p2gqa(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP333A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP333A = p2gga(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP333B(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP333B = p2ggb(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP333D(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      xx       = x       !avoid compiler warning
      qq       = qmu2    !avoid compiler warning
      jf       = nf      !avoid compiler warning

      dqcP333D = p2ggc(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP553A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP553A = p2nsma(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP553B(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP553B = p2nsb(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP553D(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP553D = p2nsmc(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP563A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP563A = (int(nf/2)-int((nf+1)/2))/real(nf)
     $*(p2nssa(x,nf)/8.D0)

      return
      end

C     =============================================
      double precision function dqcP663A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP663A = p2nssa(x,nf)/8.D0
     $+p2nsma(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP663B(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP663B = p2nsb(x,nf)/8.D0

      return
      end

C     =============================================
      double precision function dqcP663D(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP663D = p2nsmc(x,nf)/8.D0

      return
      end
