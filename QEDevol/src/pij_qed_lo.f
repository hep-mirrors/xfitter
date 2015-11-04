C     =============================================
      double precision function dqcP114R(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP114R = 5.D0/24.D0*dqcP0FFR(x,nf)

      return
      end

C     =============================================
      double precision function dqcP114S(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP114S = dqcP0FFS(x,nf)

      return
      end

C     =============================================
      double precision function dqcP114D(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      xx       = x       !avoid compiler warning
      qq       = qmu2    !avoid compiler warning
      jf       = nf      !avoid compiler warning
      dqcP114D = 5.D0/12.D0

      return
      end

C     =============================================
      double precision function dqcP124R(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP124R = 1.D0/8.D0*dqcP0FFR(x,nf)

      return
      end

C     =============================================
      double precision function dqcP124S(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP124S = dqcP0FFS(x,nf)

      return
      end

C     =============================================
      double precision function dqcP124D(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      xx       = x       !avoid compiler warning
      qq       = qmu2    !avoid compiler warning
      jf       = nf      !avoid compiler warning
      dqcP124D = 1.D0/4.D0

      return
      end

C     =============================================
      double precision function dqcP144A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP144A = 4.D0/9.D0
     $*(4.D0*int(nf/2)-int((nf+1)/2))
     $*dqcP0FGA(x,nf)

      return
      end

C     =============================================
      double precision function dqcP244A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP244A = 4.D0/9.D0
     $*(4.D0*int(nf/2)+int((nf+1)/2))
     $*dqcP0FGA(x,nf)

      return
      end

C     =============================================
      double precision function dqcP414A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP414A = 1.D0/8.D0*dqcP0GFA(x,nf)

      return
      end

C     =============================================
      double precision function dqcP424A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP424A = 5.D0/24.D0*dqcP0GFA(x,nf)

      return
      end

C     =============================================
      double precision function dqcP444D(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      xx       = x       !avoid compiler warning
      qq       = qmu2    !avoid compiler warning
      jf       = nf      !avoid compiler warning
      dqcP444D = -2.D0/9.D0
     $*(4.D0*int(nf/2)+int((nf+1)/2))

      return
      end

C     =============================================
      double precision function dqcP554R(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP554R = 5.D0/24.D0*dqcP0FFR(x,nf)

      return
      end

C     =============================================
      double precision function dqcP554S(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP554S = dqcP0FFS(x,nf)

      return
      end

C     =============================================
      double precision function dqcP554D(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      xx       = x       !avoid compiler warning
      qq       = qmu2    !avoid compiler warning
      jf       = nf      !avoid compiler warning
      dqcP554D = 5.D0/12.D0

      return
      end

C     =============================================
      double precision function dqcP564R(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP564R = 1.D0/8.D0*dqcP0FFR(x,nf)

      return
      end

C     =============================================
      double precision function dqcP564S(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP564S = dqcP0FFS(x,nf)

      return
      end

C     =============================================
      double precision function dqcP564D(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      xx       = x       !avoid compiler warning
      qq       = qmu2    !avoid compiler warning
      jf       = nf      !avoid compiler warning
      dqcP564D = 1.D0/4.D0

      return
      end

C     =============================================
      double precision function dqcP774R(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP774R = 1.D0/12.D0*dqcP0FFR(x,nf)

      return
      end

C     =============================================
      double precision function dqcP774S(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP774S = dqcP0FFS(x,nf)

      return
      end

C     =============================================
      double precision function dqcP774D(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      xx       = x       !avoid compiler warning
      qq       = qmu2    !avoid compiler warning
      jf       = nf      !avoid compiler warning
      dqcP774D = 1.D0/6.D0

      return
      end

C     =============================================
      double precision function dqcP884R(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP884R = 1.D0/3.D0*dqcP0FFR(x,nf)

      return
      end

C     =============================================
      double precision function dqcP884S(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcP884S = dqcP0FFS(x,nf)

      return
      end

C     =============================================
      double precision function dqcP884D(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      xx       = x       !avoid compiler warning
      qq       = qmu2    !avoid compiler warning
      jf       = nf      !avoid compiler warning
      dqcP884D = 2.D0/3.D0

      return
      end
