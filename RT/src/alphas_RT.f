
      FUNCTION ALPHA(T)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)



      COMMON/DYLAMB/XLAM,S0

      QS = XLAM**2*EXP(T)
cv      print*,'inside ', qs, xlam, t 
      alpha =ASFUNC(qs, nf, ierr) 
cv      print*,'xlam', qs, xlam, t, nf, ierr, asfunc(qs,nf,ierr) 
cv      stop
      RETURN
      END
