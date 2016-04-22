      Subroutine FILKMM
C                                                   -=-=- filkmm
C                       Given the Mixing parameters for NGN generations,
C                       calculate the KM matrix in the Cabbibo-KM-Maiani-
C                       -Wolfenstein-Chau-Keung...etc scheme.
C                       Used in SetEwk.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NGN = 3, NANG= NGN*(NGN-1)/2, NPHS=(NGN-1)*(NGN-2)/2)

      COMMON / KMATRX / VKM (NGN, NGN)
      COMMON / MIXPAR / CMX(NANG), DMX(NPHS)

C     (how to handle complex numbers due to phase factors efficiently,
C      bearing in mind that most applications do not need this info??)

C     For the moment, put the absolute values of the matrix elements in by hand

      DIMENSION AK (3, 3)
c CKM values from Particle data group
c W.-M. Yao et al., J. of Phys. G 33, 1 (2006)
      DATA (AK(I,1), I=1,3) / 0.97383, 0.2272, 0.00396 /
      DATA (AK(I,2), I=1,3) / 0.2271, 0.97296, 0.04221 /
      DATA (AK(I,3), I=1,3) / 0.00814, 0.04161, 0.9991 /
C           for debugging convenience ....
C      DATA (AK(I,1), I=1,3) / 1d0, 0d0, 0d0 /
C      DATA (AK(I,2), I=1,3) / 0d0, 1d0, 0d0 /
C      DATA (AK(I,3), I=1,3) / 0d0, 0d0, 1d0 /

      DO 5 I = 1, 3
      DO 6 J = 1, 3
        VKM (I,J) = AK (I,J)
    6 CONTINUE
    5 CONTINUE

      Return
C                        ****************************
      END

C                                                          =-=-= EwCpl0
