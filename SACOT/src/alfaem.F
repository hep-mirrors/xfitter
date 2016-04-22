      FUNCTION AlfaEm(Amu)
C                                                   -=-=- alfaem
C --------------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C                                                     Adapted from PYTHIA
C...Calculate real part of photon vacuum polarization.
C...For leptons simplify by using asymptotic (Q^2 >> m^2) expressions.
C...For hadrons use parametrization of H. Burkhardt et al.
C...See R. Kleiss et al, CERN 89-08, vol. 3, pp. 129-131.

      Pi = 3.1415927
      Alf0 = 1./137.04
      Alpi = Alf0/(3D0*Pi)
C                                                         AlfZ = 1./128.8
      Am2  =  Amu **2
C                                     approximate order alpha correction :
      IF(Amu .LT. 1.4D-3) Then
        Del = 0D0
      ElseIf(Amu .LT. 0.3D0) Then
        Del = Alpi*(13.4916D0 +LOG(Am2)) +0.00835D0*LOG(1D0 +Am2)
      ElseIf(Amu .LT. 3D0) Then
        Del = Alpi*(16.3200D0 + 2D0*LOG(Am2)) + 
     &  0.00238D0*LOG(1D0 + 3.927D0*Am2)
      ElseIf(Amu .LT. 1D2) Then
        Del = Alpi*(13.4955D0 + 3D0*LOG(Am2)) + 0.00165D0 + 
     &  0.00299D0*LOG(1D0 + Am2)
      Else
        Del = Alpi*(13.4955D0 + 3D0*LOG(Am2)) + 0.00221D0 + 
     &  0.00293D0*LOG(1D0 + Am2)
      EndIF
C...                                                     running alpha_em.
      AlfaEm = Alf0/(1D0-Del)

      RETURN
C                     *************************
      End
C                                                          =-=-= EwCpl1
