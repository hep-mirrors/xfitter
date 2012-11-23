
*************************************************************************  

c ... The z-dependent parts, i.e, L1**k
c ... are local pieces from the plus-distribution
c ... The z-independent part is the real delta-contribution
      double precision function C2nsdelta(z,Q2,nf,iord)
      implicit none
      double precision z,Q2
      integer nf,iord

      double precision C2nsdel(1:3)
      real*8 C2NN2C,C2NP3C

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      double precision myalphas ! Interface to alpha_s/(4 pi)
      integer i1
     
      double precision alps
      double precision CF

      double precision L1,pi,zeta2

c      C2nsdelta = 0d0
      C2nsdelta = 1d0 ! leading order contribution

      CF = 4d0/3d0

c ... initialize Wilson coefficients at O(alpha_s^i1)
      do i1 = 1,3
         C2nsdel(i1) = 0d0
      enddo

c ... alpha_s/(4 pi)(Q2) ; Note for the time being we set mur^2 = muf^2 = Q2
      alps = myalphas(Q2,iord+1,Isetalps)

      if (iord .ge. 1) then  ! O(alpha_s^1)
         L1 = dlog(1-z)
         pi = dacos(-1d0)
         zeta2 = pi**2/6d0
         C2nsdel(1) = CF * (-(9d0+4d0*zeta2)+2d0*L1**2-3d0*L1) ! Eq. (4.3) in [1]
      endif

      if (iord .ge. 2) then ! O(alpha_s^2)
         C2nsdel(2) = dble(C2NN2C(z, nf))  ! Eq. (4.8) in [1]
      endif

      if (iord .ge. 3) then ! O(alpha_s^3)
         C2nsdel(3) = dble(C2NP3C(z, nf))  ! Eq. (4.11) in [1]
      endif

      do i1 = 1,iord
         C2nsdelta = C2nsdelta + alps**i1 * C2nsdel(i1)
      enddo

      end

*************************************************************************  

      double precision function C2nsplus(z,Q2,nf,iord)
      implicit none
      double precision z,Q2
      integer nf,iord

      double precision C2nspl(1:3)
      real*8 C2NS2B,C2NS3B

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      double precision myalphas ! Interface to alpha_s/(4 pi)
      integer i1
     
      double precision alps
      double precision CF

      double precision x,x1,L1

      C2nsplus = 0d0

      CF = 4d0/3d0

c ... initialize Wilson coefficients at O(alpha_s^i1)
      do i1 = 1,3
         C2nspl(i1) = 0d0
      enddo

c ... alpha_s/(4 pi)(Q2) ; Note for the time being we set mur^2 = muf^2 = Q2
      alps = myalphas(Q2,iord+1,Isetalps)

      if (iord .ge. 1) then  ! O(alpha_s^1)
         x = z
         x1 = 1d0 - x
         L1 = dlog(x1)
         C2nspl(1) = CF * (4d0*L1 - 3d0)/x1 ! Eq. (4.3) in [1]
      endif

      if (iord .ge. 2) then ! O(alpha_s^2)
         C2nspl(2) = dble(C2NS2B(z, nf)) ! Eq. (4.8) in [1]
      endif

      if (iord .ge. 3) then ! O(alpha_s^3)
         C2nspl(3) = dble(C2NS3B(z, nf)) ! Eq. (4.11) in [1]
      endif

      do i1 = 1,iord
         C2nsplus = C2nsplus + alps**i1 * C2nspl(i1)
      enddo

      end

*************************************************************************  

      double precision function C2nsregular(z,Q2,nf,iord)
      implicit none
      double precision z,Q2
      integer nf,iord

      double precision C2nsr(1:3)
      real*8 C2NN2A,C2NP3A

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      double precision myalphas ! Interface to alpha_s/(4 pi)
      integer i1
     
      double precision alps
      double precision CF

      double precision x,x1,L0,L1

      C2nsregular = 0d0

      CF = 4d0/3d0

c ... initialize Wilson coefficients at O(alpha_s^i1)
      do i1 = 1,3
         C2nsr(i1) = 0d0
      enddo

c ... alpha_s/(4 pi)(Q2) ; Note for the time being we set mur^2 = muf^2 = Q2
c ... Note the iord+1! For example NLO F2 has iord = 1 but uses
c ... normally the two-loop alpha_s
      alps = myalphas(Q2,iord+1,Isetalps)

      if (iord .ge. 1) then  ! O(alpha_s^1)
c ... Eq. (4.3) in [1]; regular part
         x = z
         x1 = 1d0 - x
         L0 = dlog(x)
         L1 = dlog(x1)
         C2nsr(1) = CF * (-2d0*(1d0+x)*(L1-L0)-4d0/x1*L0+6d0+4d0*x)
      endif

      if (iord .ge. 2) then ! O(alpha_s^2)
         C2nsr(2) = dble(C2NN2A (z, nf)) ! from xc2ns2p.f
      endif

      if (iord .ge. 3) then ! O(alpha_s^3)
         C2nsr(3) = dble(C2NP3A (z, nf)) ! from xc2ns3p.f
      endif

      do i1 = 1,iord
         C2nsregular = C2nsregular + alps**i1 * C2nsr(i1)
      enddo

      end

*************************************************************************  

      double precision function C2psdelta(z,Q2,nf,iord)
      implicit none
      double precision z,Q2
      integer nf,iord

      double precision C2psdel(1:3)
      real*8 C2S3C

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      double precision myalphas ! Interface to alpha_s/(4 pi)
      integer i1
     
      double precision alps
      double precision CF

      C2psdelta = 0d0

      CF = 4d0/3d0

c ... initialize Wilson coefficients at O(alpha_s^i1)
      do i1 = 1,3
         C2psdel(i1) = 0d0
      enddo

c ... alpha_s/(4 pi)(Q2) ; Note for the time being we set mur^2 = muf^2 = Q2
      alps = myalphas(Q2,iord+1,Isetalps)

      if (iord .ge. 1) then  ! O(alpha_s^1)
         C2psdel(1) = 0d0 ! Eq. (4.2) in [1]
      endif

      if (iord .ge. 2) then ! O(alpha_s^2)
         C2psdel(2) = 0d0  ! Eq. (4.9) in [1]
      endif

      if (iord .ge. 3) then ! O(alpha_s^3)
         C2psdel(3) = dble(C2S3C(z,nf))  ! Eq. (4.12) in [1]; xc2sg3p.f
      endif

      do i1 = 1,iord
         C2psdelta = C2psdelta + alps**i1 * C2psdel(i1)
      enddo

      end

*************************************************************************  

      double precision function C2psregular(z,Q2,nf,iord)
      implicit none
      double precision z,Q2
      integer nf,iord

      double precision C2psr(1:3)
      real*8 C2S2A,C2S3A

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      double precision myalphas ! Interface to alpha_s/(4 pi)
      integer i1
     
      double precision alps

      C2psregular = 0d0

c ... initialize Wilson coefficients at O(alpha_s^i1)
      do i1 = 1,3
         C2psr(i1) = 0d0
      enddo

c ... alpha_s/(4 pi)(Q2) ; Note for the time being we set mur^2 = muf^2 = Q2
      alps = myalphas(Q2,iord+1,Isetalps)

      if (iord .ge. 1) then  ! O(alpha_s^1)
         C2psr(1) = 0d0 ! Eq. (4.2) in [1]
      endif

      if (iord .ge. 2) then ! O(alpha_s^2)
         C2psr(2) = dble(C2S2A (z, nf)) ! from xc2sg2p.f
      endif

      if (iord .ge. 3) then ! O(alpha_s^3)
         C2psr(3) = dble(C2S3A (z, nf)) ! from xc2sg3p.f
      endif

      do i1 = 1,iord
         C2psregular = C2psregular + alps**i1 * C2psr(i1)
      enddo

      end


























c .......................................................................

      double precision function CLgregular(z,Q2,nf,iord)
      implicit none
      double precision z,Q2
      integer nf,iord

      double precision CLgr(1:3)
      real*8 CLG2A,CLG3A

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      double precision myalphas ! Interface to alpha_s/(4 pi)
      integer i1
     
      double precision alps

      CLgregular = 0d0

c ... initialize Wilson coefficients at O(alpha_s^i1)
      do i1 = 1,3
         CLgr(i1) = 0d0
      enddo

c ... alpha_s/(4 pi)(Q2) ; Note for the time being we set mur^2 = muf^2 = Q2
c ... Note the iord+1! For example NLO F2 has iord = 1 but uses
c ... normally the two-loop alpha_s
      alps = myalphas(Q2,iord+1,Isetalps)

      if (iord .ge. 1) then  ! O(alpha_s^1)
         CLgr(1) = 8d0*nf*z*(1d0-z) ! Eq. (3) in [2]
      endif

      if (iord .ge. 2) then ! O(alpha_s^2)
         CLgr(2) = dble(CLG2A (z, nf)) ! from xclsg2p.f
      endif

      if (iord .ge. 3) then ! O(alpha_s^3)
         CLgr(3) = dble(CLG3A (z, nf)) ! from xclsg3p.f
      endif

      do i1 = 1,iord
         CLgregular = CLgregular + alps**i1 * CLgr(i1)
      enddo

      end

*
* ---------------------------------------------------------------------
*
*
* ..This is the regular non-singlet piece for the odd-moment (CC) FL, 
*    corresponding to CLNSP-CLNSN in WvN's program. The 8 numerical 
*    coefficients the NF^0 piece are fitted to his results. The NF part
*    is the same as in CLNN2A.
*




















c .......................................................................

      double precision function CLnsdelta(z,Q2,nf,iord)
      implicit none
      double precision z,Q2
      integer nf,iord

      double precision CLnsdel(1:3)
      real*8 CLNP2C,CLNP3C

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      double precision myalphas ! Interface to alpha_s/(4 pi)
      integer i1
     
      double precision alps
      double precision CF

      CLnsdelta = 0d0

      CF = 4d0/3d0

c ... initialize Wilson coefficients at O(alpha_s^i1)
      do i1 = 1,3
         CLnsdel(i1) = 0d0
      enddo

c ... alpha_s/(4 pi)(Q2) ; Note for the time being we set mur^2 = muf^2 = Q2
      alps = myalphas(Q2,iord+1,Isetalps)

      if (iord .ge. 1) then  ! O(alpha_s^1)
         CLnsdel(1) = 0d0 ! Eq. (3) in [2]
      endif

      if (iord .ge. 2) then ! O(alpha_s^2)
         CLnsdel(2) = dble(CLNP2C(z))  ! does not depend on z; Eq. (4) in [2]
      endif

      if (iord .ge. 3) then ! O(alpha_s^3)
         CLnsdel(3) = dble(CLNP3C(z,nf))  ! does not depend on z; Eq. (8) in [2]
      endif

      do i1 = 1,iord
         CLnsdelta = CLnsdelta + alps**i1 * CLnsdel(i1)
      enddo

      end

c .......................................................................

      double precision function CLnsplus(z,Q2,nf,iord)
      implicit none
      double precision z,Q2
      integer nf,iord

      double precision CLnspl(1:3)

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      double precision myalphas ! Interface to alpha_s/(4 pi)
      integer i1
     
      double precision alps
      double precision CF

      CLnsplus = 0d0

      CF = 4d0/3d0

c ... initialize Wilson coefficients at O(alpha_s^i1)
      do i1 = 1,3
         CLnspl(i1) = 0d0
      enddo

c ... alpha_s/(4 pi)(Q2) ; Note for the time being we set mur^2 = muf^2 = Q2
      alps = myalphas(Q2,iord+1,Isetalps)

      if (iord .ge. 1) then  ! O(alpha_s^1)
         CLnspl(1) = 0d0 ! Eq. (3) in [2]
      endif

      if (iord .ge. 2) then ! O(alpha_s^2)
         CLnspl(2) = 0d0 ! Eq. (4) in [2]
      endif

      if (iord .ge. 3) then ! O(alpha_s^3)
         CLnspl(3) = 0d0 ! Eq. (8) in [2]
      endif

      do i1 = 1,iord
         CLnsplus = CLnsplus + alps**i1 * CLnspl(i1)
      enddo

      end

c .......................................................................

      double precision function CLnsregular(z,Q2,nf,iord)
      implicit none
      double precision z,Q2
      integer nf,iord

      double precision CLnsr(1:3)
      real*8 CLNP2A,CLNP3A

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      double precision myalphas ! Interface to alpha_s/(4 pi)
      integer i1
     
      double precision alps
      double precision CF

      CLnsregular = 0d0

      CF = 4d0/3d0

c ... initialize Wilson coefficients at O(alpha_s^i1)
      do i1 = 1,3
         CLnsr(i1) = 0d0
      enddo

c ... alpha_s/(4 pi)(Q2) ; Note for the time being we set mur^2 = muf^2 = Q2
      alps = myalphas(Q2,iord+1,Isetalps)

      if (iord .ge. 1) then  ! O(alpha_s^1)
         CLnsr(1) = 4d0*CF*z ! Eq. (3) in [2]
      endif

      if (iord .ge. 2) then ! O(alpha_s^2)
         CLnsr(2) = dble(CLNP2A (z, nf)) ! from xclns2p.f
      endif

      if (iord .ge. 3) then ! O(alpha_s^3)
         CLnsr(3) = dble(CLNP3A (z, nf)) ! from xclns3p.f
      endif

      do i1 = 1,iord
         CLnsregular = CLnsregular + alps**i1 * CLnsr(i1)
      enddo

      end

c .......................................................................

      double precision function CLpsregular(z,Q2,nf,iord)
      implicit none
      double precision z,Q2
      integer nf,iord

      double precision CLpsr(1:3)
      real*8 CLS2A,CLS3A

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      double precision myalphas ! Interface to alpha_s/(4 pi)
      integer i1
     
      double precision alps

      CLpsregular = 0d0

c ... initialize Wilson coefficients at O(alpha_s^i1)
      do i1 = 1,3
         CLpsr(i1) = 0d0
      enddo

c ... alpha_s/(4 pi)(Q2) ; Note for the time being we set mur^2 = muf^2 = Q2
      alps = myalphas(Q2,iord+1,Isetalps)

      if (iord .ge. 1) then  ! O(alpha_s^1)
         CLpsr(1) = 0d0 ! Eq. (3) in [2]
      endif

      if (iord .ge. 2) then ! O(alpha_s^2)
         CLpsr(2) = dble(CLS2A (z, nf)) ! from xclsg2p.f
      endif

      if (iord .ge. 3) then ! O(alpha_s^3)
         CLpsr(3) = dble(CLS3A (z, nf)) ! from xclsg3p.f
      endif

      do i1 = 1,iord
         CLpsregular = CLpsregular + alps**i1 * CLpsr(i1)
      enddo

      end











c ... wilson.f TS/IS
c ...03.07.2011

c .......................................................................

      double precision function C2gregular(z,Q2,nf,iord)
      implicit none
      double precision z,Q2
      integer nf,iord

      double precision C2gr(1:3)
      real*8 C2G2A,C2G3A

      integer IsetPDF,Isetalps
      double precision A1,Z1
      common / sfs_switches / IsetPDF,Isetalps,A1,Z1

      double precision myalphas ! Interface to alpha_s/(4 pi)
      integer i1
     
      double precision alps
      double precision x,x1,L0,L1

      C2gregular = 0d0

c ... initialize Wilson coefficients at O(alpha_s^i1)
      do i1 = 1,3
         C2gr(i1) = 0d0
      enddo

c ... alpha_s/(4 pi)(Q2) ; Note for the time being we set mur^2 = muf^2 = Q2
      alps = myalphas(Q2,iord+1,Isetalps)
      
      if (iord .ge. 1) then  ! O(alpha_s^1)
         x = z
         x1 = 1d0 - x
         L0 = dlog(x)
         L1 = dlog(x1)
         C2gr(1) = nf*((2d0-4d0*x*x1)*(L1-L0)-2d0+16d0*x*x1) ! Eq. (4.4) in [1]
      endif

      if (iord .ge. 2) then ! O(alpha_s^2)
         C2gr(2) = dble(C2G2A (z, nf)) ! from xc2sg2p.f
      endif

      if (iord .ge. 3) then ! O(alpha_s^3)
         C2gr(3) = dble(C2G3A (z, nf)) ! from xc2sg3p.f
      endif

      do i1 = 1,iord
         C2gregular = C2gregular + alps**i1 * C2gr(i1)
      enddo

      end

c .......................................................................

c ... Eq. (4) in [2]
c ... This approx. is supposed to be an improvement over the approx. CLNN2A
c ... NOTE: CLNP2A has to be used together with CLNP2C !!!
       FUNCTION CLNP2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       CLNP2A = 
     #          - 37.338 + 89.53 * Y + 128./9. * Y * DL1**2 
     #          - 46.50 * Y * DL1 - 84.094 * DL * DL1
     #          + 33.82 * Y**2 + Y * DL * (32.90 + 18.41 * DL)
     #          - 128./9. * DL 
     #          + NF * 16./27. * ( 6.* Y*DL1 - 12.* Y*DL - 25.* Y + 6.)

       RETURN
       END

c .......................................................................

c ... See Eq. (4) in [2]
c ... This approx. is supposed to be an improvement over the approx. CLNN2C
c ... NOTE: CLNP2C has to be used together with CLNP2A !!!
       FUNCTION CLNP2C (Y)
       IMPLICIT REAL*8 (A-Z)
*
       CLNP2C = -0.012
*
       RETURN
       END
*
* =================================================================av==
*
* ..File: xc2sg3p.f    F_2^PS  and  F_2^G
*
*
* ..Parametrizations of the 3-loop MS(bar) pure-singlet and gluon coef-
*    ficient functions for the electromagnetic structure function F_2
*    at  mu_r = mu_f = Q. The expansion parameter is  alpha_s/(4 pi).
*
*  ..The relative accuracy of these parametrizations, as well as of
*    the convolution results, is one part in thousand or better.
*
* ..Reference: J. Vermaseren, A. Vogt and S. Moch 
*              hep-ph/0504242 = Nucl. Phys. B724 (2005) 3
* 
* =====================================================================
*
*
* ..The pure-singlet coefficient function, regular piece
*
       FUNCTION C2S3A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       DIMENSION FL(0:6), FLS(0:6)  !*** Extended to deal with NF=0  (FIO 2 Aug 2011) 
       INTEGER NF
       DATA FL  /0.0d0,-1.d0, 0.5d0, 0.d0, 0.5d0, 0.2d0, 0.5d0 /
       DATA FLS /0.0d0, 1.d0, 0.1d0, 0.d0, 0.1d0, 0.018181818d0, 0.1d0 /
*
       DL  = LOG (Y)
       Y1  = 1.-Y
       DL1 = LOG (Y1)
       D9  = 1./9.D0
       D81 = D9*D9
*
       C2S31 = ( 856.*D81 * DL1**4 - 6032.*D81 * DL1**3 + 130.57* DL1**2
     ,         - 542.0 * DL1 + 8501. - 4714.* Y + 61.50 * Y**2 ) * Y1
     ,         + DL*DL1 * (8831.* DL + 4162.* Y1) - 15.44 * Y*DL**5     
     ,         + 3333.* Y*DL**2 + 1615.* DL + 1208.* DL**2 
     ,         - 333.73 * DL**3 + 4244.*D81 * DL**4 - 40.*D9 * DL**5 
     ,         - 2731.82 * Y1/Y - 414.262 * DL/Y
       C2S32 = ( - 64.*D81 * DL1**3 + 208.*D81 * DL1**2 + 23.09 * DL1
     ,         - 220.27 + 59.80 * Y - 177.6 * Y**2) * Y1  
     ,         -  DL*DL1 * (160.3 * DL + 135.4 * Y1) - 24.14 * Y*DL**3 
     ,         - 215.4 * Y*DL**2 - 209.8 * DL - 90.38 * DL**2 
     ,         - 3568./243.* DL**3 - 184.*D81 * DL**4 + 40.2426 * Y1/Y
       C2S3F = ( ( 126.42 - 50.29 * Y - 50.15 * Y**2) * Y1 - 26.717 
     ,         - 320.*D81 * DL**2 * (DL+5.D0) + 59.59 * DL 
     ,         - Y*DL**2 * (101.8 + 34.79 * DL + 3.070 * DL**2) 
     ,         - 9.075 * Y*Y1*DL1 ) * Y
       C2S3A = NF * ( C2S31 + (FLS(NF)-FL(NF)) * C2S3F + NF * C2S32 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The (truncated) 'local' piece due to the FL11 contribution
*
       FUNCTION C2S3C (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
       DIMENSION FL(0:6), FLS(0:6)  !*** Extended to deal with NF=0  (FIO 2 Aug 2011) 
       DATA FL  /0.0d0,-1.d0, 0.5d0, 0.d0, 0.5d0, 0.2d0, 0.5d0 /
       DATA FLS /0.0d0,  1.d0, 0.1d0, 0.d0, 0.1d0, 0.018181818d0,0.1d0/
*
       FL11 = FL(NF)
       C2S3C = - (FLS(NF)-FL(NF)) * NF * 11.8880
*
       RETURN
       END
*
* =================================================================av==
*
* ..File: xc2ns3p.f    F2_NS
*
*
* ..Parametrization of the 3-loop MS(bar) non-singlet coefficient
*    functions for the structure function F_2 in electromagnetic DIS.
*    at  mu_r = mu_f = Q.  The expansion parameter is  alpha_s/(4 pi).
*
* ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to
*    the kernel superscripts [2], [3], and [1] in that equation.
*
*  ..The relative accuracy of these parametrizations, as well as of
*    the convolution results, is one part in thousand or better.
*
* ..References: [1] S. Moch, J. Vermaseren and A. Vogt, hep-ph/0209100
*               [2] J. Vermaseren, A. Vogt and S. Moch, hep-ph/0504242
*
*
* =====================================================================
*
*
* ..The regular piece. The rational end-point coefficients are exact, 
*    the rest has been fitted for x between 10^-6 and 1 - 10^-6. 
*
       FUNCTION C2NP3A (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
       DIMENSION FL(0:6)  !*** Extended to deal with NF=0  (FIO 2 Aug 2011) 
       DATA FL /0.0d0,-1.d0, 0.5d0, 0.d0, 0.5d0, 0.2d0, 0.5d0 /
*
       FL11 = FL(NF)
*
       Y1  = 1.D0 -Y
       DL  = LOG (Y)
       DL1 = LOG (Y1)
       D27  = 1./27.D0
       D243 = 1./243.D0
*
C-TS   Regular part of Eq.(4.11) in [2]. 
c...   Note - some of the fractions appear in different form such as: 8796/243 = 2932/81  
       C2NP3A = 
     ,            - 4926. + 7725.* Y + 57256.* Y**2 + 12898.* Y**3
     ,            - 32.*D27 * DL**5 - 8796.*D243 * DL**4 - 309.1 * DL**3
     ,            - 899.6 * DL**2 - 775.8 * DL + 4.719 * Y*DL**5
     ,            - 512.*D27 * DL1**5 + 6336.*D27 * DL1**4
     ,            - 3368.* DL1**3 - 2978.* DL1**2 + 18832.* DL1
     ,            - 56000.* (1.-Y)*DL1**2 - DL*DL1 * (6158. + 1836.*DL)
     ,        + NF * ( 831.6 - 6752.* Y - 2778.* Y**2
     ,            + 728.* D243 * DL**4 + 12224.* D243 * DL**3
     ,            + 187.3 * DL**2 + 275.6 * DL + 4.102 * Y*DL**4
     ,            - 1920.* D243 * DL1**4 + 153.5 * DL1**3 
     ,            - 828.7 * DL1**2 - 501.1 * DL1 + 171.0 * (1.-Y)*DL1**4
     ,            + DL*DL1 * (4365. + 716.2 * DL - 5983.* DL1) )
     ,        + NF**2 * ( 129.2 * Y + 102.5 * Y**2 - 368.* D243 * DL**3
     ,            - 1984.* D243 * DL**2 - 8.042 * DL
     ,            - 192.* D243 * DL1**3 + 18.21 * DL1**2 - 19.09 * DL1
     ,            + DL*DL1 * ( - 96.07 - 12.46 * DL + 85.88 * DL1) )
     ,        + FL11*NF * ( ( 126.42 - 50.29 * Y - 50.15 * Y**2) * Y1 
     ,           - 26.717 - 960.*D243 * DL**2 * (DL+5.D0) + 59.59 * DL
     ,           - Y*DL**2 * (101.8 + 34.79 * DL + 3.070 * DL**2)
     ,           - 9.075 * Y*Y1*DL1 ) * Y
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The exact singular piece (irrational coefficients truncated)
*
       FUNCTION C2NS3B (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
       DM  = 1./(1.-Y)
       D81 = 1./81.D0
*
C-TS    Dk part of Eq.(4.11) in Ref. [2], replacing Dk by DL1^k
c       Note - this part has to be convoluted as in first term on RHS of
c       Eq.(4.7)
       C2NS3B = 
     ,            + 1536.*D81 * DL1**5 - 16320.* D81 * DL1**4
     ,            + 5.01099E+2 * DL1**3 + 1.17154E+3 * DL1**2 
     ,            - 7.32845E+3 * DL1 + 4.44276E+3
     ,        + NF * ( 640.* D81 * DL1**4 - 6592.* D81 * DL1**3
     ,            + 220.573 * DL1**2 + 294.906 * DL1 - 729.359 )
     ,        + NF**2 * ( 64.* D81 * DL1**3 - 464.* D81 * DL1**2
     ,            + 7.67505 * DL1 + 1.00830 )
*
       C2NS3B = DM * C2NS3B
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The 'local' piece.  The coefficients of delta(1-x) have been 
*    slightly shifted with respect to their (truncated) exact values.  
*
       FUNCTION C2NP3C (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
       DIMENSION FL(0:6)  !*** Extended to deal with NF=0  (FIO 2 Aug 2011) 
       DATA FL /0.0d0,-1.d0, 0.5d0, 0.d0, 0.5d0, 0.2d0, 0.5d0 /
*
       FL11 = FL(NF)
*
       DL1 = LOG (1.-Y)
       D81 = 1./81.D0
       D3  = 1./3.D0
*
C-TS   delta function part and part of the Dk (see 2nd term on RHS of Eq.(4.7) ) 
c      of Eq.(4.11) in Ref. [2]
       C2NP3C = 
     ,            + 256.*D81 * DL1**6 - 3264.*D81 * DL1**5
     ,            + 1.252745E+2 * DL1**4 + 3.905133E+2 * DL1**3 
     ,            - 3.664225E+3 * DL1**2 + 4.44276E+3  * DL1
     ,            - 9195.48 + 25.10 
     ,        + NF * ( 128.* D81 * DL1**5 - 1648.* D81 * DL1**4
     ,            + 220.573 * D3 * DL1**3 + 147.453 * DL1**2
     ,            - 729.359 * DL1 + 2575.074 - 0.387 )
     ,        + NF**2 * ( 16.* D81 * DL1**4 - 464.* D81*D3 * DL1**3
     ,            + 7.67505 * 5.D-1 * DL1**2 + 1.0083 * DL1 - 103.2521
     ,            + 0.0155 )
     ,        - FL11*NF * 11.8880
*
       RETURN
       END
















































*
* ---------------------------------------------------------------------
*
*
* ..The gluon coefficient function
*
       FUNCTION C2G3A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       DIMENSION FLG(0:6)  !*** Extended to deal with NF=0  (FIO 2 Aug 2011) 
       INTEGER NF
       DATA FLG /0.0d0,1.d0,0.1d0,0.d0,0.1d0,0.018181818d0,0.1d0/
*
       YI  = 1./Y
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
       D9  = 1./9.D0
       D81 = D9*D9
*
       C2G31 = 
     ,           966.*D81 * DL1**5 - 935.5*D9 * DL1**4 + 89.31 * DL1**3 
     ,         + 979.2 * DL1**2 - 2405. * DL1 + 1372.* (1.-Y)* DL1**4
     ,         - 15729. - 310510.* Y + 331570.* Y**2 - 244150.* Y*DL**2
     ,         - 253.3* Y*DL**5
     ,         + DL*DL1 * (138230. - 237010.* DL) - 11860.* DL 
     ,         - 700.8 * DL**2 - 1440.* DL**3 + 2480.5*D81 * DL**4
     ,         - 134.*D9 * DL**5 - 6362.54 * YI - 932.089 * DL*YI
       C2G32 = 
     ,           131.*D81 * DL1**4 - 14.72 * DL1**3 + 3.607 * DL1**2
     ,         - 226.1 * DL1 + 4.762 - 190.0 * Y - 818.4 * Y**2
     ,         - 4019.* Y*DL**2 - DL*DL1 * (791.5 + 4646 * DL)
     ,         + 739.0 * DL + 418.0 * DL**2 + 104.3 * DL**3 
     ,         + 809.*D81 * DL**4 + 12.*D9 * DL**5 + 84.423 * YI
       C2G3F =   3.211 * DL1**2 + 19.04 * Y*DL1 + 0.623 * (1.-Y)*DL1**3 
     ,         - 64.47 * Y + 121.6 * Y**2 - 45.82 * Y**3 - Y*DL*DL1 
     ,         * ( 31.68 + 37.24 * DL) - Y*DL * (82.40 + 16.08 * DL)
     ,         + Y*DL**3 * (520.*D81 + 11.27 * Y) + 60.*D81 * Y*DL**4

       C2G3A = NF * ( C2G31 + NF * (C2G32 + FLG(NF) * C2G3F) )
*
       RETURN
       END
