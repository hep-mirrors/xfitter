C Math constants:
      double precision pi
! [--- WS 2016-01-14
      ! 21 digits --- double precision standard
      parameter (pi = 3.14159265358979323846d0)
      ! 35 digits --- gnuC standard for quad M_PIq
      ! parameter (pi = 3.1415926535897932384626433832795029d0)
      
      ! --- the Euler constant gamma_E
      double precision gammaEuler
      parameter (gammaEuler = 0.577215664901532860607d0)
      ! parameter (gammaEuler = 0.5772156649015328606065120900824024d0)
! ---]

C Electric charges:
      double precision euq,edq, e2u,e2d 
      parameter (euq =  2.d0/3.)
      parameter (edq = -1.d0/3.)
      parameter (e2u = euq*euq)
      parameter (e2d = edq*edq)      

C Electroweak common blocks
      ! Physical constants defined in ewparam.txt
      ! \c gf = \f$ G_{\rm Fermi}/(\hbar c)^3 [GeV^{-2}] \approx 1.166\cdot 10^{-5}\f$
      ! \c convfac = \f$ (\hbar c)^2 [GeV^2 pbarn] \approx 0.389\cdot 10^{9}\f$
      double precision gf, convfac
      double precision alphaem,sin2thw,cos2thw
      double precision Mz, Mw, Mh
      double precision wz, ww, wh, wtp
      double precision Vud, Vus, Vub, Vcd, Vcs, Vcb, Vtd, Vts, Vtb
      double precision men, mel, mmn, mmo, mtn, mta, mup, mdn,
     + mch, mst, mtp, mbt

      common/constants/ gf, convfac
      common/ew_couplings/alphaem,sin2thw,cos2thw
      common/boson_masses/ Mz, Mw, Mh
      common/widths/wz, ww, wh, wtp
      common/ckm_matrix/ Vud, Vus, Vub, Vcd, Vcs, Vcb, Vtd, Vts, Vtb
      common/fermion_masses/men, mel, mmn, mmo, mtn, mta, mup, mdn,
     + mch, mst, mtp, mbt
