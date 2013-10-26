C  pi :
      double precision pi
      parameter (pi = 3.141592653589793D0)

C Electric charges:
      double precision euq,edq, e2u,e2d 
      parameter (euq =  2./3.)
      parameter (edq = -1./3.)
      parameter (e2u = euq*euq)
      parameter (e2d = edq*edq)      

C Electroweak common blocks
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
      
