******************************************************************************
**  hplog: a subroutine for the evaluation of harmonic polylogarithms
**         Version 1.1         26/10/2004
**  described in: 
**  T.Gehrmann and E.Remiddi: Numerical Evaluation of the Harmonic 
**                            Polylogarithms up to Weight 4
**                            (hep-ph/0107173; Comp.Phys.Comm. 141 (2001) 296)
**  the harmonic polylogarithms are defined in: 
**  E.Remiddi and J.Vermaseren: Harmonic Polylogarithms
**                            (hep-ph/9905237; Int.J.Mod.Phys. A15 (2000) 725)
**  email:
**  Thomas.Gehrmann@cern.ch and Ettore.Remiddi@bo.infn.it
**
**  changes with respect to version 1.0 (12/07/2001):
**  improved Chebyshev expansions by factoring out leading behaviour
**  to improve the relative acurracy for very small values of the arguments
******************************************************************************
      subroutine hplog(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                      Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
****** 
** x is the argument of the 1dHPL's (1 dimensional Harmonic PolyLogarithms) 
**   to be evaluated; 
** nw is the maximum weight of the required 1dHPL's; 
**    the maximum allowed value of nw of this implementation is 4; 
** Hc1,Hc2,Hc3,Hc4 are the complex*16 values of the 1dHPL; 
**    they must all be supplied in the arguments even if some of them 
**    are not to be evaluated; 
** Hr1,Hr2,Hr3,Hr4 are the double precision real parts of 
**    Hc1,Hc2,Hc3,Hc4; 
** Hi1,Hi2,Hi3,Hi4 are the double precision immaginary parts of 
**    Hc1,Hc2,Hc3,Hc4 divided by pi=3.114159.... 
** n1,n2 is the required range of indices, the allowed ranges are 
**    (0,1), (-1,0), (-1,1) ; 
****** 
      implicit double precision (a-h,o-z) 
      complex*16 Hc1,Hc2,Hc3,Hc4 
      dimension Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), 
     $          Hc4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2), 
     $          Hr4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
      common /fillred/infilldim,infill(3) 
      parameter (r2   = 1.4142135623730950488d0) 
** check on the weight nw 
      if ( (nw.lt.1).or.(nw.gt.4) ) then 
        print*, ' illegal call of eval1dhpl with second argument', 
     $          ' (the weight) = ',nw 
        print*, ' the allowed values of the weight are 1,2,3,4 ' 
        stop
      endif
** check on the range n1:n2 
      if ( (n1.eq.-1).and.(n2.eq.0) ) then 
        infilldim =  2 
        infill(1) =  0 
        infill(2) = -1  
      elseif ( (n1.eq.0).and.(n2.eq.1) ) then 
        infilldim =  2 
        infill(1) =  0 
        infill(2) =  1  
      elseif ( (n1.eq.-1).and.(n2.eq.1) ) then 
        infilldim =  3 
        infill(1) =  0 
        infill(2) = -1  
        infill(3) =  1  
      else 
        print*, ' illegal call of eval1dhpl with the two last ', 
     $          'arguments = (',n1,',',n2,')' 
        print*, ' the allowed values are (-1,0), (0,1), (-1,1) ' 
        stop 
      endif 
** setting the immaginary parts equal to zero 
      call setzero(nw,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** looking at the range of the argument 
*      r2 = sqrt(2.d0) 
      r2m1 = r2 - 1 
      r2p1 = r2 + 1 
      if ( ( x.gt.-r2m1 ).and.( x.le.r2m1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplat0 ' 
        call eval1dhplat0(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( x.eq.1d0 ) then
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplin1 ' 
        call eval1dhplin1(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( ( x.gt.r2m1 ).and.( x.le.r2p1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplat1 ' 
        call eval1dhplat1(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( ( x.gt.r2p1 ) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplatinf ' 
        call eval1dhplatinf(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( ( x.le.-r2p1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplatminf ' 
        call eval1dhplatminf(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( x.eq.-1d0 ) then
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplinm1 ' 
        call eval1dhplinm1(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( ( x.gt.-r2p1 ).and.( x.le.-r2m1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplatm1 ' 
        call eval1dhplatm1(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      endif 
** 
      end 
************************************************************************ 
      subroutine eval1dhplat0(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the 0-range  -(r2-1) < y <= (r2-1) 
** by direct series expansion (Bernoulli-accelerated) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** evaluate the irreducible 1dHPL's first 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
      call fillirr1dhplat0(y,nw,HY1,HY2,HY3,HY4,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      return 
      end 
************************************************************************ 
      subroutine eval1dhplin1(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's for y=1 (explicit values are tabulated)
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
      parameter (pi   = 3.14159265358979324d0) 
** evaluate the irreducible 1dHPL's first 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
      call fillirr1dhplin1(y,nw,HY1,HY2,HY3,HY4,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      if (n2.eq.0) return
** correct the ill-defined entries
      HY2(1,0) = - HY2(0,1) 
      Hi2(1,0) = 0d0 
      H2(1,0) = dcmplx(HY2(1,0),Hi2(1,0)*pi)
      if ( nw.eq.2 ) return
      HY3(1,0,0) = HY3(0,0,1) 
      Hi3(1,0,0) = 0d0 
      H3(1,0,0) = dcmplx(HY3(1,0,0),Hi3(1,0,0)*pi)
      if ( nw.eq.3 ) return
      HY4(1,0,0,0) = -HY4(0,0,0,1) 
      Hi4(1,0,0,0) = 0d0 
      H4(1,0,0,0) = dcmplx(HY4(1,0,0,0),Hi4(1,0,0,0)*pi)
      return 
      end 
************************************************************************ 
      subroutine eval1dhplat1(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the 1-range  (r2-1) < y <= (r2+1) 
** evaluating first the H(..,r=(1-y)/(1+y)) by calling eval1dhplat0(r)  
** and then expressing H(..,y=(1-r)/(1+r)) in terms of H(..,r) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension HR1(-1:1),HR2(-1:1,-1:1),HR3(-1:1,-1:1,-1:1), 
     $          HR4(-1:1,-1:1,-1:1,-1:1) 
** the nw = 1 case 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
** the nw > 1 case 
      r = (1.d0-y)/(1.d0+y) 
*      print*,' eval1dhplat1: y = ',y,', r = ',r 
** the whole (-1,1) range is in general needed for any pair (n1,n2)
      call fillirr1dhplat0(r,nw,HR1,HR2,HR3,HR4,-1,1) 
** fillirr1dhplat1 takes care automatically of all the immaginary 
** parts as well as of the jump across y=1 
      call fillirr1dhplat1(r,nw,HR1,HR2,HR3,HR4, 
     $                          HY1,HY2,HY3,HY4, 
     $                          Hi1,Hi2,Hi3,Hi4,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      return 
      end 
************************************************************************ 
      subroutine eval1dhplatinf(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the inf-range  (r2+1) < abs(y) 
** evaluating first the H(..,x=1/y) by calling eval1dhplat0(x)  
** and then expressing H(..,y=1/x) in terms of H(..,x) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension HX1(n1:n2),HX2(n1:n2,n1:n2),HX3(n1:n2,n1:n2,n1:n2), 
     $          HX4(n1:n2,n1:n2,n1:n2,n1:n2) 
      parameter (pi   = 3.14159265358979324d0) 
** the nw = 1 case 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
** the nw > 1 case 
      x = 1.d0/y 
*      print*,' eval1dhplatinf: y = ',y,', x = ',x 
      call fillirr1dhplat0(x,nw,HX1,HX2,HX3,HX4,n1,n2) 
** fillirr1dhplatinf takes care automatically of all the immaginary 
** parts as well as of the jump across y=1 
      call fillirr1dhplatinf(x,nw,HX1,HX2,HX3,HX4, 
     $                            HY1,HY2,HY3,HY4, 
     $                            Hi1,Hi2,Hi3,Hi4,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      return 
      end 
************************************************************************ 
      subroutine eval1dhplinm1(y,nw,H1,H2,H3,H4, 
     $                           HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's for y=-1 (explicit values are tabulated)
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      complex*16 G1,G2,G3,G4  
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension G1(-n2:-n1),G2(-n2:-n1,-n2:-n1),
     $          G3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension GY1(-n2:-n1),GY2(-n2:-n1,-n2:-n1),
     $          GY3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension Gi1(-n2:-n1),Gi2(-n2:-n1,-n2:-n1),
     $          Gi3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      common /fillred/infilldim,infill(3) 
      dimension istorfill(3)
      dimension nphase(-1:1) 
      data nphase/-1,1,-1/ 
      parameter (pi   = 3.14159265358979324d0) 
*      print*,' eval1dhplatm1: y = ',y 
      if (infilldim.eq.2) then
         do i=1,2
            istorfill(i) = infill(i)
            infill(i) = -istorfill(i)
         enddo
      endif
** evaluate H(...,-y) 
      call setzero(nw,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      Gi1(0) = -1
      call eval1dhplin1(-y,nw,G1,G2,G3,G4, 
     $                        GY1,GY2,GY3,GY4,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      if (infilldim.eq.2) then
         do i=1,2
            infill(i) = istorfill(i)
         enddo
      endif
** fill the arrays H's 
      do k1=n1,n2 
        nph1 = nphase(k1) 
        HY1(k1) =   nph1*GY1(-k1) 
        Hi1(k1) = - nph1*Gi1(-k1) 
        H1(k1)  = dcmplx(HY1(k1),Hi1(k1)*pi) 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            nph2 = nph1*nphase(k2) 
            HY2(k1,k2) =   nph2*GY2(-k1,-k2) 
            Hi2(k1,k2) = - nph2*Gi2(-k1,-k2) 
            H2(k1,k2)  = dcmplx(HY2(k1,k2),Hi2(k1,k2)*pi) 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                nph3 = nph2*nphase(k3) 
                HY3(k1,k2,k3) =   nph3*GY3(-k1,-k2,-k3) 
                Hi3(k1,k2,k3) = - nph3*Gi3(-k1,-k2,-k3) 
                H3(k1,k2,k3)  = dcmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi) 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    nph4 = nph3*nphase(k4) 
                    HY4(k1,k2,k3,k4) =   nph4*GY4(-k1,-k2,-k3,-k4) 
                    Hi4(k1,k2,k3,k4) = - nph4*Gi4(-k1,-k2,-k3,-k4) 
                    H4(k1,k2,k3,k4)  = 
     $                    dcmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi) 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      if (n1.eq.0) return
** correct the ill-defined entries
      HY2(-1,0) = - HY2(0,-1) 
      Hi2(-1,0) = Hi1(0)*HY1(-1)
      H2(-1,0) = dcmplx(HY2(-1,0),Hi2(-1,0)*pi)
      if ( nw.eq.2 ) return
      HY3(-1,0,0) = HY1(-1)*HY2(0,0)+HY3(0,0,-1) 
      Hi3(-1,0,0) = HY1(-1)*Hi2(0,0)-HY2(0,-1)*Hi1(0)
      H3(-1,0,0) = dcmplx(HY3(-1,0,0),Hi3(-1,0,0)*pi)
      if ( nw.eq.3 ) return
      HY4(-1,0,0,0) = -HY2(0,-1)*HY2(0,0)-HY4(0,0,0,-1) 
      Hi4(-1,0,0,0) = HY1(-1)*Hi3(0,0,0)+HY3(0,0,-1)*Hi1(0)
      H4(-1,0,0,0) = dcmplx(HY4(-1,0,0,0),Hi4(-1,0,0,0)*pi)
      return 
      end 
************************************************************************ 
      subroutine eval1dhplatm1(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the (-1)-range  -(r2+1) < y <= -(r2-1) 
** evaluating first the H(..,-y) by calling eval1dhplat1(-y), 
** and then expressing H(..,y) in terms of H(..,-y) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4  
      complex*16 G1,G2,G3,G4  
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension G1(-n2:-n1),G2(-n2:-n1,-n2:-n1),
     $          G3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension GY1(-n2:-n1),GY2(-n2:-n1,-n2:-n1),
     $          GY3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension Gi1(-n2:-n1),Gi2(-n2:-n1,-n2:-n1),
     $          Gi3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
** 
      common /fillred/infilldim,infill(3) 
      dimension istorfill(3)
      dimension nphase(-1:1) 
      data nphase/-1,1,-1/ 
      parameter (pi   = 3.14159265358979324d0) 
*      print*,' eval1dhplatm1: y = ',y 
      if (infilldim.eq.2) then
         do i=1,2
            istorfill(i) = infill(i)
            infill(i) = -istorfill(i)
         enddo
      endif
** evaluate H(...,-y) 
      call setzero(nw,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      Gi1(0) = -1
      call eval1dhplat1(-y,nw,G1,G2,G3,G4, 
     $                        GY1,GY2,GY3,GY4,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      if (infilldim.eq.2) then
         do i=1,2
            infill(i) = istorfill(i)
         enddo
      endif
** fill the arrays H's 
      do k1=n1,n2 
        nph1 = nphase(k1) 
        HY1(k1) =   nph1*GY1(-k1) 
        Hi1(k1) = - nph1*Gi1(-k1) 
        H1(k1)  = dcmplx(HY1(k1),Hi1(k1)*pi) 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            nph2 = nph1*nphase(k2) 
            HY2(k1,k2) =   nph2*GY2(-k1,-k2) 
            Hi2(k1,k2) = - nph2*Gi2(-k1,-k2) 
            H2(k1,k2)  = dcmplx(HY2(k1,k2),Hi2(k1,k2)*pi) 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                nph3 = nph2*nphase(k3) 
                HY3(k1,k2,k3) =   nph3*GY3(-k1,-k2,-k3) 
                Hi3(k1,k2,k3) = - nph3*Gi3(-k1,-k2,-k3) 
                H3(k1,k2,k3)  = dcmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi) 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    nph4 = nph3*nphase(k4) 
                    HY4(k1,k2,k3,k4) =   nph4*GY4(-k1,-k2,-k3,-k4) 
                    Hi4(k1,k2,k3,k4) = - nph4*Gi4(-k1,-k2,-k3,-k4) 
                    H4(k1,k2,k3,k4)  = 
     $                    dcmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi) 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      return 
      end 


      subroutine eval1dhplatminf(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the (-1)-range y  <= -(r2+1) 
** evaluating first the H(..,-y) by calling eval1dhplatinf(-y), 
** and then expressing H(..,y) in terms of H(..,-y) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4  
      complex*16 G1,G2,G3,G4  
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension G1(-n2:-n1),G2(-n2:-n1,-n2:-n1),
     $          G3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension GY1(-n2:-n1),GY2(-n2:-n1,-n2:-n1),
     $          GY3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension Gi1(-n2:-n1),Gi2(-n2:-n1,-n2:-n1),
     $          Gi3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
** 
      common /fillred/infilldim,infill(3) 
      dimension istorfill(3)
      dimension nphase(-1:1) 
      data nphase/-1,1,-1/ 
      parameter (pi   = 3.14159265358979324d0) 
*      print*,' eval1dhplatm1: y = ',y 
      if (infilldim.eq.2) then
         do i=1,2
            istorfill(i) = infill(i)
            infill(i) = -istorfill(i)
         enddo
      endif
** evaluate H(...,-y) 
      call setzero(nw,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      Gi1(0) = -1
      call eval1dhplatinf(-y,nw,G1,G2,G3,G4, 
     $                        GY1,GY2,GY3,GY4,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      if (infilldim.eq.2) then
         do i=1,2
            infill(i) = istorfill(i)
         enddo
      endif
** fill the arrays H's 
      do k1=n1,n2 
        nph1 = nphase(k1) 
        HY1(k1) =   nph1*GY1(-k1) 
        Hi1(k1) = - nph1*Gi1(-k1) 
        H1(k1)  = dcmplx(HY1(k1),Hi1(k1)*pi) 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            nph2 = nph1*nphase(k2) 
            HY2(k1,k2) =   nph2*GY2(-k1,-k2) 
            Hi2(k1,k2) = - nph2*Gi2(-k1,-k2) 
            H2(k1,k2)  = dcmplx(HY2(k1,k2),Hi2(k1,k2)*pi) 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                nph3 = nph2*nphase(k3) 
                HY3(k1,k2,k3) =   nph3*GY3(-k1,-k2,-k3) 
                Hi3(k1,k2,k3) = - nph3*Gi3(-k1,-k2,-k3) 
                H3(k1,k2,k3)  = dcmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi) 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    nph4 = nph3*nphase(k4) 
                    HY4(k1,k2,k3,k4) =   nph4*GY4(-k1,-k2,-k3,-k4) 
                    Hi4(k1,k2,k3,k4) = - nph4*Gi4(-k1,-k2,-k3,-k4) 
                    H4(k1,k2,k3,k4)  = 
     $                    dcmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi) 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      return 
      end 
************************************************************************ 
      subroutine setzero(nw,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** initializes with 0 the elements of the arrays 
      implicit double precision (a-h,o-z) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
      do k1=n1,n2 
        Hi1(k1) = 0.d0 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            Hi2(k1,k2) = 0.d0 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                Hi3(k1,k2,k3) = 0.d0 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    Hi4(k1,k2,k3,k4) = 0.d0 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      return 
      end 
************************************************************************ 
      subroutine fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
* fills the reducible 1dhpl from the irreducible set
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
      common /fillred/infilldim,infill(3) 
      parameter (pinv = 0.318309886183790672d0) 
      parameter (pi   = 3.14159265358979324d0) 
** combining real and immaginary into the complex value 
      do k1=n1,n2 
      do k2=n1,n2 
        H2(k1,k2) = dcmplx(HY2(k1,k2),Hi2(k1,k2)*pi) 
        if ( nw.gt.2 ) then 
          do k3=n1,n2 
            H3(k1,k2,k3) = dcmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi) 
            if ( nw.gt.3 ) then 
              do k4=n1,n2 
                H4(k1,k2,k3,k4) = 
     $               dcmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi) 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      enddo 
** evaluating the reduced HPL's 
** iflag = 0 to suppress auxiliary printings of FILLREDHPLx 
      iflag = 0 
      do ia =  1,infilldim 
      do ib = ia,infilldim 
        call FILLREDHPL2(iflag,H1,H2,n1,n2,infill(ia),infill(ib)) 
        if ( nw.gt.2 ) then 
          do ic = ib,infilldim 
            call FILLREDHPL3(iflag,H1,H2,H3,n1,n2, 
     $                          infill(ia),infill(ib),infill(ic)) 
            if ( nw.gt.3 ) then 
              do id = ic,infilldim 
                call FILLREDHPL4(iflag,H1,H2,H3,H4,n1,n2, 
     $               infill(ia),infill(ib),infill(ic),infill(id)) 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      enddo 
** extractin real and immaginary parts from the complex value 
      do k1=n1,n2 
      do k2=n1,n2 
        HY2(k1,k2) =  dble(H2(k1,k2)) 
        Hi2(k1,k2) = dimag(H2(k1,k2))*pinv 
        if ( nw.gt.2 ) then 
          do k3=n1,n2 
            HY3(k1,k2,k3) =  dble(H3(k1,k2,k3)) 
            Hi3(k1,k2,k3) = dimag(H3(k1,k2,k3))*pinv 
            if ( nw.gt.3 ) then 
              do k4=n1,n2 
                HY4(k1,k2,k3,k4) =  dble(H4(k1,k2,k3,k4)) 
                Hi4(k1,k2,k3,k4) = dimag(H4(k1,k2,k3,k4))*pinv 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      enddo 
      return 
      end 
************************************************************************ 
      subroutine FILLREDHPL2(iflag,H1,H2,i1,i2,na,nb) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2 
      dimension H1(i1:i2),H2(i1:i2,i1:i2) 
*23456789012345678901234567890123456789012345678901234567890123456789012 
* must be called with ordered indices na <= nb 
*      print*,' FILLREDHPL2, iflag =',iflag 
      if ( na.eq.nb ) then 
        H2(na,na) = 1.d0/2*( H1(na) )**2 
      else 
        H2(nb,na) = + H1(na)*H1(nb) - H2(na,nb) 
        if ( iflag.eq.1 ) then 
          call printer2(na,nb) 
        endif 
      endif 
      return 
      end 
************************************************************************ 
      subroutine FILLREDHPL3(iflag,H1,H2,H3,i1,i2,ia,ib,ic) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3 
      dimension H1(i1:i2),H2(i1:i2,i1:i2),H3(i1:i2,i1:i2,i1:i2) 
* must be called with "properly ordered" indices 
* note in particular the remapping of, say, (ia,ib,ic) into 
* (na,na,nb) of ReducerTest.out 
      na = ia 
      if ( (ia.eq.ib).and.(ib.eq.ic) ) then 
* case (na,na,na) 
        H3(na,na,na) = 1.d0/6*( H1(na) )**3 
* ic cannot be anymore equal to ia 
      else if ( ic.eq.ia ) then 
        print*,' FILLREDHPL3, error 1, called with arguments ' 
        print*,'               ',ia,ib,ic 
        stop 
      else if ( ia.eq.ib ) then 
* case (na,na,nb) 
        nb = ic 
        if ( iflag.eq.1 ) then 
          call printer3(na,na,nb) 
        endif 
        H3(na,nb,na) = + H1(na)*H2(na,nb) - 2*H3(na,na,nb) 
        H3(nb,na,na) = + 1.d0/2*H1(na)*H1(na)*H1(nb) 
     $                 - H1(na)*H2(na,nb) + H3(na,na,nb) 
* ib cannot be anymore equal to ia 
      else if ( ib.eq.ia ) then 
        print*,' FILLREDHPL3, error 2, called with arguments ' 
        print*,'               ',ia,ib,ic 
        stop 
      else if ( ib.eq.ic ) then 
* case (na,nb,nb) 
        nb = ib 
        if ( iflag.eq.1 ) then 
          call printer3(na,nb,nb) 
        endif 
        H3(nb,na,nb) = + H1(nb)*H2(na,nb) - 2*H3(na,nb,nb) 
        H3(nb,nb,na) = + 1.d0/2*H1(na)*H1(nb)*H1(nb) 
     $                 - H1(nb)*H2(na,nb) + H3(na,nb,nb) 
* no need to protect against ic.eq.ib 
* when arriving here all indices are different 
      else 
* case (na,nb,nc)    all indices are different 
        nb = ib 
        nc = ic 
        if ( iflag.eq.1 ) then 
          call printer3(na,nb,nc) 
          call printer3(na,nc,nb) 
        endif 
        H3(nb,na,nc) = + H1(nb)*H2(na,nc) 
     $                 - H3(na,nb,nc) - H3(na,nc,nb) 
        H3(nb,nc,na) = + H1(na)*H2(nb,nc) 
     $                 - H1(nb)*H2(na,nc) + H3(na,nc,nb) 
        H3(nc,na,nb) = + H1(nc)*H2(na,nb) 
     $                 - H3(na,nb,nc) - H3(na,nc,nb) 
        H3(nc,nb,na) = + H1(na)*H1(nb)*H1(nc) - H1(na)*H2(nb,nc) 
     $                 - H1(nc)*H2(na,nb) + H3(na,nb,nc) 
      endif 
*23456789012345678901234567890123456789012345678901234567890123456789012 
      return 
      end 
************************************************************************ 
      subroutine FILLREDHPL4(iflag,H1,H2,H3,H4,i1,i2,ia,ib,ic,id) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(i1:i2),H2(i1:i2,i1:i2),H3(i1:i2,i1:i2,i1:i2) 
      dimension H4(i1:i2,i1:i2,i1:i2,i1:i2) 
*23456789012345678901234567890123456789012345678901234567890123456789012 
* must be called with "properly ordered" indices 
* note in particular the remapping of, say, (ia,ib,ic) into 
* (na,na,nb) of ReducerTest.out 
      na = ia 
      if ( (ia.eq.ib).and.(ib.eq.ic).and.(ic.eq.id) ) then 
* case (na,na,na,na) 
        H4(na,na,na,na) = 1.d0/24*( H1(na) )**4 
* id cannot be anymore equal to ia 
      else if ( id.eq.ia ) then 
        print*,' FILLREDHPL4, error 1, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( (ia.eq.ib).and.(ib.eq.ic) ) then 
* case (na,na,na,nb) 
        nb = id 
        H4(na,na,nb,na) = + H1(na)*H3(na,na,nb) - 3*H4(na,na,na,nb) 
        H4(na,nb,na,na) = + 1.d0/2*H1(na)*H1(na)*H2(na,nb) 
     $                    - 2*H1(na)*H3(na,na,nb) + 3*H4(na,na,na,nb) 
        H4(nb,na,na,na) = + 1.d0/6*H1(na)*H1(na)*H1(na)*H1(nb) 
     $                    - 1.d0/2*H1(na)*H1(na)*H2(na,nb) 
     $                    + H1(na)*H3(na,na,nb) - H4(na,na,na,nb) 
        if ( iflag.eq.1 ) then 
          call printer4(na,na,na,nb) 
        endif 
* ic cannot be anymore equal to ia 
      else if ( ic.eq.ia ) then 
        print*,' FILLREDHPL4, error 2, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( (ia.eq.ib).and.(ic.eq.id) ) then 
* case (na,na,nb,nb) 
        nb = ic 
        H4(na,nb,na,nb) = + 1.d0/2*H2(na,nb)*H2(na,nb) 
     $                    - 2*H4(na,na,nb,nb) 
        H4(na,nb,nb,na) = + H1(na)*H3(na,nb,nb) 
     $                    - 1.d0/2*H2(na,nb)*H2(na,nb) 
        H4(nb,na,na,nb) = + H1(nb)*H3(na,na,nb) 
     $                    - 1.d0/2*H2(na,nb)*H2(na,nb) 
        H4(nb,na,nb,na) = + H1(na)*H1(nb)*H2(na,nb) 
     $                    - 2*H1(na)*H3(na,nb,nb) 
     $                    - 2*H1(nb)*H3(na,na,nb) 
     $                    + 1.d0/2*H2(na,nb)*H2(na,nb) 
     $                    + 2*H4(na,na,nb,nb) 
        H4(nb,nb,na,na) = + 1.d0/4*H1(na)*H1(na)*H1(nb)*H1(nb) 
     $                    - H1(na)*H1(nb)*H2(na,nb) 
     $                    + H1(na)*H3(na,nb,nb) 
     $                    + H1(nb)*H3(na,na,nb) - H4(na,na,nb,nb) 
        if ( iflag.eq.1 ) then 
          call printer4(na,na,nb,nb) 
        endif 
      else if ( ia.eq.ib ) then 
* case (na,na,nb,nc) 
        nb = ic 
        nc = id 
        H4(na,nb,nc,na) = + H1(na)*H3(na,nb,nc) - 2*H4(na,na,nb,nc) 
     $                    - H4(na,nb,na,nc) 
        H4(na,nc,na,nb) = + H2(na,nb)*H2(na,nc) - 2*H4(na,na,nb,nc) 
     $                    - 2*H4(na,na,nc,nb) - H4(na,nb,na,nc) 
        H4(na,nc,nb,na) = + H1(na)*H3(na,nc,nb) - H2(na,nb)*H2(na,nc) 
     $                    + 2*H4(na,na,nb,nc) + H4(na,nb,na,nc) 
        H4(nb,na,na,nc) = + H1(nb)*H3(na,na,nc) - H4(na,na,nb,nc) 
     $                    - H4(na,na,nc,nb) - H4(na,nb,na,nc) 
        H4(nb,na,nc,na) = + H1(na)*H1(nb)*H2(na,nc) 
     $                    - H1(na)*H3(na,nb,nc) - H1(na)*H3(na,nc,nb) 
     $                    - 2*H1(nb)*H3(na,na,nc) + 2*H4(na,na,nb,nc) 
     $                    + 2*H4(na,na,nc,nb) + H4(na,nb,na,nc) 
        H4(nb,nc,na,na) = + 1.d0/2*H1(na)*H1(na)*H2(nb,nc) 
     $                    - H1(na)*H1(nb)*H2(na,nc) 
     $                    + H1(na)*H3(na,nc,nb) + H1(nb)*H3(na,na,nc) 
     $                    - H4(na,na,nc,nb) 
        H4(nc,na,na,nb) = + H1(nc)*H3(na,na,nb) - H2(na,nb)*H2(na,nc) 
     $                    + H4(na,na,nb,nc) + H4(na,na,nc,nb) 
     $                    + H4(na,nb,na,nc) 
        H4(nc,na,nb,na) = + H1(na)*H1(nc)*H2(na,nb) 
     $                    - H1(na)*H3(na,nb,nc) - H1(na)*H3(na,nc,nb) 
     $                    - 2*H1(nc)*H3(na,na,nb) + H2(na,nb)*H2(na,nc) 
     $                    - H4(na,nb,na,nc) 
        H4(nc,nb,na,na) = + 1.d0/2*H1(na)*H1(na)*H1(nb)*H1(nc) 
     $                    - 1.d0/2*H1(na)*H1(na)*H2(nb,nc) 
     $                    - H1(na)*H1(nc)*H2(na,nb) 
     $                    + H1(na)*H3(na,nb,nc) + H1(nc)*H3(na,na,nb) 
     $                    - H4(na,na,nb,nc)  
        if ( iflag.eq.1 ) then 
          call printer4(na,na,nb,nc) 
          call printer4(na,na,nc,nb) 
          call printer4(na,nb,na,nc) 
        endif 
* ib cannot be anymore equal to ia 
      else if ( ib.eq.ia ) then 
        print*,' FILLREDHPL4, error 3, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( (ib.eq.ic).and.(ic.eq.id) ) then 
* case (na,nb,nb,nb) 
        nb = ib 
        H4(nb,na,nb,nb) = + H1(nb)*H3(na,nb,nb) - 3*H4(na,nb,nb,nb) 
        H4(nb,nb,na,nb) = + 1.d0/2*H1(nb)*H1(nb)*H2(na,nb) 
     $                    - 2*H1(nb)*H3(na,nb,nb) + 3*H4(na,nb,nb,nb) 
        H4(nb,nb,nb,na) = + 1.d0/6*H1(na)*H1(nb)*H1(nb)*H1(nb) 
     $                    - 1.d0/2*H1(nb)*H1(nb)*H2(na,nb) 
     $                    + H1(nb)*H3(na,nb,nb) - H4(na,nb,nb,nb) 
        if ( iflag.eq.1 ) then 
          call printer4(na,nb,nb,nb) 
        endif 
* id cannot be anymore equal to ib 
      else if ( id.eq.ib ) then 
        print*,' FILLREDHPL4, error 4, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( ib.eq.ic ) then 
* case (na,nb,nb,nc) 
        nb = ib 
        nc = id 
        H4(nb,na,nb,nc) = + H1(nb)*H3(na,nb,nc) 
     $                    - 2*H4(na,nb,nb,nc) - H4(na,nb,nc,nb) 
        H4(nb,na,nc,nb) = + H1(nb)*H3(na,nc,nb) - H4(na,nb,nc,nb) 
     $                    - 2*H4(na,nc,nb,nb) 
        H4(nb,nb,na,nc) = + 1.d0/2*H1(nb)*H1(nb)*H2(na,nc) 
     $                    - H1(nb)*H3(na,nb,nc) - H1(nb)*H3(na,nc,nb) 
     $                    + H4(na,nb,nb,nc) + H4(na,nb,nc,nb) 
     $                    + H4(na,nc,nb,nb) 
        H4(nb,nb,nc,na) = + H1(na)*H3(nb,nb,nc) 
     $                    - 1.d0/2*H1(nb)*H1(nb)*H2(na,nc) 
     $                    + H1(nb)*H3(na,nc,nb) - H4(na,nc,nb,nb) 
        H4(nb,nc,na,nb) = - H1(nb)*H3(na,nb,nc) - H1(nb)*H3(na,nc,nb) 
     $                    + H2(na,nb)*H2(nb,nc) + H4(na,nb,nc,nb) 
     $                    + 2*H4(na,nc,nb,nb) 
        H4(nb,nc,nb,na) = + H1(na)*H1(nb)*H2(nb,nc) 
     $                    - 2*H1(na)*H3(nb,nb,nc) 
     $                    + H1(nb)*H3(na,nb,nc) 
     $                    - H2(na,nb)*H2(nb,nc) - H4(na,nb,nc,nb) 
        H4(nc,na,nb,nb) = + H1(nc)*H3(na,nb,nb) - H4(na,nb,nb,nc) 
     $                    - H4(na,nb,nc,nb) - H4(na,nc,nb,nb) 
        H4(nc,nb,na,nb) = + H1(nb)*H1(nc)*H2(na,nb) 
     $                    - 2*H1(nc)*H3(na,nb,nb) 
     $                    - H2(na,nb)*H2(nb,nc) + 2*H4(na,nb,nb,nc) 
     $                    + H4(na,nb,nc,nb) 
        H4(nc,nb,nb,na) = + 1.d0/2*H1(na)*H1(nb)*H1(nb)*H1(nc) 
     $                    - H1(na)*H1(nb)*H2(nb,nc) 
     $                    + H1(na)*H3(nb,nb,nc) 
     $                    - H1(nb)*H1(nc)*H2(na,nb) 
     $                    + H1(nc)*H3(na,nb,nb) + H2(na,nb)*H2(nb,nc) 
     $                    - H4(na,nb,nb,nc) 
        if ( iflag.eq.1 ) then 
          call printer4(na,nb,nb,nc) 
          call printer4(na,nb,nc,nb) 
          call printer4(na,nc,nb,nb) 
        endif 
* ic cannot be anymore equal to ib 
      else if ( ic.eq.ib ) then 
        print*,' FILLREDHPL4, error 5, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( ic.eq.id ) then 
* case (na,nb,nc,nc) 
        nb = ib 
        nc = ic 
        H4(nb,na,nc,nc) = + H1(nb)*H3(na,nc,nc) - H4(na,nb,nc,nc) 
     $                    - H4(na,nc,nb,nc) - H4(na,nc,nc,nb) 
        H4(nb,nc,na,nc) = - 2*H1(nb)*H3(na,nc,nc) + H2(na,nc)*H2(nb,nc) 
     $                    + H4(na,nc,nb,nc) + 2*H4(na,nc,nc,nb) 
        H4(nb,nc,nc,na) = + H1(na)*H3(nb,nc,nc) + H1(nb)*H3(na,nc,nc) 
     $                    - H2(na,nc)*H2(nb,nc) - H4(na,nc,nc,nb) 
        H4(nc,na,nb,nc) = + H1(nc)*H3(na,nb,nc) - 2*H4(na,nb,nc,nc) 
     $                    - H4(na,nc,nb,nc) 
        H4(nc,na,nc,nb) = + H1(nc)*H3(na,nc,nb) - H4(na,nc,nb,nc) 
     $                    - 2*H4(na,nc,nc,nb) 
        H4(nc,nb,na,nc) = + H1(nb)*H1(nc)*H2(na,nc) 
     $                    - H1(nc)*H3(na,nb,nc) - H1(nc)*H3(na,nc,nb) 
     $                    - H2(na,nc)*H2(nb,nc) + 2*H4(na,nb,nc,nc) 
     $                    + H4(na,nc,nb,nc) 
        H4(nc,nb,nc,na) = + H1(na)*H1(nc)*H2(nb,nc) 
     $                    - 2*H1(na)*H3(nb,nc,nc) 
     $                    - H1(nb)*H1(nc)*H2(na,nc) 
     $                    + H1(nc)*H3(na,nc,nb) + H2(na,nc)*H2(nb,nc) 
     $                    - H4(na,nc,nb,nc) 
        H4(nc,nc,na,nb) = + 1.d0/2*H1(nc)*H1(nc)*H2(na,nb) 
     $                    - H1(nc)*H3(na,nb,nc) - H1(nc)*H3(na,nc,nb) 
     $                    + H4(na,nb,nc,nc) + H4(na,nc,nb,nc) 
     $                    + H4(na,nc,nc,nb) 
        H4(nc,nc,nb,na) = + 1.d0/2*H1(na)*H1(nb)*H1(nc)*H1(nc) 
     $                    - H1(na)*H1(nc)*H2(nb,nc) 
     $                    + H1(na)*H3(nb,nc,nc) 
     $                    - 1.d0/2*H1(nc)*H1(nc)*H2(na,nb) 
     $                    + H1(nc)*H3(na,nb,nc) - H4(na,nb,nc,nc) 
        if ( iflag.eq.1 ) then 
          call printer4(na,nb,nc,nc) 
          call printer4(na,nc,nb,nc) 
          call printer4(na,nc,nc,nb) 
        endif 
* no need to protect against id.eq.ic 
* when arriving here all indices are different 
      else 
* case (na,nb,nc,nd) all indices are different 
        nb = ib 
        nc = ic 
        nd = id 
        H4(nb,na,nc,nd) = + H1(nb)*H3(na,nc,nd) - H4(na,nb,nc,nd) 
     $                    - H4(na,nc,nb,nd) - H4(na,nc,nd,nb) 
        H4(nb,na,nd,nc) = + H1(nb)*H3(na,nd,nc) - H4(na,nb,nd,nc) 
     $                    - H4(na,nd,nb,nc) - H4(na,nd,nc,nb) 
        H4(nb,nc,na,nd) = - H1(nb)*H3(na,nc,nd) - H1(nb)*H3(na,nd,nc) 
     $                    + H2(na,nd)*H2(nb,nc) + H4(na,nc,nb,nd) 
     $                    + H4(na,nc,nd,nb) + H4(na,nd,nc,nb) 
        H4(nb,nc,nd,na) = + H1(na)*H3(nb,nc,nd) + H1(nb)*H3(na,nd,nc) 
     $                    - H2(na,nd)*H2(nb,nc) - H4(na,nd,nc,nb) 
        H4(nb,nd,na,nc) = - H1(nb)*H3(na,nc,nd) - H1(nb)*H3(na,nd,nc) 
     $                    + H2(na,nc)*H2(nb,nd) + H4(na,nc,nd,nb) 
     $                    + H4(na,nd,nb,nc) + H4(na,nd,nc,nb) 
        H4(nb,nd,nc,na) = + H1(na)*H3(nb,nd,nc) + H1(nb)*H3(na,nc,nd) 
     $                    - H2(na,nc)*H2(nb,nd) - H4(na,nc,nd,nb) 
        H4(nc,na,nb,nd) = + H1(nc)*H3(na,nb,nd) - H4(na,nb,nc,nd) 
     $                    - H4(na,nb,nd,nc) - H4(na,nc,nb,nd) 
        H4(nc,na,nd,nb) = + H1(nc)*H3(na,nd,nb) - H4(na,nc,nd,nb) 
     $                    - H4(na,nd,nb,nc) - H4(na,nd,nc,nb) 
        H4(nc,nb,na,nd) = + H1(nb)*H1(nc)*H2(na,nd) 
     $                    - H1(nc)*H3(na,nb,nd) - H1(nc)*H3(na,nd,nb) 
     $                    - H2(na,nd)*H2(nb,nc) + H4(na,nb,nc,nd) 
     $                    + H4(na,nb,nd,nc) + H4(na,nd,nb,nc) 
        H4(nc,nb,nd,na) = + H1(na)*H1(nc)*H2(nb,nd) 
     $                    - H1(na)*H3(nb,nc,nd) - H1(na)*H3(nb,nd,nc) 
     $                    - H1(nb)*H1(nc)*H2(na,nd) 
     $                    + H1(nc)*H3(na,nd,nb) + H2(na,nd)*H2(nb,nc) 
     $                    - H4(na,nd,nb,nc) 
        H4(nc,nd,na,nb) = - H1(nc)*H3(na,nb,nd) - H1(nc)*H3(na,nd,nb) 
     $                    + H2(na,nb)*H2(nc,nd) + H4(na,nb,nd,nc) 
     $                    + H4(na,nd,nb,nc) + H4(na,nd,nc,nb) 
        H4(nc,nd,nb,na) = + H1(na)*H1(nb)*H2(nc,nd) 
     $                    - H1(na)*H1(nc)*H2(nb,nd) 
     $                    + H1(na)*H3(nb,nd,nc) + H1(nc)*H3(na,nb,nd) 
     $                    - H2(na,nb)*H2(nc,nd) - H4(na,nb,nd,nc) 
        H4(nd,na,nb,nc) = + H1(nd)*H3(na,nb,nc) - H4(na,nb,nc,nd) 
     $                    - H4(na,nb,nd,nc) - H4(na,nd,nb,nc) 
        H4(nd,na,nc,nb) = + H1(nd)*H3(na,nc,nb) - H4(na,nc,nb,nd) 
     $                    - H4(na,nc,nd,nb) - H4(na,nd,nc,nb) 
        H4(nd,nb,na,nc) = + H1(nb)*H1(nd)*H2(na,nc) 
     $                    - H1(nd)*H3(na,nb,nc) - H1(nd)*H3(na,nc,nb) 
     $                    - H2(na,nc)*H2(nb,nd) + H4(na,nb,nc,nd) 
     $                    + H4(na,nb,nd,nc) + H4(na,nc,nb,nd) 
        H4(nd,nb,nc,na) = + H1(na)*H1(nd)*H2(nb,nc) 
     $                    - H1(na)*H3(nb,nc,nd) - H1(na)*H3(nb,nd,nc) 
     $                    - H1(nb)*H1(nd)*H2(na,nc) 
     $                    + H1(nd)*H3(na,nc,nb) + H2(na,nc)*H2(nb,nd) 
     $                    - H4(na,nc,nb,nd) 
        H4(nd,nc,na,nb) = + H1(nc)*H1(nd)*H2(na,nb) 
     $                    - H1(nd)*H3(na,nb,nc) - H1(nd)*H3(na,nc,nb) 
     $                    - H2(na,nb)*H2(nc,nd) + H4(na,nb,nc,nd) 
     $                    + H4(na,nc,nb,nd) + H4(na,nc,nd,nb) 
        H4(nd,nc,nb,na) = + H1(na)*H1(nb)*H1(nc)*H1(nd) 
     $                    - H1(na)*H1(nb)*H2(nc,nd) 
     $                    - H1(na)*H1(nd)*H2(nb,nc) 
     $                    + H1(na)*H3(nb,nc,nd) 
     $                    - H1(nc)*H1(nd)*H2(na,nb) 
     $                    + H1(nd)*H3(na,nb,nc) 
     $                    + H2(na,nb)*H2(nc,nd) - H4(na,nb,nc,nd) 
        if ( iflag.eq.1 ) then 
          call printer4(na,nb,nc,nd) 
          call printer4(na,nb,nd,nc) 
          call printer4(na,nc,nb,nd) 
          call printer4(na,nc,nb,nd) 
          call printer4(na,nd,nb,nc) 
          call printer4(na,nd,nc,nb) 
        endif 
      endif 
*23456789012345678901234567890123456789012345678901234567890123456789012 
      return 
      end 
************************************************************************ 
      subroutine printer2(na,nb) 

      write(11,'(''g [H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('',y)] = H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('',y) ; '')') 

      write(12,'(''id H('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('',y) = H[('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('',y)] ; '')') 

      return 
      end 
*** 
      subroutine printer3(na,nb,nc) 

      write(11,'(''g [H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('','',$)') 
      call subprint(11,nc) 
      write(11,'('',y)] = H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('','',$)') 
      call subprint(11,nc) 
      write(11,'('',y) ; '')') 

      write(12,'(''id H('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('','',$)') 
      call subprint(12,nc) 
      write(12,'('',y) = H[('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('',y)] ; '')') 

      return 
      end 
*** 
      subroutine printer4(na,nb,nc,nd) 

      write(11,'(''g [H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('','',$)') 
      call subprint(11,nc) 
      write(11,'('','',$)') 
      call subprint(11,nd) 
      write(11,'('',y)] = H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('','',$)') 
      call subprint(11,nc) 
      write(11,'('','',$)') 
      call subprint(11,nd) 
      write(11,'('',y) ; '')') 

      write(12,'(''id H('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('','',$)') 
      call subprint(12,nc) 
      write(12,'('','',$)') 
      call subprint(12,nd) 
      write(12,'('',y) = H[('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('','',$)') 
      call subprint(12,nc) 
      write(12,'('','',$)') 
      call subprint(12,nd) 
      write(12,'('',y)] ; '')') 

      return 
      end 
*** 
      subroutine subprint(n,na) 
      if ( na.lt.0 ) then 
        write (n,102) na 
      else 
        write (n,101) na 
      endif 
      return 
  101 format(i1,$) 
  102 format(i2,$) 
      end 

************************************************************************
** the following routines contain th set of routines evaluating 
** irreducible 1dhpl's for various values of the arguments 
************************************************************************ 
      subroutine fillh1(y,H1,HY1,Hi1,n1,n2) 
** fillh1 evaluates the 1dhpl's of weight 1 
      implicit double precision (a-h,o-z) 
      complex*16 H1 
      dimension H1(n1:n2) 
      dimension HY1(n1:n2) 
      dimension Hi1(n1:n2) 
      parameter (pi   = 3.14159265358979324d0) 
      if ( n1.eq.-1) then 
        if ( y.ge.-1.d0 ) then 
          HY1(-1) = log(1.d0+y) 
          Hi1(-1) = 0.d0 
        elseif ( y.lt.-1.d0 ) then 
          HY1(-1) = log(-1.d0-y) 
          Hi1(-1) = 1.d0 
        endif 
        H1(-1) = dcmplx(HY1(-1),pi*Hi1(-1)) 
      endif 
      if ( y.ge.0.d0 ) then 
        HY1(0) = log(y) 
*        Hi1(0) = 0.d0 
      elseif ( y.lt.0.d0 ) then 
        HY1(0) = log(-y) 
        Hi1(0) = 1.d0 
      endif 
      H1(0) = dcmplx(HY1(0),pi*Hi1(0)) 
      if ( n2.eq.1 ) then 
        if ( y.ge.1.d0 ) then 
          HY1(1) = - log(-1.d0+y) 
          Hi1(1) = 1.d0 
        elseif ( y.lt.1.d0 ) then 
          HY1(1) = - log(1.d0-y) 
          Hi1(1) = 0.d0 
        endif 
        H1(1) = dcmplx(HY1(1),pi*Hi1(1)) 
      endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplat0(y,nw,HY1,HY2,HY3,HY4,n1,n2) 
** evaluate the HPL from their power series expansions
** fillirr1dhplat0 is called by eval1dhplat0; 
** it is guaranteed that nw is in the range 1:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
** 
** for y < 0 DOES NOT evaluates the immaginary part of H(0,y) = log(y) 
      implicit double precision (a-h,o-z) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
** evaluating the required 1dHPL of weight 1 
      if ( n1.eq.-1) then 
** 1+y = (1+ep)/(1-ep), ep = y/(2+y) 
** log(1+y) = log((1+y)/(1-y)) = 2*ep*(1+ep^2/3+ep^4/5+.....) 
** at y= -(r2-1) = - 0.4142135624, ep = - 0.26120387496 
** ep2 = 0.068227464296, ep2^13 = 6.9 x 10^(-16) 
         ep = y/(2.d0+y) 
         e2 = ep*ep 
*         v = log(1.d0+y) 
         v = 2*ep*(1+e2*(1.d0/ 3+e2*(1.d0/ 5+e2*(1.d0/ 7+e2*(1.d0/ 9 
     $              +e2*(1.d0/11+e2*(1.d0/13+e2*(1.d0/15+e2*(1.d0/17    
     $              +e2*(1.d0/19+e2*(1.d0/21+e2*(1.d0/23+e2*(1.d0/25    
     $              )))))))))))))
         HY1(-1) = v 
      endif 
      if (y.ge.0d0) then 
         HY1(0) = log(y) 
      else 
         HY1(0) = log(-y) 
** the immaginary part is evaluated in the calling routine eval1dhplat0 
**       Hi1(0) = 1d0 
      endif 
      if ( n2.eq.1) then 
** 1-y = (1-ep)/(1+ep), ep = y/(2-y) 
         ep = y/(2.d0-y) 
         e2 = ep*ep 
*         u = - log(1.d0-y) 
         u = 2*ep*(1+e2*(1.d0/ 3+e2*(1.d0/ 5+e2*(1.d0/ 7+e2*(1.d0/ 9 
     $              +e2*(1.d0/11+e2*(1.d0/13+e2*(1.d0/15+e2*(1.d0/17    
     $              +e2*(1.d0/19+e2*(1.d0/21+e2*(1.d0/23+e2*(1.d0/25    
     $              )))))))))))))
         HY1(1) = u 
      endif 
      if ( nw.eq.1 ) return 
** from now on nw > 1 
** evaluating the Cebyshev polynomials for the expansions 
      ep = y 
      if ( n2.eq.1) then 
        tu01 = 20d0/11d0*u 
        tu02 = 2d0*tu01*tu01 - 1d0 
        tu03 = 2d0*tu01*tu02 - tu01 
        tu04 = 2d0*tu01*tu03 - tu02 
        tu05 = 2d0*tu01*tu04 - tu03 
        tu06 = 2d0*tu01*tu05 - tu04 
        tu07 = 2d0*tu01*tu06 - tu05 
        tu08 = 2d0*tu01*tu07 - tu06 
        tu09 = 2d0*tu01*tu08 - tu07 
        tu10 = 2d0*tu01*tu09 - tu08 
        tu11 = 2d0*tu01*tu10 - tu09 
        tu12 = 2d0*tu01*tu11 - tu10 
        u01 = u
        u02 = u01*u01
        u03 = u01*u02
      endif 
      if ( n1.eq.-1 ) then 
        tv01 = 20d0/11d0*v 
        tv02 = 2d0*tv01*tv01 - 1d0 
        tv03 = 2d0*tv01*tv02 - tv01 
        tv04 = 2d0*tv01*tv03 - tv02 
        tv05 = 2d0*tv01*tv04 - tv03 
        tv06 = 2d0*tv01*tv05 - tv04 
        tv07 = 2d0*tv01*tv06 - tv05 
        tv08 = 2d0*tv01*tv07 - tv06 
        tv09 = 2d0*tv01*tv08 - tv07 
        tv10 = 2d0*tv01*tv09 - tv08 
        tv11 = 2d0*tv01*tv10 - tv09 
        tv12 = 2d0*tv01*tv11 - tv10 
        v01 = v
        v02 = v01*v01
        v03 = v01*v02
      endif 
** evaluating the expansions 
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,1) = u01*(
     $  - 1.3750000000000000d-01*tu01
     $  + 4.1887406497219320d-03*tu02
     $  - 3.1529487739207918d-06*tu04
     $  + 4.0387979720925612d-09*tu06
     $  - 5.9161728033102941d-12*tu08
     $  + 9.2089553676892591d-15*tu10
     $  + 1.0041918976432192d+00)
**    it would be wrong to write 
**    if ( nw.eq.2 ) return 
**    because the (n1.eq.-1).and.(n2.eq.1) case is not yet complete 
      if ( nw.gt.2 ) then 
      HY3(0,0,1) = u01*(
     $  - 2.0733063311847688d-01*tu01
     $  + 1.1909822093713566d-02*tu02
     $  - 3.5978922715974371d-04*tu03
     $  + 1.4651506635246111d-06*tu04
     $  + 2.5264978015885375d-07*tu05
     $  - 2.9129600146853522d-09*tu06
     $  - 3.1200969105278629d-10*tu07
     $  + 5.5608760227994999d-12*tu08
     $  + 4.4683018298608872d-13*tu09
     $  - 1.0377741049678742d-14*tu10
     $  - 6.8492668900504254d-16*tu11
     $  + 1.0119083540245187d+00)
      HY3(0,1,1) = u02*(
     $  - 4.5833333333333333d-02*tu01
     $  + 1.5702519995611332d-03*tu02
     $  - 1.3132234112581665d-06*tu04
     $  + 1.7663819193128645d-09*tu06
     $  - 2.6615094984230225d-12*tu08
     $  + 4.2197085610629510d-15*tu10
     $  + 2.5157156699202004d-01)
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,1) = u01*(
     $  - 2.4309920936894915d-01*tu01
     $  + 1.7710499311295207d-02*tu02
     $  - 8.2489740705617284d-04*tu03
     $  + 2.1971556431143249d-05*tu04
     $  - 9.6292059718056785d-08*tu05
     $  - 1.3396044855444690d-08*tu06
     $  + 1.9835004778127507d-10*tu07
     $  + 1.4796431466871027d-11*tu08
     $  - 3.8471818563928138d-13*tu09
     $  - 1.8899753406004423d-14*tu10
     $  + 7.2332583052852315d-16*tu11
     $  + 2.5542041547393113d-17*tu12
     $  + 1.0176885143440038d+00)
      HY4(0,0,1,1) = u02*(
     $  - 3.8496955551090814d-02*tu01
     $  + 2.7602283087853423d-03*tu02
     $  - 1.0070794072930197d-04*tu03
     $  + 7.6339942564322393d-07*tu04
     $  + 7.7318248534321445d-08*tu05
     $  - 1.4676109918929733d-09*tu06
     $  - 9.8845338864317649d-11*tu07
     $  + 2.7576552760575231d-12*tu08
     $  + 1.4399601068889699d-13*tu09
     $  - 5.1074353206574952d-15*tu10
     $  - 2.2287715896178803d-16*tu11
     $  + 1.2775946343898593d-01)
      HY4(0,1,1,1) = u03*(
     $  - 1.1458333333333333d-02*tu01
     $  + 4.1863378899004153d-04*tu02
     $  - 3.7509447620205836d-07*tu04
     $  + 5.2322893016315673d-10*tu06
     $  - 8.0632112665315483d-13*tu08
     $  + 1.2980885707671547d-15*tu10
     $  + 5.5974564963058350d-02)
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (0,1) or (-1,1) endif 
************ 
** (n1,n2) = (-1,0) or (-1,1) 
      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,-1) = v01*(
     $  + 1.3750000000000000d-01*tv01
     $  + 4.1887406497219320d-03*tv02
     $  - 3.1529487739207918d-06*tv04
     $  + 4.0387979720925612d-09*tv06
     $  - 5.9161728033102941d-12*tv08
     $  + 9.2089553676892591d-15*tv10
     $  + 1.0041918976432192d+00)
      if ( nw.gt.2 ) then 
      HY3(0,0,-1) = v01*(
     $  + 2.0733063311847688d-01*tv01
     $  + 1.1909822093713566d-02*tv02
     $  + 3.5978922715974371d-04*tv03
     $  + 1.4651506635246111d-06*tv04
     $  - 2.5264978015885375d-07*tv05
     $  - 2.9129600146853522d-09*tv06
     $  + 3.1200969105278629d-10*tv07
     $  + 5.5608760227994999d-12*tv08
     $  - 4.4683018298608872d-13*tv09
     $  - 1.0377741049678742d-14*tv10
     $  + 6.8492668900504254d-16*tv11
     $  + 1.0119083540245187d+00)
      HY3(0,-1,-1) = v02*(
     $  + 4.5833333333333333d-02*tv01
     $  + 1.5702519995611332d-03*tv02
     $  - 1.3132234112581665d-06*tv04
     $  + 1.7663819193128645d-09*tv06
     $  - 2.6615094984230225d-12*tv08
     $  + 4.2197085610629510d-15*tv10
     $  + 2.5157156699202004d-01)
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,-1) = v01*(
     $  + 2.4309920936894915d-01*tv01
     $  + 1.7710499311295207d-02*tv02
     $  + 8.2489740705617284d-04*tv03
     $  + 2.1971556431143249d-05*tv04
     $  + 9.6292059718056785d-08*tv05
     $  - 1.3396044855444690d-08*tv06
     $  - 1.9835004778127507d-10*tv07
     $  + 1.4796431466871027d-11*tv08
     $  + 3.8471818563928138d-13*tv09
     $  - 1.8899753406004423d-14*tv10
     $  - 7.2332583052852315d-16*tv11
     $  + 2.5542041547393113d-17*tv12
     $  + 1.0176885143440038d+00)
      HY4(0,0,-1,-1) = v02*(
     $  + 3.8496955551090814d-02*tv01
     $  + 2.7602283087853423d-03*tv02
     $  + 1.0070794072930197d-04*tv03
     $  + 7.6339942564322393d-07*tv04
     $  - 7.7318248534321445d-08*tv05
     $  - 1.4676109918929733d-09*tv06
     $  + 9.8845338864317649d-11*tv07
     $  + 2.7576552760575231d-12*tv08
     $  - 1.4399601068889699d-13*tv09
     $  - 5.1074353206574952d-15*tv10
     $  + 2.2287715896178803d-16*tv11
     $  + 1.2775946343898593d-01)
      HY4(0,-1,-1,-1) = v03*(
     $  + 1.1458333333333333d-02*tv01
     $  + 4.1863378899004153d-04*tv02
     $  - 3.7509447620205836d-07*tv04
     $  + 5.2322893016315673d-10*tv06
     $  - 8.0632112665315483d-13*tv08
     $  + 1.2980885707671547d-15*tv10
     $  + 5.5974564963058350d-02)
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,0) or (-1,1) endif 
** (n1,n2) = (-1,1) -- completion 
      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HY2(-1,1) = u01*(
     $  - 1.0634379058775795d-01*tu01
     $  + 3.9945736734466125d-03*tu02
     $  - 3.7507406061620431d-05*tu03
     $  - 2.6446877835643132d-06*tu04
     $  + 6.3655681886837401d-08*tu05
     $  + 2.8057793781419732d-09*tu06
     $  - 1.1163462030090596d-10*tu07
     $  - 3.1053577733469452d-12*tu08
     $  + 1.9419310124261954d-13*tu09
     $  + 3.0929819630518251d-15*tu10
     $  - 3.3172586119612366d-16*tu11
     $  + 6.9714440173006331d-01)
     $  - 6.9314718055994530d-01*HY1(-1)
      if ( nw.gt.2 ) then 
      HY3(0,-1,1) = u01*(
     $  - 1.4936500293026433d-01*tu01
     $  + 9.1294996590286675d-03*tu02
     $  - 3.1364411079828567d-04*tu03
     $  + 3.3241884419265611d-06*tu04
     $  + 1.7149504451979047d-07*tu05
     $  - 4.9072777567462259d-09*tu06
     $  - 1.6287492110060996d-10*tu07
     $  + 7.7880352395550587d-12*tu08
     $  + 1.7144761374259933d-13*tu09
     $  - 1.2581978726037913d-14*tu10
     $  - 1.7794760437458505d-16*tu11
     $  + 2.0410050886562810d-17*tu12
     $  + 7.0227335111545365d-01)
     $  - 6.9314718055994530d-01*HY2(0,-1)
      HY3(0,1,-1) = v01*(
     $  - 1.4936500293026433d-01*tv01
     $  - 9.1294996590286675d-03*tv02
     $  - 3.1364411079828567d-04*tv03
     $  - 3.3241884419265611d-06*tv04
     $  + 1.7149504451979047d-07*tv05
     $  + 4.9072777567462259d-09*tv06
     $  - 1.6287492110060996d-10*tv07
     $  - 7.7880352395550587d-12*tv08
     $  + 1.7144761374259933d-13*tv09
     $  + 1.2581978726037913d-14*tv10
     $  - 1.7794760437458505d-16*tv11
     $  - 2.0410050886562810d-17*tv12
     $  - 7.0227335111545365d-01)
     $  + 6.9314718055994530d-01*HY2(0,1)
      HY3(-1,-1,1) = u01*(
     $  - 1.3057686804989822d-01*tu01
     $  + 8.4523830620987958d-03*tu02
     $  - 3.1975649807736323d-04*tu03
     $  + 4.6539175216815103d-06*tu04
     $  + 1.5659574963291279d-07*tu05
     $  - 7.2562180394741360d-09*tu06
     $  - 9.6461310543268369d-11*tu07
     $  + 1.1508160166123477d-11*tu08
     $  - 8.9467904347194841d-15*tu09
     $  - 1.7738601968192296d-14*tu10
     $  + 2.3555578892527793d-16*tu11
     $  + 2.6064475394479255d-17*tu12
     $  + 5.9068824834184565d-01)
     $  - 5.8224052646501250d-01*HY1(-1)
     $  - 3.4657359027997265d-01*HY1(-1)*HY1(-1)
      HY3(-1,1,1) = u01*(
     $  + 1.3339975063287861d-01*tu01
     $  - 1.1138781596983798d-02*tu02
     $  + 4.2451528270924965d-04*tu03
     $  - 3.1992340565454975d-06*tu04
     $  - 3.2466452956593668d-07*tu05
     $  + 6.5128138606318713d-09*tu06
     $  + 3.7579027495485995d-10*tu07
     $  - 1.2538868973543979d-11*tu08
     $  - 4.5112950725555316d-13*tu09
     $  + 2.3153072699638146d-14*tu10
     $  + 5.0451421432830974d-16*tu11
     $  - 4.1321853554982899d-17*tu12
     $  - 2.5136208279665204d-01)
     $  + 2.4022650695910071d-01*HY1(-1)
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,-1,1) = u01*(
     $  - 1.7139866732942541d-01*tu01
     $  + 1.2827932191861988d-02*tu02
     $  - 6.2674225522140709d-04*tu03
     $  + 1.8556075778365746d-05*tu04
     $  - 1.8020791893375308d-07*tu05
     $  - 8.5359006198704024d-09*tu06
     $  + 2.4439311365943823d-10*tu07
     $  + 7.2062003807073264d-12*tu08
     $  - 3.7410898512861779d-13*tu09
     $  - 6.6835951661525073d-15*tu10
     $  + 5.9500622408716246d-16*tu11
     $  + 7.0595654813291542d-01)
     $  - 6.9314718055994530d-01*HY3(0,0,-1)
      HY4(0,0,1,-1) = v01*(
     $  - 1.7139866732942541d-01*tv01
     $  - 1.2827932191861988d-02*tv02
     $  - 6.2674225522140709d-04*tv03
     $  - 1.8556075778365746d-05*tv04
     $  - 1.8020791893375308d-07*tv05
     $  + 8.5359006198704024d-09*tv06
     $  + 2.4439311365943823d-10*tv07
     $  - 7.2062003807073264d-12*tv08
     $  - 3.7410898512861779d-13*tv09
     $  + 6.6835951661525073d-15*tv10
     $  + 5.9500622408716246d-16*tv11
     $  - 7.0595654813291542d-01)
     $  + 6.9314718055994530d-01*HY3(0,0,1)
      HY4(0,-1,0,1) = u01*(
     $  - 2.0402091563423462d-01*tu01
     $  + 1.5326796026789534d-02*tu02
     $  - 7.5141752370299855d-04*tu03
     $  + 2.2195862858607303d-05*tu04
     $  - 1.9717216111888252d-07*tu05
     $  - 1.1758559785586910d-08*tu06
     $  + 3.3692089706460695d-10*tu07
     $  + 9.7746481899849247d-12*tu08
     $  - 5.6212130740537968d-13*tu09
     $  - 7.4090460514738114d-15*tu10
     $  + 9.1690209320756988d-16*tu11
     $  + 8.3777162181970230d-01)
     $  - 8.2246703342411321d-01*HY2(0,-1)
      HY4(0,-1,-1,1) = u01*(
     $  - 1.4659170690032865d-01*tu01
     $  + 1.1272254950417420d-02*tu02
     $  - 5.7548355241851144d-04*tu03
     $  + 1.8491926476063684d-05*tu04
     $  - 2.4426966611375586d-07*tu05
     $  - 7.0735008480283508d-09*tu06
     $  + 3.2062539799800555d-10*tu07
     $  + 3.3956059929569370d-12*tu08
     $  - 4.4853380323145144d-13*tu09
     $  + 1.2349995594214108d-15*tu10
     $  + 6.2581127147090835d-16*tu11
     $  + 5.9349428241205865d-01)
     $  - 5.8224052646501250d-01*HY2(0,-1)
     $  - 6.9314718055994530d-01*HY3(0,-1,-1)
      HY4(0,-1,1,-1) = v01*(
     $  - 1.6348235416997523d-01*tv01
     $  - 1.4075851558225201d-02*tv02
     $  - 7.0421561036652880d-04*tv03
     $  - 1.8676050687724687d-05*tv04
     $  - 6.3797094702302232d-08*tv05
     $  + 9.8853199509129313d-09*tv06
     $  + 5.9389227320253100d-11*tv07
     $  - 1.0707718572496392d-11*tv08
     $  - 7.4475662801536059d-14*tv09
     $  + 1.4186776038855675d-14*tv10
     $  + 1.0406065160682246d-16*tv11
     $  - 2.0606613212542347d-17*tv12
     $  - 4.9451017952969702d-01)
     $  + 4.8045301391820142d-01*HY2(0,-1)
     $  + 6.9314718055994530d-01*HY3(0,-1,1)
      HY4(0,1,-1,-1) = v01*(
     $  - 1.0118780392493172d-01*tv01
     $  - 1.0875902950044192d-02*tv02
     $  - 6.9860578351752513d-04*tv03
     $  - 2.5530533667092047d-05*tv04
     $  - 2.8880336859297735d-07*tv05
     $  + 1.5974880084273227d-08*tv06
     $  + 5.0308622711054577d-10*tv07
     $  - 1.5805265017025195d-11*tv08
     $  - 8.6523820415066307d-13*tv09
     $  + 1.6459560669863922d-14*tv10
     $  + 1.4691445597298165d-15*tv11
     $  - 2.5107686338477598d-01)
     $  + 2.4022650695910071d-01*HY2(0,1)
      HY4(0,-1,1,1) = u01*(
     $  + 1.0118780392493172d-01*tu01
     $  - 1.0875902950044192d-02*tu02
     $  + 6.9860578351752513d-04*tu03
     $  - 2.5530533667092047d-05*tu04
     $  + 2.8880336859297735d-07*tu05
     $  + 1.5974880084273227d-08*tu06
     $  - 5.0308622711054577d-10*tu07
     $  - 1.5805265017025195d-11*tu08
     $  + 8.6523820415066307d-13*tu09
     $  + 1.6459560669863922d-14*tu10
     $  - 1.4691445597298165d-15*tu11
     $  - 2.5107686338477598d-01)
     $  + 2.4022650695910071d-01*HY2(0,-1)
      HY4(0,1,-1,1) = u01*(
     $  + 1.6348235416997523d-01*tu01
     $  - 1.4075851558225201d-02*tu02
     $  + 7.0421561036652880d-04*tu03
     $  - 1.8676050687724687d-05*tu04
     $  + 6.3797094702302232d-08*tu05
     $  + 9.8853199509129313d-09*tu06
     $  - 5.9389227320253100d-11*tu07
     $  - 1.0707718572496392d-11*tu08
     $  + 7.4475662801536059d-14*tu09
     $  + 1.4186776038855675d-14*tu10
     $  - 1.0406065160682246d-16*tu11
     $  - 2.0606613212542347d-17*tu12
     $  - 4.9451017952969702d-01)
     $  + 4.8045301391820142d-01*HY2(0,1)
     $  - 6.9314718055994530d-01*HY3(0,1,-1)
      HY4(0,1,1,-1) = v01*(
     $  + 1.4659170690032865d-01*tv01
     $  + 1.1272254950417420d-02*tv02
     $  + 5.7548355241851144d-04*tv03
     $  + 1.8491926476063684d-05*tv04
     $  + 2.4426966611375586d-07*tv05
     $  - 7.0735008480283508d-09*tv06
     $  - 3.2062539799800555d-10*tv07
     $  + 3.3956059929569370d-12*tv08
     $  + 4.4853380323145144d-13*tv09
     $  + 1.2349995594214108d-15*tv10
     $  - 6.2581127147090835d-16*tv11
     $  + 5.9349428241205865d-01)
     $  - 5.8224052646501250d-01*HY2(0,1)
     $  + 6.9314718055994530d-01*HY3(0,1,1)
      HY4(-1,-1,-1,1) = u01*(
     $  - 1.3704186675693694d-01*tu01
     $  + 1.0738664799045256d-02*tu02
     $  - 5.6405487546831413d-04*tu03
     $  + 1.8975289545570792d-05*tu04
     $  - 2.8136699841840261d-07*tu05
     $  - 6.9483354697236112d-09*tu06
     $  + 3.9086778571652884d-10*tu07
     $  + 1.5013470800132572d-12*tu08
     $  - 5.4857187012000446d-13*tu09
     $  + 7.0089201006248443d-15*tu10
     $  + 7.1985359614000227d-16*tu11
     $  - 2.2459200848103789d-17*tu12
     $  + 5.4793287616771010d-01)
     $  - 5.3721319360804020d-01*HY1(-1)
     $  - 2.9112026323250625d-01*HY1(-1)*HY1(-1)
     $  - 1.1552453009332421d-01*HY1(-1)*HY1(-1)*HY1(-1)
      HY4(-1,-1,1,1) = u01*(
     $  + 1.0577793416858150d-01*tu01
     $  - 1.0477240350312842d-02*tu02
     $  + 6.6271618168294058d-04*tu03
     $  - 2.5384291818094501d-05*tu04
     $  + 3.7532043652092099d-07*tu05
     $  + 1.3804832669451887d-08*tu06
     $  - 6.6481886939724744d-10*tu07
     $  - 8.3615840903091915d-12*tu08
     $  + 1.1154536519042777d-12*tu09
     $  - 2.8020941621207856d-15*tu10
     $  - 1.7701739081331563d-15*tu11
     $  + 2.7769076888919123d-17*tu12
     $  - 3.1927721734213725d-01)
     $  + 3.0882537509683393d-01*HY1(-1)
     $  + 1.2011325347955035d-01*HY1(-1)*HY1(-1)
      HY4(-1,1,1,1) = u01*(
     $  - 3.2833473259986549d-02*tu01
     $  + 8.5187654858627660d-03*tu02
     $  - 7.6913672477984726d-04*tu03
     $  + 3.0885336618882649d-05*tu04
     $  - 2.3861711695799581d-07*tu05
     $  - 2.5246408915769234d-08*tu06
     $  + 5.1217849936005379d-10*tu07
     $  + 3.0326173877457694d-11*tu08
     $  - 1.0166074598017197d-12*tu09
     $  - 3.7296649908045336d-14*tu10
     $  + 1.9148798730971194d-15*tu11
     $  + 4.2483769310186512d-17*tu12
     $  + 6.3991963537293034d-02)
     $  - 5.5504108664821579d-02*HY1(-1)
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,1) -- completion endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplat1(r,nw,HR1,HR2,HR3,HR4, 
     $                                HY1,HY2,HY3,HY4, 
     $                                Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates the HPL for r2m1 < y < r2p1
** fillirr1dhplat1 is called by eval1dhplat1 after calling 
** fillirr1dhplat0 with argument r=(1-y)/(1+y) 
** it is guaranteed that nw is in the range 2:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
      implicit double precision (a-h,o-z) 
      dimension HR1(-1:1),HR2(-1:1,-1:1),HR3(-1:1,-1:1,-1:1), 
     $          HR4(-1:1,-1:1,-1:1,-1:1) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,1) = 
     $  + 1.6449340668482264d+00
     $  + 6.9314718055994530d-01*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)
     $  + HR1( -1)*HR1(0)
     $  - HR1( -1)*HR1(1)
     $  + HR1(0) *HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)
     $  + HR2( -1,1)
     $  - HR2(0, -1)
     $  - HR2(0,1) 
      if (r.lt.0d0) then 
      Hi2(0,1) = 
     $  - HR1( -1) 
     $  - HR1(1) 
      endif 
      if ( nw.gt.2 ) then 
      HY3(0,0,1) = 
     $  + 1.2020569031595942d+00
     $  - 1.6449340668482264d+00*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + HR1( -1)*HR2(0,-1)
     $  + HR1( -1)*HR2(0,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(1)*HR1(1)
     $  - 1.6449340668482264d+00*HR1(1)
     $  - 3.4657359027997265d-01*HR1(1)*HR1(1)
     $  - HR1(1) *HR2(-1,1)
     $  + HR1(1) *HR2(0,-1)
     $  + HR1(1) *HR2(0,1)
     $  - HR3( -1,-1,1)
     $  + HR3( -1,1,1)
     $  - HR3(0, -1,-1)
     $  - HR3(0, -1,1)
     $  - HR3(0,1, -1)
     $  - HR3(0,1,1) 
      HY3(0,1,1) = 
     $  + 1.2020569031595942d+00
     $  - 2.4022650695910071d-01*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + HR1( -1)*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(0)*HR1(1)
     $  - HR1(0) *HR2(-1,1)
     $  + HR1(0) *HR2(0,-1)
     $  + HR1(0) *HR2(0,1)
     $  - 2.4022650695910071d-01*HR1(1)
     $  - 6.9314718055994530d-01*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(0,1)
     $  - HR3( -1,-1,1)
     $  - HR3(0, -1,-1)
     $  - HR3(0,0, -1)
     $  - HR3(0,0,1) 
     $  - HR3(0,1, -1)
      if (r.lt.0d0) then 
      HY3(0,1,1) = HY3(0,1,1) 
     $  + 4.9348022005446793d+00*HR1(-1)
     $  + 4.9348022005446793d+00*HR1(1)
      Hi3(0,0,1) = 
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1) 
     $  + HR1( -1)*HR1(1) 
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1) 
      Hi3(0,1,1) = 
     $  + 6.9314718055994530d-01*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)
     $  + HR1( -1)*HR1(0)
     $  - HR1( -1)*HR1(1)
     $  + HR1(0) *HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)
     $  + HR2( -1,1)
     $  - HR2(0, -1)
     $  - HR2(0,1) 
      endif
      endif
      if ( nw.gt.3 ) then 
      HY4(0,0,0,1) = 
     $  + 1.0823232337111381d+00
     $  - 1.2020569031595942d+00*HR1(-1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  + 1.6449340668482264d+00*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(1)*HR2(0,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,-1,1)
     $  + HR1( -1)*HR3(0,1,-1)
     $  + HR1( -1)*HR3(0,1,1)
     $  + 1.6666666666666666d-01*HR1(0)*HR1(1)*HR1(1)*HR1(1)
     $  - 1.2020569031595942d+00*HR1(1)
     $  + 8.2246703342411321d-01*HR1(1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(0,1)
     $  + HR1(1) *HR3(-1,-1,1)
     $  - HR1(1) *HR3(-1,1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,-1,1)
     $  + HR1(1) *HR3(0,1,-1)
     $  + HR1(1) *HR3(0,1,1)
     $  + HR4( -1,-1,-1,1)
     $  - HR4( -1,-1,1,1)
     $  + HR4( -1,1,1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0, -1,1,1)
     $  - HR4(0,1, -1,-1)
     $  - HR4(0,1, -1,1)
     $  - HR4(0,1,1, -1)
     $  - HR4(0,1,1,1) 
      HY4(0,0,1,1) = 
     $  + 2.7058080842778454d-01
     $  - 1.2020569031595942d+00*HR1(-1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR2(0,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(-1,1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  + HR1( -1)*HR3(0,0,1)
     $  + HR1( -1)*HR3(0,1,-1)
     $  + 2.5000000000000000d-01*HR1(0)*HR1(0)*HR1(1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(0)*HR1(1)*HR1(1)
     $  + HR1(0) *HR1(1)*HR2(-1,1)
     $  - HR1(0) *HR1(1)*HR2(0,-1)
     $  - HR1(0) *HR1(1)*HR2(0,1)
     $  + HR1(0) *HR3(-1,-1,1)
     $  - HR1(0) *HR3(-1,1,1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,-1,1)
     $  + HR1(0) *HR3(0,1,-1)
     $  + HR1(0) *HR3(0,1,1)
     $  - 1.2020569031595942d+00*HR1(1)
     $  + 1.2011325347955035d-01*HR1(1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(0,1)
     $  + HR1(1) *HR3(-1,-1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,0,-1)
     $  + HR1(1) *HR3(0,0,1)
     $  + HR1(1) *HR3(0,1,-1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,1,1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR3(0,1,-1)
     $  + 6.9314718055994530d-01*HR3(0,1,1)
     $  + 2.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  - HR4( -1,-1,1,1)
     $  - 2.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0, -1,1)
     $  - HR4(0,0,1, -1)
     $  - HR4(0,0,1,1) 
     $  - 2.0000000000000000d+00*HR4(0,1,-1,-1)
     $  - HR4(0,1, -1,1)
     $  - HR4(0,1,1, -1)
      HY4(0,1,1,1) = 
     $  + 1.0823232337111381d+00
     $  + 5.5504108664821579d-02*HR1(-1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + 1.6666666666666666d-01*HR1(0)*HR1(0)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(0)*HR1(0)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(0,1)
     $  + 2.4022650695910071d-01*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(0)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(0,1)
     $  + HR1(0) *HR3(-1,-1,1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,0,-1)
     $  + HR1(0) *HR3(0,0,1)
     $  + HR1(0) *HR3(0,1,-1)
     $  + 5.5504108664821579d-02*HR1(1)
     $  + 2.4022650695910071d-01*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR2(0,-1)
     $  - 2.4022650695910071d-01*HR2(0,1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,0,-1)
     $  + 6.9314718055994530d-01*HR3(0,0,1)
     $  + 6.9314718055994530d-01*HR3(0,1,-1)
     $  + HR4( -1,-1,-1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0,0, -1)
     $  - HR4(0,0,0,1) 
     $  - HR4(0,0,1, -1)
     $  - HR4(0,1, -1,-1)
      if (r.lt.0d0) then 
      HY4(0,0,1,1) = HY4(0,0,1,1) 
     $  - 2.4674011002723396d+00*HR1(-1)*HR1(-1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR1(1)
     $  - 2.4674011002723396d+00*HR1(1)*HR1(1)
      HY4(0,1,1,1) = HY4(0,1,1,1) 
     $  - 3.4205442319285582d+00*HR1(-1)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(-1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR1(0)
     $  + 4.9348022005446793d+00*HR1(-1)*HR1(1)
     $  - 4.9348022005446793d+00*HR1(0)*HR1(1)
     $  - 3.4205442319285582d+00*HR1(1)
     $  - 4.9348022005446793d+00*HR2(-1,1)
     $  + 4.9348022005446793d+00*HR2(0,-1)
     $  + 4.9348022005446793d+00*HR2(0,1)
      Hi4(0,0,0,1) = 
     $  - 1.666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  - 5.000000000000000d-01*HR1(-1)*HR1(1)*HR1(1) 
     $  - 1.666666666666666d-01*HR1(1)*HR1(1)*HR1(1) 
      Hi4(0,0,1,1) = 
     $  - 3.465735902799726d-01*HR1(-1)*HR1(-1) 
     $  + 1.666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0) 
     $  + 5.000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  - HR1( -1)*HR1(0)*HR1(1) 
     $  - 6.931471805599453d-01*HR1(-1)*HR1(1) 
     $  + 5.000000000000000d-01*HR1(-1)*HR1(1)*HR1(1) 
     $  + HR1( -1)*HR2(0,-1) 
     $  + HR1( -1)*HR2(0,1) 
     $  - 5.000000000000000d-01*HR1(0)*HR1(1)*HR1(1) 
     $  - 3.465735902799726d-01*HR1(1)*HR1(1) 
     $  - HR1(1) *HR2(-1,1) 
     $  + HR1(1) *HR2(0,-1) 
     $  + HR1(1) *HR2(0,1) 
     $  - HR3( -1,-1,1) 
     $  + HR3( -1,1,1) 
     $  - HR3(0, -1,-1) 
     $  - HR3(0, -1,1) 
     $  - HR3(0,1, -1) 
     $  - HR3(0,1,1) 
      Hi4(0,1,1,1) = 
     $  + 1.404707559889125d+00*HR1(-1) 
     $  + 3.465735902799726d-01*HR1(-1)*HR1(-1) 
     $  - 1.666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  + 5.000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0) 
     $  - 5.000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  - 6.931471805599453d-01*HR1(-1)*HR1(0) 
     $  - 5.000000000000000d-01*HR1(-1)*HR1(0)*HR1(0) 
     $  + HR1( -1)*HR1(0)*HR1(1) 
     $  + 6.931471805599453d-01*HR1(-1)*HR1(1) 
     $  + HR1( -1)*HR2(-1,1) 
     $  - 5.000000000000000d-01*HR1(0)*HR1(0)*HR1(1) 
     $  - 6.931471805599453d-01*HR1(0)*HR1(1) 
     $  - HR1(0) *HR2(-1,1) 
     $  + HR1(0) *HR2(0,-1) 
     $  + HR1(0) *HR2(0,1) 
     $  + 1.404707559889125d+00*HR1(1) 
     $  - 6.931471805599453d-01*HR2(-1,1) 
     $  + 6.931471805599453d-01*HR2(0,-1) 
     $  + 6.931471805599453d-01*HR2(0,1) 
     $  - HR3( -1,-1,1) 
     $  - HR3(0, -1,-1) 
     $  - HR3(0,0, -1) 
     $  - HR3(0,0,1)  
     $  - HR3(0,1, -1) 
      endif 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (0,1) or (-1,1) endif 
************ 
** (n1,n2) = (-1,0) or (-1,1) 
      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
       HY2(0,-1) = 
     $  + 8.2246703342411321d-01
     $  - 6.9314718055994530d-01*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)
     $  + HR1( -1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(1)
     $  - HR2( -1,1)
      if ( nw.gt.2 ) then 
      HY3(0,0,-1) = 
     $  + 9.0154267736969571d-01
     $  - 8.2246703342411321d-01*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 8.2246703342411321d-01*HR1(1)
     $  + 3.4657359027997265d-01*HR1(1)*HR1(1)
     $  + HR1(1) *HR2(-1,1)
     $  + HR3( -1,-1,1)
     $  - HR3( -1,1,1)
      HY3(0,-1,-1) = 
     $  + 1.5025711289494928d-01
     $  - 2.4022650695910071d-01*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR1(1)
     $  - 6.9314718055994530d-01*HR2(-1,1)
     $  - HR3( -1,-1,1)
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,-1) = 
     $  + 9.4703282949724591d-01
     $  - 9.0154267736969571d-01*HR1(-1)
     $  + 4.1123351671205660d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  - 9.0154267736969571d-01*HR1(1)
     $  + 4.1123351671205660d-01*HR1(1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  - HR1(1) *HR3(-1,-1,1)
     $  + HR1(1) *HR3(-1,1,1)
     $  - HR4( -1,-1,-1,1)
     $  + HR4( -1,-1,1,1)
     $  - HR4( -1,1,1,1)
      HY4(0,0,-1,-1) = 
     $  + 8.7785671568655302d-02
     $  - 1.5025711289494928d-01*HR1(-1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(-1,1,1)
     $  - 1.5025711289494928d-01*HR1(1)
     $  + 1.2011325347955035d-01*HR1(1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  + HR1(1) *HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,1,1)
     $  + 2.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  - HR4( -1,-1,1,1)
      HY4(0,-1,-1,-1) = 
     $  + 2.3752366322618485d-02
     $  - 5.5504108664821579d-02*HR1(-1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - 5.5504108664821579d-02*HR1(1)
     $  - 2.4022650695910071d-01*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - HR4( -1,-1,-1,1)
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,0) or (-1,1) endif 
** (n1,n2) = (-1,1) -- completion 
      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HY2(-1,1) = 
     $  + 5.8224052646501250d-01
     $  + 6.9314718055994530d-01*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)
     $  + HR1( -1)*HR1(0)
     $  - HR2(0, -1)
      if (r.lt.0d0) then 
      Hi2(-1,1) = 
     $  - HR1( -1) 
      endif 
      if ( nw.gt.2 ) then 
      HY3(0,-1,1) = 
     $  + 2.4307035167006157d-01
     $  - 5.8224052646501250d-01*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR2(-1,1)
     $  + HR1( -1)*HR2(0,-1)
     $  + HR1(0) *HR2(-1,1)
     $  - 5.8224052646501250d-01*HR1(1)
     $  + HR1(1) *HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(-1,1)
     $  + HR3( -1,-1,1)
     $  - HR3(0, -1,-1)
     $  - HR3(0, -1,1)
      HY3(0,1,-1) = 
     $  + 5.0821521280468485d-01
     $  + 1.0626935403832139d+00*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR2(-1,1)
     $  - HR1( -1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR1(1)
     $  + 1.0626935403832139d+00*HR1(1)
     $  - HR1(1) *HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR2(0,1)
     $  + HR3( -1,-1,1)
     $  + 2.0000000000000000d+00*HR3(0,-1,-1)
     $  + HR3(0, -1,1)
     $  + HR3(0,1, -1)
      HY3(-1,-1,1) = 
     $  + 9.4753004230127705d-02
     $  - 5.8224052646501250d-01*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + HR1( -1)*HR2(0,-1)
     $  - HR3(0, -1,-1)
      HY3(-1,1,1) = 
     $  + 5.3721319360804020d-01
     $  - 2.4022650695910071d-01*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + HR1(0) *HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(0,-1)
     $  - HR3(0, -1,-1)
     $  - HR3(0,0, -1)
      if (r.lt.0d0) then 
      HY3(-1,1,1) = HY3(-1,1,1) 
     $  + 4.9348022005446793d+00*HR1(-1)
      Hi3(0,-1,1) = 
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1) 
     $  + HR1( -1)*HR1(1) 
     $  - HR2( -1,1) 
      Hi3(0,1,-1) = 
     $  - 6.9314718055994530d-01*HR1(-1) 
     $  - 6.9314718055994530d-01*HR1(1) 
      Hi3(-1,-1,1) = 
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1) 
      Hi3(-1,1,1) = 
     $  + 6.9314718055994530d-01*HR1(-1) 
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1) 
     $  + HR1( -1)*HR1(0) 
     $  - HR2(0, -1) 
      endif 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,-1,1) = 
     $  + 1.1787599965050932d-01
     $  - 2.4307035167006157d-01*HR1(-1)
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR3(-1,1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,-1,1)
     $  - HR1(0) *HR1(1)*HR2(-1,1)
     $  - HR1(0) *HR3(-1,-1,1)
     $  + HR1(0) *HR3(-1,1,1)
     $  - 2.4307035167006157d-01*HR1(1)
     $  + 2.9112026323250625d-01*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  - HR1(1) *HR3(-1,-1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(-1,1,1)
     $  - 2.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  + HR4( -1,-1,1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0, -1,1,1)
      HY4(0,0,1,-1) = 
     $  + 1.7284527823898438d-01
     $  - 5.0821521280468485d-01*HR1(-1)
     $  - 5.3134677019160696d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(-1,1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR3(-1,1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR3(0,-1,1)
     $  - HR1( -1)*HR3(0,1,-1)
     $  - 3.4657359027997265d-01*HR1(0)*HR1(1)*HR1(1)
     $  - 5.0821521280468485d-01*HR1(1)
     $  - 5.3134677019160696d-01*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,1)
     $  - HR1(1) *HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,-1,-1)
     $  - HR1(1) *HR3(0,-1,1)
     $  - HR1(1) *HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(-1,1,1)
     $  - 6.9314718055994530d-01*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR3(0,1,1)
     $  - 2.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  + HR4( -1,-1,1,1)
     $  + 3.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR4(0,-1,-1,1)
     $  + 2.0000000000000000d+00*HR4(0,-1,1,-1)
     $  + HR4(0, -1,1,1)
     $  + 2.0000000000000000d+00*HR4(0,1,-1,-1)
     $  + HR4(0,1, -1,1)
     $  + HR4(0,1,1, -1)
      HY4(0,-1,0,1) = 
     $  + 2.0293560632083841d-01
     $  - 3.8889584616810632d-01*HR1(-1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  - HR1( -1)*HR1(0)*HR2(-1,1)
     $  + 1.6449340668482264d+00*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(1)*HR2(0,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR3(-1,1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,1,-1)
     $  + HR1(0) *HR1(1)*HR2(-1,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR3(-1,1,1)
     $  - 3.8889584616810632d-01*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR3(-1,-1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,1,-1)
     $  - 1.6449340668482264d+00*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR2(-1,1)*HR2(-1,1)
     $  + HR2( -1,1)*HR2(0,-1)
     $  + HR2( -1,1)*HR2(0,1)
     $  + 1.3862943611198906d+00*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR3(-1,1,1)
     $  + 4.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR4(-1,-1,1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0,1, -1,-1)
     $  - HR4(0,1, -1,1)
      HY4(0,-1,-1,1) = 
     $  + 3.4159126166513913d-02
     $  - 9.4753004230127705d-02*HR1(-1)
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR2(-1,1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1(0) *HR3(-1,-1,1)
     $  - 9.4753004230127705d-02*HR1(1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  - 5.8224052646501250d-01*HR2(-1,1)
     $  + HR2( -1,1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + HR4( -1,-1,-1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
      HY4(0,-1,1,-1) = 
     $  + 5.4653052738263652d-02
     $  - 2.1407237086670622d-01*HR1(-1)
     $  - 5.3134677019160696d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR2(-1,1)
     $  - 2.1407237086670622d-01*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,-1,-1)
     $  + 1.0626935403832139d+00*HR2(-1,1)
     $  - HR2( -1,1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR3(0,-1,1)
     $  + HR4( -1,-1,-1,1)
     $  + 3.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR4(0,-1,-1,1)
     $  + HR4(0, -1,1,-1)
      HY4(0,1,-1,-1) = 
     $  + 1.1412342741606084d-01
     $  + 4.7533770109129867d-01*HR1(-1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(0)
     $  - 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + 2.4022650695910071d-01*HR1(0)*HR1(1)
     $  + 4.7533770109129867d-01*HR1(1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + 2.4022650695910071d-01*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR2(0,-1)
     $  - 2.4022650695910071d-01*HR2(0,1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR3(0,1,-1)
     $  + HR4( -1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0,1, -1,-1)
      HY4(0,-1,1,1) = 
     $  + 9.3097125991768577d-02
     $  - 5.3721319360804020d-01*HR1(-1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  + HR1( -1)*HR1(0)*HR2(-1,1)
     $  - HR1( -1)*HR1(0)*HR2(0,-1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(-1,1)
     $  - HR1(0) *HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(-1,1)
     $  - HR1(0) *HR3(-1,-1,1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,-1,1)
     $  - 5.3721319360804020d-01*HR1(1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,0,-1)
     $  - 2.4022650695910071d-01*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,1)
     $  - HR4( -1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0, -1,1)
      HY4(0,1,-1,1) = 
     $  + 1.9355535381306524d-01
     $  + 1.4780047665430420d+00*HR1(-1)
     $  - 2.9112026323250625d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR1(0)
     $  + HR1( -1)*HR1(0)*HR2(-1,1)
     $  + HR1( -1)*HR1(0)*HR2(0,-1)
     $  - 5.8224052646501250d-01*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,0,-1)
     $  + 5.8224052646501250d-01*HR1(0)*HR1(1)
     $  + HR1(0) *HR1(1)*HR2(0,-1)
     $  - HR1(0) *HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR3(0,-1,-1)
     $  - HR1(0) *HR3(0,-1,1)
     $  - HR1(0) *HR3(0,1,-1)
     $  + 1.4780047665430420d+00*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,0,-1)
     $  + 5.8224052646501250d-01*HR2(-1,1)
     $  - HR2( -1,1)*HR2(0,-1)
     $  - 5.8224052646501250d-01*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR2(0,-1)*HR2(0,-1)
     $  + HR2(0, -1)*HR2(0,1)
     $  - 5.8224052646501250d-01*HR2(0,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR3(0,1,-1)
     $  - HR4( -1,-1,-1,1)
     $  + 4.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR4(0,-1,-1,1)
     $  - HR4(0, -1,0,1)
     $  + HR4(0, -1,1,-1)
     $  + 2.0000000000000000d+00*HR4(0,0,-1,-1)
     $  + HR4(0,1, -1,-1)
      HY4(0,1,1,-1) = 
     $  + 4.3369237704895519d-01
     $  - 1.1073038989294665d+00*HR1(-1)
     $  + 5.3134677019160696d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR1(0)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  + 1.0626935403832139d+00*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  - 3.4657359027997265d-01*HR1(0)*HR1(0)*HR1(1)
     $  - 1.0626935403832139d+00*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR2(0,1)
     $  - 1.1073038989294665d+00*HR1(1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,0,-1)
     $  - 1.0626935403832139d+00*HR2(-1,1)
     $  + HR2( -1,1)*HR2(0,-1)
     $  + 1.0626935403832139d+00*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR2(0,-1)*HR2(0,-1)
     $  - HR2(0, -1)*HR2(0,1)
     $  + 1.0626935403832139d+00*HR2(0,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR3(0,0,-1)
     $  - 6.9314718055994530d-01*HR3(0,0,1)
     $  - 6.9314718055994530d-01*HR3(0,1,-1)
     $  - HR4( -1,-1,-1,1)
     $  - HR4(0, -1,-1,1)
     $  + HR4(0, -1,0,1)
     $  + HR4(0,0, -1,1)
     $  + HR4(0,0,1, -1)
     $  + HR4(0,1, -1,-1)
      HY4(-1,-1,-1,1) = 
     $  + 1.4134237214990008d-02
     $  - 9.4753004230127705d-02*HR1(-1)
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  - HR4(0, -1,-1,-1)
      HY4(-1,-1,1,1) = 
     $  + 4.0758239159309251d-02
     $  - 5.3721319360804020d-01*HR1(-1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - HR1( -1)*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  - HR4(0,0, -1,-1)
      HY4(-1,1,1,1) = 
     $  + 5.1747906167389938d-01
     $  + 5.5504108664821579d-02*HR1(-1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(0,-1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,0,-1)
     $  - 2.4022650695910071d-01*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,0,-1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0,0, -1)
      if (r.lt.0d0) then 
      HY4(0,-1,1,1) = HY4(0,-1,1,1) 
     $  - 2.4674011002723396d+00*HR1(-1)*HR1(-1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR1(1)
     $  + 4.9348022005446793d+00*HR2(-1,1)
      HY4(0,1,1,-1) = HY4(0,1,1,-1) 
     $  + 3.4205442319285582d+00*HR1(-1)
     $  + 3.4205442319285582d+00*HR1(1)
      HY4(-1,-1,1,1) = HY4(-1,-1,1,1) 
     $  - 2.4674011002723396d+00*HR1(-1)*HR1(-1) 
      HY4(-1,1,1,1) = HY4(-1,1,1,1)
     $  - 3.4205442319285582d+00*HR1(-1)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(-1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR1(0)
     $  + 4.9348022005446793d+00*HR2(0,-1)
      Hi4(0,0,-1,1) = 
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1) 
     $  + HR1(1) *HR2(-1,1) 
     $  + HR3( -1,-1,1) 
     $  - HR3( -1,1,1) 
      Hi4(0,0,1,-1) = 
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1) 
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1) 
     $  + 3.4657359027997265d-01*HR1(1)*HR1(1) 
      Hi4(0,-1,0,1) = 
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  + HR1( -1)*HR2(-1,1) 
     $  - HR1(1) *HR2(-1,1) 
     $  - 2.0000000000000000d+00*HR3(-1,-1,1) 
     $  + 2.0000000000000000d+00*HR3(-1,1,1) 
      Hi4(0,-1,-1,1) = 
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  + HR1( -1)*HR2(-1,1) 
     $  - HR3( -1,-1,1) 
      Hi4(0,-1,1,-1) = 
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1) 
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1) 
     $  - 6.9314718055994530d-01*HR2(-1,1) 
      Hi4(0,1,-1,-1) =  
     $  - 2.4022650695910071d-01*HR1(-1) 
     $  - 2.4022650695910071d-01*HR1(1) 
      Hi4(0,-1,1,1) = 
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR2(-1,1)
     $  + HR1( -1)*HR2(0,-1)
     $  + HR1(0) *HR2(-1,1)
     $  + HR1(1) *HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(-1,1)
     $  + HR3( -1,-1,1)
     $  - HR3(0, -1,-1)
     $  - HR3(0, -1,1)
      Hi4(0,1,-1,1) = 
     $  - 5.8224052646501250d-01*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR2(-1,1)
     $  - HR1( -1)*HR2(0,-1)
     $  - 5.8224052646501250d-01*HR1(1)
     $  - HR1(1) *HR2(0,-1)
     $  + HR3( -1,-1,1)
     $  + 2.0000000000000000d+00*HR3(0,-1,-1)
     $  + HR3(0, -1,1)
     $  + HR3(0,1, -1)
      Hi4(0,1,1,-1) = 
     $  + 1.0626935403832139d+00*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(0)*HR1(1)
     $  + 1.0626935403832139d+00*HR1(1)
     $  + 6.9314718055994530d-01*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR2(0,1)
      Hi4(-1,-1,-1,1) = 
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
      Hi4(-1,-1,1,1) = 
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + HR1( -1)*HR2(0,-1)
     $  - HR3(0, -1,-1)
      Hi4(-1,1,1,1) = 
     $  + 1.4047075598891257d+00*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + HR1(0) *HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(0,-1)
     $  - HR3(0, -1,-1)
     $  - HR3(0,0, -1)
      endif  
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,1) -- completion endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplatinf(x,nw,HX1,HX2,HX3,HX4, 
     $                                HY1,HY2,HY3,HY4, 
     $                                Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates the HPL for y > r2p1
** fillirr1dhplatinf is called by eval1dhplatinf after calling 
** fillirr1dhplat0 with argument r=1/y 
** it is guaranteed that nw is in the range 2:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
      implicit double precision (a-h,o-z) 
      dimension HX1(n1:n2),HX2(n1:n2,n1:n2),HX3(n1:n2,n1:n2,n1:n2), 
     $          HX4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,1) = 
     $  + 3.2898681336964528d+00 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  - HX2(0,1) 
      Hi2(0,1) = 
     $  - HX1(0) 
      if ( nw.gt.2 ) then 
      HY3(0,0,1) = 
     $  - 3.2898681336964528d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3(0,0,1) 
      HY3(0,1,1) = 
     $  + 1.2020569031595942d+00 
     $  + 4.9348022005446793d+00*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX2(0,1) 
     $  + HX3(0,0,1) 
     $  - HX3(0,1,1) 
      Hi3(0,0,1) = 
     $  + 5.000000000000000d-01*HX1(0)*HX1(0) 
      Hi3(0,1,1) = 
     $  + 1.6449340668482264d+00 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  - HX2(0,1) 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,1) = 
     $  + 2.1646464674222763d+00 
     $  + 1.6449340668482264d+00*HX1(0)*HX1(0) 
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - HX4(0,0,0,1) 
      HY4(0,0,1,1) = 
     $  + 2.1646464674222763d+00 
     $  - 1.2020569031595942d+00*HX1(0) 
     $  - 2.4674011002723396d+00*HX1(0)*HX1(0) 
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX3(0,0,1) 
     $  - 2.0000000000000000d+00*HX4(0,0,0,1) 
     $  + HX4(0,0,1,1) 
      HY4(0,1,1,1) = 
     $  - 5.1410353601279064d+00 
     $  + 2.4674011002723396d+00*HX1(0)*HX1(0) 
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1) 
     $  + HX1(0) *HX3(0,0,1) 
     $  - HX1(0) *HX3(0,1,1) 
     $  + 4.9348022005446793d+00*HX2(0,1) 
     $  - HX4(0,0,0,1) 
     $  + HX4(0,0,1,1) 
     $  - HX4(0,1,1,1) 
      Hi4(0,0,0,1) = 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
      Hi4(0,0,1,1) = 
     $  - 1.2020569031595942d+00 
     $  - 1.6449340668482264d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3(0,0,1) 
      Hi4(0,1,1,1) = 
     $  + 1.6449340668482264d+00*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX2(0,1) 
     $  + HX3(0,0,1) 
     $  - HX3(0,1,1) 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (0,1) or (-1,1) endif 
************ 
** (n1,n2) = (-1,0) or (-1,1) 
      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,-1) = 
     $  + 1.6449340668482264d+00 
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  - HX2(0, -1) 
      if ( nw.gt.2 ) then 
      HY3(0,0,-1) = 
     $  - 1.6449340668482264d+00*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3(0,0, -1) 
      HY3(0,-1,-1) = 
     $  + 1.2020569031595942d+00 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,-1) 
     $  - HX3(0, -1,-1) 
     $  - HX3(0,0, -1) 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,-1) = 
     $  + 1.8940656589944918d+00 
     $  + 8.2246703342411321d-01*HX1(0)*HX1(0) 
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - HX4(0,0,0, -1) 
      HY4(0,0,-1,-1) = 
     $  - 1.8940656589944918d+00 
     $  - 1.2020569031595942d+00*HX1(0) 
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX3(0,0,-1) 
     $  + HX4(0,0, -1,-1) 
     $  + 2.0000000000000000d+00*HX4(0,0,0,-1) 
      HY4(0,-1,-1,-1) = 
     $  + 1.0823232337111381d+00 
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1) 
     $  + HX1(0) *HX3(0,-1,-1) 
     $  + HX1(0) *HX3(0,0,-1) 
     $  - HX4(0, -1,-1,-1) 
     $  - HX4(0,0, -1,-1) 
     $  - HX4(0,0,0, -1) 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,0) or (-1,1) endif 
** (n1,n2) = (-1,1) -- completion 
      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HY2(-1,1) = 
     $  + 2.4674011002723396d+00 
     $  + HX1( -1)*HX1(0) 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  + HX2( -1,1) 
     $  - HX2(0, -1) 
     $  - HX2(0,1) 
      Hi2(-1,1) = 
     $  - 6.9314718055994530d-01 
     $  + HX1( -1) 
     $  - HX1(0) 
      if ( nw.gt.2 ) then 
      HY3(0,-1,1) = 
     $  - 2.5190015545588625d+00 
     $  - 2.4674011002723396d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX2(0,-1) 
     $  - HX3(0, -1,1) 
     $  + 2.0000000000000000d+00*HX3(0,0,-1) 
     $  + HX3(0,0,1) 
      HY3(0,1,-1) = 
     $  + 4.3220869092982539d+00 
     $  + 2.4674011002723396d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,1) 
     $  - HX3(0,0, -1) 
     $  - 2.0000000000000000d+00*HX3(0,0,1) 
     $  - HX3(0,1, -1) 
      HY3(-1,-1,1) = 
     $  - 2.7620719062289241d+00 
     $  + 2.4674011002723396d+00*HX1(-1) 
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0) 
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0) 
     $  - HX1( -1)*HX2(0,-1) 
     $  - HX1( -1)*HX2(0,1) 
     $  - 2.4674011002723396d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3( -1,-1,1) 
     $  + HX3(0, -1,-1) 
     $  + HX3(0,0, -1) 
     $  + HX3(0,0,1) 
     $  + HX3(0,1, -1) 
      HY3(-1,1,1) = 
     $  + 2.7620719062289241d+00 
     $  - 4.9348022005446793d+00*HX1(-1) 
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0) 
     $  + 4.9348022005446793d+00*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(-1,1) 
     $  - HX1(0) *HX2(0,-1) 
     $  - HX1(0) *HX2(0,1) 
     $  + HX3( -1,1,1) 
     $  - HX3(0, -1,1) 
     $  + HX3(0,0, -1) 
     $  + HX3(0,0,1) 
     $  - HX3(0,1,1) 
      Hi3(0,-1,1) = 
     $  + 8.2246703342411321d-01 
     $  + 6.9314718055994530d-01*HX1(0) 
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  - HX2(0, -1) 
      Hi3(0,1,-1) = 
     $  - 6.9314718055994530d-01*HX1(0) 
      Hi3(-1,-1,1) = 
     $  + 2.4022650695910071d-01 
     $  - 6.9314718055994530d-01*HX1(-1) 
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1) 
     $  - HX1( -1)*HX1(0) 
     $  + 6.9314718055994530d-01*HX1(0) 
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0) 
      Hi3(-1,1,1) = 
     $  + 1.8851605738073271d+00 
     $  + HX1( -1)*HX1(0) 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  + HX2( -1,1) 
     $  - HX2(0, -1) 
     $  - HX2(0,1) 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,-1,1) = 
     $  + 3.9234217222028759d+00
     $  + 2.5190015545588625d+00*HX1(0)
     $  + 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX3(0,0,-1)
     $  + HX4(0,0, -1,1)
     $  - 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  - HX4(0,0,0,1) 
      HY4(0,0,1,-1) = 
     $  - 4.1940025306306604d+00
     $  - 4.3220869092982539d+00*HX1(0)
     $  - 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1(0) *HX3(0,0,1)
     $  + HX4(0,0,0, -1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,1)
     $  + HX4(0,0,1, -1)
      HY4(0,-1,0,1) = 
     $  + 9.4703282949724591d-01
     $  + 1.8030853547393914d+00*HX1(0)
     $  + 1.6449340668482264d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,0,-1)
     $  - 3.2898681336964528d+00*HX2(0,-1)
     $  + HX4(0, -1,0,1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  - HX4(0,0,0,1) 
      HY4(0,-1,-1,1) = 
     $  + 2.5209599327464717d+00
     $  + 2.7620719062289241d+00*HX1(0)
     $  + 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - HX1(0) *HX3(0,-1,-1)
     $  - HX1(0) *HX3(0,0,-1)
     $  - 2.4674011002723396d+00*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX2(0,-1)*HX2(0,-1)
     $  - HX4(0, -1,-1,1)
     $  + HX4(0, -1,0,1)
     $  + HX4(0,0, -1,1)
     $  - HX4(0,0,0,1) 
      HY4(0,-1,1,-1) = 
     $  - 8.5266539820739622d+00
     $  - 5.5241438124578482d+00*HX1(0)
     $  - 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  + HX1(0) *HX3(0,-1,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,0,-1)
     $  - HX1(0) *HX3(0,0,1)
     $  + 2.4674011002723396d+00*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX2(0,-1)*HX2(0,-1)
     $  - HX4(0, -1,0,1)
     $  - HX4(0, -1,1,-1)
     $  + 2.0000000000000000d+00*HX4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,-1,1)
     $  + 4.0000000000000000d+00*HX4(0,0,0,-1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,1)
     $  + HX4(0,0,1, -1)
      HY4(0,1,-1,-1) = 
     $  + 5.8027584430066521d+00
     $  + 2.7620719062289241d+00*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  + HX1(0) *HX3(0,0,-1)
     $  + 2.0000000000000000d+00*HX1(0)*HX3(0,0,1)
     $  + HX1(0) *HX3(0,1,-1)
     $  - HX4(0,0, -1,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,0,-1)
     $  - 3.0000000000000000d+00*HX4(0,0,0,1)
     $  - 2.0000000000000000d+00*HX4(0,0,1,-1)
     $  - HX4(0,1, -1,-1)
      HY4(0,-1,1,1) = 
     $  + 6.2689427375197987d-01
     $  - 2.7620719062289241d+00*HX1(0)
     $  - 2.4674011002723396d+00*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - HX1(0) *HX3(0,-1,1)
     $  + 2.0000000000000000d+00*HX1(0)*HX3(0,0,-1)
     $  + HX1(0) *HX3(0,0,1)
     $  + 4.9348022005446793d+00*HX2(0,-1)
     $  - HX4(0, -1,1,1)
     $  + 2.0000000000000000d+00*HX4(0,0,-1,1)
     $  - 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,0,1)
     $  + HX4(0,0,1,1) 
      HY4(0,1,-1,1) = 
     $  - 4.3326514514433017d+00
     $  - 1.3169446513992682d+00*HX1(0)
     $  - 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  - HX1(0) *HX3(0,0,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,0,1)
     $  - HX1(0) *HX3(0,1,-1)
     $  + HX2(0, -1)*HX2(0,1)
     $  - 2.4674011002723396d+00*HX2(0,1)
     $  + 5.0000000000000000d-01*HX2(0,1)*HX2(0,1)
     $  - HX4(0, -1,0,1)
     $  - 3.0000000000000000d+00*HX4(0,0,-1,1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  + 4.0000000000000000d+00*HX4(0,0,0,1)
     $  - 2.0000000000000000d+00*HX4(0,0,1,1)
     $  - HX4(0,1, -1,1)
      HY4(0,1,1,-1) = 
     $  - 1.5001934240460787d-01
     $  + 4.0790165576281924d+00*HX1(0)
     $  + 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  - HX1(0) *HX3(0,0,1)
     $  + HX1(0) *HX3(0,1,1)
     $  - HX2(0, -1)*HX2(0,1)
     $  + 2.4674011002723396d+00*HX2(0,1)
     $  - 5.0000000000000000d-01*HX2(0,1)*HX2(0,1)
     $  + HX4(0, -1,0,1)
     $  + 2.0000000000000000d+00*HX4(0,0,-1,1)
     $  - HX4(0,0,0, -1)
     $  + HX4(0,0,1, -1)
     $  - HX4(0,1,1, -1)
      HY4(-1,-1,-1,1) = 
     $  + 2.4278628067547031d+00
     $  - 2.7620719062289241d+00*HX1(-1)
     $  + 1.2337005501361698d+00*HX1(-1)*HX1(-1)
     $  + 1.6666666666666666d-01*HX1(-1)*HX1(-1)*HX1(-1)*HX1(0)
     $  - 2.5000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX2(0,1)
     $  - 2.4674011002723396d+00*HX1(-1)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1( -1)*HX3(0,-1,-1)
     $  + HX1( -1)*HX3(0,0,-1)
     $  + HX1( -1)*HX3(0,0,1)
     $  + HX1( -1)*HX3(0,1,-1)
     $  + 2.7620719062289241d+00*HX1(0)
     $  + 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX4( -1,-1,-1,1)
     $  - HX4(0, -1,-1,-1)
     $  - HX4(0,0, -1,-1)
     $  - HX4(0,0,0, -1)
     $  - HX4(0,0,0,1) 
     $  - HX4(0,0,1, -1)
     $  - HX4(0,1, -1,-1)
      HY4(-1,-1,1,1) = 
     $  + 2.0293560632083841d+00
     $  + 2.7620719062289241d+00*HX1(-1)
     $  - 2.4674011002723396d+00*HX1(-1)*HX1(-1)
     $  + 2.5000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)
     $  + 4.9348022005446793d+00*HX1(-1)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1( -1)*HX1(0)*HX2(0,-1)
     $  - HX1( -1)*HX1(0)*HX2(0,1)
     $  - HX1( -1)*HX3(0,-1,1)
     $  + HX1( -1)*HX3(0,0,-1)
     $  + HX1( -1)*HX3(0,0,1)
     $  - HX1( -1)*HX3(0,1,1)
     $  - 2.7620719062289241d+00*HX1(0)
     $  - 2.4674011002723396d+00*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX3(-1,-1,1)
     $  + HX1(0) *HX3(0,-1,-1)
     $  + HX1(0) *HX3(0,0,-1)
     $  + HX1(0) *HX3(0,0,1)
     $  + HX1(0) *HX3(0,1,-1)
     $  + HX4( -1,-1,1,1)
     $  + HX4(0, -1,-1,1)
     $  + HX4(0, -1,1,-1)
     $  - HX4(0,0, -1,-1)
     $  + HX4(0,0, -1,1)
     $  - 2.0000000000000000d+00*HX4(0,0,0,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,0,1)
     $  - HX4(0,0,1, -1)
     $  + HX4(0,0,1,1) 
     $  + HX4(0,1, -1,1)
     $  + HX4(0,1,1, -1)
      HY4(-1,1,1,1) = 
     $  - 6.4865749331714713d+00
     $  - 4.9348022005446793d+00*HX1(-1)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  + 2.4674011002723396d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(-1,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  + HX1(0) *HX3(-1,1,1)
     $  - HX1(0) *HX3(0,-1,1)
     $  + HX1(0) *HX3(0,0,-1)
     $  + HX1(0) *HX3(0,0,1)
     $  - HX1(0) *HX3(0,1,1)
     $  - 4.9348022005446793d+00*HX2(-1,1)
     $  + 4.9348022005446793d+00*HX2(0,-1)
     $  + 4.9348022005446793d+00*HX2(0,1)
     $  + HX4( -1,1,1,1)
     $  - HX4(0, -1,1,1)
     $  + HX4(0,0, -1,1)
     $  - HX4(0,0,0, -1)
     $  - HX4(0,0,0,1) 
     $  + HX4(0,0,1,1) 
     $  - HX4(0,1,1,1) 
      Hi4(0,0,-1,1) = 
     $  - 9.0154267736969571d-01 
     $  - 8.2246703342411321d-01*HX1(0) 
     $  - 3.4657359027997265d-01*HX1(0)*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3(0,0, -1) 
      Hi4(0,0,1,-1) = 
     $  + 3.4657359027997265d-01*HX1(0)*HX1(0) 
      Hi4(0,-1,0,1) = 
     $  + 1.8030853547393914d+00 
     $  + 8.2246703342411321d-01*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,-1) 
     $  - 2.0000000000000000d+00*HX3(0,0,-1) 
      Hi4(0,-1,-1,1) = 
     $  + 4.8170908494321862d-01 
     $  - 2.4022650695910071d-01*HX1(0) 
     $  - 3.4657359027997265d-01*HX1(0)*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,-1) 
     $  + 6.9314718055994530d-01*HX2(0,-1) 
     $  - HX3(0, -1,-1) 
     $  - HX3(0,0, -1) 
      Hi4(0,-1,1,-1) = 
     $  + 5.7009070532142637d-01 
     $  + 4.8045301391820142d-01*HX1(0) 
     $  + 3.4657359027997265d-01*HX1(0)*HX1(0) 
     $  - 6.9314718055994530d-01*HX2(0,-1) 
      Hi4(0,1,-1,-1) = 
     $  - 2.4022650695910071d-01*HX1(0) 
      Hi4(0,-1,1,1) = 
     $  - 2.7620719062289241d+00 
     $  - 1.8851605738073271d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX2(0,-1) 
     $  - HX3(0, -1,1) 
     $  + 2.0000000000000000d+00*HX3(0,0,-1) 
     $  + HX3(0,0,1) 
      Hi4(0,1,-1,1) = 
     $  + 2.6736902858507163d+00 
     $  + 1.3029200473423146d+00*HX1(0) 
     $  + 3.4657359027997265d-01*HX1(0)*HX1(0) 
     $  + 1.6666666666666665d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,1) 
     $  + 6.9314718055994530d-01*HX2(0,1) 
     $  - HX3(0,0, -1) 
     $  - 2.0000000000000000d+00*HX3(0,0,1) 
     $  - HX3(0,1, -1) 
      Hi4(0,1,1,-1) =  
     $  + 1.1401814106428527d+00 
     $  + 5.8224052646501250d-01*HX1(0) 
     $  - 3.4657359027997265d-01*HX1(0)*HX1(0) 
     $  - 6.9314718055994530d-01*HX2(0,1) 
      Hi4(-1,-1,-1,1) = 
     $  - 5.5504108664821579d-02
     $  + 2.4022650695910071d-01*HX1(-1)
     $  - 3.4657359027997265d-01*HX1(-1)*HX1(-1)
     $  + 1.6666666666666666d-01*HX1(-1)*HX1(-1)*HX1(-1)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)
     $  + 6.9314718055994530d-01*HX1(-1)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0)
     $  - 2.4022650695910071d-01*HX1(0)
     $  - 3.4657359027997265d-01*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)
      Hi4(-1,-1,1,1) = 
     $  - 2.4532465311320902d+00
     $  + 1.8851605738073271d+00*HX1(-1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0)
     $  - HX1( -1)*HX2(0,-1)
     $  - HX1( -1)*HX2(0,1)
     $  - 1.8851605738073271d+00*HX1(0)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)
     $  + HX3( -1,-1,1)
     $  + HX3(0, -1,-1)
     $  + HX3(0,0, -1)
     $  + HX3(0,0,1) 
     $  + HX3(0,1, -1)
      Hi4(-1,1,1,1) = 
     $  - 5.5504108664821579d-02
     $  - 1.6449340668482264d+00*HX1(-1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0)
     $  + 1.6449340668482264d+00*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX2(-1,1)
     $  - HX1(0) *HX2(0,-1)
     $  - HX1(0) *HX2(0,1)
     $  + HX3( -1,1,1)
     $  - HX3(0, -1,1)
     $  + HX3(0,0, -1)
     $  + HX3(0,0,1) 
     $  - HX3(0,1,1) 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,1) -- completion endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplin1(y,nw,HY1,HY2,HY3,HY4,n1,n2) 
** evaluates the irreducible HPL for y =1
** it is guaranteed that nw is in the range 2:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
      implicit double precision (a-h,o-z) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,1) = 
     $  + 1.6449340668482264d+00
      if (nw.gt.2) then 
      HY3(0,0,1) = 
     $  + 1.2020569031595942d+00
      HY3(0,1,1) = 
     $  + 1.2020569031595942d+00
      endif
      if (nw.gt.3) then 
      HY4(0,0,0,1) = 
     $  + 1.0823232337111381d+00
      HY4(0,0,1,1) = 
     $  + 2.7058080842778454d-01
      HY4(0,1,1,1) = 
     $  + 1.0823232337111381d+00
      endif
      endif
** (n1,n2) = (0,1) or (-1,1) endif 
************ 
** (n1,n2) = (-1,0) or (-1,1) 
      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,-1) = 
     $  + 8.2246703342411321d-01
      if (nw.gt.2) then
      HY3(0,-1,-1) = 
     $  + 1.5025711289494928d-01
      HY3(0,0,-1) = 
     $  + 9.0154267736969571d-01
      endif
      if (nw.gt.3) then
      HY4(0,-1,-1,-1) = 
     $  + 2.3752366322618485d-02
      HY4(0,0,-1,-1) = 
     $  + 8.7785671568655302d-02
      HY4(0,0,0,-1) = 
     $  + 9.4703282949724591d-01
      endif
      endif 
** (n1,n2) = (-1,0) or (-1,1) endif 
** (n1,n2) = (-1,1) -- completion 
      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HY2(-1,1) = 
     $  + 5.8224052646501250d-01
      if (nw.gt.2) then
      HY3(0,-1,1) = 
     $  + 2.4307035167006157d-01
      HY3(0,1,-1) = 
     $  + 5.0821521280468485d-01
      HY3(-1,-1,1) = 
     $  + 9.4753004230127705d-02
      HY3(-1,1,1) = 
     $  + 5.3721319360804020d-01
      endif
      if (nw.gt.3) then 
      HY4(0,0,-1,1) = 
     $  + 1.1787599965050932d-01
      HY4(0,0,1,-1) = 
     $  + 1.7284527823898438d-01
      HY4(0,-1,0,1) = 
     $  + 2.0293560632083841d-01
      HY4(0,-1,-1,1) = 
     $  + 3.4159126166513913d-02
      HY4(0,-1,1,-1) = 
     $  + 5.4653052738263652d-02
      HY4(0,1,-1,-1) = 
     $  + 1.1412342741606084d-01
      HY4(0,-1,1,1) = 
     $  + 9.3097125991768577d-02
      HY4(0,1,-1,1) = 
     $  + 1.9355535381306524d-01
      HY4(0,1,1,-1) = 
     $  + 4.3369237704895519d-01
      HY4(-1,-1,-1,1) = 
     $  + 1.4134237214990008d-02
      HY4(-1,-1,1,1) = 
     $  + 4.0758239159309251d-02
      HY4(-1,1,1,1) = 
     $  + 5.1747906167389938d-01
      endif
      endif
** (n1,n2) = (-1,1) -- completion endif 
      return
      end
