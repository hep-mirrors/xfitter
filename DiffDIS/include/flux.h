/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

/*
  Author: Wojtek Slominski
  Date: 2011-11
  Description:
		Regge type flux
		ver. 2.02
	Changes:
		tmax added
*/
#ifndef FLUX_H_
#define FLUX_H_

#include <iostream>
#include <cstdio>
#include <cmath>

#ifndef NO_XDR
#include "xdr++.h"
#endif

/**
  \brief Regge-type flux integrated over \a t
  <!-- \section intro_sec Introduction -->
 
  \f$ \displaystyle
  \Phi(\xi) = \int\limits_{t_{\rm min}}^{t_{\rm max}} dt \, \frac{e^{bt}}{\xi^{2\alpha(t) -1}}
  \f$
  
  where \f$ \alpha(t) = \alpha_0 + \alpha' t \f$
  and
  \f$ \displaystyle
  t_{\rm max} \leq t_{\rm MAX} = -{\xi^2 \over 1-\xi} m_p^2\;.
  \f$
  
  See SetParams().
*/
// oooooooooooooooooooooooooooooooooooooooooooooooo
class Flux_t {
  // static const double mp2 = 0.93827231*0.93827231;
  double a0;
  double a1;
  double b0;
  double tmin,tmax;
  double om2a0; //--- 1 - 2*a0
  // double cJLAP;

  public:
    //================================================
    void SetParams(double _a0,  ///< \f$ \alpha_0 \f$
                   double _a1,  ///< \f$ \alpha' \f$
                   double _b,  ///< \a b
                   double _tmin=-1, ///< min. \a t value
                   double _tmax=0 ///< max. \a t value; 0 means t_MAX
                   ) {
      SetA0(_a0);
      a1 = _a1;
      b0 = _b;
      // cJLAP = (1 - exp(-b0))/b0;
      tmin = _tmin;
      tmax = _tmax;
    }

    //================================================
    Flux_t() {}

    //================================================
    Flux_t(double _a0, double _a1, double _b, double _tmin=-1, double _tmax=0) {
      SetParams(_a0, _a1, _b, _tmin, _tmax);
    }

    //================================================
    void SetA0(double _a0) {
      a0 = _a0;
      om2a0 = 1 - 2*a0;
    }

    //================================================
    void Settmin(double _tmin) {
      tmin = _tmin;
    }

    //================================================
    void Settrng(double _tmin, double _tmax=0) {
      tmin = _tmin;
      tmax = _tmax;
    }

    // //================================================
    // double f0(double xP) {
      // return cJLAP*pow(xP, om2a0);
    // }

    //================================================
    double f(double xP) {
      const double mp2 = 0.93827231*0.93827231;
      double t0 = tmax ? tmax : -mp2*xP*xP/(1-xP);
      double b = b0 - 2*a1*log(xP);
      return (exp(b*t0) - exp(b*tmin))*pow(xP, om2a0)/b;
    }

    //================================================
    void Save(FILE* df) {
      fwrite(&a0, sizeof(double), 1, df);
      fwrite(&b0, sizeof(double), 1, df);
      fwrite(&a1, sizeof(double), 1, df);
      fwrite(&tmin, sizeof(double), 1, df);
    }

    //================================================
    void Load(FILE* df) {
      fread(&a0, sizeof(double), 1, df);
      fread(&b0, sizeof(double), 1, df);
      fread(&a1, sizeof(double), 1, df);
      fread(&tmin, sizeof(double), 1, df);
      //cout <<"a0 = "<< a0 << endl;
      SetParams(a0, a1, b0, tmin);
    }

    //==========================================================
		friend ostream &operator<<(ostream &ostr, const Flux_t& f) {
			ostr << "alpha0 = " << f.a0
					 << ", alpha1 = " << f.a1
					 << ", B = " << f.b0
					 << ", " << f.tmin << " < t < " << f.tmax;
			return ostr << endl;
		}		
		
#ifndef NO_XDR
  #define XDR_R(a) if(!xdr.RWdouble(&a)) return 0;
    //==========================================================
    bool RW(XDRio_t& xdr) {
      XDR_R(a0);
      XDR_R(b0);
      XDR_R(a1);
      XDR_R(tmin);
      if(xdr.isReading()) SetParams(a0, a1, b0, tmin);
      return 1;
    }

  #undef XDR_R
#endif

};

// #ifdef _MAIN_
// const double Flux_t::mp2 = 0.93827231*0.93827231;
// #endif

#endif
