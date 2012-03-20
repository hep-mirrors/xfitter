/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef PHYSPARAMS_H_
#define PHYSPARAMS_H_

#include <iostream>
#include <cmath>
#include <cstring>
using namespace std;
#include "types.h"
#include "basic_defs.h"

//-----------------------------------------------
class PhysParams_t {
    bool OK;
  public:
    enum {VFNS, FFNS=3};
    int QCDorder; //--- QCD order
    int nFixedFlavors; //--- 0 means VFNS
    bool isPhoton;
    real_type m[MAX_FLAVORS+1];
    real_type m2[MAX_FLAVORS+1];
    // real_type alpha0, QQ0;
    // private:
      // void SetTranScales(double mc=-1, double mb=-1, double mt=-1);
  public:
    bool isOK() {return OK;}
    
    PhysParams_t() {
      OK = false; 
      isPhoton = false;
      QCDorder = nFixedFlavors = -1;
      memset(m, 0, sizeof(m));
    }
    
    // =================================================
    void set(int _nFixedFlavors, int QCDord,
        double mc=MASS_C, double mb=MASS_B, double mt=MASS_T
        // , double aqcd=ALPHAS_Z0, double Q0=MASS_Z0
    ) {
      if(QCDord >= 0 && QCDord <= MAX_QCD_ORDER) QCDorder = QCDord;
      else throw Fiasco("Illegal QCD order: %d", QCDord);
      nFixedFlavors = _nFixedFlavors;
      // alpha0 = aqcd;
      // QQ0 = Q0*Q0;
      // SetTranScales(mc, mb, mt);
      DBG_SHOW(mc)
      DBG_SHOW(mb)
      DBG_SHOW(mt)
      if(mc >= 0) {m[c_quark] = mc; }
      if(mb >= 0) {m[b_quark] = mb; }
      if(mt >= 0) {m[t_quark] = mt; }
      for(int j=1; j <= MAX_FLAVORS; j++) m2[j] = m[j]*m[j];
      OK = true;
    }
    
    // =================================================
    PhysParams_t(int _nFixedFlavors, int ord,
        double mc=MASS_C, double mb=MASS_B, double mt=MASS_T
        // , double aqcd=ALPHAS_Z0, double Q0=MASS_Z0
        ) {
          isPhoton = false;
          set(_nFixedFlavors, ord, mc, mb, mt);
    }
    
    // void FixAlphaQCD(double alpha, double QQ, double scale=1);
};

#endif
