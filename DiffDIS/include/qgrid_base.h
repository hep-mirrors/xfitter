/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef QGRID_BASE_H_
#define QGRID_BASE_H_

#include <cmath>
#include <cstring> // for memcpy

using namespace std;

#include "types.h"
#include "basic_defs.h"
#include "tblipol.h"

#ifndef NO_XDR
  #include "xdr++.h"
#endif

class pdf_base_t;

//#pragma pack(push,4)

/// Q-grid allocation + file io
// ooooooooooooooooooooooooooooooooooooooooooooooooo
class qgrid_base_t {

  // --- DATA
  
  protected:
    bool OK;
    int npt, it0, it_c, it_b, it_t;
    real_type QQstart, QQlo, QQhi, tStep, t0, tlo, thi;
    int *nFlavors;
    real_type *tval;

  // --- METHODS
  
    // ===================
    void Alloc(int n) {
      if(n == npt) return;
      Free();
      if(n > 0) {
        nFlavors = new int[n];
        tval = new real_type[n];
        npt = n;
      }
      // OK = true;
    }

  public:
    // ===================
    void Free() {
      if(!npt) return;
      delete[] nFlavors; 
      delete[] tval;
      nFlavors = NULL;
      tval = NULL;
      npt = 0;
      OK = false;
    }

    bool isOK() const {return OK;}
    real_type tVal(real_type QQ) const { return log(QQ);}
    real_type qqVal(real_type t) const { return exp(t);}

    // ===============================
    qgrid_base_t() : OK(false) {
      npt = 0;
      nFlavors = NULL;
      tval = NULL;
    }

    // ===============================
    ~qgrid_base_t() {
      delete[] nFlavors;
      delete[] tval;
    }
    
    void Show(ostream& out = cout, bool full=false);
    int GetIndQ0() const {return it0;}
    int GetFixedInd(int spec) const {
      switch(spec) {
        case 0: return it0;
        case c_quark: return it_c;
        case b_quark: return it_b;
        case t_quark: return it_t;
        default: throw Fiasco("Illegal arg. of GetFixedInd: %d", spec);
      }
    }
    int nIniFlavors() const {return nFlavors[it0];}
    int GetNFlavors(int iq) const {return nFlavors[iq];}
    real_type getQQ(int iq) const {return qqVal(tval[iq]);}
    real_type Gett(int iq) const {return tval[iq];}
    real_type operator[] (int iq) const {return qqVal(tval[iq]);}
    int GetN() const {return npt;}

    #define COPY(a) a = qgrd.a;
    //========================
    qgrid_base_t& operator= (qgrid_base_t& qgrd) {
      if(this == &qgrd) return *this;
      Alloc(qgrd.npt);
      // if(npt) {delete[] x; delete[] d;}
      // memcpy(this, &xgrd, sizeof(xgrid_base_t));
      COPY(it0)
      COPY(it_c)
      COPY(it_b)
      COPY(it_t)
      COPY(QQstart)
      COPY(QQlo)
      COPY(QQhi)
      COPY(tStep)
      COPY(t0)
      COPY(tlo)
      COPY(thi)
      if(npt) {
        memcpy(nFlavors, qgrd.nFlavors, npt*sizeof(int));
        memcpy(tval, qgrd.tval, npt*sizeof(real_type));
        OK = true;
      }
      return *this;
    }
    #undef COPY
 
    bool Save(FILE* df);
    bool Load(FILE* df);
#ifndef NO_XDR
    bool RW(XDRio_t& xdr);
#endif

  // friend void pdf_base_t::Vals(real_type xx, real_type QQ, real_type f[], int xpow_out);
  friend class pdf_base_t;
};
//#pragma pack(pop)

#endif
