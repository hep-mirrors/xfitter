/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef XGRID_BASE_H_
#define XGRID_BASE_H_

// #include <cstdio>
#include "genut.h"
#include "types.h"
#include <cstring> // for memcpy

#ifndef NO_XDR
  #include "xdr++.h"
#endif

#include "dbgtools.h"

namespace CWeights {
  void SaveWeightsBin(const char* fname);
  void LoadWeightsBin(const char* fname);
}

/// x-grid allocation + file io
// oooooooooooooooooooooooooooooooooooooooooooooooo
//#pragma pack(push,4)
class xgrid_base_t {

  // --- DATA
  
  protected:
    bool OK;
    int npt, smode;
    real_type xlo, xcrit;
    real_type *x, *d;

  // --- METHODS
  
    // ===================
    void Alloc(int n) {
      if(n == npt) return;
      Free();
      if(n > 0) {
        x = new real_type[n];
        d = new real_type[n];
        npt = n;
      }
      // OK = true;
    }

  public:
    // virtual void fill(real_type _xlo, int _nx, int _smode, real_type _xcrit=0) {throw Fiasco("Cannot fill xgrid_base_t");}
    
    // ===================
    void Free() {
      if(!npt) return;
      delete[] x; 
      delete[] d;
      x = NULL;
      d = NULL;
      npt = 0;
      OK = false;
    }

    //enum {LIN, LOGLIN, TANH=5};
    xgrid_base_t() : OK(false), npt(0), x(0), d(0) {}
    ~xgrid_base_t() {delete[] x; delete[] d;}
    bool isOK() const {return OK;}
    real_type operator[] (int j) const {return x[j];}
    real_type at(int j) const {
      if(j<0 || j >= npt) throw out_of_range("Xgrid");
      return x[j];
    }
    const real_type* GetX() const {return x;}
    const real_type* GetDX() const {return d;}
    int GetN() const {return npt;}
    
    #define COPY(a) a = xgrd.a;
    //========================
    xgrid_base_t& operator= (xgrid_base_t& xgrd) {
      DBG_SHOW(xgrd.npt)
      if(this == &xgrd) return *this;
      Alloc(xgrd.npt);
      // if(npt) {delete[] x; delete[] d;}
      // memcpy(this, &xgrd, sizeof(xgrid_base_t));
      COPY(smode)
      COPY(xlo)
      COPY(xcrit)
      if(npt) {
        memcpy(x, xgrd.x, npt*sizeof(real_type));
        memcpy(d, xgrd.d, npt*sizeof(real_type));
        OK = true;
      }
      return *this;
    }
    #undef COPY
    
    bool Load(FILE* df);
    bool Save(FILE* df);
    
  #ifndef NO_XDR
    bool RW(XDRio_t& xdr);
  #endif
  
  friend void CWeights::SaveWeightsBin(const char* fname);
  friend void CWeights::LoadWeightsBin(const char* fname);
};
//#pragma pack(pop)

#endif
