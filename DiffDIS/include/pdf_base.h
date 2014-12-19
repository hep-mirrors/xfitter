/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

/*
  Name: pdf_base.h
  Author: W. Slominski
  Date: 2011-11-15
  Description:
    Base class for pdf table over t = ln(Q^2) cross x grid
*/

#ifndef PDF_BASE_H_
#define PDF_BASE_H_

#include <iostream>

// #include "qgrid_base.h"
// #include "xgrid_base.h"
#include "basic_defs.h"
// #ifdef USE_FULL_GRIDS
  #include "qgrid.h"
  #include "xgrid.h"
// #endif

#include "dbgtools.h"

#define Distr_t RealPtr2_t

/// @cond NIEWAZNE
#define QUAL_NAME(n) n
/// @endcond

/**
\brief Base PDF class. Use \c pdf_t
\internal
Arrays for pdf values over the 2-dim grid in t = log(QQ) and x.

PDF(x,QQ) = fi(x,QQ) / (M_Norm_fac * x^M_Norm_xpow)

One array Nq × Nx for each PDF – in general 2+2*Nflav arrays for quarks, anti-quarks, gluon and f_+

The constructor
does not allocate any memory for pdfs. This can only be done, once Q- and X- grid sizes are known.

The only method to fill the grids and pdf values is to read from an XDR stream.
\endinternal
*/
// oooooooooooooooooooooooooooooooooooooooooo
class pdf_base_t {
  constexpr static const double DBLstamp=1234567/4.0;
  constexpr static const int INTstamp=12345;
  
  // --- DATA
  // ------------
  
protected:  
  int M_nFlavors;
  int M_pn_min;
  // int M_nAllocDistr; //!--- # alocated Distr_t arrays
  int M_nX; //--- = Xgrid.npt
  int M_nQ; //--- = Qgrid.npt
  const int M_nQextra; //--- = 1 for R-K
  bool M_FlavorSymm; //--- if pdf[-k] = pdf[k]
public:
  enum PDFmode_t {mode_Quark, mode_PM, mode_SNS};
protected:  
  PDFmode_t M_PDFmode, M_PDFcurmode;
  real_type M_Norm_fac;
  int M_Norm_xpow;
  
  // --- used by io
  // string Label;
  // bool M_used_distr_mem[N_FLAV_INDICES], *M_used_distr;
  
public:
  // #ifdef USE_FULL_GRIDS
    qgrid_t Qgrid;
    xgrid_t Xgrid;
  // #else
    // qgrid_base_t Qgrid;
    // xgrid_base_t Xgrid;
  // #endif
  // xgrid_base_t* Xgrid;
  // #define X_GRD (*Xgrid)
  #define X_GRD Xgrid
  
  Distr_t fi_a[N_FLAV_INDICES], *fi;
  
  // --- METHODS
  // ------------
  
  /*
    fi[iFlav][it][ix] is the PDF value
    iFlav = 0,...,MAX_FLAVORS+1; for M_FlavorSymm == true
    iFlav = -MAX_FLAVORS,...,MAX_FLAVORS+1; for M_FlavorSymm == false
      MAX_FLAVORS+1 = q_sum
    it = 0,...,Nq;  it = Nq is the working space for R-K
    ix = 0,...,Nx-1
  */

protected:
  //===================================================================
  void QUAL_NAME(AllocWksp_1)(Distr_t* f, int pn);
  void QUAL_NAME(AllocWksp)(Distr_t* f);
  void QUAL_NAME(DeAllocWksp_1)(Distr_t* f, int pn);
  void QUAL_NAME(DeAllocWksp)(Distr_t* f);
  void QUAL_NAME(ClearAll)();
  
public:
  // bool isOK() const {return M_nAllocDistr;}
  // bool isA(int pn) const {return M_used_distr[pn];}
  bool isA(int pn) const {return fi[pn] != NULL;}
  // int GetNpAll() {return M_nAllocDistr;} //--- who needs this?
  int GetMaxFlavor() const {return M_nFlavors;}

  // ================================================================
  void SetFlavSymm(bool flav_symm=true) {
    M_FlavorSymm = flav_symm;
    M_pn_min = M_FlavorSymm? 0 : -M_nFlavors;
  }

  /**
    @brief Normalization of the internal tables, fi(x,QQ)
    
    fi(x,QQ) = M_Norm_fac * x^M_Norm_xpow * PDF(x,QQ)
    
    i.e.
    PDF(x,QQ) = fi(x,QQ) / (M_Norm_fac * x^M_Norm_xpow)
  */
  // ================================================================
  void SetNorm(real_type c, int p) {
    M_Norm_fac = c;
    M_Norm_xpow = p;
  }

  // ================================================================
  // pdf_base_t(int nf=5, bool flav_symm=false) : Xgrid(M_Xgrid_base), M_nFlavors(nf) {
  pdf_base_t(int nf=MAX_FLAVORS, bool flav_symm=false) : M_nFlavors(nf) 
    , M_nQextra(1)
  {
    SetFlavSymm(flav_symm);
    SetNorm(1,0);
    M_PDFcurmode = M_PDFmode = mode_Quark;
    // M_nAllocDistr = 0;
    M_nQ = 0;
    M_nX = 0;
    // memset(M_used_distr_mem, 0, sizeof(M_used_distr_mem));
    // M_used_distr = M_used_distr_mem + MAX_FLAVORS;
    
    // --- set pointers to all Distr_t arrays to NULL
    memset(fi_a, 0, sizeof(fi_a));
    fi = fi_a + MAX_FLAVORS;
    // Xgrid = new xgrid_base_t;
    // AllocWksp();
  }
  
  // =======================
  ~pdf_base_t() {DeAllocWksp(fi);}

  // ========================================================
  void QUAL_NAME(Init)() {
    if(!X_GRD.isOK()) throw Fiasco("Xgrid not initialised");
    if(!Qgrid.isOK()) throw Fiasco("Qgrid not initialised");
    AllocWksp(fi);
  }

  // ================================================================
  virtual void Vals(real_type xx, real_type QQ, real_type f[], int xpow_out=0);
  virtual real_type Val(real_type xx, real_type QQ, int pn, int xpow_out=0);

  /**
    Save PDFs(x,QQ) * x^xpow_out
  */
  void QUAL_NAME(Save)(const char* dname, int xpow_out=1);
  
  //==========================================================
  #ifndef NO_XDR
    bool QUAL_NAME(RW1)(XDRio_t& xdr, int pn);
    bool QUAL_NAME(RW)(XDRio_t& xdr);  
  #endif
   
  // void QUAL_NAME(ToPM)(int mt);

};
#undef QUAL_NAME

#endif
