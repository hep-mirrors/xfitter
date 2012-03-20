/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

/* \file hdpdf.h
  \author W. Slominski
  \version 2.12
  \date 2011-11-13
*/

#ifndef ZDPDF_HDR_
#define ZDPDF_HDR_

#define ZDPDF_VER 2.12

#include <iostream>
#include "pdf_base.h"

#ifndef NO_XDR
  #include "GrdHeader.h"
  //#define GRD_EXT ".xdr"
#else
  #include "GrdHeaderBin.h"
  //#define GRD_EXT ".bin"
#endif

#include "flux.h"

#define N_FLAVORS 6
// #define N_PDFS 4

/**
  \brief 2-component (Pomeron + Reggeon) diffractive PDFs.
  
  <!-- \section intro_sec Introduction -->
  Integrated over \a t.
 
  \f$ \def\Pom{{I\!\!P}} \def\Reg{{I\!\!R}}
  f^{D(3)}(\xi,\beta,Q^2) = \Phi_\Pom(\xi)\, f_\Pom(\beta,Q^2)
  + N_\Reg \Phi_\Reg(\xi)\, f_\Reg(\beta,Q^2) 
  \f$
  
  where \f$ \def\Reg{{I\!\!R}} N_\Reg \f$ = \c Reggeon_factor.

  This is the main class for ZEUS diffractive PDFs 2009.

  \internal
  \todo N_FLAVORS = 5 or 6 as a parameter.  
  \endinternal
*/
// oooooooooooooooooooooooooooooooooo
class hdpdf_t : public pdf_base_t {
  double Reggeon_factor;
  int PionOrder;
  // bool GridsLoaded;
  int Verbose;
  char CommentChar;
  GrdHeader_t header;
  string Label;
  int ErrInd;

  // void ReadGrid(XDRio_t&, int pn);
  void LoadGrids_f(const char* fn);
  // void LoadGrids();
  
  //==========================================================
  void LoadGrids() {
    require(sizeof(real_type) == sizeof(double), "real_type must be double");
    string gname(Label);
    if(ErrInd) {
      char ts[8];
      if(ErrInd > 0) sprintf(ts, "_P%02d", ErrInd);
      else sprintf(ts, "_M%02d", -ErrInd);
      gname += ts;
    }
    gname += ".xdr";
    LoadGrids_f(gname.c_str());
    // fi[3] = fi[2] = fi[1];
    // GridsLoaded = true;
    if(Verbose) cout << "Grids loaded.\n" << endl;
    PionOrder = header.PiOrd;
  }

  public:
    // static int QQrelax;
    Flux_t Pflux;
    Flux_t Rflux;

    // ==========================================
    hdpdf_t(const string& Lbl, ///< grid file name w/o extension
      int ierr=0 ///< error index
      ) : pdf_base_t(5, true) {
      Verbose = 1;
      CommentChar = '#';
      // GridsLoaded = 0;
      Label = Lbl;
      ErrInd = ierr;
      LoadGrids();
      //QQrelax = 0;
    }

    void SetOpts(const string& opt);
    void ShowInfo();
    // void ClipQQ(bool y=1) {QQrelax = y;}
    double GetReggeon_factor() {return Reggeon_factor;}
    int GetPionOrder() {return PionOrder;}

    // ===================================
    real_type Pom1(real_type zP, real_type QQ, int pn, int xpow=1) {
      if(pn == 6) return 0;
      if(pn == 2) pn = 1;
      else if(pn == 3) pn = 1;
      return Val(zP, QQ, pn, xpow);
    }
    
    // ===================================
    void Pom(real_type zP, real_type QQ, real_type f[], int xpow=1) {
      Vals(zP, QQ, f, xpow);
      f[3] = f[2] = f[1];
      if(N_FLAVORS >= 6) f[6] = 0;
    }
    
    // ===================================
    real_type fxP1(real_type xP, real_type zP, real_type QQ, int pn, int xpow=1) {
      return Pflux.f(xP)*Pom1(zP, QQ, pn, xpow);
    }
    
    // ===================================
    void fxP(real_type xP, real_type zP, real_type QQ, real_type f[], int xpow=1);
    
    void pi0xf(real_type x, real_type QQ, real_type f[]);
    void fxR(real_type xP, real_type zP, real_type QQ, real_type f[], int xpow=1);
    void GetDPDF3(real_type xP, real_type zP, real_type QQ, real_type f[], int xpow=1);
};

#endif
