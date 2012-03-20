/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

/*
  Author: Wojtek Slominski
  Date: 2009-12-17 22:22
  Description:
    non-class access to hdpdf_t.
*/

#include "gZDPDF.h"

//#define gZDPDF_DEBUG

static hdpdf_t* gZDPDF;
static hdpdf_t* gZDPDFarr00[1+2*MAX_N_ERR];
static hdpdf_t** gZDPDFarr = gZDPDFarr00 +MAX_N_ERR;
static int CurGridIndex;
static int DeltaMode;
#define _INLINE_ inline

//=======================================================================
hdpdf_t* CurGrid() {return gZDPDFarr[CurGridIndex];}

//=======================================================================
_INLINE_ void ZeusDpdf3Pom(double xP, double zP, double QQ, double f[7], int xpow) {
  #ifdef gZDPDF_DEBUG
  static int FirstCall=1;
  if(FirstCall) {
    FirstCall = 0;
    cout << "DeltaMode = " << DeltaMode << endl;
    cout << "ZeusDpdf3Pom: "
    << xP <<" "<< zP <<" "<< QQ <<" "<< xpow << endl;
  }
  #endif
  if(CurGridIndex && DeltaMode) {
    double fm[N_FLAVORS+1];
    gZDPDFarr[CurGridIndex]->fxP(xP, zP, QQ, f, xpow);
    gZDPDFarr[-CurGridIndex]->fxP(xP, zP, QQ, fm, xpow);
    for(int pn = 0; pn <= N_FLAVORS; pn++) f[pn] = 0.5*(f[pn] - fm[pn]);
  }
  else gZDPDF->fxP(xP, zP, QQ, f, xpow);
}

//=======================================================================
_INLINE_ void ZeusDpdf3Reg(double xP, double zP, double QQ, double f[7], int xpow) {
  #ifdef gZDPDF_DEBUG
  static int FirstCall=1;
  if(FirstCall) {
    FirstCall = 0;
    cout << "DeltaMode = " << DeltaMode << endl;
    cout << "ZeusDpdf3Reg: "
    << xP <<" "<< zP <<" "<< QQ <<" "<< xpow << endl;
  }
  #endif
  if(CurGridIndex && DeltaMode) {
    double fm[N_FLAVORS+1];
    gZDPDFarr[CurGridIndex]->fxR(xP, zP, QQ, f, xpow);
    gZDPDFarr[-CurGridIndex]->fxR(xP, zP, QQ, fm, xpow);
    for(int pn = 0; pn <= N_FLAVORS; pn++) f[pn] = 0.5*(f[pn] - fm[pn]);
  }
  else gZDPDF->fxR(xP, zP, QQ, f, xpow);
}

//=======================================================================
_INLINE_ void ZeusDpdf3(double xP, double zP, double QQ, double f[7], int xpow) {
  #ifdef gZDPDF_DEBUG
  static int FirstCall=1;
  if(FirstCall) {
    FirstCall = 0;
    cout << "DeltaMode = " << DeltaMode << endl;
    cout << "ZeusDpdf3: "
    << xP <<" "<< zP <<" "<< QQ <<" "<< xpow << endl;
  }
  #endif
  if(CurGridIndex && DeltaMode) {
    double fm[N_FLAVORS+1];
    gZDPDFarr[CurGridIndex]->GetDPDF3(xP, zP, QQ, f, xpow);
    gZDPDFarr[-CurGridIndex]->GetDPDF3(xP, zP, QQ, fm, xpow);
    for(int pn = 0; pn <= N_FLAVORS; pn++) f[pn] = 0.5*(f[pn] - fm[pn]);
  }
  else gZDPDF->GetDPDF3(xP, zP, QQ, f, xpow);
}

//=======================================================================
void LoadGrid(const string& Label, int n) {
  if(abs(n) > MAX_N_ERR) {
    cerr << "Grid index "<< n <<" outside the allowed range!" << endl;
    exit(2);
  }
  if(gZDPDFarr[n]) delete gZDPDFarr[n];
  gZDPDF = gZDPDFarr[CurGridIndex=n] = new hdpdf_t(Label, n);
  if(!n) gZDPDF->ShowInfo();
}

//=======================================================================
void Set_tmin(double t) {
  //--- must be called after loading ALL grids
  for(int n=-MAX_N_ERR; n <= MAX_N_ERR; n++) if(gZDPDFarr[n]) {
    gZDPDFarr[n]->Pflux.Settmin(t);
    gZDPDFarr[n]->Rflux.Settmin(t);
  }
}

//=======================================================================
_INLINE_ void SelectGrid(int n) {gZDPDF = gZDPDFarr[CurGridIndex=n];}

//=======================================================================
_INLINE_ void SetDeltaMode(int dm) {DeltaMode = dm;}

//=======================================================================
//---   FORTRAN   ---

string F2Cstring(char *tx, int len) {
  string ss(tx,len);
  while(ss[--len] == ' ');
  ss.erase(len+1);
  return ss;
}

#ifdef __cplusplus
extern "C" {
#endif

  void loadgrid_(char* Label, int* n, int len) {
    int n1 = *n;
    LoadGrid(F2Cstring(Label, len), n1);
  }

  void selectgrid_(int* n) {SelectGrid(*n);}

  void setdeltamode_(int* n) {SetDeltaMode(*n);}

  void dpdf3p_(double* xP, double* zP, double* QQ, double* f, int* xpow/*, int flen*/) {
    ZeusDpdf3Pom(*xP, *zP, *QQ, f, *xpow);
  }

  void dpdf3r_(double* xP, double* zP, double* QQ, double* f, int* xpow/*, int flen*/) {
    ZeusDpdf3Reg(*xP, *zP, *QQ, f, *xpow);
  }

  void dpdf3_(double* xP, double* zP, double* QQ, double* f, int* xpow/*, int flen*/) {
    ZeusDpdf3(*xP, *zP, *QQ, f, *xpow);
  }

  void settmin_(double* t) {
    //--- must be called after loading ALL grids
    Set_tmin(*t);
  }

#ifdef __cplusplus
}
#endif
