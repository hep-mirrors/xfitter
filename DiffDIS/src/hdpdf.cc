/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

/*
  Author: W. Slominski
  Date: 2011-11
*/

#include <cstdio>
#include <cmath>
#include "hdpdf.h"
#include "emsg.h"
#include "grvpi.h"
// #include "grvpi.h"

//#define Z_TEST_

// int hdpdf_t::QQrelax;

//==========================================================
void hdpdf_t::SetOpts(const string& opt) {
  int p;
  if((p=opt.find('c')) != string::npos) {
    if(++p < opt.length()) CommentChar = opt[p];
  }
  if((p=opt.find('v')) != string::npos) Verbose = 1;
}

//==========================================================
void hdpdf_t::ShowInfo() {
  cout << CommentChar <<" Grid label: " << header.lbl << endl;
  //cout << "Reggeon_factor "<< Reggeon_factor <<endl;
  //Qgrid.show();
  cout << CommentChar << " Q^2 range = "<< Qgrid[0] <<" -- "<< Qgrid[Qgrid.GetN()-1]
    <<" GeV^2" <<endl;
  int i0 = Qgrid.GetFixedInd(c_quark);
  if(i0 >= 0) cout << CommentChar << " m_c = "<< sqrt(Qgrid[i0])<<" GeV"<<endl;
  i0 = Qgrid.GetFixedInd(b_quark);
  if(i0 >= 0) cout << CommentChar << " m_b = "<< sqrt(Qgrid[i0])<<" GeV"<<endl;
  i0 = Qgrid.GetFixedInd(t_quark);
  if(i0 >= 0) cout << CommentChar << " m_t = "<< sqrt(Qgrid[i0])<<" GeV"<<endl;
  //cout << "nx "<< Xgrid.nx <<endl;
  cout << CommentChar << " min. x = "<< Xgrid[0] <<endl;
  cout << CommentChar << " N distr. "<< header.Ndistr <<endl;
  cout << CommentChar << "----------------------------"<< endl;
}

//================================================================
void hdpdf_t::pi0xf(real_type x, real_type QQ, real_type f[]) {
//--- f[0:N_FLAVORS]
  int k;
  //--- pion
  double xF[N_FLAVORS+1]; //--- g, uval=dval, usea=dsea, s,c,b
  if(PionOrder) GRVpi_HO(x, QQ, xF); else  GRVpi_LO(x, QQ, xF);
  f[0] = xF[0];
  f[1] = f[2] = xF[1]/2 + xF[2];
  for(k=3; k <= 5; k++) f[k] = xF[k];
  if(N_FLAVORS >= 6) f[6] = 0;
}

//================================================================
void hdpdf_t::fxP(real_type xP, real_type zP, real_type QQ,
    real_type f[], int xpow) {
  real_type flux = Pflux.f(xP);
  Pom(zP, QQ, f, xpow);
  #ifdef Z_TEST_
    for(int pn = 0; pn <= N_FLAVORS; pn++) cout << f[pn] << endl;
  #endif
  for(int pn = 0; pn <= N_FLAVORS; pn++) f[pn] *= flux;
}

//================================================================
void hdpdf_t::fxR(real_type xP, real_type zP, real_type QQ,
    real_type f[], int xpow) {
  real_type flux = Rflux.f(xP);
  real_type cx = pow(zP, xpow-1);
  pi0xf(zP, QQ, f);
  for(int pn = 0; pn <= N_FLAVORS; pn++) f[pn] *= Reggeon_factor*cx*flux;
}

//================================================================
void hdpdf_t::GetDPDF3(real_type xP, real_type zP, real_type QQ,
    real_type f[], int xpow) {
  real_type fP[N_FLAVORS+1], fR[N_FLAVORS+1];
  fxP(xP, zP, QQ, fP, xpow);
  fxR(xP, zP, QQ, fR, xpow);
  for(int pn = 0; pn <= N_FLAVORS; pn++) f[pn] = fP[pn] + fR[pn];
}

//==========================================================
void hdpdf_t::LoadGrids_f(const char* fn) {
  #define LG_ERR_MSG "LoadGrid error"
  #define XDR_R(a) require(xdr.RWdouble(&a),LG_ERR_MSG);
  #define XDR_I(a) require(xdr.RWint(&a), LG_ERR_MSG);

  const char* Magic="XDRG"; 
  // require(sizeof(real_type) == sizeof(double), "real_type must be double");
  FILE* gf=fopen(fn, "rb");
  require(gf, "Cannot open '%s'.", fn);
  printf("\n==> Reading Grid '%s'...\n", fn);
  char mg[4];
  fread(mg, 4, 1, gf);
  require(!memcmp(mg, Magic, 4),"Bad grid file.");

  //require(XDRcreate(&xdrs, gf, XDR_DECODE), "XDR creation failed.");
  XDRio_t xdr(gf, XDRio_t::R);
  require(header.RW(xdr), "Header read failed");
  // if(Verbose) 
  printf("Grid version %d.%d\nGrid label: '%s'\n", header.Ver, header.SubVer, header.lbl);
  //cout << "Grid label: " << header.lbl << endl;
  require(Pflux.RW(xdr), LG_ERR_MSG);
  require(Rflux.RW(xdr), LG_ERR_MSG);
  XDR_R(Reggeon_factor)
  
  require(RW(xdr), LG_ERR_MSG);
  fclose(gf);
  SetNorm(2*M_PI, header.x_exp); //--- must be AFTER reading the grid values
  // GridsLoaded = true;
  // if(Verbose) 
  // cout << "Grids loaded." << endl;
  // PionOrder = header.PiOrd;
}
