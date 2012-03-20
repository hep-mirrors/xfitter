/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#include <iostream>
#include <cstdio>
//#include <iomanip>
//#include <cstdlib>
#include "emsg.h"

using namespace std;

#include "qgrid_base.h"
  
//==============================================
void qgrid_base_t::Show(ostream& out, bool full) {
  char buf[128];
  out << "Q^2 grid: " << npt <<" points, Q^2 range = [" << qqVal(tval[0]) <<", "<< qqVal(tval[npt-1]) <<"]"<< endl;
  out << "Q0^2 = "<< qqVal(tval[it0]) << endl;
  if(full) {
    out << "Q^2 factor = "<< qqVal(tStep) << endl;
    sprintf(buf, "Fixed points at %d %d %d %d", it0, it_c, it_b, it_t); out << buf << endl;
    sprintf(buf, "%3s %12s %8s  %4s", "ind", "Q^2 ", "t ", "Nf"); out << buf << endl;
    int ind;
    for(ind = 0; ind < npt; ind++) {
      // sprintf(buf, "%3d %10.3f %6.3f %2d", ind, qqVal(tval[ind]), tval[ind], nFlavors[ind]);
      sprintf(buf, "%3d %13.4f %9.4f  %d", ind, qqVal(tval[ind]), tval[ind], nFlavors[ind]);
      out << buf << endl;
    }
  }
#ifdef TEST
	out << "Init_q: " << npt <<", "<< tlo << " - " << thi << endl;
#endif
  out << "--------------------------------------" << endl;
}

//=================================================
bool qgrid_base_t::Save(FILE* df) {
  //FILE* df = fopen(fname, "ab");
  int nb= sizeof(qgrid_base_t);
  fwrite(&nb, sizeof(int) , 1, df);
  fwrite(this, nb, 1, df);
  fwrite(nFlavors, sizeof(nFlavors[0]), npt, df);
  fwrite(tval, sizeof(tval[0]), npt, df);
  //fclose(df);
  return 1;
}

//=================================================
bool qgrid_base_t::Load(FILE* df) {
  int nb;
  fread(&nb, sizeof(int) , 1, df);
  require(nb == sizeof(qgrid_base_t), "Bad sizeof(qgrid_base_t)");
  fread(this, sizeof(qgrid_base_t), 1, df);
  //cout << "FILE LINE: " << __FILE__ <<": "<< __LINE__ << endl;
  //cout << "npt: " << npt << hex << npt << endl;
  nFlavors = new int[npt];
  fread(nFlavors, sizeof(nFlavors[0]), npt, df);
  tval = new real_type[npt];
  fread(tval, sizeof(tval[0]), npt, df);
  return OK=1;
}

#ifndef NO_XDR
//=================================================
bool qgrid_base_t::RW(XDRio_t& xdr) {
  #define XDR_R(a) if(!xdr.RWdouble(&a)) return 0;
  #define XDR_I(a) if(!xdr.RWint(&a)) return 0;

  OK = 0;

  int N = npt;
  XDR_I(N)
  XDR_I(it0)
  XDR_I(it_c)
  XDR_I(it_b)
  XDR_I(it_t)

  XDR_R(QQstart)
  XDR_R(QQlo)
  XDR_R(QQhi)
  XDR_R(tStep)
  XDR_R(t0)
  XDR_R(tlo)
  XDR_R(thi)

  if(xdr.isReading()) {
    Alloc(N);
    // nFlavors = new int[npt];
    // tval = new real_type[npt];
  }
  if(!xdr.RWintv(nFlavors, npt)) return 0;
  if(!xdr.RWdoublev(tval, npt)) return 0;

  return OK=1;
}

#endif
