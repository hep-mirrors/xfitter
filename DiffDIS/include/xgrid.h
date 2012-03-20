/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef XGRID_H_
#define XGRID_H_

#include <iostream>
#include "xgrid_base.h"
//#include "emsg.h"

//#pragma pack(push,4)
//#pragma pack(show)

/// x-grid creation + file io
// ooooooooooooooooooooooooooooooooooooooooo
class xgrid_t : public xgrid_base_t {
  public:
    enum {mode_USER=-1, mode_LIN, mode_LOGLIN, mode_TANH=5};
    void kill() {OK = false;}
    xgrid_t() : xgrid_base_t() {}
    // ~xgrid_t() : ~xgrid_base_t() {}
    
    xgrid_t(real_type _xlo, int _nx, int _smode, real_type _xcrit=0) : xgrid_base_t() {
      fill(_xlo, _nx, _smode, _xcrit);
    }
    void fill(real_type _xlo, int _nx, int _smode, real_type _xcrit=0);
    void fill(int _nx, real_type *xvals, real_type *dvals=NULL);
    //real_type operator[](const int i) const {return x[i];}
    
  //==============================================
  void Show(ostream& out=cout, bool full=false) {
    out << "X grid: " << npt <<" points, scl_mode = ";
    switch(smode) {
      case mode_USER: out << "USER"; break;
      case mode_LIN: out << "LIN"; break;
      case mode_LOGLIN: out << "LOGLIN"; break;
      case mode_TANH: out << "TANH"; break;
      default: out << "-unknown-"; break;
    }
    out << endl;
    out << "par = "<< xcrit << endl;
    out << "x_min = "<< x[0] << endl;
    if(full) for(int ind = 0; ind < npt; ind++) out << x[ind] << endl;
    out << "----------------------" << endl;
  }

};
//#pragma pack(pop)
//#pragma pack(show)

#endif
