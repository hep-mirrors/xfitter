/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef QGRID_H_
#define QGRID_H_

#include "qgrid_base.h"
// #include <iostream>
//#include <iomanip>
//#include <cstdlib>

using namespace std;

//#include "emsg.h"
#include "physparams.h"
// #include "tblipol.h"

//#pragma pack(push,4)

/// Q-grid creation + file io
// oooooooooooooooooooooooooooooooooooooooooooo
class qgrid_t : public qgrid_base_t {
    int SetNode(real_type t);
  public:
    void kill() {OK = false;}
    void fill(real_type qqlo, real_type qqhi, real_type qqfac, real_type qq0, const PhysParams_t& ph);
    
    qgrid_t() : qgrid_base_t() {}
    qgrid_t(real_type qqlo, real_type qqhi, real_type qqfac, real_type qq0, const PhysParams_t& ph) {
      fill(qqlo, qqhi, qqfac, qq0, ph);
    }
};
//#pragma pack(pop)

#endif
