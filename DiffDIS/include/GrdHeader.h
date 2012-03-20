/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef _GRD_HDR
#define _GRD_HDR
#include <cstdio>
#include <string.h> //--- memset, memcmp
//#include <unistd.h>
#include "xdr++.h"

using namespace std;

//#define HDR_VER 300

class GrdHeader_t {
  //FILE* GridFile;
  public:
    unsigned char Ver;
    unsigned char SubVer;
    //int ver;
    int QCDClib_ver;
    char lbl[64];
    char info[128];
    int Ndistr;
    int x_exp;
    int PiOrd;
    //int reserved[11];
    
    GrdHeader_t() {memset(this, 0, sizeof(GrdHeader_t));}
    bool RW(XDRio_t& xdr);
};

#endif
