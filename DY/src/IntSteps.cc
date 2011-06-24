#include <cmath>
#include <iostream>
#include <string>
#include <stdlib.h>

#include "IntSteps.h"

IntSteps::~IntSteps()
{
  delete[] _mr;
  delete[] _yr;
  delete[] _msteps;
  delete[] _ysteps;
}

IntSteps::IntSteps(const IntSteps& ist)
{
  _boz = ist._boz;
  _mr = new double[2];
  _yr = new double[2];
  _mr[0] = ist._mr[0]; _mr[1] = ist._mr[1];
  _yr[0] = ist._yr[0]; _yr[1] = ist._yr[1];

  _nms = ist._nms;
  _nys = ist._nys;
  _msteps = new double[_nms+1];
  for ( int ims=0; ims<_nms+1; ims++) _msteps[ims]=ist._msteps[ims];
  _ysteps = new double[_nys+1];
  for ( int iys=0; iys<_nys+1; iys++) _ysteps[iys]=ist._ysteps[iys];

  _leptPtCut = ist._leptPtCut;
}

IntSteps::IntSteps(const std::string boz, const double *m_range, 
          const double *y_range, const double lept_pt_cut)
{
  using namespace std;
  _boz = boz;
  _mr = new double[2];
  _yr = new double[2];
  _mr[0] = m_range[0]; _mr[1] = m_range[1];
  _yr[0] = y_range[0]; _yr[1] = y_range[1];
  _leptPtCut = lept_pt_cut;

  if ( string("W") == _boz ){
    _makeMstepsW();
  } else if ( string("Z") == _boz ){
    _makeMstepsZ();
  }

  _makeYsteps();
}

void IntSteps::_makeMstepsZ()
{
  double m_low(_mr[0]), m_up(_mr[1]);
  if ( m_up < 2.*_leptPtCut ) {
    std::cout << "DY calc: lepton pT cut excludes mass range. Exit"<<std::endl;
    exit(1);
  }
  if ( m_low < 2.*_leptPtCut ) m_low = 2*_leptPtCut;

  double d(0.1);
  double a(m_low), b(m_low+d);
  _nms=0;
  while (1) {
    if ( a<76.) d = 5.;
    else if ( a>=76. && a<89. ) d = 1.;
    else if ( a>=89. && a<94. ) d = 0.1;
    else if ( a>=94. && a<101. ) d = 1.;
    else if ( a>=101. ) d = 5.*exp(a/2000.);
    
    _nms++;
    if ( b>=m_up ) break;

    a= b;
    b+= d;
  }

  _msteps = new double[_nms+1];
  for (int ims=0; ims<_nms+1; ims++) _msteps[ims] = 0.;

  int ims(0);
  d = 0.1; a = m_low; b= m_low+d;
  while (1) {
    if ( a<76.) d = 5.;
    else if ( a>=76. && a<89. ) d = 1.;
    else if ( a>=89. && a<94. ) d = 0.1;
    else if ( a>=94. && a<101. ) d = 1.;
    else if ( a>=101. ) d = 5.*exp(a/2000.);
    
    _msteps[ims]=a;
    if ( b>=m_up ) {
      _msteps[ims+1] = m_up;
      break;
    }

    a= b;
    b+= d;
    ims++;
  }
  //for (int ims=0; ims<_nms+1; ims++) cout << ims << " " << _msteps[ims] << endl;

}



void IntSteps::_makeMstepsW()
{
  double m_low(_mr[0]), m_up(_mr[1]);
  if ( m_low < 2.*_leptPtCut ) m_low = 2*_leptPtCut;

  double d(0.1);
  double a(m_low), b(m_low+d);
  _nms=0;
  while (1) {
    if ( a<65.) d = 5.;
    else if ( a>=65. && a<78. ) d = 1.;
    else if ( a>=78. && a<83. ) d = 0.1;
    else if ( a>=83. && a<90. ) d = 1.;
    else if ( a>=90. ) d = 5.*exp(a/2000.);
    
    _nms++;
    if ( b>=m_up ) break;

    a= b;
    b+= d;
  }

  _msteps = new double[_nms+1];
  for (int ims=0; ims<_nms+1; ims++) _msteps[ims] = 0.;

  int ims(0);
  d = 0.1; a = m_low; b= m_low+d;
  while (1) {
    if ( a<65.) d = 5.;
    else if ( a>=65. && a<78. ) d = 1.;
    else if ( a>=78. && a<83. ) d = 0.1;
    else if ( a>=83. && a<90. ) d = 1.;
    else if ( a>=90. ) d = 5.*exp(a/2000.);
    
    _msteps[ims]=a;
    if ( b>=m_up ) {
      _msteps[ims+1] = m_up;
      break;
    }

    a= b;
    b+= d;
    ims++;
  }
  //for (int ims=0; ims<_nms+1; ims++) cout << ims << " " << _msteps[ims] << endl;

}

void IntSteps::_makeYsteps()
{
  double y_low(_yr[0]), y_up(_yr[1]);

  double d(0.01);
  double a(y_low), b(y_low+d);
  _nys=0;
  while (1) {
  
    if ( a<-6.) d = 0.5;
    else if ( a>=-6. && a<-2.1 ) d = 0.2;
    else if ( a>=-2.1 && a<2.1 ) d = 0.2;
    else if ( a>=2.1 && a<6. ) d = 0.2;
    else if ( a>=6. ) d = 0.5;
  
    _nys++;
    if ( b>=y_up ) break;

    a= b;
    b+= d;
  }

  _ysteps = new double[_nys+1];
  for (int iys=0; iys<_nys+1; iys++) _ysteps[iys] = 0.;

  int iys(0);
  d = 0.01; a = y_low; b= y_low+d;
  while (1) {
    if ( a<-6.) d = 0.5;
    else if ( a>=-6. && a<-2.1 ) d = 0.2;
    else if ( a>=-2.1 && a<2.1 ) d = 0.2;
    else if ( a>=2.1 && a<6. ) d = 0.2;
    else if ( a>=6. ) d = 0.5;
  
    _ysteps[iys]=a;
    if ( b>=y_up ) {
      _ysteps[iys+1] = y_up;
      break;
    }

    a= b;
    b+= d;
    iys++;
  }
 // for (int iys=0; iys<_nys+1; iys++) cout << iys << " " << _ysteps[iys] << endl;

}

