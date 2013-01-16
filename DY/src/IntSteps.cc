#include <cmath>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector>

#include "IntSteps.h"

using namespace std;

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
  _etar = new double[2];
  _mr[0] = ist._mr[0]; _mr[1] = ist._mr[1];
  _yr[0] = ist._yr[0]; _yr[1] = ist._yr[1];
  _etar[0] = ist._etar[0]; _etar[1] = ist._etar[1];

  _var_name = ist._var_name;
  _nbins = ist._nbins;
  _bins = new double[_nbins+1];
  for (int ib=0;ib<_nbins+1;ib++){
    _bins[ib] = ist._bins[ib];
  }
  _ssb = new int[_nbins];
  for (int ib=0;ib<_nbins;ib++){
    _ssb[ib] = ist._ssb[ib];
  }

  _nms = ist._nms;
  _nys = ist._nys;
  _msteps = new double[_nms+1];
  for ( int ims=0; ims<_nms+1; ims++) _msteps[ims]=ist._msteps[ims];
  _ysteps = new double[_nys+1];
  for ( int iys=0; iys<_nys+1; iys++) _ysteps[iys]=ist._ysteps[iys];

  _leptPtCut = ist._leptPtCut;
}

IntSteps::IntSteps(const std::string boz, const double *ranges, 
    const string var_name, const int n_bins, const double *var_bins)
{
  using namespace std;
  _boz = boz;
  _mr = new double[2];
  _yr = new double[2];
  _etar = new double[2];
  _mr[0] = ranges[0]; _mr[1] = ranges[1];
  _yr[0] = ranges[2]; _yr[1] = ranges[3];
  _etar[0] = ranges[4]; _etar[1] = ranges[5];
  _leptPtCut = ranges[6];

  _var_name = var_name;
  _nbins = n_bins;
  _bins = new double[n_bins+1];
  _ssb = new int[_nbins];

  if ( string("eta") == var_name ) {
    for (int ib=0;ib<n_bins+1;ib++){
      _bins[ib] = tanh(*((double*)var_bins+ib));
    }
  } else if ( string("y") == var_name ){
    for (int ib=0;ib<n_bins+1;ib++){
      _bins[ib] = *((double*)var_bins+ib);
    }
  } else {
    cout << "IntSteps: binning in selected variables is not supported. \n\
             Exit. " << endl;
    exit(1);
  }

  if ( string("W") == _boz ){
    _makeMstepsW();
  } else if ( string("Z") == _boz ){
    _makeMstepsZ();
  }

  if ( string("Z") == _boz && string("y") == var_name) {
    _makeYstepsBinned();
  } else {
    _makeYsteps();
  }
}

void IntSteps::_makeMstepsZ()
{
  double m_low(_mr[0]), m_up(_mr[1]);
  if ( m_up < 2.*_leptPtCut ) {
    std::cout << "IntSteps: lepton pT cut excludes mass range. \n\
                  Exit"<<std::endl;
    exit(1);
  }
  if ( m_low < 2.*_leptPtCut ) m_low = 2*_leptPtCut;

  double d(0.1);
  double a(m_low);
  vector <double> va;
  while ( a<m_up) {
    if ( a<76.) d = 5.;
    else if ( a>=76. && a<89. ) d = 1.;
    else if ( a>=89. && a<94. ) d = 0.1;
    else if ( a>=94. && a<111. ) d = 1.;
    else if ( a>=111. ) d = 30.*exp(a/2000.);

    va.push_back(a);
    a+=d;
  }
  va.push_back(m_up);

  _nms = va.size()-1;
  _msteps = new double[_nms+1];
  for ( size_t ims=0;ims<va.size();ims++){
    _msteps[ims] = va.at(ims);
  }

  //for (int ims=0; ims<_nms+1; ims++) cout << ims << " " << _msteps[ims] << endl;

}


void IntSteps::_makeMstepsW()
{
  double m_low(_mr[0]), m_up(_mr[1]);
  if ( m_up < 2.*_leptPtCut ) {
    std::cout << "DY calc: lepton pT cut excludes mass range. Exit"<<std::endl;
    exit(1);
  }
  if ( m_low < 2.*_leptPtCut ) m_low = 2*_leptPtCut;

  double d(0);
  double a(m_low);
  vector<double> va;
  while ( a<m_up) {
    if ( a<70.) d = 5.;
    else if ( a>=70. && a<78. ) d = 1.;
    else if ( a>=78. && a<83. ) d = 0.1;
    else if ( a>=83. && a<100. ) d = 1.;
    else if ( a>=100. ) d = 30.*exp(a/2000.);
    
    va.push_back(a);
    a+=d;
  }
  va.push_back(m_up);

  _nms = va.size()-1;
  _msteps = new double[_nms+1];
  for ( size_t ims=0;ims<va.size();ims++){
    _msteps[ims] = va.at(ims);
  }
  //for (int ims=0; ims<_nms+1; ims++) cout << ims << " " << _msteps[ims] << endl;
}


void IntSteps::_makeYsteps()
{
  double y_low(_yr[0]), y_up(_yr[1]);

  double d(0);
  double a(y_low);
  vector<double> va;
  while ( a < y_up ) {
    if ( a<-5.5) d = 0.5;
    else if ( a>=-5.5 && a<-2.1 ) d = 0.2;
    else if ( a>=-2.1 && a<2.1 ) d = 0.2;
    else if ( a>=2.1 && a<5.5 ) d = 0.2;
    else if ( a>=5.5 ) d = 0.5;
    
    va.push_back(a);
    a+=d;
  }
  va.push_back(y_up);

  _nys = va.size()-1;
  _ysteps = new double[_nys+1];
  for ( size_t iys=0;iys<va.size();iys++){
    _ysteps[iys] = va.at(iys);
  }

  //for (int iys=0; iys<_nys+1; iys++) std::cout << iys << " " << _ysteps[iys] << std::endl;

}


void IntSteps::_makeYstepsBinned()
{
  const double d0(0.04);
  double d(0.);
  int ns(0), nsib(0); // number of all steps, number of steps in a bin
  double a(_bins[0]);
  vector<double> va;
  for ( int ib=0; ib<_nbins; ib++){
    _ssb[ib] = ns+=nsib;
//    cout << ns << " " << _ssb[ib] << endl;
    nsib = ceil((_bins[ib+1]-_bins[ib])/d0);
    d = (_bins[ib+1]-_bins[ib])/nsib;
    for ( int isib = 0; isib<nsib; isib++){
      va.push_back(a);
      a+=d;
    }  
  }
  va.push_back(_bins[_nbins]);

  _nys = va.size()-1;
  _ysteps = new double[_nys+1];
  for ( size_t iys=0;iys<va.size();iys++){
    _ysteps[iys] = va.at(iys);
  }

  //for (int iys=0; iys<_nys+1; iys++) std::cout << iys << " " << _ysteps[iys] << std::endl;

}

