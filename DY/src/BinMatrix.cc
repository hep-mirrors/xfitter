#include <iostream>
#include <string>
#include <cmath>

#include "PhysPar.h"
#include "BinMatrix.h"
#include "IntSteps.h"
#include <stdlib.h>

using namespace std;

BinMatrix::~BinMatrix()
{
  for (int ims=0;ims<2*_nms+1;ims++){
    for (int iys=0;iys<2*_nys+1;iys++){
      for (int ib=0; ib<_nbins; ib++){
        delete[] BM[ims][iys][ib] ;
      }
      delete[] BM[ims][iys] ;
    }
    delete[] BM[ims] ;
  }
  delete[] BM;

}

void BinMatrix::setBins(const string &var_name, const int n_bins, 
           const double *var_bins)
{
  _var_name = var_name;
  _nbins = n_bins;
  _bins = new double[n_bins+1];

  if ( string("eta") == var_name && 
       ( string("W") == _boz ) ) {
    for (int ib=0;ib<n_bins+1;ib++){
      _bins[ib] = tanh(*((double*)var_bins+ib));
    }
  } else {
    cout << "Binning other than in lepton eta is not supported so far. \n\
             Exit. " << endl;
    exit(1);
  }
    
  if ( string("W") == _boz ){
    BuildBM_W();
  } else if ( string("Z") == _boz ){
    BuildBM_Z();
  }

}

void BinMatrix::BuildBM_Z()
{
}

void BinMatrix::BuildBM_W()
{
  using namespace PhysPar;

  int nmp=2*_nms+1, nyp=2*_nys+1;

  // create and initialize with 0 the bin matrix
  BM = new double***[nmp];
  for (int imp=0;imp<nmp;imp++){ // loop on mass steps
    BM[imp] = new double**[nyp];
    for (int iyp=0;iyp<2*nyp+1;iyp++){ // loop on rapidity steps
      BM[imp][iyp] = new double*[_nbins];
      for (int ib=0; ib<_nbins; ib++){ // loop on eta bins
        BM[imp][iyp][ib] = new double[2];
	for ( int id = 0; id<2; id++){ // loop on direction
          BM[imp][iyp][ib][id] = 0.;
	}
      }
    }
  }

  // assign theory coefficients
  for (int imp=0; imp<nmp; imp++){
    double m = (_msteps[imp/2]+_msteps[(imp+1)/2])/2.;
    double m2 = m*m;
    double X = conhc * pow(gfermi*sqrt(2.)*cw2*m2z/M_PI/_beam_en,2)*m*m2*M_PI/ 48.
                / (pow(m2-m2w,2)+pow(ww*mw,2));
    for(int iyp=0;iyp<nyp;iyp++){
      double y = (_ysteps[iyp/2]+_ysteps[(iyp+1)/2])/2.;
      for(int ib=0;ib<_nbins;ib++){
        // using +- index to specify different directions
	// add 1 since we start with 0 and -0 == +0
        BM[imp][iyp][ib][0] = X*CosthAnIntW(m,y, ib+1)/2.;
        BM[imp][iyp][ib][1] = X*CosthAnIntW(m,-y, -(ib+1))/2.;
	/*
        BM[imp][iyp][ib][0] = X*ai_c(m,y, _bins[ib], _bins[ib+1])/2.;
        BM[imp][iyp][ib][1] = X*ai_c(m,-y, -_bins[ib+1], -_bins[ib])/2.;
	*/
      }
    }
  }
  
}

// integrate in costheta analytically
double BinMatrix::CosthAnIntW(const double &m,const double &y,const int ibp1){
            //const double &c_low_eta_lab, const double &c_up_eta_lab){

  double c_low_eta_lab(0.), c_up_eta_lab(0.);
  if ( ibp1 > 0 ) {
    c_low_eta_lab = _bins[ibp1-1];
    c_up_eta_lab = _bins[ibp1];
  } else if ( ibp1 < 0 ) {
    c_low_eta_lab = -_bins[-ibp1];
    c_up_eta_lab = -_bins[-ibp1-1];
  }

  // calculate actual costh range applying lorentz transformation
  // to LAB cos_theta cuts
  double g = cosh(y); // gamma, lorentz factor
  double b = tanh(y); // beta, velocity

  // for eta defined cuts
  double c_low_eta_cms = costhLT ( b, c_low_eta_lab);
  double c_up_eta_cms = costhLT ( b, c_up_eta_lab);

  // for pt defined ones
  double c_low_pt_cms = - sqrt(1.-4.*pow(_leptPtCut/m,2));
  double c_up_pt_cms  = sqrt(1.-4.*pow(_leptPtCut/m,2));
  double c_low_pt_lab = costhLT ( -b, c_low_pt_cms);
  double c_up_pt_lab = costhLT ( -b, c_up_pt_cms);
  //cout << c_up_eta_lab << "\t" << c_up_eta << "\t" << b << endl;

  double c_low_cms(-1.), c_low_lab(-1.), c_up_cms(1.), c_up_lab(1.);
  if ( c_low_eta_cms < c_low_pt_cms ){
    c_low_cms = c_low_pt_cms;
    c_low_lab = c_low_pt_lab;
  } else {
    c_low_cms = c_low_eta_cms;
    c_low_lab = c_low_eta_lab;
  }

  if ( c_up_eta_cms > c_up_pt_cms ){
    c_up_cms = c_up_pt_cms;
    c_up_lab = c_up_pt_lab;
  } else {
    c_up_cms = c_up_eta_cms;
    c_up_lab = c_up_eta_lab;
  }
  if ( c_up_cms < c_low_cms ) return 0.;
  //cout << c_low_cms << "\t" << c_up_cms << endl;

  double ai = - pow(1.-c_up_cms,3)/3. + pow(1.-c_low_cms,3)/3.;
  return ai;
}

// transforms cosine theta to a frame defined by b and g
// arguments are betta, cosine theta
double BinMatrix::costhLT(const double &b, const double &c0){
  double c = (c0-b)/(1-b*c0);
  return c;
}

