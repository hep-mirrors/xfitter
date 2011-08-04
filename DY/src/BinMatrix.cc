#include <iostream>
#include <string>
#include <cmath>

#include "PhysPar.h"
#include "BinMatrix.h"
#include "IntSteps.h"
#include <stdlib.h>

using namespace std;
using namespace PhysPar;

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

BinMatrix::BinMatrix(const double *be, const IntSteps* ist)
               :IntSteps(*ist)
{
  _beam_en = *be;

  if ( string("eta") == _var_name ) {
    if ( string("W") == _boz ){
      BuildBM_W_eta();
    } else if ( string("Z") == _boz ){
      BuildBM_Z_eta();
    }
  } else if ( string("y") == _var_name && 
       string("Z") == _boz ) {
    BuildBM_Z_y();
  } else {
    cout << "BinMatrix: binning is not supported. \n\
             Exit. " << endl;
    exit(1);
  }
    
}

void BinMatrix::BuildBM_Z_eta()
{
  int nmp=2*_nms+1, nyp=2*_nys+1;

  // create and initialize with 0 the bin matrix
  BM = new double***[2*_nms+1];
  for (int ims=0;ims<2*_nms+1;ims++){ // mass steps
    BM[ims] = new double**[2*_nys+1];
    for (int iys=0;iys<2*_nys+1;iys++){ // boz rapidity steps
      BM[ims][iys] = new double*[_nbins];
      for (int ieb=0; ieb<_nbins; ieb++){ // eta bins
        BM[ims][iys][ieb] = new double[12];
        for (int icdf = 0; icdf<12; icdf++){  // loop on component, direction, flavour
          BM[ims][iys][ieb][icdf] = 0.;
        }
      }
    }
  }

  for (int imp=0; imp<nmp; imp++){
    double m = (_msteps[imp/2]+_msteps[(imp+1)/2])/2.;
    double m2 = m*m;
    const double denZ = pow(m2-m2z,2)+pow(wz*mz,2);
    double X[3] = {0.};
    X[0] = conhc*M_PI*pow(alpha,2)/(3.*m*pow(_beam_en,2));
    X[1] = X[0] * (-2.)*m2*(m2-m2z)/(sw2*cw2*denZ);
    X[2] = X[0] * 3.*pow(m2,2)/(sw2*cw2*denZ*3.*sw2*cw2);
    for(int iyp=0;iyp<nyp;iyp++){
      double y = (_ysteps[iyp/2]+_ysteps[(iyp+1)/2])/2.;

      double x1(0.), x2(0.);
      x1 = m/_beam_en*exp(y);
      x2 = m/_beam_en*exp(-y);
      if ( x1<0. || x1 > 1. || x2<0. || x2 > 1.) continue;
      
      for(int ib=0;ib<_nbins;ib++){
        for(int icdf = 0; icdf<12; icdf++){
	  double AI = CosthAnIntZ(m,y,icdf, _bins[ib], _bins[ib+1]);
	  BM[imp][iyp][ib][icdf] = X[icdf/4]*AI;
	}
      }
    }
  }
}

/**
 * Builds bin matrix for Z process with rapidity binning.
 */
void BinMatrix::BuildBM_Z_y()
{
  int nmp=2*_nms+1, nyp=2*_nys+1;

  // create and initialize with 0 the bin matrix
  BM = new double***[2*_nms+1];
  for (int ims=0;ims<2*_nms+1;ims++){ // mass steps
    BM[ims] = new double**[2*_nys+1];
    for (int iys=0;iys<2*_nys+1;iys++){ // boz rapidity steps
      BM[ims][iys] = new double*[1];
      for (int ieb=0; ieb<1; ieb++){ // eta bin
        BM[ims][iys][ieb] = new double[12];
        for (int icdf = 0; icdf<12; icdf++){  // loop on component, direction, flavour
          BM[ims][iys][ieb][icdf] = 0.;
        }
      }
    }
  }

  for (int imp=0; imp<nmp; imp++){
    double m = (_msteps[imp/2]+_msteps[(imp+1)/2])/2.;
    double m2 = m*m;
    const double denZ = pow(m2-m2z,2)+pow(wz*mz,2);
    double X[3] = {0.};
    X[0] = conhc*M_PI*pow(alpha,2)/(3.*m*pow(_beam_en,2));
    X[1] = X[0] * (-2.)*m2*(m2-m2z)/(sw2*cw2*denZ);
    X[2] = X[0] * 3.*pow(m2,2)/(sw2*cw2*denZ*3.*sw2*cw2);
    for(int iyp=0;iyp<nyp;iyp++){
      double y = (_ysteps[iyp/2]+_ysteps[(iyp+1)/2])/2.;

      double x1(0.), x2(0.);
      x1 = m/_beam_en*exp(y);
      x2 = m/_beam_en*exp(-y);
      if ( x1<0. || x1 > 1. || x2<0. || x2 > 1.) continue;
      
      for(int ieb=0;ieb<1;ieb++){
        for(int icdf = 0; icdf<12; icdf++){
	  double AI = CosthAnIntZ(m,y,icdf, tanh(_etar[0]), tanh(_etar[1]) );
	  BM[imp][iyp][ieb][icdf] = X[icdf/4]*AI;
	}
      }
    }
  }
}

void BinMatrix::BuildBM_W_eta()
{
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
        BM[imp][iyp][ib][0] = X*CosthAnIntW(m,y, _bins[ib], _bins[ib+1])/2.;
        BM[imp][iyp][ib][1] = X*CosthAnIntW(m,-y, -_bins[ib+1], -_bins[ib])/2.;
      }
    }
  }
  
}

// integrate in costheta analytically
double BinMatrix::CosthAnIntW(const double &m,const double &y,
            const double &c_low_eta_lab, const double &c_up_eta_lab){

  // calculate actual costh range applying lorentz transformation
  // to LAB cos_theta cuts
  double b = tanh(y); // beta, velocity

  // for eta defined cuts
  double c_low_eta_cms = costhLT ( b, c_low_eta_lab);
  double c_up_eta_cms = costhLT ( b, c_up_eta_lab);

  // for pt defined ones
  double c_low_pt_cms = - sqrt(1.-4.*pow(_leptPtCut/m,2));
  double c_up_pt_cms  = sqrt(1.-4.*pow(_leptPtCut/m,2));
  double c_low_pt_lab = costhLT ( -b, c_low_pt_cms);
  double c_up_pt_lab = costhLT ( -b, c_up_pt_cms);

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

  double ai = - pow(1.-c_up_cms,3)/3. + pow(1.-c_low_cms,3)/3.;
  return ai;
}

// integrate in costheta analytically
double BinMatrix::CosthAnIntZ(const double &m,const double &ry,
    const int cdf, 
    const double &rc_low_eta_lab, const double &rc_up_eta_lab){

  /*
  double c_low_eta_lab(0.), c_up_eta_lab(0.);
  c_low_eta_lab = _bins[ibp1-1];
  c_up_eta_lab = _bins[ibp1];
  */

  int comp = cdf/4; // component: 0 - gamma, 1 - gZ, 2 - Z
  double dir = pow(-1.,cdf%4/2); // 0, 1 -> 1 (forw) ; 2,3 -> -1 (bw);
  int flav = cdf%4%2; // 0 -> 0 (d), 1 -> 1(u), 2 -> 0(d), 3 -> 1(u)
  // calculate actual costh range applying lorentz transformation
  // to LAB cos_theta cuts
  double c_low_eta_lab(rc_low_eta_lab), c_up_eta_lab(rc_up_eta_lab);
  double y =dir*ry;
  if ( 0 > dir ){
    double ctmp = c_up_eta_lab;
    c_up_eta_lab = -c_low_eta_lab;
    c_low_eta_lab = -ctmp;
  }
  double b = tanh(y); // beta, velocity

  // for eta defined cuts
  double c_low_eta_cms = costhLT ( b, c_low_eta_lab);
  double c_up_eta_cms = costhLT ( b, c_up_eta_lab);

  // for pt defined ones
  double c_low_pt_cms = - sqrt(1.-4.*pow(_leptPtCut/m,2));
  double c_up_pt_cms  = sqrt(1.-4.*pow(_leptPtCut/m,2));
  double c_low_pt_lab = costhLT ( -b, c_low_pt_cms);
  double c_up_pt_lab = costhLT ( -b, c_up_pt_cms);

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

  double a1(0.), a2(0.), ai(0.);
  a1 = c_up_cms + pow(c_up_cms,3)/3. - ( c_low_cms + pow(c_low_cms,3)/3.);
  a2 = (pow(c_up_cms,2)-pow(c_low_cms,2))/2.;
  if ( idlept < 0 ) a2*=-1.; // positron
  if ( 0 == comp ){
    ai = pow(chq[flav],2)*a1;
  } else if ( 1 == comp ){
    ai = chq[flav]*((-0.25+sw2)*vq[flav]/4.*a1+2.*(-0.25)*I3q[flav]/2.*a2);
  } else if ( 2 == comp ){
    ai = (pow(-0.25,2)+pow(-0.25+sw2,2))
       * ( pow(vq[flav]/4.,2) + pow(I3q[flav]/2.,2)) * a1
       + 8.*(-0.25+sw2)*(-0.25)*vq[flav]/4.*I3q[flav]/2. * a2;
  }
  //cout << flav<< " " << dir << " " <<a1 << endl;
  return ai;
}


// transforms cosine theta to a frame defined by b and g
// arguments are betta, cosine theta
double BinMatrix::costhLT(const double &b, const double &c0){
  double c = (c0-b)/(1-b*c0);
  return c;
}

