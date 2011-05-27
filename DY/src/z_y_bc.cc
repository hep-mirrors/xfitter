#include <algorithm>
#include <vector>
#include <math.h>
#include <iostream>
#include <gsl/gsl_integration.h>

#include <time.h>

#include "PhysPar.h"
#include "dy_integ_utils.h"

#define MAX_CHAR 1000
#define NALLOC   100000
//#define PRECISION	1.e-3

typedef double**** double4d;

using namespace std;
using namespace PhysPar;

extern "C"{
void fpdfxq_(int *iset, const double *x, const double *q2, double *pdfs, int *ichk);
}

extern "C" {
  int z_set_ybc__(void *mass_range, void *y_range, double *pel_pt_cut, 
      int *pneb, void *eb);
  int z_get_ybc_xs__(void *xs);
  int z_free_ybc__();
}

extern "C" {
  void dywm_check__();
}

namespace IntZ_YBC
{

int integrate_M(double *, double , double );
int iy(int , double *, double *);
int get_pdfconv(double &, double , double , const double &, double *);
double ai_c(const double &,const double &,
            const double &, const double &);
void set_m_bins();
void set_y_bins();
void fill_bin_matrix();

// eta bins data
int nmb, nyb, neb; // number of integration bins
//int nmp, nyp, nep; // number of simson integration points for mass, y, eta
double *mbins;
double *ybins;
double *cbins;
double *cb_sig_em;
double mr[2], yr[2], el_pt_cut;
// bin matrix to store the coefficients for the bins
// for up-quark coming from forward and backward protons
double4d BM;

int integrate_M(double *pars, double m_low, double m_up){
  // check the pt_cut
  if ( m_up < 2.*el_pt_cut ) {
    for (int ib = 0; ib<neb; ib++)  {
      cb_sig_em[ib]=0.;
    }
    return 0;
  }

  double iya[2*neb], iyb[2*neb], iym[2*neb];
  for (int ib=0;ib<2*neb;ib++){
    iya[ib]=0.; iyb[ib]=0.; iym[ib]=0.;
  }

  iy(0,pars,iya);
  for (int imb=0; imb<nmb; imb++){
    double a = mbins[imb], b = mbins[imb+1];
    //cout << "simpson   " << a << " " << b-a << " " << cb_sig_em[0] << endl;
    iy(2*(imb+1),pars,iyb);
    iy(2*imb+1, pars, iym);
    for (int ib = 0; ib<neb; ib++){
      cb_sig_em[ib]+=(b-a)/6.*(iya[ib]+4.*iym[ib]+iyb[ib]);
      iya[ib]=iyb[ib];
    }
  }

  return 1;
}

int get_pdfconv(double &m, double y, double dir, const double &scale, 
     double &xfxcD, double &xfxcU)
{
  double x1 = m/ebeam*exp(dir*y);
  double x2 = m/ebeam*exp(-dir*y);

  double xq1[13]={0.};
  double xq2[13]={0.};
  int   ichk = 0;
  int   iset = 1;

  double q2 = (scale) * (scale);

  if ( (x1< 1.e-6)||(x2<1.e-6)||(x1>1.-1.e-6)||(x2>1.-1.e-6)) {
    xfxcD=0.;
    xfxcU=0.;
    return 0;
  }
  else{
  // Access QCDNUM pdfs:
    fpdfxq_(&iset,&x1,&q2,xq1,&ichk);
    fpdfxq_(&iset,&x2,&q2,xq2,&ichk);
  }

  // apply beam composition
  double xqt[13] = {0.};
  if ( -1 == ih1 ) {
    if ( 0 < dir ){
      for (int iq=0; iq<13; iq++)  xqt[iq]=xq1[iq];
      for (int iq=0; iq<13; iq++)  xq1[iq]=xqt[12-iq];
    } else {
      for (int iq=0; iq<13; iq++)  xqt[iq]=xq2[iq];
      for (int iq=0; iq<13; iq++)  xq2[iq]=xqt[12-iq];
    }
  }
  if ( -1 == ih2 ) {
    if ( 0 < dir ){
      for (int iq=0; iq<13; iq++)  xqt[iq]=xq2[iq];
      for (int iq=0; iq<13; iq++)  xq2[iq]=xqt[12-iq];
    } else {
      for (int iq=0; iq<13; iq++)  xqt[iq]=xq1[iq];
      for (int iq=0; iq<13; iq++)  xq1[iq]=xqt[12-iq];
    }
  }
  // numbering:
  // \bar{t  b  c  s  u  d} g  d  u  s  c  b  t
  //      0  1  2  3  4  5  6  7  8  9 10 11 12
  xfxcD = 0.; xfxcU = 0.;
  for ( int iflav = 1; iflav<6; iflav++ ) {
    int ifl       = 6+iflav  ;
    int iflbar    = 6-iflav ;
    if ( 1 == iflav%2 ) 
      xfxcD += xq1[ifl]*xq2[iflbar]/(x1*x2);
    else if ( 0 == iflav%2 )
      xfxcU += xq1[ifl]*xq2[iflbar]/(x1*x2);
  }
  return 1;
}

int iy(int imp, double *pars, double *iyval){
  double y_low     = *(double*) pars;
  double y_up      = *((double*)pars+1);
  double y_cent    = (y_low+y_up)/2.;
  double m = (mbins[imp/2]+mbins[(imp+1)/2])/2.;

  const double scale = mz;
  double ica[2*neb];
  for (int ib=0;ib<2*neb;ib++){
    ica[ib]=0.;
  }
  double xfxc[4] = {0.}; // pdf convolutions for dir=-1,+1 and fl=d,u
  //  dir \ flav  |   d |   u
  //  1           |   0 |   1
  //  -1          |   2 |   3
  get_pdfconv(m, y_cent, 1., scale, xfxc[0], xfxc[1]);
  get_pdfconv(m, y_cent, -1., scale, xfxc[2], xfxc[3]);
  for (int ib=0;ib<neb;ib++){
    for (int icdf = 0; icdf<12; icdf ++){
      ica[ib] += xfxc[icdf%4]*BM[imp][0][ib][icdf];
    }
    iyval[ib] = ica[ib];
  }

  return 1;
}

double ai_c(const double &m, const double &ry, const int cdf,
            const double &rc_low_eta_lab, const double &rc_up_eta_lab){

  int comp = cdf/4; // component: 0 - gamma, 1 - gZ, 2 - Z
  double dir = pow(-1.,cdf%4/2); // 0, 1 -> 1 (forw) ; 2,3 -> -1 (bw);
  int flav = cdf%4%2; // 0 -> 0 (d), 1 -> 1(u), 2 -> 0(d), 3 -> 1(u)
  // calculate actual costh range applying lorentz transformation
  // to LAB cos_theta cuts
  double y =dir*ry;
  double c_low_eta_lab(rc_low_eta_lab), c_up_eta_lab(rc_up_eta_lab);
  if ( 0 > dir ){
    double ctmp = c_up_eta_lab;
    c_up_eta_lab = -c_low_eta_lab;
    c_low_eta_lab = -ctmp;
  }
  double g = cosh(y); // gamma, lorentz factor
  double b = tanh(y); // beta, velocity

  // for eta defined cuts
  double c_low_eta_cms = Utils::costh_LT ( b, g, c_low_eta_lab);
  double c_up_eta_cms = Utils::costh_LT ( b, g, c_up_eta_lab);

  // for pt defined ones
  double c_low_pt_cms = - sqrt(1.-4.*pow(el_pt_cut/m,2));
  double c_up_pt_cms  = sqrt(1.-4.*pow(el_pt_cut/m,2));
  double c_low_pt_lab = Utils::costh_LT ( -b,g, c_low_pt_cms);
  double c_up_pt_lab = Utils::costh_LT ( -b,g, c_up_pt_cms);

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

// calculates number of mass integration points
void set_m_bins(){
  double m_low(mr[0]), m_up(mr[1]);
  if ( m_low < 2.*el_pt_cut ) m_low = 2*el_pt_cut;

  double d(0.1);
  double a(m_low), b(m_low+d);
  nmb=0;
  while (1) {
    if ( a<76.) d = 5.;
    else if ( a>=76. && a<89. ) d = 1.;
    else if ( a>=89. && a<94. ) d = 0.1;
    else if ( a>=94. && a<101. ) d = 1.;
    else if ( a>=101. ) d = 5.*exp(a/2000.);
    
    nmb++;
    if ( b>=m_up ) break;

    a= b;
    b+= d;
  }

  mbins = new double[nmb+1];
  for (int imb=0; imb<nmb+1; imb++) mbins[imb] = 0.;

  int imb(0);
  d = 0.1; a = m_low; b= m_low+d;
  while (1) {
    if ( a<76.) d = 5.;
    else if ( a>=76. && a<89. ) d = 1.;
    else if ( a>=89. && a<94. ) d = 0.1;
    else if ( a>=94. && a<101. ) d = 1.;
    else if ( a>=101. ) d = 5.*exp(a/2000.);
    
    mbins[imb]=a;
    if ( b>=m_up ) {
      mbins[imb+1] = m_up;
      break;
    }

    a= b;
    b+= d;
    imb++;
  }
  //for (int imb=0; imb<nmb+1; imb++) cout << imb << " " << mbins[imb] << endl;

}

// calculates number of rapidity integration points
void set_y_bins(){
  double y_low(yr[0]), y_up(yr[1]);

  if ( y_low > y_up ) {
    cout << "z_integ_ybins: rapidity range is wrong: " << y_low << "\t" << y_up << endl;
    exit(1);
  }
  nyb = 1;
  double d = (y_up-y_low)/nyb;
  ybins = new double[nyb+1];
  for (int iyb=0; iyb<nyb+1; iyb++) {
    ybins[iyb]=y_low+d*iyb;
  }
  //for (int iyb=0; iyb<nyb+1; iyb++) cout << iyb << " " << ybins[iyb] << endl;

}


// fill the bin matrix
void fill_bin_matrix(){
  int nmp=2*nmb+1, nyp=2*nyb+1;
  for (int imp=0; imp<nmp; imp++){
    double m = (mbins[imp/2]+mbins[(imp+1)/2])/2.;
    double m2 = m*m;
    alpha = 1./M_PI*gfermi*sqrt(2.)*sw2*cw2*m2z;
    const double denZ = pow(m2-m2z,2)+pow(wz*mz,2);
    double X[3] = {0.};
    X[0] = conhc*M_PI*pow(alpha,2)/(3.*m*pow(ebeam,2));
    X[1] = X[0] * (-2.)*m2*(m2-m2z)/(sw2*cw2*denZ);
    X[2] = X[0] * 3.*pow(m2,2)/(sw2*cw2*denZ*3.*sw2*cw2);
    for(int iyp=0;iyp<nyp;iyp++){
      double y = (ybins[iyp/2]+ybins[(iyp+1)/2])/2.;
      for(int ieb=0;ieb<neb;ieb++){
        for(int icdf = 0; icdf<12; icdf++){
	  // first 4 are for photon component
	  // second 4 are for interference component
	  // third 4 are for Z component
	  int idf = icdf%4;
	  int flav = idf%2; // 0 -> 0 (d), 1 -> 1(u), 2 -> 0(d), 3 -> 1(u)
	  double dir = pow(-1,idf/2); // 0, 1 -> 1 (forw) ; 2,3 -> -1 (bw);
	  double AI = ai_c(m,y,icdf, cbins[ieb], cbins[ieb+1]);
	  BM[imp][iyp][ieb][icdf] = X[icdf/4]*AI;
	}
      }
    }
  }
}
} // IntZ_YBC namespace


// set cos binning corresponding to eta one
// and fill the bins matrix
int z_set_ybins__(void *mass_range, void *y_range, double *pel_pt_cut, int *pneb, void *eb ){
  using namespace IntZ_YBC;
  neb = *pneb;
  
  cbins = new double[neb+1];
  cb_sig_em = new double[neb];
  for (int ib=0;ib<neb+1;ib++){
    cbins[ib] = tanh(*((double*)eb+ib));
  }
  mr[0] = *((double*)mass_range);
  mr[1] = *((double*)mass_range+1);
  yr[0] = *((double*)y_range);
  yr[1] = *((double*)y_range+1);
  el_pt_cut = *pel_pt_cut;

  for (int ib=0; ib<neb; ib++){
    cb_sig_em[ib] = 0.;
  }
  
  set_m_bins();
  set_y_bins();
  BM = new double***[2*nmb+1];
  for (int imb=0;imb<2*nmb+1;imb++){
    BM[imb] = new double**[2*nyb+1];
    for (int iyb=0;iyb<2*nyb+1;iyb++){
      BM[imb][iyb] = new double*[neb];
      for (int ieb=0; ieb<neb; ieb++){
        BM[imb][iyb][ieb] = new double[12];
        for (int icdf = 0; icdf<12; icdf++){  // loop on component, direction, flavour
          BM[imb][iyb][ieb][icdf] = 0.;
        }
      }
    }
  }

  fill_bin_matrix();
  return 1;
}

int z_get_ybins_xs__(void *xs){
//clock_t start, finish;
//start = clock();

  double ypt_pars[2];
  ypt_pars[0] = IntZ_YBC::yr[0]; ypt_pars[1] = IntZ_YBC::yr[1];

  for (int ib=0; ib<IntZ_YBC::neb; ib++){
    IntZ_YBC::cb_sig_em[ib] = 0.;
  }
  
  int integ_res = 
    IntZ_YBC::integrate_M(ypt_pars, IntZ_YBC::mr[0], IntZ_YBC::mr[1]);

//finish = clock();
//cout << double(finish - start)/CLOCKS_PER_SEC<< endl;

  for (int ib=0;ib<IntZ_YBC::neb;ib++){
    *((double*)xs+ib) = IntZ_YBC::cb_sig_em[ib];
  }

  return 1;
}

int z_free_ybins__(){
  using namespace IntZ_YBC;
  for (int imb=0;imb<2*nmb+1;imb++){
    for (int iyb=0;iyb<2*nyb+1;iyb++){
      for (int ieb=0; ieb<neb; ieb++){
        delete[] BM[imb][iyb][ieb] ;
      }
      delete[] BM[imb][iyb] ;
    }
    delete[] BM[imb] ;
  }
  delete[] BM;

  return 1;
}
