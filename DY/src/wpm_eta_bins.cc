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

typedef double*** double3d;

using namespace std;
using namespace PhysPar;

extern "C"{
void fpdfxq_(int *iset, const double *x, const double *q2, double *pdfs, int *ichk);
}

extern "C" {
  int w_set_etabins__(void *mass_range, void *y_range, double *pel_pt_cut, 
      int *pneb, void *eb );
  int w_get_etabins_xs__(void *bs_wm, void *bs_wp);
  int w_free_etabins__();
}

namespace IntW_EB
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
double *cb_sig_wm;
double *cb_sig_wp;
double mr[2], yr[2], el_pt_cut;
double3d BMF, BMB; // bin matrix to store the coefficients for the bins
                   // for up-quark coming from forward and backward protons

int integrate_M(double *pars, double m_low, double m_up){
  // check the pt_cut
  if ( m_up < 2.*el_pt_cut ) {
    for (int ib = 0; ib<neb; ib++)  {
      cb_sig_wm[ib]=0.;
      cb_sig_wp[ib]=0.;
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
    //cout << "simpson   " << a << " " << b-a << " " << cb_sig_wm[0] << endl;
    iy(2*(imb+1),pars,iyb);
    iy(2*imb+1, pars, iym);
    for (int ib = 0; ib<neb; ib++){
      cb_sig_wm[ib]+=(b-a)/6.*(iya[ib]+4.*iym[ib]+iyb[ib]);
      cb_sig_wp[ib]+=(b-a)/6.*(iya[neb+ib]+4.*iym[neb+ib]+iyb[neb+ib]);
      iya[ib]=iyb[ib];
      iya[neb+ib]=iyb[neb+ib];
    }
  }

  return 1;
}

int get_pdfconv(double &m, double y, double dir, const double &scale, double *Vxf)
{
  double x1 = m/ebeam*exp(dir*y);
  double x2 = m/ebeam*exp(-dir*y);

  double xq1[13]={0.};
  double xq2[13]={0.};
  int   ichk = 0;
  int   iset = 1;

  double q2 = (scale) * (scale);

  if ( (x1< 1.e-6)||(x2<1.e-6)||(x1>1.-1.e-6)||(x2>1.-1.e-6)) {
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
  Vxf[0] = (V2[0]*xq1[4]*xq2[7]
          + V2[1]*xq1[2]*xq2[9]
          + V2[2]*xq1[4]*xq2[9]
          + V2[3]*xq1[2]*xq2[7])/(x1*x2)*2.;

  Vxf[1] = (V2[0]*xq1[8] *xq2[5]
          + V2[1]*xq1[10]*xq2[3]
          + V2[2]*xq1[8] *xq2[3]
          + V2[3]*xq1[10]*xq2[5])/(x1*x2)*2.;
  
  return 1;
}

int iy(int imp, double *pars, double *iyval){
  double y_low     = *(double*) pars;
  double y_up      = *((double*)pars+1);
  double m = (mbins[imp/2]+mbins[(imp+1)/2])/2.;

  const double scale = mw;
  double ica[2*neb], icb[2*neb], icm[2*neb];
  for (int ib=0;ib<2*neb;ib++){
    ica[ib]=0.; icb[ib]=0.; icm[ib]=0.;
  }
  double pcaF[2]={0., 0.}; // pdf convolution
  double pcaB[2]={0., 0.}; // pdf convolution
  get_pdfconv(m, y_low, 1., scale, pcaF);
  get_pdfconv(m, y_low, -1., scale, pcaB);
  for (int ib=0;ib<neb;ib++){
    ica[ib] = pcaF[0]*BMF[imp][0][ib]+pcaB[0]*BMB[imp][0][ib];
    ica[neb+ib] = pcaF[1]*BMF[imp][0][ib]+pcaB[1]*BMB[imp][0][ib];
    iyval[ib] = 0.;
    iyval[neb+ib] = 0.;
  }
  //exit(1);

  for (int iyb=0;iyb<nyb;iyb++){
    double ya = ybins[iyb];
    double yb = ybins[iyb+1];
    double ym = (ya+yb)/2.;
    double pcbF[2]={0.,0.}, pcmF[2]={0.,0.};
    double pcbB[2]={0.,0.}, pcmB[2]={0.,0.};
    get_pdfconv(m, yb, 1., scale, pcbF);
    get_pdfconv(m, ym, 1., scale, pcmF);
    get_pdfconv(m, yb, -1.,scale, pcbB);
    get_pdfconv(m, ym, -1.,scale, pcmB);
    for(int ieb=0;ieb<neb;ieb++){
      icb[ieb]     = pcbF[0]*BMF[imp][2*(iyb+1)][ieb] + pcbB[0]*BMB[imp][2*(iyb+1)][ieb];
      icb[neb+ieb] = pcbF[1]*BMF[imp][2*(iyb+1)][ieb] + pcbB[1]*BMB[imp][2*(iyb+1)][ieb];
      icm[ieb]     = pcmF[0]*BMF[imp][2*iyb+1][ieb]   + pcmB[0]*BMB[imp][2*iyb+1][ieb];
      icm[neb+ieb] = pcmF[1]*BMF[imp][2*iyb+1][ieb]   + pcmB[1]*BMB[imp][2*iyb+1][ieb];
      iyval[ieb] += (yb-ya)/6.*(ica[ieb]+4*icm[ieb]+icb[ieb]);
      iyval[neb+ieb] += (yb-ya)/6.*(ica[neb+ieb]+4*icm[neb+ieb]+icb[neb+ieb]);
      ica[ieb] = icb[ieb];
      ica[neb+ieb] = icb[neb+ieb];
    }
  }
      
  return 1;
}

double ai_c(const double &m, const double &y, 
            const double &c_low_eta_lab, const double &c_up_eta_lab){

  // calculate actual costh range applying lorentz transformation
  // to LAB cos_theta cuts
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

// calculates number of mass integration points
void set_m_bins(){
  double m_low(mr[0]), m_up(mr[1]);
  if ( m_low < 2.*el_pt_cut ) m_low = 2*el_pt_cut;

  double d(0.1);
  double a(m_low), b(m_low+d);
  nmb=0;
  while (1) {
    if ( a<65.) d = 5.;
    else if ( a>=65. && a<78. ) d = 1.;
    else if ( a>=78. && a<83. ) d = 0.1;
    else if ( a>=83. && a<90. ) d = 1.;
    else if ( a>=90. ) d = 5.*exp(a/2000.);
    
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
    if ( a<65.) d = 5.;
    else if ( a>=65. && a<78. ) d = 1.;
    else if ( a>=78. && a<83. ) d = 0.1;
    else if ( a>=83. && a<90. ) d = 1.;
    else if ( a>=90. ) d = 5.*exp(a/2000.);
    
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

  double d(0.01);
  double a(y_low), b(y_low+d);
  nyb=0;
  while (1) {
  
    if ( a<-6.) d = 0.5;
    else if ( a>=-6. && a<-2.1 ) d = 0.2;
    else if ( a>=-2.1 && a<2.1 ) d = 0.2;
    else if ( a>=2.1 && a<6. ) d = 0.2;
    else if ( a>=6. ) d = 0.5;
  
    nyb++;
    if ( b>=y_up ) break;

    a= b;
    b+= d;
  }

  ybins = new double[nyb+1];
  for (int iyb=0; iyb<nyb+1; iyb++) ybins[iyb] = 0.;

  int iyb(0);
  d = 0.01; a = y_low; b= y_low+d;
  while (1) {
    if ( a<-6.) d = 0.5;
    else if ( a>=-6. && a<-2.1 ) d = 0.2;
    else if ( a>=-2.1 && a<2.1 ) d = 0.2;
    else if ( a>=2.1 && a<6. ) d = 0.2;
    else if ( a>=6. ) d = 0.5;
  
    ybins[iyb]=a;
    if ( b>=y_up ) {
      ybins[iyb+1] = y_up;
      break;
    }

    a= b;
    b+= d;
    iyb++;
  }
 // for (int iyb=0; iyb<nyb+1; iyb++) cout << iyb << " " << ybins[iyb] << endl;

}


// fill the bin matrix
void fill_bin_matrix(){
  int nmp=2*nmb+1, nyp=2*nyb+1;
  for (int imp=0; imp<nmp; imp++){
    double m = (mbins[imp/2]+mbins[(imp+1)/2])/2.;
    double m2 = m*m;
    double X = conhc * pow(gfermi*sqrt(2.)*cw2*m2z/M_PI/ebeam,2)*m*m2*M_PI/ 48.
                / (pow(m2-m2w,2)+pow(ww*mw,2));
    for(int iyp=0;iyp<nyp;iyp++){
      double y = (ybins[iyp/2]+ybins[(iyp+1)/2])/2.;
      for(int ieb=0;ieb<neb;ieb++){
        BMF[imp][iyp][ieb] = X*ai_c(m,y, cbins[ieb], cbins[ieb+1])/2.;
        BMB[imp][iyp][ieb] = X*ai_c(m,-y, -cbins[ieb+1], -cbins[ieb])/2.;
      }
    }
  }
}
} // IntW_EB namespace


// set cos binning corresponding to eta one
// and fill the bins matrix
int w_set_etabins__( void *mass_range, void *y_range, double *pel_pt_cut,
                     int *pneb, void *eb){
  using namespace IntW_EB;
  neb = *pneb;
  
  cbins = new double[neb+1];
  cb_sig_wm = new double[neb];
  cb_sig_wp = new double[neb];
  for (int ib=0;ib<neb+1;ib++){
    cbins[ib] = tanh(*((double*)eb+ib));
  }
  mr[0] = *((double*)mass_range);
  mr[1] = *((double*)mass_range+1);
  yr[0] = *((double*)y_range);
  yr[1] = *((double*)y_range+1);
  el_pt_cut = *pel_pt_cut;

  for (int ib=0; ib<neb; ib++){
    cb_sig_wm[ib] = 0.;
    cb_sig_wp[ib] = 0.;
  }
  
  set_m_bins();
  set_y_bins();
  BMF = new double**[2*nmb+1];
  BMB = new double**[2*nmb+1];
  for (int imb=0;imb<2*nmb+1;imb++){
    BMF[imb] = new double*[2*nyb+1];
    BMB[imb] = new double*[2*nyb+1];
    for (int iyb=0;iyb<2*nyb+1;iyb++){
      BMF[imb][iyb] = new double[neb];
      BMB[imb][iyb] = new double[neb];
      for (int ieb=0; ieb<neb; ieb++){
        BMF[imb][iyb][ieb] = 0.;
        BMB[imb][iyb][ieb] = 0.;
      }
    }
  }

  fill_bin_matrix();
  return 1;
}

int w_get_etabins_xs__(void *bs_wm, void *bs_wp){
//clock_t start, finish;
//start = clock();

  double ypt_pars[2];
  ypt_pars[0] = IntW_EB::yr[0]; ypt_pars[1] = IntW_EB::yr[1];

  for (int ib=0; ib<IntW_EB::neb; ib++){
    IntW_EB::cb_sig_wm[ib] = 0.;
    IntW_EB::cb_sig_wp[ib] = 0.;
  }
  
  int integ_res = 
    IntW_EB::integrate_M(ypt_pars, IntW_EB::mr[0], IntW_EB::mr[1]);

//finish = clock();
//cout << double(finish - start)/CLOCKS_PER_SEC<< endl;

  for (int ib=0;ib<IntW_EB::neb;ib++){
    *((double*)bs_wm+ib) = IntW_EB::cb_sig_wm[ib];
    *((double*)bs_wp+ib) = IntW_EB::cb_sig_wp[ib];
  }

  return 1;
}

int w_free_etabins__(){
  using namespace IntW_EB;
  for (int imb=0;imb<2*nmb+1;imb++){
    for (int iyb=0;iyb<2*nyb+1;iyb++){
      delete[] BMF[imb][iyb];
      delete[] BMB[imb][iyb];
    }
    delete[] BMF[imb];
    delete[] BMB[imb];
  }
  delete[] BMF, BMB;
}
