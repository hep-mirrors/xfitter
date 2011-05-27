#include <algorithm>
#include <vector>
#include <math.h>
#include <iostream>
#include "PhysPar.h"
#include <gsl/gsl_integration.h>

#define MAX_CHAR 1000
#define NALLOC   100000
//#define PRECISION	1.e-3

using namespace std;
using namespace PhysPar;

extern "C"{
void fpdfxq_(int *iset, const double *x, const double *q2, double *pdfs, int *ichk);
}

extern "C" {
  int wmenu_integ__(double *integral, double *integral_err, void *mass_range, 
                    void *y_range, void *eta_range, double *el_pt_cut);
}

extern "C" {
  void dywm_check__();
}

static double CurrentY = 0;

namespace IntegDY_Wm
{

double intgd_m(double , void *);
//double intgd_y(double , void *);
double dsig_Wm(double , void *);
double get_PhiWm(const double &, const double &, const double &);
double ai_c(const double &,const double &,const double &,
            const double &, const double &);
double costh_LT(const double &, const double &, const double &);

double INT_PREC(1.e-3);

double intgd_m(double m, void *ycpt_pars){
//cout << "m: " << m << endl;
  if ( m < 2.* (*((double*)ycpt_pars+4)) ) return 0.;
  gsl_integration_workspace *wsp = 
    gsl_integration_workspace_alloc(NALLOC);

  // prepare y-range for integration
  double y_low = *(double*)ycpt_pars;
  double y_up  = *((double*)ycpt_pars+1);
  //cout << y_low << "\t" << y_up << endl;
  double mcpt_pars[4];
  mcpt_pars[0] = m;
  mcpt_pars[1] = *((double*)ycpt_pars+2); // costh low
  mcpt_pars[2] = *((double*)ycpt_pars+3); // costh up
  mcpt_pars[3] = *((double*)ycpt_pars+4); // electron pt cut
  // prepare gsl_function
  double inty, inty_err;
  gsl_function F;
  F.function = &dsig_Wm;
  F.params = &mcpt_pars;

  gsl_integration_qag(&F, y_low, y_up, 0, INT_PREC, NALLOC, 
     GSL_INTEG_GAUSS41, wsp, &inty, &inty_err);

  gsl_integration_workspace_free(wsp);
 //cout << "gnuplot " << m << "\t" << inty << endl;
 //exit(1);
  
  return inty;
}

double get_PhiWm(const double &x1, const double &x2, const double &scale)
{
  double xq1[13];
  double xq2[13];
  int   ichk = 0;
  int   iset = 1;

  double q2 = (scale) * (scale);

  if ( (x1< 1.e-6)||(x2<1.e-6)||(x1>1.-1.e-6)||(x2>1.-1.e-6)) {
    return 0.;
  }
  else{
  // Access QCDNUM pdfs:
    fpdfxq_(&iset,&x1,&q2,xq1,&ichk);
    fpdfxq_(&iset,&x2,&q2,xq2,&ichk);
  }


/*
cout << x1<< "\t" << x2 << endl;
for (int i = 0 ; i < 13; i++){
cout << xq1[i] << "  " ;
}
cout << endl;
for (int i = 0 ; i < 13; i++){
cout << xq2[i] << "  " ;
}
cout << endl;
*/


  // numbering:
  // \bar{t  b  c  s  u  d} g  d  u  s  c  b  t
  //      0  1  2  3  4  5  6  7  8  9 10 11 12
  double ret = (V2[0]*xq1[4]*xq2[7]
              + V2[1]*xq1[2]*xq2[9]
              + V2[2]*xq1[4]*xq2[9]
              + V2[3]*xq1[2]*xq2[7])/(x1*x2)*2.;

  /*
  for (int i=0; i<=12; i++){
    cout <<  pdfs1[i] << " ";
  }
  cout << "\n";
  cout << "IFlav = " << iflav << endl;
  cout << "Check: "<< xq1 <<  "\t" << xq2 << "\t" << xqbar1<<"\t"<<xqbar2<<endl;
  */
  return ret;
}

double dsig_Wm(double y, void *mcpt_pars){
  double m         = *(double*) mcpt_pars;
  double c_low_eta_lab = *((double*)mcpt_pars+1); // cos_theta ranges in LAB frame
  double c_up_eta_lab  = *((double*)mcpt_pars+2);
  double el_pt_cut = *((double*)mcpt_pars+3);
//cout << m << "\t" << y <<"\t" <<  endl;
  double m2 = m*m;
  double x1 = m/ebeam*exp(y);
  double x2 = m/ebeam*exp(-y);
  if ( (x1< 1.e-6)||(x2<1.e-6)||(x1>1.-1.e-6)||(x2>1.-1.e-6)) {
    return 0.;
  }
//x1=0.1;x2=0.1;
  //cout << x1 << "\t" << x2 << endl;


  double phi_Wm(0.);
  phi_Wm = get_PhiWm(x1,x2,mw);

  // in g-fermi scheme, as in MCFM
  double ds_dmdydcosth = conhc
              * pow(gfermi*sqrt(2.)*cw2*m2z/M_PI/ebeam,2)*m*m2*M_PI/ 48.
              * phi_Wm / (pow(m2-m2w,2)+pow(ww*mw,2))
	      * ( ai_c(m,y, el_pt_cut, c_low_eta_lab, c_up_eta_lab) 
	        + ai_c(m,y, el_pt_cut, -c_up_eta_lab, -c_low_eta_lab) )/2.;

  //cout << "gnuplot  "  << pt << "\t" << dsig << endl;
  //cout << "gnuplot  " << m << " " << y << " " << pt << "\t" << phi_Wm << endl;
  //cout << "gnuplot  " << m << " " << y << " " << pt << "\t" << ds_dmdydcosth << endl;
  return ds_dmdydcosth;
}

double ai_c(const double &m, const double &y, const double &el_pt_cut, 
            const double &c_low_eta_lab, const double &c_up_eta_lab){

  // calculate actual costh range applying lorentz transformation
  // to LAB cos_theta cuts
  double g = cosh(y); // gamma, lorentz factor
  double b = tanh(y); // beta, velocity

  // for eta defined cuts
  double c_low_eta_cms = costh_LT ( b, g, c_low_eta_lab);
  double c_up_eta_cms = costh_LT ( b, g, c_up_eta_lab);
  //cout << costh_LT ( b, g, c_low_eta_lab)<< "\t" << costh_LT ( b, g, c_up_eta_lab) << endl ;
  //exit(1);

  // for pt defined ones
  double c_low_pt_cms = - sqrt(1.-4.*pow(el_pt_cut/m,2));
  double c_up_pt_cms  = sqrt(1.-4.*pow(el_pt_cut/m,2));
  double c_low_pt_lab = costh_LT ( -b,g, c_low_pt_cms);
  double c_up_pt_lab = costh_LT ( -b,g, c_up_pt_cms);
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
// arguments are betta, gamma and cosine theta
double costh_LT(const double &b, const double &g, const double &c0){
  double c = (c0-b)/(1-b*c0);
  return c;
}
} // IntegDY_Wm namespace

int wmenu_integ__( double *integral, double *integral_err, void *mass_range, 
                   void *y_range, void *eta_range, double *el_pt_cut){

  //int proc_id(0);
  double m_low = (*(double*)mass_range);
  double m_up  = (*((double*)mass_range+1));
  double y_low = (*(double*)y_range);
  double y_up  = (*((double*)y_range+1));
  double c_low = tanh(*(double*)eta_range);
  double c_up  = tanh(*((double*)eta_range+1));

  //  cout << m_low << "\t" << m_up <<endl;
  //  cout << y_low << "\t" << y_up <<endl;
  //  cout << c_low << "\t" << c_up <<endl;

  gsl_integration_workspace *wsp = 
    gsl_integration_workspace_alloc(NALLOC);
  
  double ycpt_pars[5];
  ycpt_pars[0] = y_low; ycpt_pars[1] = y_up;
  ycpt_pars[2] = c_low; ycpt_pars[3] = c_up;
  ycpt_pars[4] = *el_pt_cut;

  //triple integral, integral_err;
  gsl_function gF;
  gF.function = &(IntegDY_Wm::intgd_m);
  gF.params = &ycpt_pars;
  
  int integ_res = 
    gsl_integration_qags(&gF, m_low, m_up, 0, IntegDY_Wm::INT_PREC, 
    NALLOC, wsp, integral, integral_err);

  gsl_integration_workspace_free(wsp);

  return integ_res;
}

void dywm_check__() {
  double mass = 80.;
  for ( int i = 0; i<=0; i++){
    mass += i / 10.;
    double val = IntegDY_Wm::dsig_Wm(0.0 , &mass);
    cout << "mass="<<mass<<" Double diff = " << val << "\n";
  }

  double integral(0); double err(0);
  double mass_range[2]; double y_range[2], eta_range[2];
  mass_range[0] = 1.; mass_range[1] = 7000.;
  y_range[0] = -10.; y_range[1] = 10.;
  eta_range[0] = -10.; eta_range[1] = 10.;
  double pt_cut = 0.;

  int res = wmenu_integ__(&integral,&err,mass_range,y_range,eta_range,&pt_cut);
  /*
  double yrange[2]; 
  yrange[0] = -10.;
  yrange[1] =  10.;
  double amint = intgd_mm(mass,yrange);
  cout << " mass int=" << amint << "\n";  
*/
  cout << " Integral ="<< integral <<" +- "<<err<<"\n";
}

