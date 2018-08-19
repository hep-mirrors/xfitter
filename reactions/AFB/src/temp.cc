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
  int zee_integ__(double *integral, double *integral_err, void *mass_range, void *y_range);
}

extern "C" {
  void dy_check__();
}

static double CurrentY = 0;

namespace IntegDY_Z
{

double intgd_m(double , void *);
double dsig_gammaZ(double , void *);
double get_Fqqbar(const double &, const double &, const double &, const int );

double INT_PREC(1.e-3);

double intgd_m(double m, void *y_range){
//cout << "m: " << m << endl;
  gsl_integration_workspace *wsp = 
    gsl_integration_workspace_alloc(NALLOC);

  // prepare y-range for integration
  double y_low = (*(double*)y_range);
  double y_up = (*((double*)y_range+1));
  //  cout << "here" << y_low << "\t" << y_up << endl;
  // prepare gsl_function
  double inty, inty_err;
  gsl_function F;
  F.function = &dsig_gammaZ;
  F.params = &m;

  gsl_integration_qags(&F, y_low, y_up, 0, INT_PREC, NALLOC, wsp, &inty, &inty_err);

  gsl_integration_workspace_free(wsp);
//  cout << "gnuplot " << m << "\t" << inty << endl;
  
  return inty;
}

double get_Fqqbar(const double &x1, const double &x2, const double &scale, const int iflav, const int iget)
{
//  short int bar_sign = ( iflav>2 )? 1: -1;
//  if ( iflav>2 ) bar_sign = 1;
  static double pdfs1[13];
  static double pdfs2[13];
  int   ichk = 0;
  int   iset = 1;

  double q2 = (scale) * (scale);


  if ( (x1< 1.e-6)||(x2<1.e-6)||(x1>1.-1.e-6)||(x2>1.-1.e-6)) {
    return 0.;
  }
  else{
  // Access QCDNUM pdfs:
    if (iget == 0) {
      fpdfxq_(&iset,&x1,&q2,pdfs1,&ichk);
      fpdfxq_(&iset,&x2,&q2,pdfs2,&ichk);
    }
  }

  int ifl       = 6+iflav  ;
  int iflbar    = 6-iflav ;

  double xq1    = pdfs1[ifl];  //LHAPDF::xfx(x1, scale, iflav);
  double xqbar2 = pdfs2[iflbar]; //LHAPDF::xfx(x2, scale, bar_sign*iflav);
  double xq2    = pdfs2[ifl]; //lLHAPDF::xfx(x2, scale, iflav);
  double xqbar1 = pdfs1[iflbar]; //LHAPDF::xfx(x1, scale, bar_sign*iflav);
  double ret = (xq1*xqbar2+xqbar1*xq2)/(x1*x2);

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

// http://cdsweb.cern.ch/record/1232420/files/CERN-THESIS-2010-006.pdf version
double dsig_gammaZ(double y, void *pm)
{
  const double m = *(double*) pm;
//cout << m << "\t" << y << endl;
  const double m2 = m*m;
  const double x1 = m/ebeam*exp(y);
  const double x2 = m/ebeam*exp(-y);
  if ( (x1< 1.e-6)||(x2<1.e-6)||(x1>1.-1.e-6)||(x2>1.-1.e-6)) {
    return 0.;
  }
  //cout << x1 << "\t" << x2 << endl;

  double phi_gamma(0.), phi_Z(0.), phi_intf(0.);
  for ( int iflav = 0; iflav<5; iflav++ ) {
    double Fqqbar = get_Fqqbar(x1,x2,mz, iflav+1,iflav);  // iflav = 0: read new pdfs from QCDNUM

    phi_gamma += pow(chq[iflav%2],2) * Fqqbar;
    phi_intf  += vq[iflav%2]*chq[iflav%2]*Fqqbar;
    phi_Z     += (pow(vq[iflav%2],2)+1.) * Fqqbar;
  }

  //const double kZ = 1./4./sw2/cw2;
  alpha = 1./M_PI*gfermi*sqrt(2.)*sw2*cw2*m2z;
  //  cout << 1./alpha << endl;
  const double denZ = pow(m2-m2z,2)+pow(wz*mz,2);
  double dsig = conhc*8.*M_PI*pow(alpha,2)/(9.*m*pow(ebeam,2))*(
    phi_gamma 
  + m2*(m2-m2z)*(1.-4.*sw2)/(8.*sw2*cw2*denZ)*phi_intf
  + 3*pow(m2,2)*(1.+pow(1.-4.*sw2,2))/(16.*sw2*cw2*denZ*48.*sw2*cw2)*phi_Z
  );

  return dsig;
}

} // IntegDY_Z namespace

int zee_integ__( double *integral, double *integral_err, void *mass_range, void *y_range){

  //int proc_id(0);
  double m_low = (*(double*)mass_range);
  double m_up  = (*((double*)mass_range+1));
  double y_low = (*(double*)y_range);
  double y_up  = (*((double*)y_range+1));

  //  cout << m_low << "\t" << m_up <<endl;
  //  cout << y_low << "\t" << y_up <<endl;

  gsl_integration_workspace *wsp = 
    gsl_integration_workspace_alloc(NALLOC);
  

  //double integral, integral_err;
  gsl_function gF;
  gF.function = &(IntegDY_Z::intgd_m);
  gF.params = y_range;
  
  int integ_res = 
    gsl_integration_qags(&gF, m_low, m_up, 0, IntegDY_Z::INT_PREC, NALLOC, wsp, integral, integral_err);

  gsl_integration_workspace_free(wsp);

  return integ_res;
}

void dy_check__() {
  double mass = 90.;
  for ( int i = 0; i<=0; i++){
    mass += i / 10.;
    double val = IntegDY_Z::dsig_gammaZ(0.0 , &mass);
    cout << "mass="<<mass<<" Double diff = " << val << "\n";
  }

  double integral(0); double err(0);
  double mass_range[2]; double y_range[2];
  mass_range[0] = 50.; mass_range[1] = 110.;
  y_range[0] = -10.; y_range[1] = 10.;

  int res = zee_integ__(&integral,&err,mass_range,y_range);
  /*
  double yrange[2]; 
  yrange[0] = -10.;
  yrange[1] =  10.;
  double amint = intgd_m(mass,yrange);
  cout << " mass int=" << amint << "\n";  
*/
  cout << " Integral ="<< integral <<" +- "<<err<<"\n";
}
