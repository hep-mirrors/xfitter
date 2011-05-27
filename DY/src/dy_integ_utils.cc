#include "dy_integ_utils.h"
#include "iostream"
#include "math.h"
#include "PhysPar.h"

using namespace PhysPar;

/**
 * External config routines
 **/
extern "C" {
  int set_beams__(int *p1, int *p2, double *cme);
  int set_lepton__(int *leptid);
}

int set_beams__(int *ip1, int *ip2, double *cme){
  ih1 = *ip1;
  ih2 = *ip2;

  ebeam = *cme;
  return 1;
}

int set_lepton__(int *leptid){
  idlept = *leptid;
}

/**
 * namespace for internal usage, mostly not used.
 */
namespace Utils
{
// calculates numerical integral using simpson rule
int simpson( double (*F)(double, void*), void *pars, const double xmin, 
            const double xmax, const int N, double &integ, double &integ_err){

  double d = (xmax-xmin)/N;
  integ=0;
  double rel_sum2(0);
  double Fa = F(xmin,pars), Fb(0.);
//    std::cout << "simpson " << Fa <<std::endl;
//    exit(1);
  for (int ib=0; ib<N; ib++){
    double a = xmin + d*ib;
    double b = xmin + d*(ib+1);

    Fb = F(b,pars);
    integ+=(b-a)/6.*(Fa+4.*F((a+b)/2.,pars)+Fb);
    Fa = Fb;
    /*
    rel_sum2+=pow((6.*(b-a)/81.*(F(a,pars)-2.*F((a+b)/2.,pars)+F(b,pars)))
             /((b-a)/6.*(F(a,pars)+4.*F((a+b)/2.,pars)+F(b,pars))),2);
	     */
  }
  //integ_err = integ*sqrt(rel_sum2);

  return 1;
}

// simpson with binning adjusted for W-mass integration
int simpsonW( double (*F)(double, void*), void *pars, const double xmin, 
            const double xmax, const int N, double &integ, double &integ_err){

  integ=0;
  double d(0.1);
  double a(xmin), b(xmin+d);
  double Fa = F(xmin,pars), Fb(0.);
  int n(0);
  while (1) {
  //break;
    //std::cout << "simpson   " << a << " " << d << " " <<  integ <<std::endl;
    if ( a<65.) d = 5.;
    else if ( a>=65. && a<78. ) d = 1.;
    else if ( a>=78. && a<83. ) d = 0.1;
    else if ( a>=83. && a<90. ) d = 1.;
    else if ( a>=90. ) d = 5.*exp(a/2000.);
    
    if ( b>=xmax ){
      integ+=(xmax-a)/6.*(Fa+4.*F((a+xmax)/2.,pars)+F(xmax,pars));
      break;
    }

    Fb = F(b,pars);
    integ+=(b-a)/6.*(Fa+4.*F((a+b)/2.,pars)+Fb);
    Fa = Fb;

    a= b;
    b+= d;
  }

  return 1;
}

// simpson with binning adjusted for Z-mass integration
int simpsonZ( double (*F)(double, void*), void *pars, const double xmin, 
            const double xmax, const int N, double &integ, double &integ_err){

  integ=0;
  double d(0.1);
  double a(xmin), b(xmin+d);
  double Fa = F(xmin,pars), Fb(0.);
  while (1) {
    if ( a<75.) d = 5.;
    else if ( a>=75. && a<88. ) d = 1.;
    else if ( a>=88. && a<94. ) d = 0.1;
    else if ( a>=94. && a<100. ) d = 1.;
    else if ( a>=100. ) d = 5.*exp(a/2000.);
    
    if ( b>=xmax ){
      integ+=(xmax-a)/6.*(Fa+4.*F((a+xmax)/2.,pars)+F(xmax,pars));
      break;
    }

    Fb = F(b,pars);
    integ+=(b-a)/6.*(Fa+4.*F((a+b)/2.,pars)+Fb);
    Fa = Fb;

    a= b;
    b+= d;
  }

  return 1;
}

// transforms cosine theta to a frame defined by b and g
// arguments are betta, gamma and cosine theta
double costh_LT(const double &b, const double &g, const double &c0){
  double c = (c0-b)/(1-b*c0);
  return c;
}

} // namespace utils

