/**
* @file konst.h
* @author dwilhelm
* @date	26.12.2012
* @brief Definition of constants*
*/

#ifndef KONST_h
#define KONST_h

//#include <config.h>
#include <math.h>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <cstdlib>
typedef unsigned int  uint;

using std::cout;
using std::endl;



//! @brief Print GSL_ERROR.
//! @param status return status of a gsl function
using namespace std;
inline void GSL_ERROR_CODE(int status){
	switch(status)
	{
	case GSL_SUCCESS:	break;
	case GSL_CONTINUE: 	break;
	case GSL_EBADFUNC: cout << "GSL_EBADFUNC" << endl;	break;
	case GSL_EBADLEN: cout << "GSL_EBADLEN" << endl;	break;
	case GSL_EBADTOL: cout << "GSL_EBADTOL" << endl;	break;
	case GSL_ECACHE: cout << "GSL_ECACHE" << endl;	break;
	case GSL_EDOM: cout << "GSL_EDOM" << endl;	break;
	case GSL_EFACTOR: cout << "GSL_EFACTOR" << endl;	break;
	case GSL_EFAILED: cout << "GSL_EFAILED" << endl;	break;
	case GSL_EFAULT: cout << "GSL_EFAULT" << endl;	break;
	case GSL_EINVAL: cout << "GSL_EINVAL" << endl;	break;
	case GSL_ELOSS: cout << "GSL_ELOSS" << endl;	break;
	case GSL_ENOMEM: cout << "GSL_ENOMEM" << endl;	break;
	case GSL_ENOPROG: cout << "GSL_ENOPROG" << endl;	break;
	case GSL_ENOPROGJ: cout << "GSL_ENOPROGJ" << endl;	break;
	case GSL_ENOTSQR: cout << "GSL_ENOTSQR" << endl;	break;
	case GSL_EOF: cout << "GSL_EOF" << endl;	break;
	case GSL_EOVRFLW: cout << "GSL_EOVRFLW" << endl;	break;
	case GSL_ERANGE: cout << "GSL_ERANGE" << endl;	break;
	case GSL_ERUNAWAY: cout << "GSL_ERUNAWAY" << endl;	break;
	case GSL_ESANITY: cout << "GSL_ESANITY" << endl;	break;
	case GSL_ETABLE: cout << "GSL_ETABLE" << endl;	break;
	case GSL_ETOL: cout << "GSL_ETOL" << endl;	break;
	case GSL_ETOLF: cout << "GSL_ETOLF" << endl;	break;
	case GSL_ETOLG: cout << "GSL_ETOLG" << endl;	break;
	case GSL_ETOLX: cout << "GSL_ETOLX" << endl;	break;
	case GSL_EUNDRFLW: cout << "GSL_EUNDRFLW" << endl;	break;
	case GSL_EUNIMPL: cout << "GSL_EUNIMPL" << endl;	break;
	case GSL_EUNSUP: cout << "GSL_EUNSUP" << endl;	break;
	case GSL_EZERODIV: cout << "GSL_EZERODIV" << endl;	break;
	case GSL_EMAXITER: cout << "the maximum number of subdivisions was exceeded." << endl;	break;
	case GSL_EROUND: cout << "cannot reach tolerance because of roundoff error, or roundoff error was detected in the extrapolation table." << endl;	break;
	case GSL_ESING: cout << "a non-integrable singularity or other bad integrand behavior was found in the integration interval." << endl;	break;
	case GSL_EDIVERGE: cout << "the integral is divergent, or too slowly convergent to be integrated numerically." << endl;	break;
	default: cout << "ERROR_CODE: " << status << endl;	break;
	}
}

#ifndef PI
#define PI (M_PI)
#endif

#ifndef PI2
#define PI2 (M_PI*M_PI)
#endif

#ifndef PI3
#define PI3 (M_PI*M_PI*M_PI)
#endif

#ifndef PI4
#define PI4 (M_PI*M_PI*M_PI*M_PI)
#endif


#ifndef gamma_E
#define gamma_E 0.577215664901532	//!<@def gamma_E Eulergamma \f$\gamma_E\f$
#endif

static const double b_0= 2*exp(-gamma_E);//!< b_0 blabla

//define schemes
#define Gmu_scheme
//Gmu scheme (MCFM default): Input MZ=91.1876 MW=80.398, GF=1.16639e-5 -> Output alpha(MZ)=1/132.3384323, xw=0.2226459
#ifdef MCFM_scheme
static double ALPHA_EM = (1./132.3384323);
static const double ThetaW = 0.491392051938186802;
static const double G_F = 1.16639e-05;
#endif

//Gmu scheme (PDG 2012): Input MZ=91.1876 MW=80.385, GF=1.1663787e-05 -> Output alpha(MZ)=1/132.2332298, xw=0.2228972
#ifdef Gmu_scheme
static double ALPHA_EM = (1./132.2332298);
static const double ThetaW = 0.491694018114958109;
static const double G_F = 1.1663787e-05;
#endif

//alpha(MZ) scheme (PDG 2012): Input MZ=91.1876, xw=0.2312, alpha(MZ)=1/128.89 -> Output MW=79.9544195, GF=1.16612e-05
#ifdef alphamz_scheme
static double ALPHA_EM = (1./128.89);
static const double G_F = 1.16612e-05;
static const double ThetaW = 0.501604054189120463;
#endif

// M_W=80.385 M_Z=91.1876 => ThetaW
//static const double ThetaW	= 0.4916940451746219; //(xw=0.222897)
static const double SThetaW 	= sin(ThetaW);
static const double SThetaW2 	= gsl_pow_2(SThetaW);
static const double CThetaW 	= cos(ThetaW);
static const double CThetaW2 	= gsl_pow_2(CThetaW);

#ifndef T_F
#define T_F 0.5
#endif

#ifndef N_G
#define N_G 8.0
#endif

#ifndef N_c
#define N_c 3.0
#endif

#ifndef T_A
#define T_A 3.0
#endif

#ifndef C_A
#define C_A 3.0
#endif

#ifndef C_F
#define C_F (4.0/3.0)
#endif

#ifndef zeta_3
#define zeta_3 1.2020569031595942
#endif
#ifndef zeta_5
#define zeta_5 1.03692775514337
#endif

//static double Lambda_QCD;
//static double alpha_s_order;
//static double ALPHA_EM = (1./127.916);
//static const double G_F = 1.1620765257672414e-5;

static const double Gamma_0F = (4.0*C_F);
static const double Gamma_0 = (4.0);
static const double Gamma_0A = (4.0*C_A);
static const double Gamma__0 = 1.0/(Gamma_0);
static const double Gamma__0F = 1.0/(Gamma_0F);
static const double Gamma__0A = 1.0/(Gamma_0A);
inline double Gamma_1_0(int const n_f){return (((67.0-3.0*PI2)*C_A-20.0*T_F*n_f)/9.0);}
inline double Gamma_1(int const n_f){return Gamma_1_0(n_f)*Gamma_0;}
inline double Gamma_1F(int const n_f){return Gamma_1_0(n_f)*Gamma_0F;}
inline double Gamma_1A(int const n_f){return Gamma_1_0(n_f)*Gamma_0A;}
inline double Gamma_2_0(int const n_f){return  ((
		C_A*C_A*(245.0/6.0-134.0/27.0*PI2+11.0*PI4/45.0+22.0/3.0*zeta_3)
		+C_A*T_F*n_f*(-418.0/27.0+40.0*PI2/27.0-56.0/3.0*zeta_3)
		+C_F*T_F*n_f*(-55.0/3.0+16.0*zeta_3)
		-16.0/27.0*T_F*T_F*n_f*n_f));}
inline double Gamma_2(int const n_f){return  (Gamma_2_0(n_f)*Gamma_0);}
inline double Gamma_2F(int const n_f){return  (Gamma_2_0(n_f)*Gamma_0F);}
inline double Gamma_2A(int const n_f){return  (Gamma_2_0(n_f)*Gamma_0A);}
inline double Gamma_3F(int const n_f){return  ( (n_f==5) ? 1553. : ( (n_f==4) ? 4313. : 7849. ) );}
inline double Gamma_3(int const n_f){return  Gamma_3F(n_f)/C_F;}
inline double Gamma_3_0(int const n_f){return  Gamma_3(n_f)*Gamma__0;}



static const double gamma_0q=-3*C_F;
static const double gamma__0q=1.0/gamma_0q;
inline double gamma_1q(int const n_f){return  (	C_F*C_F*(-3.0/2.0+2.0*PI2-24.0*zeta_3)
								+C_F*C_A*(-961.0/54.0-11.0/6.0*PI2+26.0*zeta_3)
								+C_F*T_F*n_f*(130.0/27.0+2.0/3.0*PI2)	);}
inline double gamma_1_0q(int const n_f){return gamma_1q(n_f)*gamma__0q;}
inline double gamma_2q(int const n_f){return  (C_F*n_f*n_f*pow(T_F,2)*(13.262002743484224 - (40*pow(M_PI,2))/27. - (32*zeta_3)/27.) +
		   pow(C_F,2)*n_f*T_F*(109.37037037037037 - (26*pow(M_PI,2))/9. - (28*pow(M_PI,4))/27. +
		      (512*zeta_3)/9.) + C_A*(C_F*n_f*T_F*(-23.755829903978054 + (2594*pow(M_PI,2))/243. +
		         (22*pow(M_PI,4))/45. - (1928*zeta_3)/27.) +
		      pow(C_F,2)*(-37.75 + (205*pow(M_PI,2))/9. + (247*pow(M_PI,4))/135. - (844*zeta_3)/3. -
		         (8*pow(M_PI,2)*zeta_3)/3. - 120*zeta_5)) +
		   pow(C_A,2)*C_F*(-47.7863511659808 - (7163*pow(M_PI,2))/486. - (83*pow(M_PI,4))/90. +
		      (3526*zeta_3)/9. - (44*pow(M_PI,2)*zeta_3)/9. - 136*zeta_5) +
		   pow(C_F,3)*(-14.5 - 3*pow(M_PI,2) - (8*pow(M_PI,4))/5. - 68*zeta_3 + (16*pow(M_PI,2)*zeta_3)/3. +
		      240*zeta_5));}


static const double gamma_0S = 0.0;
inline double gamma_1S(int const n_f){return  (gsl_pow_2(C_A)*(-160./27.+11.*PI2/9.+4.*zeta_3)
								+C_A*T_F*n_f*(-208./3.-4.*PI2)/9.
								-8.*C_F*T_F*n_f);}
inline double gamma_2S(int const n_f){return  gsl_pow_3(C_A)*(37045./729. + 6109./243.*PI2 - 319./135.*PI4 + (244.-40./3.*PI2)/3.*zeta_3 - 32.*zeta_5)
								+gsl_pow_2(C_A)*T_F*n_f*(-167800./729. -2396./243.*PI2 + 164./135.*PI4 + 1424./27.*zeta_3)
								+C_A*C_F*T_F*n_f*(1178./27. - 4./3.*PI2 - 16./45.*PI4 - 608./9.*zeta_3)
								+gsl_pow_2(C_F)*T_F*n_f*8.
								+C_A*gsl_pow_2(T_F)*n_f*(24520./729. + 80./81.*PI2 - 448./27.*zeta_3)
								+C_F*gsl_pow_2(T_F)*n_f*176./9.;}

inline double beta_0(int const n_f){return  ((11.0*C_A-4.0*T_F*n_f)/3.0);}
inline double beta__0(int const n_f){return  (1.0/beta_0(n_f));}
inline double beta_1(int const n_f){return  (34.0/3.0*C_A*C_A-20.0/3.0*C_A*T_F*n_f-4.0*C_F*T_F*n_f);}
inline double beta_1_0(int const n_f){return  (beta_1(n_f)/beta_0(n_f));}
inline double beta_2(int const n_f){return  (2857.0/54.0*C_A*C_A*C_A
								+(2.0*C_F*C_F-205.0/9.0*C_A*C_F-1415.0/27.0*C_A*C_A)*T_F*n_f
								+(44.0/9.0*C_F+158.0/27.0*C_A)*T_F*T_F*n_f*n_f 				);}
inline double beta_2_0(int const n_f){return  (beta_2(n_f)/beta_0(n_f));}
inline double beta_3(int const n_f){return  149753./6. + 3564.*zeta_3 - (1078361./162. + 6508./27.*zeta_3)*n_f
								+(50065./162 + 6472./81.*zeta_3)*n_f*n_f+1093./729.*gsl_pow_3(n_f);}
inline double beta_3_0(int const n_f){return  (beta_3(n_f)/beta_0(n_f));}

inline double Gb_0(int const n_f){return  (Gamma_0/beta_0(n_f));}
inline double gqb_0(int const n_f){return  (gamma_0q/beta_0(n_f));}

//static const double c = (Gamma_0/beta_0);

//inline double d_2(int const n_f){return C_A*(808.0/27.0-28.0*zeta_3)-224.0/27.0*T_F*n_f;}
inline double d_2(int const n_f){return C_A*(202.0/27.0-7.0*zeta_3)-56.0/27.0*T_F*n_f;}

static const double _MZ = 91.1778;
static const double _MU0 = 10;
static const double _Alpha0 = 0.175895;
static const double _Alpha0MZ = 0.11707;
static const double _Slhc_I = (7000.*7000.);
static const double _Slhc_II = (8000.*8000.);
static const double _STeva = (1960.*1960.);



enum PROD {PROD_PHOTON, PROD_Z, PROD_W, PROD_H, PROD_WRONG_STATEMENT};
enum COL {COL_PP, COL_PbarP, COL_PbarPbar};
enum CUTOFF {CUTOFF_NO, CUTOFF_GAUSS, CUTOFF_HARD, CUTOFF_DIPOL};
enum RES {RES_LO, RES_LL, RES_LLO, RES_MATCHED, RES_MATCHING_CORR, RES_WRONG_STATEMENT};
enum PDF_ERR {PDF_ERR_NO, PDF_ERR_HESSIAN, PDF_ERR_GAUSSIAN, PDF_ERR_WRONG_STATEMENT};


#define N_PARTONS 13	//!< @def N_PARTONS Number of partons.
#define N_FMAX 6	//!< @def N_PARTONS Number maximum flavours.


typedef int lha_parton;
static const lha_parton lha_tbar = -6;
static const lha_parton lha_bbar = -5;
static const lha_parton lha_cbar = -4;
static const lha_parton lha_sbar = -3;
static const lha_parton lha_ubar = -2;
static const lha_parton lha_dbar = -1;
static const lha_parton lha_g =  0;
static const lha_parton lha_d =  1;
static const lha_parton lha_u =  2;
static const lha_parton lha_s =  3;
static const lha_parton lha_c =  4;
static const lha_parton lha_b =  5;
static const lha_parton lha_t =  6;

typedef unsigned int e_parton;
static const e_parton	e_tbar = 0;
static const e_parton	e_bbar = 1;
static const e_parton	e_cbar = 2;
static const e_parton	e_sbar = 3;
static const e_parton	e_ubar = 4;
static const e_parton	e_dbar = 5;
static const e_parton	e_g = 6;
static const e_parton	e_d = 7;
static const e_parton	e_u = 8;
static const e_parton	e_s = 9;
static const e_parton	e_c = 10;
static const e_parton	e_b = 11;
static const e_parton	e_t = 12;


inline e_parton const bar(e_parton q){return (N_PARTONS-(q+1));}
inline lha_parton const bar(lha_parton q){return (-q);}
inline e_parton const e_(lha_parton q){return (q +N_FMAX);}
inline lha_parton const lha_(e_parton q){return (q -N_FMAX);}

inline uint const lha_uptype(lha_parton q){ return (1-abs(q)%2);}
inline uint const lha_downtype(lha_parton q){ return (abs(q)%2);}


static const int nMuErr=3;
static const double MuFac=2.;


#endif
