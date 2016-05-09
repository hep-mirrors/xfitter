/**
* @file apkernels.h
* @author dwilhelm
* @date	26.12.2012
* @brief Definition of plus-distributions and Altarelli-Parisi-Splitting-Kernels
*
*/
#ifndef APKERNELS_h
#define APKERNELS_h

#include "konst.h"

namespace AP_Kernels {

//					-------- PLUS-Distributions  --------

typedef double (*DistFunc)(double z, double x, double fz, double fx, int n);

struct DistParam{
	double a,b;
	double DistArg;
	DistFunc f;
	void* params;
};

/*inline double diracdelta(double z, void*params){
	DistParam* p = (DistParam*)params;
	return (p->f(p->DistArg,p->params))/(p->b-p->a);
}*/

inline double diracdeltaf(double z, double x, double fz, double f, int n_f){
	return ( f/(1.-x));
}


//----------- P
inline double _Pf(double z, double fz, double f){
	return ( (fz-f)/(1.-z) );
}

inline double _Pf(double x, double f){
	return ( log(1.-x)*f);
}

inline double Pf(double z, double x, double fz, double f){
	return ( _Pf(z, fz, f) +  _Pf(x, f)/(1.-x));
}

inline double Pf(double z, double x, double fz, double f, int n_f){
	return ( Pf(z,x,fz,f));
}

//----------- PP
inline double _PPf(double z, double fz, double f){
	return ( ((fz-f)*2.*log(1.-z) - log(z)*fz)/(1.-z) );
}

inline double _PPf(double x, double f){
	return ( (gsl_pow_2(log(1.-x))-PI2/6.)*f );
}

inline double PPf(double z, double x, double fz, double f){
	return ( _PPf(z, fz, f) +  _PPf(x, f)/(1.-x) );
}

inline double PPf(double z, double x, double fz, double f, int n_f){
	return ( PPf(z,x,fz,f));
}

//					-------- ALTARELLI-PARISI  --------


//----------- P_qq
inline double _P_qqf(double z, double fz, double f){
	return ( 4.*C_F*(2*_Pf(z,fz,f) - (1.+z)*fz ) );
}

inline double _P_qqf(double x, double f){
	return ( 4.*C_F*(2*_Pf(x,f) + 1.5*f) );
}

inline double P_qqf(double z, double x, double fz, double f){
	return	( _P_qqf(z, fz, f) + _P_qqf(x, f)/(1.-x) );
}

inline double P_qqf(double z, double x, double fz, double f, int n){
	return	( _P_qqf(z, fz, f) + _P_qqf(x, f)/(1.-x) );
}

//----------- P_qg
inline double P_qg(double z){
	return ( 4.*T_F*(gsl_pow_2(z)+gsl_pow_2(1-z)));
}

inline double P_qgf(double z, double x, double fz, double f){
	return ( P_qg(z)*fz );
}

inline double P_qgf(double z, double x, double fz, double f, int n){
	return ( P_qg(z)*fz );
}

//----------- P_gq
inline double P_gq(double z){
	return ( 4.*C_F*(1+gsl_pow_2(1.-z))/z );
}

inline double P_gqf(double z, double x, double fz, double f){
	return ( P_gq(z)*fz );
}

inline double P_gqf(double z, double x, double fz, double f, int n){
	return ( P_gq(z)*fz );
}

//----------- P_gg
inline double _P_ggf(double z, double fz, double f, int n_f){
	return ( 8.*C_A*( _Pf(z,fz,f) + (1./z-2+z*(1.-z))*fz ) );
}

inline double _P_ggf(double x, double f, int n_f){
	return ( 8.*C_A*_Pf(x,f) + 2.*beta_0(n_f)*f );
}

inline double P_ggf(double z, double x, double fz, double f, int n_f){
	return ( _P_ggf(z,fz,f,n_f) + _P_ggf(x,f,n_f)/(1.-x) );
}

//----------- P_gg
inline double _P_ggPrimef(double z, double fz, double f, int n_f){
	return ( 8.*C_A*( _Pf(z,fz,f) + (1./z-2+z*(1.-z))*fz ) );
}

inline double _P_ggPrimef(double x, double f, int n_f){
	return ( 8.*C_A*_Pf(x,f));
}

inline double P_ggPrimef(double z, double x, double fz, double f, int n_f){
	return ( _P_ggPrimef(z,fz,f,n_f) + _P_ggPrimef(x,f,n_f)/(1.-x) );
}


//					-------- ALTARELLI-PARISI Convolutions --------


//----------- P_qqq
inline double _P_qqqf(double z, double fz, double f){
	return ( 16.*gsl_pow_2(C_F)*(4.*_PPf(z,fz,f) + 6.*_Pf(z,fz,f) +(-5. -z + 3.*(1.+z)*log(z) - 4.*(1.+z)*log(1.-z))*fz) );
}

inline double _P_qqqf(double x, double f){
	return ( 16.*gsl_pow_2(C_F)*( 4.*_PPf(x,f) + 6.*_Pf(x,f) + 0.25*9.*f ) );
}

inline double P_qqqf(double z, double x, double fz, double f){
	return ( _P_qqqf(z,fz,f) + _P_qqqf(x,f)/(1.-x) );
}

inline double P_qqqf(double z, double x, double fz, double f, int n_f){
	return ( P_qqqf(z,x,fz,f));
}

//----------- P_qqg
inline double P_qqg(double z){
	return ( 16.*C_F*T_F*( -0.5 + 2.*z -2.*z*z*log(z) + (z*z+gsl_pow_2(1.-z))*log(gsl_pow_2(1.-z)/z) ) );
}

inline double P_qqgf(double z, double x, double fz, double f){
	return ( P_qqg(z)*fz );
}

inline double P_qqgf(double z, double x, double fz, double f, int n_f){
	return ( P_qqgf(z,x,fz,f));
}

//----------- P_qgq
inline double P_qgq(double z){
	return ( 16.*C_F*T_F*( 1. + 4./(3.*z) -z -4./3.*z*z +2.*(1.+z)*log(z) ) );
}

inline double P_qgqf(double z, double x, double fz, double f){
	return ( P_qgq(z)*fz );
}

inline double P_qgqf(double z, double x, double fz, double f, int n_f){
	return ( P_qgqf(z,x,fz,f));
}

//----------- P_qgg
inline double P_qgg(double z, int n_f){
	return ( 32.*C_A*T_F*( 2./(3.*z) + 0.5 + 4.*z - 31./6.*z*z + (z*z + gsl_pow_2(1-z))*log(1-z) + (1+4.*z)*log(z) )
			+2.*beta_0(n_f)*P_qg(z));
}

inline double P_qggf(double z, double x, double fz, double f, int n_f){
	return ( P_qgg(z, n_f)*fz );
}


//----------- P_gqq
inline double P_gqq(double z){
	return ( 16.*C_F*C_F*(2. -0.5*z  + 2./z*log(z) + (1.+gsl_pow_2(1.-z))/z*log(gsl_pow_2(1.-z)/z)) );
}

inline double P_gqqf(double z, double x, double fz, double f){
	return ( P_gqq(z)*fz );
}

inline double P_gqqf(double z, double x, double fz, double f, int n_f){
	return ( P_gqqf(z,x,fz,f));
}

//----------- P_gqg	= P_qgq
inline double P_gqg(double z, int n_f){
	return (2.*n_f*P_qgq(z));
}
//inline double _P_gqgf(double x, double f, int n_f){ return 0.0}

inline double P_gqgf(double z, double x, double fz, double f, int n_f){
	return  (2.*n_f*P_qgqf(z,x,fz,f));
}

//----------- P_ggq
inline double P_ggq(double z, int n_f){
	return ( 32.*C_A*C_F*( -31./(6.*z) + 4. + 0.5*z + 2./3.*gsl_pow_2(z) + (1.+gsl_pow_2(1-z))/z*log((1.-z)/z) - (4.+z)*log(z) )
			+2.*beta_0(n_f)*P_gq(z));
}

//inline double _P_ggqf(double x, double f, int n_f){ return 0.0}

inline double P_ggqf(double z, double x, double fz, double f, int n_f){
	return ( P_ggq(z,n_f)*fz );
}


//----------- P_ggg

inline double _P_gggf(double z, double fz, double f, int n_f){
	return ( 4.*beta_0(n_f)*_P_ggf(z,fz,f,n_f)
			+ 64.*gsl_pow_2(C_A)*(  _PPf(z,fz,f)
									+(-11./(3.*z) +3. -3.*z +11./3.*z*z
										-2.*(1.+z)*log(z)
										+(1./z-2.+z-z*z)*log(gsl_pow_2(1.-z)/z) )*fz ));
}

inline double _P_gggf(double x, double f, int n_f){
	return ( 4.*beta_0(n_f)*(_P_ggf(x,f,n_f) -beta_0(n_f)*f )+64.*gsl_pow_2(C_A)*_PPf(x,f));
}

inline double P_gggf(double z, double x, double fz, double f, int n_f){
	return ( _P_gggf(z,fz,f,n_f) + _P_gggf(x,f,n_f)/(1.-x) );
}


//					-------- D Functions --------

//----------- D_qq
inline double _D_qqf(double z, double fz, double f){
	return ( _P_qqqf(z,fz,f) + P_qgq(z)*fz);
}

inline double _D_qqf(double x, double f){
	return ( _P_qqqf(x,f));
}

inline double D_qqf(double z, double x, double fz, double f){
	return ( _D_qqf(z,fz,f) + _D_qqf(x,f)/(1.-x));
}

inline double D_qqf(double z, double x, double fz, double f, int n){
	return ( _D_qqf(z,fz,f) + _D_qqf(x,f)/(1.-x));
}


//----------- D_qg
inline double _D_qgf(double z, double fz, double f, int n_f){
	return ( (P_qqg(z) + P_qgg(z,n_f))*fz);
}

inline double _D_qgf(double x, double f){
	return ( 0.0 );
}

inline double D_qgf(double z, double x, double fz, double f, int n_f){
	return ( _D_qgf(z, fz, f,n_f));
}

//----------- D_gq
inline double _D_gqf(double z, double fz, double f, int n_f){
	return ( (P_gqq(z) + P_ggq(z,n_f))*fz);
}

inline double _D_gqf(double x, double f){
	return ( 0.0);
}

inline double D_gqf(double z, double x, double fz, double f, int n_f){
	return ( P_gqqf(z,x,fz,f) + P_ggqf(z,x,fz,f,n_f));
}

//----------- D_gg
inline double _D_ggf(double z, double fz, double f, int n_f){
	return (( P_gqg(z,n_f)*fz + _P_gggf(z,fz,f,n_f)));
}

inline double _D_ggf(double x, double f, int n_f){
	return (  _P_gggf(x,f,n_f));
}

inline double D_ggf(double z, double x, double fz, double f, int n_f){
	return (P_gqgf(z,x,fz,f,n_f) + P_gggf(z,x,fz,f,n_f));
}

}  // namespace AP_Kernels





#endif
