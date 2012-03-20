/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef _FF_STR_FNS_H_
#define _FF_STR_FNS_H_

#include "mathut.h"
#include "ifun.h"
// #define MASS_Z0 91.1876

#define ALPHAS_KNOWN

typedef void (*PDFget_t)(double x, double QQ, double f[]);
//--- f[-6:6]

class FFStrFns_t {
	PDFget_t pd;
	double xf__[13], *xf;
	static double eq2[];
	double eq2sum[7];
	static const double cF;
	double GAUSS_ACC;
	int QCDord;
	static double cur_x, cur_fqx, cur_QQ;//, m2crit;
	int cur_flav;
	
public:
	#ifdef ALPHAS_KNOWN
		FFStrFns_t (PDFget_t _pd=NULL, int ord=0);
	#else
		void Initialize(PDFget_t _pd, int ord, double mc, double mb, double mt, double Lambda4);
	#endif
	
//=======================================
void SetPDF(PDFget_t _pd) {
  pd = _pd;
}

//=======================================
inline void Eval(double x, double QQ) {
  pd(x,QQ,xf);
}

//===============================
static real_type C2F(real_type x) {
//--- Furmanski, Pertonzio App. I, divided by cF
//--- (Gluck Reya PRD28, Eq. 3.7) B_q(x) = 2*cF*C2F(x)
//--- cf. C2Fcnv.mws
  double xb=1-x;
  return (1+x*x)/xb*(log(xb/x)-3./4) + (9+5*x)/4.0;
}

//===============================
static real_type C2G(real_type x) {
//--- Furmanski, Pertonzio App. I, divided by NF
//--- (Gluck Reya PRD28, Eq. 3.7) B_G(x) = 2*NF*C2G(x)
  double xb=1-x;
  return (x*x+xb*xb)*log(xb/x)-1.0+8.0*x*xb;
}

//===============================
static real_type CLF(real_type x) {
//--- Furmanski, Pertonzio App. I, CL = C2 - C1, divided by cF
  return 2*x;
}

//===============================
static real_type CLG(real_type x) {
//--- Furmanski, Pertonzio App. I, CL = C2 - C1, divided by NF
  return 4*x*(1-x);
}

#define CSID(nn) nn
#define NEW_GF(nn) \
	class nn##_t : public ifun_t { \
		FFStrFns_t *ap; \
		double f(double x); \
		public: \
			nn##_t(FFStrFns_t* a, double acc_=1e-5) : ifun_t(acc_) {ap = a;} \
	} *nn; \
	friend class nn##_t;

	// class C2G_fG_t : public ifun_t {
		// FFStrFns_t *ap;
		// double f(double x);
		// public:
			// C2G_fG_t(FFStrFns_t* a, double acc_=1e-5) : ifun_t(acc_) {ap = a;}
	// } *C2G_fG;
	// friend class C2G_fG_t;

	NEW_GF(C2G_fG)
	NEW_GF(CLG_fG)
	NEW_GF(C2F_fqsum)
	NEW_GF(CLF_fqsum)


//******************************************************************
double fqsumLight(real_type x) {
  return F2LO(x,cur_QQ,cur_flav);
}

double F2LO(double x, double QQ, int Nlight=3);
// static real_type C2G_fG(real_type x1);
// real_type CLG_fG(real_type x1);
// real_type C2F_fqsum(real_type x1);
// real_type CLF_fqsum(real_type x1);
	
  // void SetPDF(PDFget_t _pd);
  double F2(double x, double QQ, int Nf);
  double FL(double x, double QQ, int Nf);
	
};

#endif
