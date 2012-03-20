/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

/*--------------------------------------------------------
  F2, FL
----------------------------------------------------------*/

#include "FFstrFns.h"

double FFStrFns_t::eq2[]={0, 1./9, 4./9, 1./9, 4./9, 1./9, 4./9};
const double FFStrFns_t::cF = 4./3;
double FFStrFns_t::cur_x, FFStrFns_t::cur_fqx, FFStrFns_t::cur_QQ;

#include "aqcd.h"

#define Gauss Gauss16

#ifdef ALPHAS_KNOWN
extern
#endif
AlphaS *alphas;
//bool alphas_via_RG;

#ifndef ALPHAS_KNOWN
//=======================================
void Initialize(PDFget_t _pd, int ord, double mc, double mb, double mt, double Lambda4) {
  pd = _pd;
  QCDord = ord;
  //if(alphas_via_RG) alphas = new AlphaS_RG;
  //else
  //alphas = new AlphaS_RG;
  alphas = new AlphaS_Lam;
  cout << "Fixing Lambdas: Lambda[4] = " << Lambda4 << endl;
  alphas->Initialize(Lambda4, QCDord, mc*mc, mb*mb, mt*mt);
  //alphas->Initialize(0.118,MASS_Z0*MASS_Z0, QCDord, mc*mc, mb*mb, mt*mt);
  //cout << "alpha(M_Z) = " << alphas.val(MASS_Z0*MASS_Z0, order) << endl;
  int k;
  for(k=1; k <= 6; k++) eq2sum[k] = eq2sum[k-1] + eq2[k];
}
#else
//=======================================
FFStrFns_t::FFStrFns_t(PDFget_t _pd, int ord) {
	xf = xf__ +6;
	// GAUSS_ACC = 1e-5;
  pd = _pd;
  QCDord = ord;
  int k;
  for(k=1; k <= 6; k++) eq2sum[k] = eq2sum[k-1] + eq2[k];
	C2G_fG = new C2G_fG_t(this);
}
#endif

// //=======================================
// inline double BHfactor(double QQ, double m2, double Lam2) {
  // if(QQ <= m2) return 1;
  // double tau = log(log(QQ/Lam2)/
         // log(m2/Lam2));
  // return tau > 1 ? 0 : (1 - tau);
// }




//*************************************************************************
double FFStrFns_t::F2LO(double x, double QQ, int Nlight) {
  /*
    LO contr. to F2gamma from quarks + anti-quarks
  */
  int ii;
  double val=0;
  Eval(x, QQ);
  //--- xf = x * parton_density
  for(ii=1; ii <= Nlight; ii++) val += eq2[ii]*(xf[-ii]+xf[ii]);
  return val;
}

double FFStrFns_t::C2G_fG_t::f(double x1) {
  ap->pd(ap->cur_x/x1, ap->cur_QQ, ap->xf);
  //cout << x1 <<", "<<cur_x<<", "<<xf[g] <<endl;
  return C2G(x1)*ap->xf[0];
}

//******************************************************************
real_type FFStrFns_t::CLG_fG_t::f(real_type x1) {
  ap->pd(ap->cur_x/x1, ap->cur_QQ, ap->xf);
  //cout << x1 <<", "<<cur_x<<", "<<xf[g] <<endl;
  return CLG(x1)*ap->xf[0];
}

//******************************************************************
real_type FFStrFns_t::C2F_fqsum_t::f(real_type x1) {
  return C2F(x1)*(ap->fqsumLight(ap->cur_x/x1) - ap->cur_fqx);
}

//******************************************************************
real_type FFStrFns_t::CLF_fqsum_t::f(real_type x1) {
  return CLF(x1)*ap->fqsumLight(ap->cur_x/x1);
}

#ifndef doopa
// typedef double (*RealFcn1)(double);
//================================================
double FFStrFns_t::F2(double x, double QQ, int Nf) {
  if(x > 0.9999) return 0;
  //F2mode = F2typ;
  double val;

  val = F2LO(x, QQ, cur_flav = Nf);
  #ifdef TEST_F2
    F2contr.LO = val;
  #endif
  if(QCDord < 1) return val;
  
  //---  NLO
    //double vgrid[3];
    double val1, valG;
	  cur_x = x;
	  cur_QQ = QQ;
    //--- gamma* + gluon -> q + qbar (single massless flavor)
    valG = C2G_fG->Integral(x, 1.0);
    //cout << x<<"  g = " << valG <<endl;
    //--- light quarks
    cur_fqx = fqsumLight(x);
    val1 = eq2sum[Nf]*valG +
           cF*(
						 C2F_fqsum->Integral(x, 1.0)
						 - cur_fqx*Gauss(C2F, 0.0, x, GAUSS_ACC)
					 );
    //cout << "q = " << val1 <<endl;
    //2*NLF2[1] + F0coef[gluon]*NLF2[0];
  return val + alphas->val(QQ)/(2*M_PI)*val1;
}

//================================================
double FFStrFns_t::FL(double x, double QQ, int Nf) {
  if(x > 0.9999) return 0;
  //F2mode = F2typ;
  cur_flav = Nf;
  if(QCDord < 1) return 0;

  //---  NLO
    //double vgrid[3];
    double val1, valG;
	  cur_x = x;
	  cur_QQ = QQ;
    //--- gamma* + gluon -> q + qbar (single massless flavor)
    valG = CLG_fG->Integral(x, 1.0);
    //cout << x<<"  g = " << valG <<endl;
    //--- light quarks
    val1 = eq2sum[Nf]*valG +
           cF*CLF_fqsum->Integral(x, 1.0);
    //cout << "qL = " << val1 <<endl;
  return alphas->val(QQ)/(2*M_PI)*val1;
}
#else
double FFStrFns_t::F2(double x, double QQ, int Nf) {return 0;}
double FFStrFns_t::FL(double x, double QQ, int Nf) {return 0;}

#endif
