/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

/*
  \file aqcd.h
  \details QCD alpha for VFNS
  \author W. Slominski
  \date 2011-09-09
  \version 2.05
*/

#ifndef AQCD_H_
#define AQCD_H_
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#ifdef USE_GSL
  #include "gsl/gsl_sf_lambert.h"
  #define LambertW1 gsl_sf_lambert_Wm1
#else
  double LambertW1(double z);
#endif
// #include "fiasco.h"
// #include "basic_defs.h"
#include "physparams.h"
#include "dbgtools.h"

/// Pure virtual base class for running alpha_QCD
// oooooooooooooooooooooooooooooooooooooooooooo
class AlphaS {
  protected:
    void bas_init(int QCD_ord, double Qc2, double Qb2, double Qt2);

    double beta[2][1+MAX_FLAVORS],
           brat[1+MAX_FLAVORS],TranScale[1+MAX_FLAVORS];
    double scale;
    int QCDord;

  public:
    //===============================================
    int Afl(double QQ) {
      int nf;
      for(nf=4; nf <= MAX_FLAVORS; nf++) if(QQ < TranScale[nf]) break;
      return nf-1;
    }

    virtual double val(double QQ) = 0;
    virtual void Initialize(double a0, double QQ0, int QCDord, double Qc2, double Qb2, double Qt2) = 0;
    virtual void Initialize(double Lam4, int QCDord, double Qc2, double Qb2, double Qt2) = 0;
    
    // #ifdef PHYSPARAMS_H_
    #if 1
      virtual void Initialize(double Lam4, const PhysParams_t& Phys)
      = 0;
      virtual void Initialize(double a0, double QQ0, const PhysParams_t& Phys)
      = 0;
    #endif
    
    /// Set scale \c s such that \c val(mu^2) returns s*alpha_s(mu^2).
    virtual void SetScale(double s) = 0;
};

  //----------   RG (beta fcn) mode  -------------------------
/// Running alpha_s from RG equation.
class AlphaS_RG : public AlphaS {
  double etaTran[1+MAX_FLAVORS];

  //===============================================
  double alphainv(double QQ, double eta0, double QQ0, int nf) {
    //--- returns 4*pi/alpha_s
    DBG_SHOW(sqrt(QQ0))
    DBG_SHOW(sqrt(QQ))
    DBG_SHOW(nf)
    DBG_SHOW(4*M_PI/eta0)
    double etaLO = eta0 + beta[0][nf]*log(QQ/QQ0);
    if(!QCDord) return etaLO;
    double r10 = brat[nf];
    //printf("%g/%g -> %g\n",QQ, QQ0, -r10*(1+LambertW(-exp(-1-etaLO/r10)/r10)));
    // double z = -exp(-1-etaLO/r10)*(1+eta0/r10);
    // DBG_SHOW(z)
    // return -r10*(1 + LambertW1(z) );
    return -r10*(1 + LambertW1(-exp(-1-etaLO/r10)*(1+eta0/r10)) );
  }

  public:
    //===============================================
    AlphaS_RG() : AlphaS() {
      QCDord = -1;
      scale = 4*M_PI;
    }

    void Initialize(double Lam4, int QCDord, double Qc2, double Qb2, double Qt2)
      {throw "Initialize(Lam4, ...) - not implemented for AlphaS_RG.";}
      
    void Initialize(double a0, double QQ0, int QCDord, double Qc2, double Qb2, double Qt2);
                    
    //===============================================
    // #ifdef PHYSPARAMS_H_
    #if 1
      void Initialize(double Lam4, const PhysParams_t& Phys) {
        Initialize(Lam4, Phys.QCDorder, Phys.m2[c_quark], Phys.m2[b_quark], Phys.m2[t_quark]);
      }
      void Initialize(double a0, double QQ0, const PhysParams_t& Phys) {
        // SHOW(a0)
        // cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
        Initialize(a0, QQ0, Phys.QCDorder, Phys.m2[c_quark], Phys.m2[b_quark], Phys.m2[t_quark]);
      }
    #endif
    

    //===============================================
    void SetScale(double s) { scale = 4*M_PI*s; }

    //===============================================
    double val(double QQ) {
      #ifdef SAFE_MODE
      // #error "SAFE_MODE"
      // SHOW(QCDord)
      if(QCDord < 0) throw "alpha_s not initialized.";
      #endif
      // for(int n=4; n <= MAX_FLAVORS; n++) cout <<n <<": " << sqrt(TranScale[n]) <<"  "<< 4*M_PI/etaTran[n] << endl;
      // cout <<"mu = "<< setprecision(10) << sqrt(QQ) << endl;
      int nf = Afl(QQ);
      DBG_SHOW(scale/(4*M_PI))
      DBG_SHOW(nf)
      return scale/alphainv(QQ, etaTran[nf], TranScale[nf], nf);
    }

};

  //----------   Lambda mode  -------------------------
/// Running alpha_s via Lambda_QCD[Nflavors]
class AlphaS_Lam : public AlphaS {
  bool vLambda;
  double Lambda2[MAX_FLAVORS+1], DlnLam2[MAX_FLAVORS+1];
  double m1sq, Alf0;
  int Nflav;
  int FindZero(double xini, double acc, double *x0);
  double Alpha(double mu2, int Nf, double Lam2);
  double dAlphaNext(double Lam2);
  
  void FixLambdas(double Lam4);
  void FixLambdas(double alpha0, double QQ0, int nf0);

  public:
    //===============================================
    AlphaS_Lam() : AlphaS() {
      QCDord = -1;
      scale = 2*M_PI;
    }

    //===============================================
    void SetScale(double s) { scale = 2*M_PI*s; }
    
    //===============================================
    void Initialize(double a0, double QQ0, int QCDord,
                            double Qc2, double Qb2, double Qt2) {
      bas_init(QCDord, Qc2, Qb2, Qt2);
      FixLambdas(a0, QQ0, Afl(QQ0));
    }

    //===============================================
    void Initialize(double Lam4, int QCDord,
                            double Qc2, double Qb2, double Qt2) {
      bas_init(QCDord, Qc2, Qb2, Qt2);
      FixLambdas(Lam4);
    }

    //===============================================
    // #ifdef PHYSPARAMS_H_
    #if 1
      void Initialize(double Lam4, const PhysParams_t& Phys) {
        Initialize(Lam4, Phys.QCDorder, Phys.m2[c_quark], Phys.m2[b_quark], Phys.m2[t_quark]);
      }
      void Initialize(double a0, double QQ0, const PhysParams_t& Phys) {
        Initialize(a0, QQ0, Phys.QCDorder, Phys.m2[c_quark], Phys.m2[b_quark], Phys.m2[t_quark]);
      }
    #endif
    
    //=================================
    double alpha_s_o2pi(double t, int nf);

    //=======================================
    double val(double QQ) {
      #ifdef SAFE_MODE
      if(QCDord < 0) throw "alpha_s not initialized.";
      #endif
      return scale*alpha_s_o2pi(log(QQ), Afl(QQ));
    }

    //=======================================
    double GetLambda(int n) {return sqrt(Lambda2[n]);}

};

#endif
