/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2011--2012
  \copyright Creative Commons license CC-BY-NC 3.0
  \version 1.1.3
_____________________________________________________________*/

#ifndef SOCK_DDIS_H_
#define SOCK_DDIS_H_

// #define DBGT_ON 1
#include "plug_DDIS.h"

// #include "dbgtools.h"

// oooooooooooooooooooooooooooooooooooooooo
class sock_DDIS_t : public plug_DDIS_t {
  double EcmSquared;
  
public:

  // ====================================================
  // sock_DDIS_t();
  
  // ====================================================
  void GetExpParams(const Coca_t& EDcoca) {
    EcmSquared = pow(EDcoca.GetReal("ECM"), 2);
    if(EcmSquared <= 0) EcmSquared = EDcoca.GetReal("ECM^2");
    if(EcmSquared <= 0) throw "Cannot get E_CM";
  }
  
  // ====================================================
  void SetEcmSquared(double s) {EcmSquared = s;}
  
  /// nedded by Fitter
  // ====================================================
  void SetParams(const vector<double>& fit_params) {
    SetParams(&fit_params.front()+10);
    /*
    // M_TheModel->SetParams(fit_params);
    int ind=10;
    M_theory.Pflux.SetA0(fit_params[ind++]);
    M_theory.Reggeon_factor = fit_params[ind++];
    M_theory.Rflux.SetA0(fit_params[ind++]);  
    #ifdef HI_TWISTS
      M_theory.SetHTxi0(fit_params[ind++]);
      M_theory.SetHTas0(fit_params[ind++]);
    #endif
    */
  }
  
  // ====================================================
  void SetParams(const double* fit_params) {
  
    M_theory.Pflux.SetA0(*fit_params++);
    M_theory.Reggeon_factor = *fit_params++;
    M_theory.Rflux.SetA0(*fit_params++);  
    #ifdef HI_TWISTS
      M_theory.SetHTxi0(*fit_params++);
      M_theory.SetHTas0(*fit_params++);
    #endif
  }
  
  /// needed for ChiSqr calculation
  //==========================================================
  double Value(double *args) {
    //--- args: xi==xIP, QQ, MX
    double xi = *args++;
    double QQ = *args++;
    double MX = *args++;
    DBG_SHOW(MX)
    DBG_SHOW(xi)
    DBG_SHOW(QQ)
    DBG_SHOW(EcmSquared)
    
    // M_theory.UseFitFlux = true;
    double beta = QQ/(QQ + MX*MX);
    double y = (QQ + MX*MX)/EcmSquared/xi;
    
    // cout << "Pomeron flux: " << M_theory.Pflux;
    // cout << "Reggeon flux: " << M_theory.Rflux;

    // double sig = M_theory.xSigmaRed(y, xi, beta, QQ);
    double sig = xi*M_theory.SigmaRed(y, xi, beta, QQ);
    // cout << "sig "<< sig << endl;
    #ifdef HI_TWISTS
      // DBG_HERE("Value")
      // double ht = M_theory.xSigmaRed_HiTwist(y, xi, beta, QQ);
      // cout << "ht "<< ht << endl;
      if(hiTwist) {
        sig += M_theory.xSigmaRed_HiTwist(y, xi, beta, QQ);
        // cout << "sig+HT "<< sig << endl;
      }
    #endif
    return sig;
  }

};

#endif
