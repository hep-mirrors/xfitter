
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionAFB

  @brief A wrapper class for AFB reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2018-07-16
  */


class ReactionAFB : public ReactionTheory
{
public:
  ReactionAFB()
  {
    flagInit = false;
  }

  //    ~ReactionAFB(){};
  //    ~ReactionAFB(const ReactionAFB &){};
  //    ReactionAFB & operator =(const ReactionAFB &r){return *(new ReactionAFB(r));};

public:
  virtual string getReactionName() const { return  "AFB" ;};
  virtual void initTerm(TermData *td) override final;
  virtual void compute(TermData *td, valarray<double> &val, map<string, valarray<double> > &err);
protected:
  virtual int parseOptions(){ return 0;};

private:
  // Define a structure to pass the parameters
  struct integration_params {
    double Minv;
    ReactionTheory* ptr;
  };

  static string integration_param;
  static int integration_switch;
  static int key_param;
  static size_t alloc_space;
  static size_t calls;

  static double epsabs;
  static double epsrel;
  static double PI;
  static double GeVtofb_param, alphaEM_param, stheta2W_param, MZ_param, GammaZ_param;
  static double energy_param, eta_cut_param, pT_cut_param, y_min_param, y_max_param;

  static double e_param, gsm_param, smangle_param;
  static double foton_Vu, foton_Au, foton_Vd, foton_Ad, foton_Vl, foton_Al, foton_Vnu, foton_Anu;
  static double Z_Vu, Z_Au, Z_Vd, Z_Ad, Z_Vl, Z_Al, Z_Vnu, Z_Anu;
  static double even_foton_up, even_foton_down, even_interf_up, even_interf_down, even_Z_up, even_Z_down;
  static double odd_foton_up, odd_foton_down, odd_interf_up, odd_interf_down, odd_Z_up, odd_Z_down;

  bool flagInit;

  static double *propagators (double);


  static double uubarEF_funct (double, void *);
  static double integration_uubarEF_y (double, void *);
  double integration_uubarEF (double, double, void *);

  static double uubarEB_funct (double, void *);
  static double integration_uubarEB_y (double, void *);
  double integration_uubarEB (double, double, void *);

  static double uubarOF_funct (double, void *);
  static double integration_uubarOF_y (double, void *);
  double integration_uubarOF (double, double, void *);

  static double uubarOB_funct (double, void *);
  static double integration_uubarOB_y (double, void *);
  double integration_uubarOB (double, double, void *);


  static double ubaruEF_funct (double, void *);
  static double integration_ubaruEF_y (double, void *);
  double integration_ubaruEF (double, double, void *);

  static double ubaruEB_funct (double, void *);
  static double integration_ubaruEB_y (double, void *);
  double integration_ubaruEB (double, double, void *);

  static double ubaruOF_funct (double, void *);
  static double integration_ubaruOF_y (double, void *);
  double integration_ubaruOF (double, double, void *);

  static double ubaruOB_funct (double, void *);
  static double integration_ubaruOB_y (double, void *);
  double integration_ubaruOB (double, double, void *);


  static double ddbarEF_funct (double, void *);
  static double integration_ddbarEF_y (double, void *);
  double integration_ddbarEF (double, double, void *);

  static double ddbarEB_funct (double, void *);
  static double integration_ddbarEB_y (double, void *);
  double integration_ddbarEB (double, double, void *);

  static double ddbarOF_funct (double, void *);
  static double integration_ddbarOF_y (double, void *);
  double integration_ddbarOF (double, double, void *);

  static double ddbarOB_funct (double, void *);
  static double integration_ddbarOB_y (double, void *);
  double integration_ddbarOB (double, double, void *);


  static double dbardEF_funct (double, void *);
  static double integration_dbardEF_y (double, void *);
  double integration_dbardEF (double, double, void *);

  static double dbardEB_funct (double, void *);
  static double integration_dbardEB_y (double, void *);
  double integration_dbardEB (double, double, void *);

  static double dbardOF_funct (double, void *);
  static double integration_dbardOF_y (double, void *);
  double integration_dbardOF (double, double, void *);

  static double dbardOB_funct (double, void *);
  static double integration_dbardOB_y (double, void *);
  double integration_dbardOB (double, double, void *);


  double AFB (double, double);
};
