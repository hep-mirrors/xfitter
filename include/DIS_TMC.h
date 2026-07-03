#pragma once
// DIS target mass corrections according to prescription
// H. Georgi and H. D. Politzer, Phys. Rev. D9, 416 (1974)

class TermData;
#include "ABM.h"
namespace abm { enum class SFproc; };

class DIS_TMC {
  // if needed to be used by other DIS reactions, add them here
  // NOTE: TMC are applied to F2, FL and different flavours separately
  // it is not possible to implement it at the level of BaseDIS classes
  template<class BaseDISC, abm::SFproc Proc> friend class ReactionBaseFFABM;
  friend class ReactionFFABM_DISNC;
  friend class ReactionFFABM_DISCC;
  template<class BaseDISC, abm::SFproc Proc> friend class ReactionBaseTBSF;
  friend class ReactionTBSF_DISNC;
  friend class ReactionTBSF_DISCC;
  public:
    DIS_TMC(TermData* td);
    bool getFlagL() {return _flag_l;}
    bool getFlagC() {return _flag_c;}
    bool getFlagB() {return _flag_b;}
  private:
    bool _flag_l;
    bool _flag_c;
    bool _flag_b;
    bool _flag_f2;
    bool _flag_fl;
    bool _flag_f3;
    int _integration_method;
    double _xmin;
    double _logxlogq2min;
    const double* _mpr;
    double apply_one_flavour(double& f2, double& fl, double& f3, 
      const abm::SFflav flav, const double q2, const double x, const abm::SFproc ncflag, 
      const int orderDefault, const int orderHQ, const int orderFL, const bool msbarm, 
      const int charge, const double polarity, const double cos2thw, const double mz);
    // private methods can be used by friends
    void apply(double& f2, double& fl, double& f3, 
      double& f2c, double& flc, double& f3c, double& f2b, double& flb, double& f3b, 
      const double q2, const double x, const abm::SFproc ncflag, 
      const int orderDefault, const int orderHQ, const int orderFL, const bool msbarm, 
      const int charge, const double polarity, const double cos2thw, const double mz);
};
