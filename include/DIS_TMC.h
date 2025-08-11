#pragma once
// DIS target mass corrections according to prescription
// H. Georgi and H. D. Politzer, Phys. Rev. D9, 416 (1974)

class TermData;

class DIS_TMC {
  // if needed to be used by other DIS reactions, add them here
  // NOTE: TMC are applied to F2, FL and different flavours separately
  // it is not possible to implement it at the level of BaseDIS classes
  friend class ReactionBaseDISNC;
  friend class ReactionBaseDISCC;
  friend class ReactionFFABM_DISNC;
  friend class ReactionFFABM_DISCC;
  public:
    DIS_TMC(TermData* td);
  private:
    bool _flag_l;
    bool _flag_c;
    bool _flag_b;
    int _integration_method;
    double _xmin;
    double _logxlogq2min;
    const double* _mpr;
    int _FLAG_FAST = 0;
    double apply_one_flavour(double& f2, double& fl, double& f3, 
      const bool flag_fl, const bool flag_f3, const int flag_flavour, const double q2, const double x,
      const int ncflag, const int charge, const double polarity, const double cos2thw, const double mz);
    // private methods can be used by friends
    void apply(double& f2, double& fl, double& f3, 
      double& f2c, double& flc, double& f3c, double& f2b, double& flb, double& f3b, 
      const bool flag_fl, const bool flag_f3, const double q2, const double x, const int ncflag, 
      const int charge, const double polarity, const double cos2thw, const double mz, void* ptr_reaction);
};
