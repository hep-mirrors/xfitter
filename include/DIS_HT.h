#pragma once
#include <vector>

class TermData;
namespace tk {class spline;}


// DIS higher twist spline parametrisation
class DIS_HT {
  // if needed to be used by other DIS reactions, add them here
  // NOTE: HT is applied only to F2 and FL light flavour part
  // it is not possible to implement it at the level of BaseDIS classes
  friend class ReactionBaseDISNC;
  friend class ReactionBaseDISCC;
  friend class ReactionFFABM_DISNC;
  friend class ReactionFFABM_DISCC;
  public:
    DIS_HT(TermData* td, const bool flag_cc);
    ~DIS_HT();
  private:
    // CC flag
    bool _flag_cc;
    // splines
    tk::spline* _spline_f2;
    tk::spline* _spline_ft;
    // parameters via pointers (can be fitted)
    std::vector<const double*> _ht_x;
    std::vector<const double*> _ht_2;
    std::vector<const double*> _ht_t;
    const double* _ht_alpha_2;
    const double* _ht_alpha_t;
    // flags which store info whether the above point to newly allocated (in this class) memory
    std::vector<bool> _isnew_ht_x;
    std::vector<bool> _isnew_ht_2;
    std::vector<bool> _isnew_ht_t;
    bool _isnew_ht_alpha_2;
    bool _isnew_ht_alpha_t;
    // private methods can be used by friends
    void update();
    void apply(const double q2, const double x, double& f2, double& fl);
};
