#include <vector>

class TermData;
namespace tk {class spline;}


// DIS higher twist spline parametrisation
class DIS_HT {
  // TODO do not use friends, instead make sure this class is used only by ReactionBaseDISNC
  // and ReactionBaseDISCC reactions and make the corresponding pointers there private
  friend class ReactionBaseDISNC;
  friend class ReactionBaseDISCC;
  friend class ReactionFFABM_DISNC;
  friend class ReactionFFABM_DISCC;
  public:
    DIS_HT(TermData* td);
    ~DIS_HT();
  private:
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
  private: // TODO: make public when used only by ReactionBaseDISNC, ReactionBaseDISCC
    void update();
    void apply(const double q2, const double x, double& f2, double& fl);
};
