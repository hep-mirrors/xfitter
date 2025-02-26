#include <spline.h>

class TermData;

// DIS higher twist spline parametrisation
class DIS_HT {
  private:
    bool _flag_ht;
    tk::spline _spline_f2;
    tk::spline _spline_ft;
    std::vector<const double*> _ht_x;
    std::vector<const double*> _ht_2;
    std::vector<const double*> _ht_t;
    const double* _ht_alpha_2;
    const double* _ht_alpha_t;
  public:
    void init(TermData *td);
    void update();
    void apply(const double q2, const double x, double& f2, double& fl);
};
