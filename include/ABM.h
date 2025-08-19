#pragma once

#include <vector>

namespace abm {
  enum class SFproc
  {
    nc,
    cc,
  };
  enum class SFtype
  {
    f2,
    fl,
    f3,
  };
  enum class SFflav
  {
    l,
    c,
    b,
  };
  void initgridconst();
  void pdffillgrid();
  void update_ckm_matrix(const std::vector<const double*>& ckm);
  void set_scheme_and_order(const int kschemepdf, const int kordpdf, 
                 const int msbarm, const int flord, 
                 int kordhq = -1, int kordf2 = -1, int kordfl = -1, int kordf3 = -1, int kordalps = -1,
                 const bool hqnons = false
                 );
  void set_hq_masses(const double mc, const double mb);
  void set_hq_scales(const double hqscale1, const double hqscale2);
  void set_xbmin(const double val);
  void set_xbmax(const double val);
  double calc_point_strfun(const SFproc proc, const SFtype ftype, const SFflav flav, const double q2, const double x, 
    const int order, const int orderDefault, const int orderHQ, const int orderFL, const bool msbar, 
    const int charge, const double sin2thw, const double polar, const double mz, const double* f2c=nullptr);
  double calc_point_strfun_CC(const SFtype ftype, const SFflav flav, const double q2, const double x, 
    const int order, const int orderDefault, const int orderHQ, const int orderFL, const bool msbar, 
    const int charge, const double* f2c=nullptr);
  double calc_point_strfun_NC(const SFtype ftype, const SFflav flav, const double q2, const double x, 
    const int order, const int orderDefault, const int orderHQ, const int orderFL, const bool msbar, 
    const int charge, const double sin2thw, const double polar, const double mz);
}
