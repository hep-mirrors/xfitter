#pragma once
namespace abm {
  enum class SFprc
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
  void set_input(const int& kschemepdfin, const int& kordpdfin, const double& rmass8in, const double& rmass10in, const int& msbarmin, double& hqscale1in, const double& hqscale2in, const int& flord);
  double calc_point_strfun(const SFprc prc, const SFtype ftype, const SFflav flav, const double q2, const double x, const int order, const int charge, const double polar, const double sin2thw, const double mz, const double* f2c_cc = nullptr);
}