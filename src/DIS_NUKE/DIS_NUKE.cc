#include "DIS_NUKE.h"
#include "TermData.h"
#include "xfitter_cpp_base.h"

extern "C" {
  double nuke_fast_(const double& xb, const double& q2, const int& nsf, const int& ityp, const int& kint1, const int& kord, const int& ftyp, const float& syst);
}

DIS_NUKE::DIS_NUKE(TermData* td) {
  _ftyp = td->getParamI("nuke_ftyp");
  _kint = td->getParamI("nuke_kint");
  _kord = OrderMap(td->getParamS("Order"));
  _need_f3bar = _ftyp && _kint < 0;
  // load nuclear correction tables once, otherwise they will be loaded multiple time in parallel
  double f2(1.), fl(1.), f3(1.);
  double ret = apply(0.1, 10., f2, fl, f3);
}

double DIS_NUKE::apply(const double q2, const double x, double& f2, double& fl, double& f3, double const* f3_bar_ptr/*=nullptr*/) const {
  static constexpr int ityp = 0;
  static constexpr float syst = 0.;
  double cor_f1 = nuke_fast_(x, q2, 1, ityp, _kint, _kord, _ftyp, syst);
  double cor_f2 = nuke_fast_(x, q2, 2, ityp, _kint, _kord, _ftyp, syst);
  double cor_f3 = nuke_fast_(x, q2, 3, ityp, _kint, _kord, _ftyp, syst);
  if (_kint < 0 && f3_bar_ptr) 
  {
    const double& f3_bar = *f3_bar_ptr;
    int kint_bar = -1 * _kint;
    double cor_f3_bar = nuke_fast_(x, q2, 3,  ityp, kint_bar, _kord, _ftyp, syst);
    cor_f3 = ((f3 + f3_bar) * cor_f3 - f3_bar * cor_f3_bar) / f3;
  }
  double f1 = (f2 - fl) / (2 * x);
  f1 *= cor_f1;
  f2 *= cor_f2;
  fl = f2 - 2 * x * f1;
  f3 *= cor_f3;
  return cor_f2;
}