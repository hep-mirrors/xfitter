#include "HathorFast.h"

//template<class T>
//class FriendCaller : public Vegas<T> {
//  public:
//    FriendCaller(T& obj) : Vegas<T>(obj) {}
//    double call_f(const double x[], double wgt, double res[]) {
//      return this->integ.f(x, wgt, res);
//    }
//};

//static FriendCaller<AbstractHathor>* g_caller = nullptr;

extern "C" {
  double ttbari_(double xx[10]);
  void ttbarr_(int ll, double xx[2], double aa[2], double bb[2]);
  double sd2_(double acc, double (*f)(double*), void (*r)(int, double*, double*, double*));
  //double f_wrapper(const double x[], const double* wgt, double res[]) {
  //  return g_caller->call_f(x, *wgt, res);
  //}
}

double HathorFast::getXsection(double m_, double mur_, double muf_) {
  printf("dupa\n");
  //return 1.0;

  m = m_; mur = mur_; muf = muf_;
  update();
  
  //Vegas<AbstractHathor> integrator(*this, pdfmax);
  ////integrator.setnprn(1);
  //integrator.vegas();
  
  double acc = 1e-4;
  //g_integ = &this->f;
  double wgt;
  double x[MAXDIM];
  double res[FDIMMAX];
  //this->f(x, wgt,res);
  //FriendCaller<AbstractHathor> caller(*_hathor);
  //FriendCaller<AbstractHathor> caller(*this);
  //g_caller = &caller;
  double s1 = 0.;
  //double s1=sd2_(acc,ttbari_,ttbarr_);
  
  for (int i = 0; i < pdfmax; i++) {
    integral[i] = s1;//integrator.avgi_[i];
    error[i] = 0.;//integrator.sd_[i];
    chi2a[i] = 0.;//integrator.chi2a_[i];
  }
  
  // reset pdf member
  if (pdfmax > 1)
    pdf.InitMember(0);
  
  //return integrator.avgi;
  return s1;
}
