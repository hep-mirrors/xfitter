#include "Hathor.h"

extern "C" {
  //double ttbari_(double xx[10]);
  void ttbarr_(int ll, double xx[2], double aa[2], double bb[2]);
  double sd2_(double acc, double (*f)(double*), void (*r)(int, double*, double*, double*));
  //double f_wrapper(const double x[], const double* wgt, double res[]) {
  //  return g_caller->call_f(x, *wgt, res);
  //}
}

template<class T>
class FriendCaller : public Vegas<T> {
  public:
    using Vegas<T>::Vegas;

    static FriendCaller* self;

    //static double callback(const double x[], const double* wgt, double res[]) {
    static double callback(double x[]) {
      const double* wgt;
      double res[1];
      return self->integ.f(x, *wgt, res);
      //return self->f(x, wgt, res);
    }
};

template<class T>
FriendCaller<T>* FriendCaller<T>::self = nullptr;

class HathorGenericIntegrator: public Hathor {
  void test() {
    //FriendCaller<AbstractHathor> caller(*_hathor);
    FriendCaller<AbstractHathor> caller(*this);
    FriendCaller<AbstractHathor>::self = &caller;
    double acc = 1e-4;
    sd2_(acc, FriendCaller<AbstractHathor>::callback, ttbarr_);
  }
};
