#include "Hathor.h"

template<class T>
class FriendCaller : public Vegas<T> {
  public:
    using Vegas<T>::Vegas;

    static FriendCaller* self;

    static double callback(const double x[],
                           const double* wgt,
                           double res[]) {
      return self->integ.f(x, *wgt, res);
    }
};

template<class T>
FriendCaller<T>* FriendCaller<T>::self = nullptr;

class HathorTEST {
  void test() {
    //FriendCaller<AbstractHathor> caller(*_hathor);
    FriendCaller<AbstractHathor> caller(*this);
    FriendCaller<AbstractHathor>::self = &caller;
    double acc = 1e-4;
    sd2(acc, FriendCaller<AbstractHathor>::callback, ttbarr);
  }
};
