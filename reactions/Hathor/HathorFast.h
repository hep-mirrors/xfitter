#include "Hathor.h"

class HathorFast : public Hathor {
  public:
    HathorFast(Pdf & pdf_) : Hathor(pdf_) {;};
    double getXsection(double m, double mur, double muf);
};
