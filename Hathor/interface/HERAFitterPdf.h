#ifndef Hathor_HERAFitterPdf_
#define Hathor_HERAFitterPdf_

#include "HathorPdf.h"

class HERAFitterPdf : public Pdf {

protected:

  std::string PDFname;
  int imember;

public:

  HERAFitterPdf(const std::string str="HERAFitterPdf");
  void GetPdf(double x, double muf, double h[13]);
  double GetAlphas(double mu);
  void InitMember(int i);
  int NumberPdf(void);
  std::string GetName(void);

};

#endif 
