#ifndef Hathor_xFitterPdf_
#define Hathor_xFitterPdf_

#include "HathorPdf.h"

class xFitterPdf : public Pdf {

protected:

  std::string PDFname;
  int imember;

public:

  xFitterPdf(const std::string str="xFitterPdf");
  void GetPdf(double x, double muf, double h[13]);
  double GetAlphas(double mu);
  void InitMember(int i);
  int NumberPdf(void);
  std::string GetName(void);

};

#endif 
