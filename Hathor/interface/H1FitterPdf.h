#ifndef Hathor_H1FitterPdf_
#define Hathor_H1FitterPdf_

#include "HathorPdf.h"

class H1FitterPdf : public Pdf {

protected:

  std::string PDFname;
  int imember;

public:

  H1FitterPdf(const std::string str="H1FitterPdf");
  void GetPdf(double x, double muf, double h[13]);
  double GetAlphas(double mu);
  void InitMember(int i);
  int NumberPdf(void);
  std::string GetName(void);

};

#endif 
