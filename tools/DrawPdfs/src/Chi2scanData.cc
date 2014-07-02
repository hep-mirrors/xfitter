#include "Chi2scanData.h"

#include <fstream>
#include <iostream>

Chi2scanData::Chi2scanData(string dirname)
{
  //read central chi2 scan
  string fname = dirname + "/chi2scan.txt";
  ifstream f(fname.c_str());
  if (!f.good())
    return;
  getline(f,label);
  f >> min >> delta >> chi2min;
  double val, c2;
  while (f >> val >> c2)
    chi2[val] = c2;

  //read chi2 scan for PDF variations
  int i = 1;
  while (true)
    {
      char chi2name[100];
      sprintf(chi2name, "chi2scan_%d.txt", i);
      string fname = dirname + "/" + chi2name;
      ifstream fpdf(fname.c_str());
      if (!fpdf.good())
	return;
      string lab;
      double pmin, pdelta, pchi2min;
      getline(fpdf,lab);
      fpdf >> pmin >> pdelta >> pchi2min;
      pdfmin[i] = pmin;
      pdfdelta[i] = pdelta;
      pdfchi2min[i] = pchi2min;

      double val, c2;
      while (f >> val >> c2)
	pdfchi2[i][val] = c2;
      i++;
    }
}
