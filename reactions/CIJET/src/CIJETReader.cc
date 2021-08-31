#include <string>
#include <cstring>
#include "CIJETReader.h"
#include "ReactionTheory.h"

using namespace std;

// Interface to xFitter FORTRAN routines
extern "C"
{
  void cijetinit_(const char* fname, double* mufip, double* murip, int* od, double* norm);
  void cijetxsec_(const char* fname, double* lam, double* cl, int* npt, double* pres);
  // wrappers for applgrid functionality
  double appl_fnalphas_(double* Q);
  void appl_fnpdf_(double* x, double* mf, double* pdf);
}
double appl_fnalphas_(double* Q)
{
    return alphas_wrapper_(*Q);
}
void appl_fnpdf_(double* x, double* mf, double* pdf)
{
    double q = *mf;
    std::valarray<double> pdfV(13);
    pdf_xfxq_wrapper_(*x, q, &pdfV[0]);
    for(int i = 0; i < 13; i++) pdf[i] = pdfV[i];
}

//------------------------------------------------------------------------------------

void CIJETReader::setinit()
{
	//note the character lenth in Fortran is hardwired to 100
        char ff[101];
        strncpy(ff, filepds.c_str(), 100);
	cijetinit_(ff, &muf, &mur, &Oqcd, &fut);
}

vector<double> CIJETReader::calcxsec(void)
{
        char ff[101];
        strncpy(ff, filepds.c_str(), 100);
        double invlamsq, cpl[3];
        double pres[200];
        int npt;

        invlamsq=1.0/lambda/lambda;
        cpl[0]=c1;
        cpl[1]=c2;
        cpl[2]=c3;

        cijetxsec_(ff, &invlamsq, cpl, &npt, pres);

	vector<double> xsec;
	for (int i=0; i<npt; i++) { 
//        	cout<<"gaojun "<<pres[i]<<endl;
		xsec.push_back(pres[i]);
	}
	return xsec;
}

void CIJETReader::setnorm(double norm)
{
	fut=norm; // overall normalization factor
//	cout<<"=========> CInorm set to "<<fut<<endl;
}

void CIJETReader::setorder(int od)
{
	Oqcd=od; //0 or 1 for LO or NLO
//	cout<<"=========> CIorder set to "<<Oqcd<<endl;
}

void CIJETReader::setscales(double mf, double mr)
{
	muf=mf;
	mur=mr;
//	cout<<"=========> CIscales set to "<<muf<<" "<<mur<<endl;
}

void CIJETReader::setci(double lam, vector<double> coe)
{
	lambda=lam;  // Lambda in TeV
	c1=coe[0];  // color singlet coefficient, LL
	c2=coe[1];  //LR
	c3=coe[2];  //RR
}

void CIJETReader::GetPDF_(double x, double mf, double pdf[13])
{
  double q = mf;
  std::valarray<double> pdfV(13);
  //_reactionTheory->xfx(x, q, &pdfV[0]);	//Not supported in xFitter 2.1
  pdf_xfxq_wrapper_(x, q, &pdfV[0]);
  for(int i = 0; i < 13; i++)
    pdf[i] = pdfV[i];
}

double CIJETReader::GetAs_(double mr) 
{
  double q = mr;
  //return _reactionTheory->alphaS(q);	//Not supported in xFitter 2.1
  return alphas_wrapper_(q);
}

