#include <string>
#include <cstring>
#include "CIJETReader.h"
#include "ReactionTheory.h"

using namespace std;

// Interface to xFitter FORTRAN routines
extern "C"
{
    void cijetinit_(const char* fname, double* mufip, double* murip, int* od, double* norm, int* estat);
    void cijetxsec_(const char* fname, double* lam, double* cl, int* npt, double* pres);
}

//------------------------------------------------------------------------------------

void CIJETReader::setinit()
{
    //note the character length in Fortran is hardwired to 100
    char ff[101];
    strncpy(ff, filepds.c_str(), 100);
    int estat = 0;  //Error status for checking amount of CI grids
    cijetinit_(ff, &muf, &mur, &Oqcd, &fut, &estat);
    if (estat == 21111801) {
        hf_errlog(21111801,"F: cset in cijet.f exceeds mset in cijetInclude.h. Increase mset and recompile.");
    }
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
    for (int i=0; i<npt; i++) xsec.push_back(pres[i]);

    return xsec;
}

void CIJETReader::setnorm(double norm)
{
    fut=norm; // overall normalization factor
//        cout<<"=========> CInorm set to "<<fut<<endl;
}

void CIJETReader::setorder(int od)
{
    Oqcd=od; //0 or 1 for LO or NLO
//        cout<<"=========> CIorder set to "<<Oqcd<<endl;
}

void CIJETReader::setscales(double mf, double mr)
{
    muf=mf;
    mur=mr;
//        cout<<"=========> CIscales set to "<<muf<<" "<<mur<<endl;
}

void CIJETReader::setci(double lam, vector<double> coe)
{
    lambda=lam;  // Lambda in TeV
    c1=coe[0];  // color singlet coefficient, LL
    c2=coe[1];  //LR
    c3=coe[2];  //RR
}
