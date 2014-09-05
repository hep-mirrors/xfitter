// Author: Daniel Britzger
// DESY, 20/04/2012

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLO_reader_2.1.0                                                //
//  D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch           //
//                                                                      //
//  The projects web page can be found at:                              //
//    http://projects.hepforge.org/fastnlo                              //
//                                                                      //
//  If you use this code, please cite:                                  //
//    T. Kluge, K. Rabbertz and M. Wobisch, hep-ph/0609285              //
//    D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch,        //
//       arXiv:1109.1310                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//
//  fastNLOAlphas
//  This class inherits the PDF interface from
//  fastNLOLHAPDF, while the alpha_s evolution
//  is superseeded by the Alphas.h class.
//
//////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <LHAPDF/LHAPDF.h>
//#include "fastnlotk/fastNLOReader.h"
//#include "fastnlotk/speaker.h"
#include "fastNLOLHAPDF.h"
#include "fastNLOQCDNUMAS.h"


using namespace std;


//______________________________________________________________________________
//
//
fastNLOQCDNUMAS::fastNLOQCDNUMAS(string name) : fastNLOLHAPDF(name) {
    SetPDGValues();
};
fastNLOQCDNUMAS::fastNLOQCDNUMAS(string name, string LHAPDFFile, int PDFSet = 0) : fastNLOLHAPDF(name,LHAPDFFile,PDFSet), fAlphasMz(0.1184) {
    SetPDGValues();
};
// Getters
double fastNLOQCDNUMAS::GetMz() const {
    return fMz;
}
double fastNLOQCDNUMAS::GetQMass(int pdgid) const {
    return QMass[pdgid];
}
int fastNLOQCDNUMAS::GetNFlavor(int nflavor) const {
    return nflavor;
}
int fastNLOQCDNUMAS::GetNLoop() const {
    return fnLoop;
}
double fastNLOQCDNUMAS::GetAlphasMz() const {
    return fAlphasMz;
};

void fastNLOQCDNUMAS::SetMz(double Mz) {
   fMz = Mz;
}

void fastNLOQCDNUMAS::SetNFlavor(int  nflavor) {
   fnFlavor = nflavor;
}

void fastNLOQCDNUMAS::SetNLoop(int  nloop) {
   fnLoop = nloop;
}

void fastNLOQCDNUMAS::SetQMass(int pdgid, double qmass) {
   QMass[pdgid] = qmass;
}

void fastNLOQCDNUMAS::SetAlphasMz(double AlphasMz , bool ReCalcCrossSection) {
   debug["SetAlphasMz"]<<"Setting alpha_s(Mz)="<<AlphasMz<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the alpha_s value at M_Z
   //
   fAlphasMz    = AlphasMz;             // new alpha_s value
   if (ReCalcCrossSection) CalcCrossSection();

}

void fastNLOQCDNUMAS::SetPDGValues() {
   // Initialize with PDG values
   QMass[0]  = PDG_MD;
   QMass[1]  = PDG_MU;
   QMass[2]  = PDG_MS;
   QMass[3]  = PDG_MC;
   QMass[4]  = PDG_MB;
   QMass[5]  = PDG_MT;
   fMz       = PDG_MZ;
   fAlphasMz = PDG_ASMZ;
   //Variable flavor number scheme
   fnFlavor = 0;
   //2-loop alpha_s evolution
   fnLoop = 2;
}

void fastNLOQCDNUMAS::SetLHAPDFValues() {
   FillPDFCache();
   fAlphasMz = LHAPDF::alphasPDF(fMz);
   fnLoop = LHAPDF::getOrderAlphaS() + 1;
   for (int i = 0; i < 6; i++)
      QMass[i] = LHAPDF::getQMass(i+1);
}

void fastNLOQCDNUMAS::InitEvolveAlphas() {
   //Ensure reasonable values are set
   //TODO Really neccessary?
   char filename[] = " ";
   int len_filename = strlen(filename);
   int lun = 6;
   qcinit_(&lun, filename, len_filename);


   //LHAPDF LO=0 while QCDNUM LO=1
   int iord = fnLoop;

   //TODO Set correct Array in q2. maybe fnloreader. getQScale...
   double qarr[2] = {1.0, 1000000};
   double wgt[2] =  {1.0, 1.0};
   //Length of array
   int n= 2;
   //Number of grid points
   int nqin = 140;
   //Real number generated grid points
   int nqout = 0;
   //Create Q2 Grid
   gqmake_(qarr, wgt, &n, &nqin, &nqout);
   setord_(&iord);
   double r2 = fMz * fMz;
   setalf_(&fAlphasMz, &r2);
   //Get Indices of Flavor Thresholds (currently just the Q mass)
   double Q2Mass[6];
   for (int i = 0; i < 6; i++)
      Q2Mass[i] = QMass[i]*QMass[i];

   int iqc = iqfrmq_(&Q2Mass[3]) ;
   int iqb = iqfrmq_(&Q2Mass[4]);
   int iqt = iqfrmq_(&Q2Mass[5]);

   //cout << iqc << " " << iqb << " " << iqt << endl;
   //When fNFlavor = 0 VFNS if >0 then FFNS
   //iqc,b,t are neglected if fnflavor =0
   setcbt_(&fnFlavor, &iqc, &iqb, &iqt);

}


//______________________________________________________________________________
double fastNLOQCDNUMAS::EvolveAlphas(double Q) const {
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   double mu2 = Q*Q;
   int ierr = 9876;
   //Number of really used flavors
   int nf = 9;
   double as = asfunc_(&mu2, &nf , &ierr);
   //cout << as << "  " << mu2 << " " << nf << endl;
   if (ierr > 0)
      error["EvolveAlphas"]<<"Alphas evolution failed. ierr = "<<ierr<<", Q = "<<Q<<endl;
   return as;

}

void fastNLOQCDNUMAS::CalcCrossSection() {
   InitEvolveAlphas();
   fastNLOLHAPDF::CalcCrossSection();
}

//______________________________________________________________________________
