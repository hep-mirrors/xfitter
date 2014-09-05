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

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <LHAPDF/LHAPDF.h>
//#include "fastnlotk/speaker.h"
#include "fastNLOReader.h"
#include "fastNLOLHAPDF.h"

using namespace std;



//______________________________________________________________________________


fastNLOLHAPDF::fastNLOLHAPDF(string name) : fastNLOReader(name) , fnPDFs(0) , fiPDFMember(0) , fchksum(0.) {
   info["fastNLOLHAPDF"]<<"Please initialize a PDF file using SetLHAPDFFilename( PDFFile ) and a PDF set using SetLHAPDFMember(int PDFMember)"<<std::endl;
}


//______________________________________________________________________________


fastNLOLHAPDF::fastNLOLHAPDF(string name, string LHAPDFFile, int PDFMember) : fastNLOReader(name) , fchksum(0.) {
   SetLHAPDFFilename(LHAPDFFile);
   SetLHAPDFMember(PDFMember);
   // Call additional initialization. Not necessary for LHAPDF.
   InitEvolveAlphas();
   // Everything set. Do cross sections calculation.
   CalcCrossSection();
}

   // Getters
   int fastNLOLHAPDF::GetIPDFMember() const {
      return fiPDFMember;
   };
   int fastNLOLHAPDF::GetNPDFMembers() const {
      return fnPDFs;
   };
   int fastNLOLHAPDF::GetNPDFMaxMember() const {
      return fnPDFs-1;
   };


//______________________________________________________________________________


double fastNLOLHAPDF::EvolveAlphas(double Q) const {
   //debug<<"EvolveAlphas with Q="<<Q<<endl;
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   // the alpha_s evolution is done within LHAPDF.
   //
   // WARNING: You cannot change alpha_s(Mz), but is is
   // defined with the pdf. 'alphasMz' is not used here!
   //
   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   return PDF->alphasQ(Q);
   #else
   return LHAPDF::alphasPDF(Q);
   #endif
}


//______________________________________________________________________________


bool fastNLOLHAPDF::InitPDF() {
   //
   //  Initalize some necessary LHAPDF parameters
   //  return true, if successful initialization
   //  return false, if PDF initialization failed
   //
   // LHAPDF interface:
   // security, if multiple instance with different pdfs are instantiated.
   // we always reinitialized the set PDF-set.

   if (fLHAPDFFilename == "") {
      error["InitPDF"]<<"Empty LHAPDF filename! Please define a PDF set here!\n";
      return false;
   }

   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   //Not needed in LHAPDF6 case
   return true;

   #else
   //LHAPDF::setVerbosity(LHAPDF::SILENT);
   LHAPDF::setVerbosity(LHAPDF::LOWKEY);
   // Do not use the ByName feature, destroys ease of use on the grid without LHAPDF
   //LHAPDF::initPDFSetByName(fLHAPDFFilename);
   //cout << "PDF set name " << fLHAPDFFilename << endl;
   if (fchksum == 0 || fchksum != CalcChecksum(1.)) {
      // need to reset LHAPDF.
      debug["InitPDF"]<<"Need to reset lhapdf. fchksum="<<fchksum<<"\tCalcChecksum(1.)="<<CalcChecksum(1.)<<endl;
      LHAPDF::initPDFSet(fLHAPDFFilename);
      fnPDFs = LHAPDF::numberPDF()+1; // LHAPDF counts 0-44 and returns, 44 which must be 45
      if (fnPDFs < fiPDFMember+1) {
         error["InitPDF"]<<"There are only "<<fnPDFs<<" pdf sets within this LHAPDF file. You were looking for set number "<<fiPDFMember<<std::endl;
         return false;
      }
      LHAPDF::initPDF(fiPDFMember);
   }
   fchksum = CalcChecksum(1.);
   return true;
   #endif
}


//______________________________________________________________________________



vector<double> fastNLOLHAPDF::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //
   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   vector <double> xfx;
   for (int id=-6; id<7; id++) {
      xfx.push_back(PDF->xfxQ(id, xp, muf));
   }
   return xfx;
   #else
   return LHAPDF::xfx(xp,muf);
   #endif
}


//______________________________________________________________________________


void fastNLOLHAPDF::SetLHAPDFFilename(string filename) {
   if (filename != fLHAPDFFilename) fchksum = 0;
   fLHAPDFFilename = filename;
   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   PDFSet = new LHAPDF::PDFSet(filename);
   fnPDFs = PDFSet->size();
   #else
   // Reset pdfset member to zero
   fiPDFMember = 0;
   // KR: Reactivated this. Why was it switched off?
   // --> Mass settings etc. can be read from LHAPDF after setting the filename, i.e. the set.
   InitPDF();
   #endif
}


//______________________________________________________________________________


void fastNLOLHAPDF::SetLHAPDFMember(int set) {
   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   PDF = PDFSet->mkPDF(set);
   #else
   fiPDFMember = set;
   if (fchksum == CalcChecksum(1.)) {  // nothin has changed? we set only the pdfmember
      debug["SetLHAPDFMember"]<<"Changing only pdfmember!"<<endl;
      LHAPDF::initPDF(fiPDFMember);
      fchksum = CalcChecksum(1.);
   } else  {
      debug["SetLHAPDFMember"]<<"Demanding full re-initalization of PDF."<<endl;
      fchksum = 0;
   }
   //InitPDF();
   #endif
}


//______________________________________________________________________________


void fastNLOLHAPDF::PrintPDFInformation() const {
   //
   // print out the information about the currently used LHAPDF file.
   // unfortunately there is no getter for lhapdf-filename or
   // used pdf-member-id available.
   // One must take care, that one is always using the desired pdf.
   //
   // e.g. If one has two fastNLOReader instances and one initalizes the
   // second instance with another pdf. Then also the first one is using this
   // pdf when evaluating CalcCrossSection (after a PDFCacheRefilling).
   //

   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   cout << PDFSet->description();
   #else
   printf(" ##################################################################################\n");
   printf(" #  fastNLOLHAPDF::PrintCurrentLHAPDFInformation.\n");
   printf(" #      Your currently initalized pdf is called:\n");
   LHAPDF::getDescription();
   printf(" #      Information about current PDFMember in current LHAPDF-file cannot be displayed.\n");
   printf(" #      Please use fastNLOReader::SetLHAPDFMember(int) to choose a pdf-set.\n");
   printf(" ##################################################################################\n");
   #endif
}

void fastNLOLHAPDF::SetMz(double Mz) {
   warn["SetMz"]<<"WARNING! The Z mass cannot be changed in alpha_s evolution of LHAPDF!"<<endl;
}

void fastNLOLHAPDF::SetQMass(int pdgid, double mq) {
   warn["SetQMass"]<<"WARNING! The quark masses cannot be changed in alpha_s evolution of LHAPDF!"<<endl;
}

void fastNLOLHAPDF::SetNFlavor(int nflavor) {
   warn["SetNFlavor"]<<"WARNING! The no. of active flavors cannot be changed in alpha_s evolution of LHAPDF!"<<endl;
}

void fastNLOLHAPDF::SetNLoop(int nloop) {
   warn["SetNLoop"]<<"WARNING! The no. of loops cannot be changed in alpha_s evolution of LHAPDF!"<<endl;
}

void fastNLOLHAPDF::SetAlphasMz(double AlphasMz, bool ReCalcCrossSection) {
   warn["SetAlphasMz"]<<"WARNING! alpha_s(M_Z) cannot be changed in alpha_s evolution of LHAPDF!"<<endl;
}

void fastNLOLHAPDF::InitEvolveAlphas() {
   // For LHAPDF do nothing
}

double fastNLOLHAPDF::GetQMass(int pdgid) const {
   if (pdgid < 1 || pdgid > 6 ) {
      error["GetQMass"]<<"PDG code out of quark range 1-6! Aborted\n";
      exit(1);
   }
   return LHAPDF::getQMass(pdgid);
}

int fastNLOLHAPDF::GetNLoop() const {
   return (LHAPDF::getOrderAlphaS() + 1);
}

int fastNLOLHAPDF::GetNFlavor() const {
   return (LHAPDF::getNf());
}

double fastNLOLHAPDF::GetAlphasMz(double Q) {
   return LHAPDF::alphasPDF(Q);
}
