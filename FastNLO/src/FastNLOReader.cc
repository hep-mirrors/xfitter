// Author: Daniel Britzger
// DESY, 23/07/2011

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

//#include "FastNLOReaderVersion.cc"
#ifndef FNLO_VERSION
#define FNLO_VERSION    "2.1.0"
#define FNLO_SVNREV     "XXXX"
#define FNLO_AUTHORS    "D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch"
#define FNLO_WEBPAGE    "http://projects.hepforge.org/fastnlo"
#define FNLO_AUTHORSv14 "T. Kluge, K. Rabbertz, M. Wobisch"
#define FNLO_QUOTEv14   "hep-ph/0609285"
#define FNLO_AUTHORSv2  "D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch"
#define FNLO_QUOTEv2    "arXiv:1109.1310"
#endif

#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
//#include "fastnlo/FastNLOReader.h"
#include "FastNLOReader.h"

using namespace std;

//______________________________________________________________________________


// some names for nice output
const string FastNLOReader::fContrName[20] = {
   "Fixed order calculation", "Threshold corrections", "Electroweak corrections", "Non-perturbative corrections",
   "Undefined", "Undefined", "Undefined", "Undefined", "Undefined", "Undefined", "Undefined",
   "Quark compositeness", "ADD-LED", "TeV 1-ED", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown"
};
const string FastNLOReader::fOrdName[4][4] = { { "LO",     "NLO",    "NNLO"   , "N3LO"    },
   { "1-loop", "2-loop", "3-loop" , "4-loop"  },
   { "Undef" , "Undef" , "Undef"  , "Undef"   },
   { "LO MC" , "NLO MC", "NNLO MC", "N3LO MC" }
};
const string FastNLOReader::fNSDep[6] = {"v2.0","v2.0","v2.0","v2.1","v2.2","v2.2"};
bool FastNLOReader::WelcomeOnce = false;

//______________________________________________________________________________

FastNLOReader::FastNLOReader(string filename) : PrimalScream("FastNLOReader") {
   debug["FastNLOReader"]<<"New FastNLOReader reading filename="<<filename<<endl;
   BlockB_Data          = NULL;
   BlockB_LO_Ref        = NULL;
   BlockB_NLO_Ref       = NULL;
   fUnits               = fastNLO::kPublicationUnits;
   fMuRFunc             = fastNLO::kScale1;
   fMuFFunc             = fastNLO::kScale1;
   fPDFSuccess          = false;
   fAlphasCached        = 0.;
   fPDFCached           = 0.;
   SetFilename(filename);
   if (!WelcomeOnce) PrintWelcomeMessage();
}


//______________________________________________________________________________


FastNLOReader::~FastNLOReader(void) {
   for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
      if (!BBlocksSMCalc[j].empty()) {
         for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
            delete BBlocksSMCalc[j][i];
         }
         BBlocksSMCalc.clear();
      }
   }
}


//______________________________________________________________________________



void FastNLOReader::PrintWelcomeMessage() {
   //---  Initialization for nice printing
   const string CSEPS = " ##################################################################################\n";
   const string LSEPS = " #---------------------------------------------------------------------------------\n";

   char fnlo[100];
   //sprintf(fnlo,"27[0;31mfast27[0;34mNLO\033[0m",27,0,31,27,0,34);
   sprintf(fnlo,"%c[%d;%dmfast%c[%d;%dmNLO\033[0m",27,0,31,27,0,34);
   char package_version[100]    = FNLO_VERSION;
   char svnrev[100]             = FNLO_SVNREV;
   char authors[500]            = FNLO_AUTHORS;
   char webpage[500]    = FNLO_WEBPAGE;
   char authorsv14[200] = FNLO_AUTHORSv14;
   char quotev14[200]   = FNLO_QUOTEv14;
   char authorsv2[200]  = FNLO_AUTHORS;
   char quotev2[200]    = FNLO_QUOTEv2;

   shout>>"\n";
   shout>>""<<CSEPS;
   shout<<"\n";
   shout<<" "<<fnlo<<"_reader"<<endl;
   shout<<" Version "<<package_version<<"_"<<svnrev<<endl;
   shout<<"\n";
   shout<<" C++ program to read fastNLO v2 tables and"<<endl;
   shout<<" derive QCD cross sections using PDFs, e.g. from LHAPDF"<<endl;
   shout<<"\n";
   shout>>""<<LSEPS;
   shout<<"\n";
   shout<<" Copyright © 2011,2012 "<<fnlo<<" Collaboration"<<endl;
   shout<<" "<<authors<<endl;
   shout<<"\n";
   shout>>" # This program is free software: you can redistribute it and/or modify"<<endl;
   shout>>" # it under the terms of the GNU General Public License as published by"<<endl;
   shout>>" # the Free Software Foundation, either version 3 of the License, or"<<endl;
   shout>>" # (at your option) any later version."<<endl;
   shout>>" #\n";
   shout>>" # This program is distributed in the hope that it will be useful,"<<endl;
   shout>>" # but WITHOUT ANY WARRANTY; without even the implied warranty of"<<endl;
   shout>>" # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the"<<endl;
   shout>>" # GNU General Public License for more details."<<endl;
   shout>>" #\n";
   shout>>" # You should have received a copy of the GNU General Public License"<<endl;
   shout>>" # along with this program. If not, see <http://www.gnu.org/licenses/>."<<endl;
   shout>>" #\n";
   shout>>""<<LSEPS;
   shout>>" #\n";
   shout<<" The projects web page can be found at:"<<endl;
   shout<<"   "<<webpage<<endl;
   shout<<"\n";
   shout<<" If you use this code, please cite:"<<endl;
   shout<<"   "<<authorsv14<<", "<<quotev14<<endl;
   shout<<"   "<<authorsv2<<", "<<quotev2<<endl;
   shout<<"\n";
   shout>>""<<CSEPS;
   WelcomeOnce = true;

}


//______________________________________________________________________________



void FastNLOReader::SetFilename(string filename) {
   debug["SetFilename"]<<"New filename="<<filename<<endl;
   ffilename    = filename;
   Init();
}


//______________________________________________________________________________



void FastNLOReader::Init() {
   debug["Init"]<<endl;
   ReadTable();
   //int iprint = 2;
   //PrintFastNLOTableConstants(iprint);
   InitScalevariation();
}


//______________________________________________________________________________



void FastNLOReader::InitScalevariation() {
   debug["InitScalevariation"]<<endl;
   fScaleFacMuR  = 1.;
   fScaleFacMuF  = 1.;
   fScalevar     = -1;

   if (!GetIsFlexibleScaleTable()) {
      for (int iscls=0; iscls<GetNScaleVariations(); iscls++) {
         const double muFac = BBlocksSMCalc[0][1]->ScaleFac[0][iscls];
         if (abs(muFac-1.0) < 1.e-7) {
            SetScaleVariation(iscls,true);
            break;
         }
      }
      if (fScalevar == -1) {
         error["InitScalevariation"]<<"Could not found scale variation with scale factor 1.0. Exiting.\n";
         exit(1);
      }
   } else {
      // this is a MuVar table. You can vary mu_f and mu_r independently by any factor
      // and you can choose the functional form of mu_f and mu_r as functions of
      // scale1 and scale1 (called partly scaleQ2 and scalePt).

      if (BBlocksSMCalc[0][0]->ScaleDescript[0].size() <0) {
         warn["InitScalevariation"]<<"No scaledescription available.\n";
         SetFunctionalForm(kScale1 , kMuR);
         SetFunctionalForm(kScale1 , kMuF);
         return;
      }

      // ---- DIS ---- //
      if (BBlocksSMCalc[0][0]->NPDFDim == 0) {
         SetFunctionalForm(kQuadraticMean , kMuR);
         SetFunctionalForm(kScale1 , kMuF);
      }
      // ---- HHC --- //
      else if (BBlocksSMCalc[0][0]->NPDFDim == 1) {
         SetFunctionalForm(kScale1 , kMuR);
         SetFunctionalForm(kScale1 , kMuF);
      } else {
         error<<"Unknown process.\n";
         exit(1);
      }
   }
}


//______________________________________________________________________________



double FastNLOReader::CalcMu(FastNLOReader::EMuX kMuX , double scale1, double scale2, double scalefac) {
   //
   //  Calculate the scales with the defined function and the
   //  corresponding prefactor.
   //

   if (kMuX == kMuR && fScaleFacMuR != scalefac) error<<"Sth. went wrong with the scales.\n";
   if (kMuX == kMuF && fScaleFacMuF != scalefac) error<<"Sth. went wrong with the scales.\n";

   EScaleFunctionalForm Func = (kMuX == kMuR) ? fMuRFunc : fMuFFunc;

   double mu = 0;

   if (Func == fastNLO::kScale1)            mu      = scale1;
   else if (Func == fastNLO::kScale2)            mu      = scale2;
   else if (Func == fastNLO::kQuadraticSum)      mu      = FuncMixedOver1(scale1,scale2);
   else if (Func == fastNLO::kQuadraticMean)     mu      = FuncMixedOver2(scale1,scale2);
   else if (Func == fastNLO::kQuadraticSumOver4) mu      = FuncMixedOver4(scale1,scale2);
   else if (Func == fastNLO::kLinearMean)        mu      = FuncLinearMean(scale1,scale2);
   else if (Func == fastNLO::kLinearSum)         mu      = FuncLinearSum(scale1,scale2);
   else if (Func == fastNLO::kScaleMax)          mu      = FuncMax(scale1,scale2);
   else if (Func == fastNLO::kScaleMin)          mu      = FuncMin(scale1,scale2);
   else if (Func == fastNLO::kExpProd2)          mu      = FuncExpProd2(scale1,scale2);
   else if (Func == fastNLO::kExtern)           mu      = (kMuX==kMuR) ? (*Fct_MuR)(scale1,scale2) : (*Fct_MuF)(scale1,scale2);
   else error["CalcMu"]<<"Could not identify functional form for scales calculation.\n";

   return scalefac * mu;

}


//______________________________________________________________________________
double FastNLOReader::FuncMixedOver1(double scale1 , double scale2) {
   return (sqrt((pow(scale1,2) + pow(scale2,2))  / 1.));
}

//______________________________________________________________________________
double FastNLOReader::FuncMixedOver2(double scale1 , double scale2) {
   return (sqrt((pow(scale1,2) + pow(scale2,2))  / 2.));
}

//______________________________________________________________________________
double FastNLOReader::FuncMixedOver4(double scale1 , double scale2) {
   return (sqrt((pow(scale1,2) + pow(scale2,2))  / 4.));
}

//______________________________________________________________________________
double FastNLOReader::FuncLinearMean(double scale1 , double scale2) {
   return (scale1 + scale2) / 2.;
}

//______________________________________________________________________________
double FastNLOReader::FuncLinearSum(double scale1 , double scale2) {
   return scale1 + scale2;
}

//______________________________________________________________________________
double FastNLOReader::FuncMax(double scale1 , double scale2) {
   if (scale1 > scale2) return scale1;
   else return scale2;
}

//______________________________________________________________________________
double FastNLOReader::FuncMin(double scale1 , double scale2) {
   if (scale1 < scale2) return scale1;
   else return scale2;
}

//______________________________________________________________________________
double FastNLOReader::FuncExpProd2(double scale1 , double scale2) {
   return (scale1 * exp(0.3*scale2));
}




//______________________________________________________________________________



double FastNLOReader::SetScaleVariation(int scalevar , bool FirstCall) {
   debug["SetScaleVariation"]<<"Setting to scalevar="<<scalevar<<endl;
   // ------------------------------------------------
   //   Set the scalevariation factor for determining the
   //   'theory'-error. Usually, you have tables stored with
   //   factors of 0.5, 1 and 2 times the nominal scale.
   //     corresponding to:
   //     scalevar -> scalefactor
   //        '0'   ->   1.00
   //        '1'   ->   0.50
   //        '2'   ->   2.00
   //   This method returns the scalefactor correspoding to
   //   the chosen 'scalevar'.
   // ------------------------------------------------

   if (GetIsFlexibleScaleTable()) {
      info["SetScaleVariation"]<<"This is a flexible-scale table. No Scalevariation tables available!"<<endl;
      man<<"You can choose freely (within reason) a factorization scale factor. Your Scalevar has to be '0'.\n";
      man<<"Please use SetScaleFacMuR(double) and SetScaleFacMuF(double) to set scale factors.\n";
      return 0;
   }

   // Check for maximal scale variation of all rel. and active SM calcs
   int scalevarmax = GetNScaleVariations();
   if (scalevar >= scalevarmax) {
      warn["SetScaleVariation"]<<"This table has only "<<scalevarmax<<" scale variation(s) stored!"<<endl;
      man<<"For the currently active contributions. You wanted to access the non-existing number "<<scalevar<<endl;
      man<<"Using '0' instead."<<endl;;
      fScalevar = 0;
      return B_NLO()->ScaleFac[0][0];
   }

   fScalevar     = scalevar;
   fScaleFacMuF  = B_NLO()->ScaleFac[0][fScalevar];
   info["SetScaleVariation"]
         <<"Selecting MuF table according to a multiplicative scale factor of the factorization scale of "
         <<fScaleFacMuF<<" times the nominal scale."<<endl;

   if (!BBlocksSMCalc[kThresholdCorrection].empty()) {
      bool lkth = false;
      for (unsigned int i = 0 ; i <BBlocksSMCalc[kThresholdCorrection].size() ; i++) {
         if (bUseSMCalc[kThresholdCorrection][i]) {
            lkth = true;
         }
      }
      if (lkth && abs(fScaleFacMuR-fScaleFacMuF) > DBL_MIN) {
         fScaleFacMuR = fScaleFacMuF;
         warn["SetScaleVariation."]<<"Threshold corrections do not allow variations of the renormalization scale!"<<endl;
         man<<"The scale factor for MuR has been set equal to the one for MuF = "<<fScaleFacMuF<<endl;
         man<<"Either select a different simultaneous scale variation i, if possible, via FastNLOReader::SetScaleVariation(i)"<<endl;
         man<<"or deactivate first all threshold corrections using FastNLOReader::SetContributionON(kTresholdCorrections,Id,false)."<<endl;
      }
   }

   return B_NLO()->ScaleFac[0][fScalevar];
}



//______________________________________________________________________________



void FastNLOReader::SetFunctionalForm(EScaleFunctionalForm func , FastNLOReader::EMuX MuX) {
   //
   //  For MuVar tables this method sets the functional form of
   //  the renormalization or the factorization scale.
   //     func:  Choose a pre-defined function
   //     kMuX:  is it for mu_r or for mu_f ?
   //

   if (!GetIsFlexibleScaleTable()) {
      warn<<"This is not a flexible-scale table. SetFunctionalForm cannot be used.\n";
      return;
   }

   // ---- setting scale ---- //
   if (MuX == kMuR) fMuRFunc = func;
   else fMuFFunc = func;


   // ---- cross check ---- //
   if (func == kScale2 || func == kQuadraticSum ||  func == kQuadraticMean || func == kQuadraticSumOver4
         || func == kLinearMean || func == kLinearSum  ||  func == kScaleMax|| func == kScaleMin) {
      if (BBlocksSMCalc[0][1]->ScaleNodeScale2[0].size() <= 3) {
         error<<"There is no second scale variable available in this table. Using fastNLO::kScale1 only.\n";
         SetFunctionalForm(kScale1,MuX);
      }
      for (int i=0; i<NObsBin; i++) {
         if (BBlocksSMCalc[0][1]->ScaleNodeScale2[i].size() < 4) {
            warn<<"Scale2 has only very little nodes (n="<<BBlocksSMCalc[0][0]->ScaleNodeScale2[i].size()<<") in bin "<<i<<".\n";
         }
      }
   }
   PrintScaleSettings(MuX);
}


//______________________________________________________________________________


void FastNLOReader::SetMuRFunctionalForm(EScaleFunctionalForm func) {
   SetFunctionalForm(func,kMuR);
}


//______________________________________________________________________________


void FastNLOReader::SetMuFFunctionalForm(EScaleFunctionalForm func) {
   SetFunctionalForm(func,kMuF);
}



//______________________________________________________________________________



bool FastNLOReader::SetScaleFactorsMuRMuF(double xmur, double xmuf) {
   debug["SetScaleFactorsMuRMuF"];
   //
   // Set renormalization and factorization scale factors simultaneously for scale variations in all v2 tables.
   // You have to ReFill your cache!
   // This is done automatically, but if you want to do it by yourself set ReFillCache = false.
   //
   // The function aborts the whole program if non-sensical scale factors < 1.E-6 are requested.
   // The function returns true if the requested scale factors can be used with the available table.
   // If threshold corrections are selected, only xmur / xmuf = 1. is allowed.
   //    If this is not the case, xmur and xmuf are unchanged,
   //    a warning is printed and the function returns false!
   // If it is NOT a flexibleScaleTable and there is no scalevar-table for xmuf,
   //    xmur and xmuf are unchanged, a warning is printed and the function returns false!

   // Check whether xmur and xmuf are positive and at least larger than 1.E-6
   if (xmur < 1.E-6 || xmuf < 1.E-6) {
      error<<"Selected scale factors too small ( < 1.E-6 )! Ignoring call."<<endl;
      return false;
   }

   // Check whether threshold corrections exist and are activated
   if (!B_ThC()) {
      bool lkth = false;
      //for (vector<bool>::const_iterator it = bUseSMCalc[kThresholdCorrection].begin(); it!=bUseSMCalc[kThresholdCorrection].end(); ++it) lkth+=(*it);
      for (unsigned int i = 0 ; i <BBlocksSMCalc[kThresholdCorrection].size() ; i++) {
         if (bUseSMCalc[kThresholdCorrection][i]) {
            lkth = true;
            break;
         }
      }
      if (lkth && abs(xmur-xmuf) > DBL_MIN) {
         warn["SetScaleFactorsMuRMuFSetScaleFactorsMuRMuF"]
               <<"Threshold corrections do not allow different scale factors for MuR and MuF, nothing changed!\n";
         man<<"Please do only symmetric scale variation, i.e. xmur = xmuf,\n";
         man<<"or deactivate first all threshold corrections using\n";
         man<<"FastNLOReader::SetContributionON(kTresholdCorrections,Id,false).\n";
         return false;
      }
   }

   // Deal with factorization scale first
   // Check whether corresponding xmuf variation exists in case of v2.0 table
   if (!GetIsFlexibleScaleTable()) {
      //const double xmuf0 = B_NLO()->ScaleFac[0][fScalevar];
      const int ns = GetNScaleVariations();
      debug["SetScaleFactorsMuRMuF"]<<"NScaleVarMax="<<ns<<" must be >= than B->ScaleFac[0].size()="<<B_NLO()->ScaleFac[0].size()<<endl;
      int sf = -1;
      for (int is = 0 ; is<ns ; is++) {
         if (abs(B_NLO()->ScaleFac[0][is]-xmuf) < DBL_MIN) {
            sf = is;
            break;
         }
      }
      if (sf == -1) {
         warn["SetScaleFactorsMuRMuF"]<<"Could not find table with given mu_f scale factor of "<<xmuf<<". Nothing changed."<<endl;
         return false;
      }
      // set factorization scale
      SetScaleVariation(sf);
      // Now the renormalization scale
      fScaleFacMuR = xmur;
      PrintScaleSettings();
   } else {
      fScaleFacMuR = xmur;
      fScaleFacMuF = xmuf;
      PrintScaleSettings(kMuR);
      PrintScaleSettings(kMuF);
   }
   return true;
}

//______________________________________________________________________________

void FastNLOReader::PrintScaleSettings(FastNLOReader::EMuX MuX) {
   if (!GetIsFlexibleScaleTable()) {
      info<<"fastNLO. Renormalization scale chosen to be mu_r = "<<fScaleFacMuR<<" * "<<GetScaleDescription()<<endl;
      info<<"fastNLO. Factorization scale chosen to be   mu_f = "<<fScaleFacMuF<<" * "<<GetScaleDescription()<<endl;
   } else {
      // ---- prepare printout ---- //
      static const string sname[2] = {"Renormalization","Factorization"};
      static const string smu[2] = {"mu_r","  mu_f"};
      const int isc = MuX==kMuR?0:1;
      const double sfac = MuX==kMuR?fScaleFacMuR:fScaleFacMuF;
      EScaleFunctionalForm func = MuX==kMuR?fMuRFunc:fMuFFunc;
      char fname[100];
      switch (func) {
      case kScale1:
         sprintf(fname,"%s^2",GetScaleDescription(0).c_str());
         break;
      case kScale2:
         sprintf(fname,"%s^2",GetScaleDescription(1).c_str());
         break;
      case kQuadraticSum:
         sprintf(fname,"(%s^2 + %s^2)",GetScaleDescription(0).c_str(),GetScaleDescription(1).c_str());
         break;
      case kQuadraticMean:
         sprintf(fname,"(%s^2 + %s^2)/2",GetScaleDescription(0).c_str(),GetScaleDescription(1).c_str());
         break;
      case kQuadraticSumOver4:
         sprintf(fname,"(%s^2 + %s^2)/4",GetScaleDescription(0).c_str(),GetScaleDescription(1).c_str());
         break;
      case kLinearMean:
         sprintf(fname,"((%s+%s)/2)^2",GetScaleDescription(0).c_str(),GetScaleDescription(1).c_str());
         break;
      case kLinearSum:
         sprintf(fname,"(%s+%s)^2",GetScaleDescription(0).c_str(),GetScaleDescription(1).c_str());
         break;
      case kScaleMax:
         sprintf(fname,"max(%s^2,%s^2)",GetScaleDescription(0).c_str(),GetScaleDescription(1).c_str());
         break;
      case kScaleMin:
         sprintf(fname,"min(%s^2,%s^2)",GetScaleDescription(0).c_str(),GetScaleDescription(1).c_str());
         break;
      case kExpProd2:
         sprintf(fname,"(%s*exp(0.3*%s)^2)",GetScaleDescription(0).c_str(),GetScaleDescription(1).c_str());
         break;
      case kExtern:
         sprintf(fname,"f_ext(%s,%s)",GetScaleDescription(0).c_str(),GetScaleDescription(1).c_str());
         break;
      default:
         error<<"unknown scale choice.\n";
      }
      info<<"fastNLO. "<<sname[isc]<<" scale chosen to be "<<smu[isc]<<"^2 = "
          <<sfac<<"^2 * "<<fname<<endl;
   }
}


//______________________________________________________________________________



void FastNLOReader::ReadTable(void) {
   //
   // Read in the FastNLO Table
   //

   // Check whether file exists
   FILE* fp = fopen(ffilename.c_str(), "r");
   if (fp) {
      fclose(fp);
   } else {
      error<<"Table file '"<<ffilename.c_str()<<"' not found, exiting!"<<endl;
      exit(1);
   }

   // Open stream
   ifstream* instream = new ifstream(ffilename.c_str(),ios::in);

   // Read block A1 and A2
   ReadBlockA1(instream);
   ReadBlockA2(instream);

   // Initialize lists for BlockB's
   BBlocksSMCalc.resize(10);
   BBlocksNewPhys.resize(10);

   bUseSMCalc.resize(BBlocksSMCalc.size());
   bUseNewPhys.resize(BBlocksNewPhys.size());

   // Initialize BlockB's
   FastNLOBlockB* BlockB_LO   = NULL;
   FastNLOBlockB* BlockB_NLO  = NULL;
   FastNLOBlockB* BlockB_THC1 = NULL;
   FastNLOBlockB* BlockB_THC2 = NULL;
   FastNLOBlockB* BlockB_NPC1 = NULL;
   int nblocks = Ncontrib;

   for (int i=0; i<nblocks; i++) {
      // Read block
      FastNLOBlockB* blockb = new FastNLOBlockB("ReadingBlockB", NObsBin , instream);
      blockb->SetIc(i+1);
      char nbuf[400];
      if (blockb->IDataFlag && !BlockB_Data) {   // Data
         blockb->SetName("Data");
         BlockB_Data = blockb;
      } else if (blockb->IDataFlag && BlockB_Data) { // Data, but data already initalized
         warn["ReadTable"]<<"Only one data block is allowed, skipped!"<<endl;
      } else if (blockb->IRef == 1) { // Reference table, implemented only for LO or NLO
         if (blockb->IContrFlag1==1 && blockb->IContrFlag2==1) {
            if (blockb->NScaleDep < 3)   blockb->SetName("BlockB. LO Reference. v2.0.");
            if (blockb->NScaleDep >= 3)   blockb->SetName("BlockB. LO Reference. v2.1.");
            BlockB_LO_Ref           = blockb;
         } else if (blockb->IContrFlag1==1 && blockb->IContrFlag2==2) {
            if (blockb->NScaleDep < 3)   blockb->SetName("BlockB. NLO Reference. v2.0.");
            if (blockb->NScaleDep >= 3)   blockb->SetName("BlockB. NLO Reference. v2.1.");
            BlockB_NLO_Ref  = blockb;
         } else {
            error["ReadTable"]<<"Reference tables are only implemented for fixed order, stopped!\n";
            exit(1);
         }
      } else if (blockb->IRef == 0 && !blockb->IAddMultFlag) { // Additive corrections
         if (blockb->IContrFlag1==1) {   // Fixed order
            if (blockb->IContrFlag2==1) {  // LO
               sprintf(nbuf,"BlockB. %s %s %s",fContrName[blockb->IContrFlag1-1].c_str(),fOrdName[blockb->IContrFlag1-1][blockb->Npow-ILOord].c_str(),fNSDep[blockb->NScaleDep].c_str());
               blockb->SetName(nbuf);
               BlockB_LO  = blockb;
            } else if (blockb->IContrFlag2==2) { //NLO
               sprintf(nbuf,"BlockB. %s %s %s",fContrName[blockb->IContrFlag1-1].c_str(),fOrdName[blockb->IContrFlag1-1][blockb->Npow-ILOord].c_str(),fNSDep[blockb->NScaleDep].c_str());
               blockb->SetName(nbuf);
               BlockB_NLO = blockb;
            }
         } else if (blockb->IContrFlag1==2) { // Threshold corrections
            if (blockb->IContrFlag2==1) {  // 1-loop
               sprintf(nbuf,"BlockB. %s %s %s",fContrName[blockb->IContrFlag1-1].c_str(),fOrdName[blockb->IContrFlag1-1][blockb->Npow-ILOord].c_str(),fNSDep[blockb->NScaleDep].c_str());
               blockb->SetName(nbuf);
               BlockB_THC1 = blockb;
            } else if (blockb->IContrFlag2==2) { // 2-loop
               sprintf(nbuf,"BlockB. %s %s %s",fContrName[blockb->IContrFlag1-1].c_str(),fOrdName[blockb->IContrFlag1-1][blockb->Npow-ILOord].c_str(),fNSDep[blockb->NScaleDep].c_str());
               blockb->SetName(nbuf);
               BlockB_THC2 = blockb;
            } else {
               error["ReadTable"]<<"Threshold correction implemented only up to 2-loops, exiting!\n";
               exit(1);
            }
         }
         //        else if ( blockb->IContrFlag1>=3 ){ //
         // sprintf(nbuf,"BlockB. %s. %s. %s",fNPName[blockb->IContrFlag2].c_str(),fOrdName[blockb->Npow-ILOord].c_str(),fNSDep[blockb->NScaleDep].c_str());
         // blockb->SetName(nbuf);
         // BBlocksNewPhys[blockb->IContrFlag2].push_back(blockb);
         // bUseNewPhys[blockb->IContrFlag2].push_back(true);
         // }
         else {
            error["ReadTable"]<<"Further additive corrections not yet implemented, stopped!\n";
            exit(1);
         }
      } else if (blockb->IRef == 0 && blockb->IAddMultFlag) { // Multiplicative corrections
         if (blockb->IContrFlag1==4) {   // Non-perturbative corrections
            sprintf(nbuf,"BlockB. %s %s %s",fContrName[blockb->IContrFlag1].c_str(),fOrdName[blockb->IContrFlag1-1][blockb->IContrFlag2-1].c_str(),fNSDep[blockb->NScaleDep].c_str());
            blockb->SetName(nbuf);
            BlockB_NPC1 = blockb;
         }

      } else {
         error["ReadTable"]<<"Further multiplicative corrections not yet implemented, stopped!\n";
         exit(1);
      }
   }

   // Assign NPC, switch off by default
   if (BlockB_NPC1) {
      BBlocksSMCalc[BlockB_NPC1->IContrFlag1-1].push_back(BlockB_NPC1);
      bUseSMCalc[BlockB_NPC1->IContrFlag1-1].push_back(false);
   }

   // Assign THC, switch off by default
   if (BlockB_THC2) {
      BBlocksSMCalc[BlockB_THC2->IContrFlag1-1].push_back(BlockB_THC2);
      bUseSMCalc[BlockB_THC2->IContrFlag1-1].push_back(false);
   }
   if (BlockB_THC1) {
      BBlocksSMCalc[BlockB_THC1->IContrFlag1-1].push_back(BlockB_THC1);
      bUseSMCalc[BlockB_THC1->IContrFlag1-1].push_back(false);
   }

   // Assign fixed order calculations (LO must be [0]), switch on by default
   if (BlockB_LO)  {
      BBlocksSMCalc[0].push_back(BlockB_LO);
      bUseSMCalc[0].push_back(true);
   } else {
      error["ReadTable"]<<"Could not find any LO Calculation. Exiting!\n";
      exit(1);
   }
   if (BlockB_NLO) {
      BBlocksSMCalc[0].push_back(BlockB_NLO);
      bUseSMCalc[0].push_back(true);
   }

   // Some printout?
   // PrintTableInfo();

}

//______________________________________________________________________________



void FastNLOReader::PrintTableInfo(const int iprint) const {
   debug["PrintTableInfo"]<<"iprint="<<iprint<<endl;

   //---  Initialization for nice printing
   string CSEPS = "##################################################################################\n";
   string LSEPS = "#---------------------------------------------------------------------------------\n";
   printf("\n");
   printf(" %s",CSEPS.c_str());
   printf(" # Overview on contribution types and numbers contained in table:\n");
   printf(" %s",LSEPS.c_str());
   printf(" # Number of contributions: %2i\n",Ncontrib);

   unsigned int icnt = 0;
   for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
      if (!BBlocksSMCalc[j].empty()) {
         for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
            icnt++;
            cout << " # "<< "  No.: " << icnt << ", type: " << fContrName[j] <<", Id: " << i <<
                 ", order: " << BBlocksSMCalc[j][i]->CtrbDescript[0] <<
                 ", by: " << BBlocksSMCalc[j][i]->CodeDescript[0] << endl;
         }
         if (iprint > 0) {
            for (unsigned int k = 0 ; k<BBlocksSMCalc[j][0]->CodeDescript.size(); k++) {
               printf(" # \t\t%s\n",BBlocksSMCalc[j][0]->CodeDescript[k].c_str());
            }
            //BBlocksSMCalc[j][0]->Print(0,0);
         }
      }
   }

   for (unsigned int j = 0 ; j<BBlocksNewPhys.size() ; j++) {
      if (!BBlocksNewPhys[j].empty()) {
         for (unsigned int i = 0 ; i<BBlocksNewPhys[j].size() ; i++) {
            icnt++;
            cout << " # "<< "  No.: " << icnt << ", type: " << fContrName[j] <<", Id: " << i <<
                 ", order: " << BBlocksNewPhys[j][i]->CtrbDescript[0] <<
                 ", by: " << BBlocksNewPhys[j][i]->CodeDescript[0] << endl;
         }
         printf(" #   -> SM extensions can not be evaluated by this reader! Just skipping those ...\n");
      }
   }

   if (BlockB_Data) {
      icnt++;
      cout << " # "<< "  No.: " << icnt << ", type: Data, Id: 0" <<
         ", order: " << BlockB_Data->CtrbDescript[0] <<
         ", by: " << BlockB_Data->CodeDescript[0] << endl;
      if (iprint > 0) {
         for (unsigned int k = 0 ; k<BlockB_Data->CodeDescript.size(); k++) {
            printf(" * \t\t%s\n",BlockB_Data->CodeDescript[k].c_str());
         }
      }
   }
   printf(" %s",CSEPS.c_str());

}

//______________________________________________________________________________


void FastNLOReader::SetContributionON(ESMCalculation eCalc , unsigned int Id , bool SetOn) {
   // sanity check
   if (bUseSMCalc[eCalc].empty() || BBlocksSMCalc.empty()) {
      warn["SetContributionON"]
            <<"This contribution ("<<fContrName[eCalc]<<") does not exist in this table. Cannot switch it On/Off. Ignoring call.\n";
      return;
   }
   if (bUseSMCalc[eCalc].size() < Id || BBlocksSMCalc[eCalc].size() < Id || !BBlocksSMCalc[eCalc][Id]) {
      warn["SetContributionON"]
            <<"This Id="<<Id<<" does not exist for this contribtion. Cannot switch it On/Off. Ignoring call.\n";
      return;
   }

   info<<(SetOn?"Activating":"Deactivating")
       <<" contribution '"<<fContrName[eCalc]
       <<" with Id="<<Id<<endl;

   if (!bUseSMCalc[eCalc][Id] && SetOn) {
      debug["SetContributionON"]<<"Call FillAlphasCache for contribution eCalc="<<eCalc<<"\tId="<<Id<<endl;
      if (!BBlocksSMCalc[eCalc][Id]->IAddMultFlag) {
         if (!GetIsFlexibleScaleTable())
            FillAlphasCacheInBlockBv20(BBlocksSMCalc[eCalc][Id]);
         else
            FillAlphasCacheInBlockBv21(BBlocksSMCalc[eCalc][Id]);
      }
   }
   // set the new value
   bUseSMCalc[eCalc][Id] = SetOn;
}

//______________________________________________________________________________


void FastNLOReader::ReadBlockA1(istream *table) {
   debug["ReadBlockA1"]<<endl;
   //
   //  Read in information of Block A1 from file
   //

   table->peek();
   if (table->eof()) {
      error["ReadBlockA1"]<<"Cannot read from file.\n";
      return;
   }

   int key = 0;
   *table >> key;
   if (key != tablemagicno) {
      error["ReadBlockA1"]<<"At beginning of block found "<<key<<" instead of "<<tablemagicno<<endl;
      return;
   }
   *table >> Itabversion;
   if (Itabversion < 20000) {
      error["ReadBlockA1"]<<"This reader is only compatible with FastNLO v2.0 tables and higher. Exiting.\n";
      man<<"This FastNLO-table (file) is of version "<<Itabversion*1./10000.<<endl;;
      man<<"Please download a compatible reader from the website or use the APPL_grid interface.\n";
      exit(1);
   }
   *table >> ScenName;
   //--- Ncontrib: All contributions including additive, multiplicative or data
   //---           (Then: Nadd = Ncontrib - Nmult - Ndata)
   *table >> Ncontrib;
   *table >> Nmult;
   *table >> Ndata;
   *table >> NuserString;
   *table >> NuserInt;
   *table >> NuserFloat;
   *table >> Imachine;
   key=0;
   *table >> key;
   if (key != tablemagicno) {
      error["ReadBlockA1"]<<"At end of block found "<<key<<" instead of "<<tablemagicno<<endl;
      return;
   };
   // Put magic number back
   for (int i=0; i<(int)(log10((double)key)+1); i++) {
      table->unget();
   }
   if (debug.GetSpeak()) PrintBlockA1();

}


//______________________________________________________________________________


void FastNLOReader::ReadBlockA2(istream *table) {
   debug["ReadBlockA2"]<<endl;
   //
   //  Read in information which is called Block A2 from file
   //

   table->peek();
   if (table->eof()) {
      printf("FastNLOReader::Read: Cannot read from file.\n");
      return;
   }

   int key = 0;
   *table >> key;
   if (key != tablemagicno) {
      printf("FastNLOReader::Read: At beginning of block found %d instead of %d.\n",key,tablemagicno);
      return;
   };

   *table >> Ipublunits;
   int NScDescript;
   *table >> NScDescript;
   ScDescript.resize(NScDescript);
   char buffer[257];
   table->getline(buffer,256);
   for (int i=0; i<NScDescript; i++) {
      table->getline(buffer,256);
      ScDescript[i] = buffer;
      StripWhitespace(&ScDescript[i]);
   }

   *table >> Ecms;
   *table >> ILOord;
   *table >> NObsBin;
   *table >> NDim;
   DimLabel.resize(NDim);
   table->getline(buffer,256);
   for (int i=0; i<NDim; i++) {
      table->getline(buffer,256);
      DimLabel[i] = buffer;
      StripWhitespace(&DimLabel[i]);
   }

   IDiffBin.resize(NDim);
   for (int i=0; i<NDim; i++) {
      *table >>  IDiffBin[i];
   }
   LoBin.resize(NObsBin);
   UpBin.resize(NObsBin);
   // Set rapidity index also when reading a table
   RapIndex.push_back(0);
   int irap = 0;
   for (int i=0; i<NObsBin; i++) {
      LoBin[i].resize(NDim);
      UpBin[i].resize(NDim);
      for (int j=0; j<NDim; j++) {
         *table >>  LoBin[i][j];
         if (IDiffBin[j]==2) *table >>  UpBin[i][j];
      }
      debug << "iobs1: " << i << ", LoBin i: " << LoBin[i][1] << endl;
      if (i > 0) {
         if (LoBin[i][1] != LoBin[i-1][1]) {
            debug << "iobs2: " << i << ", LoBin i-1: " << LoBin[i-1][1] << ", LoBin i: " << LoBin[i][1] << endl;
            RapIndex.push_back(i);
            irap++;
            debug << "irap: " << irap << ", RapIndex: " << RapIndex[irap] << endl;
         }
      }
   }

   BinSize.resize(NObsBin);
   for (int i=0; i<NObsBin; i++) {
      *table >> BinSize[i];
   }

   *table >> INormFlag;
   if (INormFlag>1) {
      *table >> DenomTable;
   }
   if (INormFlag>0) {
      IDivLoPointer.resize(NObsBin);
      IDivUpPointer.resize(NObsBin);
      for (int i=0; i<NObsBin; i++) {
         *table >> IDivLoPointer[i];
         *table >> IDivUpPointer[i];
      }
   }

   key=0;
   *table >> key;
   if (key != tablemagicno) {
      error["Read"]<<"At end of block found "<<key<<" instead of "<<tablemagicno<<endl;;
      return;
   };
   // Put magic number back
   for (int i=0; i<(int)(log10((double)key)+1); i++) {
      table->unget();
   }
   if (debug.GetSpeak()) PrintBlockA2();
}



//______________________________________________________________________________



void FastNLOReader::PrintFastNLOTableConstants(const int iprint) const {
   //if ( debug.GetSpeak() ) iprint=10000;
   //
   // Define different levels of detail for printing out table content
   // The minimum (iprint = 0) just gives basic scenario information
   // including the employed code with references. The additional levels
   // are mostly for debugging purposes.
   // (Partially to be implemented!)
   //
   // iprint = 0: No additional printout
   //          1: Print Block A1 & A2 (A1, A2)
   //          2: Also print basic values of Block B (B0)
   //          3: Also print x nodes of Block B for each contribution (BX)
   //          4: Also print scale nodes of Block B for each contribution (BS)
   //          5: Also print sigma tilde of Block B (not implemented yet)

   //---  Initialization for nice printing
   string CSEPS = "##################################################################################\n";
   string LSEPS = "#---------------------------------------------------------------------------------\n";

   //
   // Print basic scenario information (always)
   //
   printf("\n");
   printf(" %s",CSEPS.c_str());
   printf(" # Information on fastNLO scenario: %s\n",ScenName.data());
   printf(" %s",LSEPS.c_str());
   printf(" # Description:\n");
   for (unsigned int i=0; i<ScDescript.size(); i++) {
      printf(" #   %s\n",ScDescript[i].data());
   }
   printf(" #\n");
   printf(" # Centre-of-mass energy Ecms: % -#10.4g GeV\n",Ecms);
   printf(" #\n");
   printf(" # Tot. no. of observable bins: %3i in %1i dimensions:\n",NObsBin,NDim);
   printf(" #\n");
   printf(" # No. of contributions: %1i\n",Ncontrib);
   for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
      if (!BBlocksSMCalc.empty()) {
         for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
            if (BBlocksSMCalc[j][i]) {
               int iContr = BBlocksSMCalc[j][i]->GetIc();
               BBlocksSMCalc[j][i]->Print(iContr,iprint);
            }
         }
      }
   }
   if (BlockB_Data) {
      int iContr = BlockB_Data->GetIc();
      BlockB_Data->Print(iContr,iprint);
   }

   //
   // Print additional info for debugging including internal variables as selected by iprint
   //
   if (iprint > 0) {
      PrintBlockA1();
      PrintBlockA2();
   }
   if (iprint > 1) {
      for (unsigned int j = 0; j<BBlocksSMCalc.size(); j++) {
         if (!BBlocksSMCalc.empty()) {
            for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
               if (BBlocksSMCalc[j][i]) {
                  int iContr = BBlocksSMCalc[j][i]->GetIc();
                  BBlocksSMCalc[j][i]->Print(iContr,iprint);
               }
            }
         }
      }
      if (BlockB_LO_Ref) {
         int iContr = BlockB_LO_Ref->GetIc();
         BlockB_LO_Ref->Print(iContr,iprint);
      }
      if (BlockB_NLO_Ref) {
         int iContr = BlockB_NLO_Ref->GetIc();
         BlockB_NLO_Ref->Print(iContr,iprint);
      }
   }

   printf(" #\n");
   printf(" %s",CSEPS.c_str());
}



//______________________________________________________________________________



void FastNLOReader::PrintBlockA1() const {
   printf("\n *****************************************\n");
   printf(" * fastNLO Table: Block A1\n");
   printf(" *****************************************\n");
   printf("  A1  ISep                              %10i\n",tablemagicno);
   printf("  A1  Itabversion                       %10i\n",Itabversion);
   printf("  A1  ScenName                          %s\n",ScenName.data());
   printf("  A1  Ncontrib                          %10i\n",Ncontrib);
   printf("  A1  Nmult                             %10i\n",Nmult);
   printf("  A1  Ndata                             %10i\n",Ndata);
   printf("  A1  NuserString                       %10i\n",NuserString);
   printf("  A1  NuserInt                          %10i\n",NuserInt);
   printf("  A1  NuserFloat                        %10i\n",NuserFloat);
   printf("  A1  Imachine                          %10i\n",Imachine);
   printf(" #########################################\n");
}



//______________________________________________________________________________



void FastNLOReader::PrintBlockA2() const {
   printf("\n *****************************************\n");
   printf(" * fastNLO Table: Block A2\n");
   printf(" *****************************************\n");
   printf("  A2  ISep                              %10i\n",tablemagicno);
   printf("  A2  IpublUnits                        %10i\n",Ipublunits);
   printf("  A2  NscDescript                       %10zi\n",ScDescript.size());
   for (unsigned int i=0; i<ScDescript.size(); i++) {
      printf("  A2    ScDescript(%1i)                   %s\n",i+1,ScDescript[i].data());
   }
   printf("  A2  Ecms                              % -#10.4g\n",Ecms);
   printf("  A2  ILOord                            %10i\n",ILOord);
   printf("  A2  NobsBin                           %10i\n",NObsBin);
   printf("  A2  NDim                              %10i\n",NDim);
   for (int i=0; i<NDim; i++) {
      printf("  A2    DimLabel(%1i)                     %s\n",i+1,DimLabel[i].data());
   }
   for (int i=0; i<NDim; i++) {
      printf("  A2    IDiffBin(%1i)                     %10i\n",i+1,IDiffBin[i]);
   }
   for (int i=0; i<NObsBin; i++) {
      for (int j=0; j<NDim; j++) {
         printf("  A2      LoBin(%3i,%1i)              % #10.4g\n", i+1,j+1,LoBin[i][j]);
         if (IDiffBin[j]==2)
            printf("  A2      UpBin(%3i,%1i)              % #10.4g\n", i+1,j+1,UpBin[i][j]);
      }
   }
   for (int i=0; i<NObsBin; i++) {
      printf("  A2    BinSize(%3i)                    % -#10.4g\n", i+1,BinSize[i]);
   }
   printf("  A2  INormFlag                         %10i\n",INormFlag);

   if (INormFlag>1) {
      printf("  A2  DenomTable                        %s\n",DenomTable.data());
   }
   if (INormFlag>0) {
      for (int i=0; i<NObsBin; i++) {
         printf("  A2    IDivLoPointer(%3i)              %10i\n",i+1,IDivLoPointer[i]);
         printf("  A2    IDivUpPointer(%3i)              %10i\n",i+1,IDivUpPointer[i]);
      }
   }
   printf(" #########################################\n");
}


//______________________________________________________________________________


void FastNLOReader::PrintCrossSections() const {
   //
   // Print Cross sections in NLO, k-factors and Reference table cross sections
   //

   //   if ( XSection.empty() )    CalcCrossSection();
   //   if ( XSectionRef.empty() && XSectionRef_s1.empty() )    CalcReferenceCrossSection();}

   vector < double > xs = XSection;

   printf(" *  \n");
   printf(" *  FastNLO Cross sections for\n");
   for (unsigned int i = 0 ; i < ScDescript.size() ; i++) {
      printf(" *     %s\n",ScDescript[i].c_str());
   }
   printf(" *  at sqrt(s) = %8.2f GeV\n", Ecms);
   printf(" *  \n");
   printf(" *  This is a %s-differential table in %s", ((NDim==1)?"single":"double"),DimLabel[0].c_str());
   if (NDim==2) printf(" and in %s",DimLabel[1].c_str());
   printf(".\n");
   printf(" *\n");

   string Aunits[16]  = { "[b] --   ","","","[mb] --  ","","","[mu b] --","","","[nb] --  ","","","[pb] --  ","","","[fb] --  "};
   string Nounits[16] = { " --      ","",""," --      ","",""," --      ","",""," --      ","",""," --      ","",""," --      "};
   string* unit = fUnits==kAbsoluteUnits ? Aunits : Nounits;


   if (NDim == 2) {
      double lobindim2 = -42;
      printf(" *  - Bin - |   ---  %5s  ---        -- XS-FNLO %s -- k-factor -- |\n",DimLabel[0].c_str(),unit[Ipublunits].c_str());
      printf(" *  --------------------------------------------------------------------\n");
      for (unsigned int i=0; i<xs.size(); i++) {
         if (LoBin[i][1] != lobindim2) {
            printf(" *                  ---->  from %9.3f to %9.3f in %s  <----\n",LoBin[i][1],UpBin[i][1],DimLabel[1].c_str());
            lobindim2 = LoBin[i][1];
         }
         printf(" *   %4.0f   | %9.3f - %9.3f       % 9.4e           % 5.2f      |\n",i*1.,LoBin[i][0],UpBin[i][0],xs[i],kFactor[i]);
      }
   }

   else {
      printf("   ---  %5s  ---        - Bin -       -- XS-FNLO --  \n",DimLabel[0].c_str());
      for (unsigned int i=0; i<xs.size(); i++) {
         printf("  %9.3f - %9.3f   %3.0f         % 9.4e\n",LoBin[i][0],UpBin[i][0],i*1.,xs[i]);
      }
   }
   printf(" *  --------------------------------------------------------------------\n");
}


//______________________________________________________________________________


void FastNLOReader::PrintFastNLODemo() {
   //
   // This method prints out cross sections for different scale
   // variation tables. Though it also changes the currently stored
   // settings of this instance!
   //

   info["PrintFastNLODemo"]<<"PrintFastNLODemo is changing settings (like scale choices) of this reader."<<endl;

   // If flexible-scale table, set MuR and MuF functional forms
   if (GetIsFlexibleScaleTable()) {
      SetMuRFunctionalForm(fastNLO::kScale1);
      SetMuFFunctionalForm(fastNLO::kScale1);
      //SetMuRFunctionalForm(fastNLO::kExpProd2);
      //SetMuRFunctionalForm(fastNLO::kExpProd2);
   }

   // Check on existence of LO and NLO (Id = -1 if not existing)
   int ilo   = ContrId(fastNLO::kFixedOrder, fastNLO::kLeading);
   int inlo  = ContrId(fastNLO::kFixedOrder, fastNLO::kNextToLeading);
   if (ilo < 0 || inlo < 0) {
      error["PrintFastNLODemo"]<<"LO and/or NLO not found, nothing to be done!\n";
      return;
   }
   // Check on existence of 2-loop threshold corrections
   int ithc2 = ContrId(fastNLO::kThresholdCorrection, fastNLO::kNextToLeading);
   // Switched off by default. Don't do scale variations. Not available for the moment.
   //  if ( ithc2 > -1 ) {
   //    SetContributionON( fastNLO::kThresholdCorrection, ithc2, false, false );
   //  }

   // Pre-define desired order of scale variations
   const int nxmu = 4;
   double xmu[nxmu] = {1.0, 0.25, 0.5, 2.0};
   int   ixmu[nxmu] = { -1,   -1,  -1,  -1};
   // Get number of available scale variations and check on available scale factors,
   // in particular for MuF; set pointers
   int nscls = GetNScaleVariations();
   // With threshold corrections, allow only default scale (0)
   if (ithc2 > -1) {
      nscls = 1;
   }
   for (int iscls=0; iscls<nscls; iscls++) {
      SetScaleVariation(iscls);
      double fxmu = fScaleFacMuF;
      for (int i=0; i<nxmu; i++) {
         if (abs(xmu[i]-fxmu) < 0.000001) {
            ixmu[i] = iscls;
         }
      }
   }

   // Loop over scales
   for (int iscls=0; iscls<nxmu; iscls++) {
      // First result is with NLO, LO result via division by K factor
      if (ixmu[iscls] > -1) {
         SetScaleVariation(ixmu[iscls]);
         CalcCrossSection();

         // Second result: Include threshold corrections for NLO if available
         vector < double > kthc;
         if (ithc2 > -1) {
            vector < double > stdk = kFactor;
            SetContributionON(fastNLO::kThresholdCorrection, ithc2, true);
            CalcCrossSection();
            kthc = kFactor;
            // Threshold K factor is NLO including 2-loop vs. NLO
            for (unsigned int i=0; i<kthc.size(); i++) {
               if (abs(kFactor[i]) > DBL_MIN) {
                  kthc[i] = kFactor[i]/stdk[i];
               } else {
                  kthc[i] = -1.;
               }
            }
            SetContributionON(fastNLO::kThresholdCorrection, ithc2, false);
         }

         PrintCrossSectionsDefault(kthc);
      }
   }
}


//______________________________________________________________________________


void FastNLOReader::PrintCrossSectionsDefault(const vector <double> kthc) const {
   //
   // Print observable binnning and cross sections at
   // LO, NLO and K factors like in Fortran Reader for comparison
   //

   // Some initialization
   const string CSEP41("#########################################");
   const string DSEP41("=========================================");
   const string SSEP41("-----------------------------------------");
   const string CSEP = CSEP41 + CSEP41 + CSEP41 + CSEP41;
   const string DSEP = DSEP41 + DSEP41 + DSEP41 + DSEP41;
   const string SSEP = SSEP41 + SSEP41 + SSEP41 + SSEP41;


   // Check on existence of 2-loop threshold corrections
   //const int ithc2 = kthc.empty() ? -1 : ContrId( fastNLO::kThresholdCorrection, fastNLO::kNextToLeading);
   const int ithc2 = kthc.empty() ? -1 : ContrId(kThresholdCorrection,kNextToLeading);

   cout << DSEP << endl;
   printf(" Cross Sections\n");
   if (!GetIsFlexibleScaleTable())
      printf(" The scale chosen here are: mu_f = % #6.3f * %s, and mu_r = % #6.3f * %s \n",fScaleFacMuF,GetScaleDescription().c_str(),fScaleFacMuR,GetScaleDescription().c_str());
   cout << SSEP << endl;

   if (NDim == 2) {

      // non-perturbative corrections (just first np correction)
      const int inpc1 = ContrId(fastNLO::kNonPerturbativeCorrection, fastNLO::kLeading);
      const vector < double > knpc = inpc1>-1 ? BBlocksSMCalc[3][0]->fact : vector<double>(NObsBin);

      string header0 = "  IObs  Bin Size IODim1 ";
      string header1 = "   IODim2 ";
      string header2 = " LO cross section   NLO cross section   K NLO";
      if (ithc2>-1)header2 += "     K THC";
      if (inpc1>-1)header2 += "     K NPC";
      unsigned int NDimBins[NDim];
      printf("%s [ %-12s ] %s [  %-12s  ] %s\n",
             header0.c_str(),DimLabel[0].c_str(),header1.c_str(),DimLabel[1].c_str(),header2.c_str());
      cout << SSEP << endl;
      for (int i=0; i<NObsBin; i++) {
         for (int j=0; j<NDim; j++) {
            if (i==0)                                  NDimBins[j] = 1;
            else if (LoBin[i-1][j] < LoBin[i][j])       NDimBins[j]++;
            else if (LoBin[i][j] < LoBin[i-1][j])       NDimBins[j] = 1;
         }
         if (ithc2<0 && inpc1<0) {
            printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F",
                   i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                   NDimBins[1],LoBin[i][1],UpBin[i][1],XSection_LO[i],XSection[i],kFactor[i]);
         } else if (inpc1<0 && ithc2 != -1) {
            printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F %#9.5F",
                   i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                   NDimBins[1],LoBin[i][1],UpBin[i][1],XSection_LO[i],XSection[i],kFactor[i],kthc[i]);
         } else if (inpc1>-1 && ithc2 == -1) {
            printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F %#9.5F",
                   i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                   NDimBins[1],LoBin[i][1],UpBin[i][1],XSection_LO[i],XSection[i],kFactor[i],knpc[i]);
         } else {
            printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F %#9.5F %#9.5F",
                   i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                   NDimBins[1],LoBin[i][1],UpBin[i][1],XSection_LO[i],XSection[i],kFactor[i],kthc[i],knpc[i]);
         }
         printf("\n");
      }
   } else {
      warn["PrintCrossSectionsDefault"]<<"Print out optimized for two dimensions. No output for "<<NDim<<" dimensions."<<endl;
   }

}


//______________________________________________________________________________


void FastNLOReader::PrintCrossSectionsData() const {
   //
   // Print summary information on data table if available
   //

   if (!BlockB_Data) {
      info["PrintCrossSectionsData"]<<"No data table found. Nothing to do."<<endl;
      return;
   }

   vector < double > xs = BlockB_Data->Value;

   string CSEP41("#########################################");
   string DSEP41("=========================================");
   string SSEP41("-----------------------------------------");
   string CSEP = CSEP41 + CSEP41 + CSEP41 + CSEP41;
   string DSEP = DSEP41 + DSEP41 + DSEP41 + DSEP41;
   string SSEP = SSEP41 + SSEP41 + SSEP41 + SSEP41;

   cout << DSEP << endl;
   printf(" Measurement\n");
   //--- Additionally print out data description
   // for ( unsigned int k = 0 ; k<BlockB_Data->CodeDescript.size();k++ ) {
   //   printf( "\t\t%s\n",BlockB_Data->CodeDescript[k].c_str());
   // }
   cout << SSEP << endl;

   if (NDim == 2) {
      string header[3] = { "  IObs  Bin Size IODim1 ",
                           "   IODim2 ",
                           "   X Section QSumUnc.+ QSumUnc.- QSumCor.+ QSumCor.-\n"
                         };
      string label[2] = { "[ " + DimLabel[0] + "     ]", "[  " + DimLabel[1] + "           ]"};
      unsigned int NDimBins[NDim];
      printf("%s %s %s %s %s",header[0].c_str(),label[0].c_str(),header[1].c_str(),label[1].c_str(),header[2].c_str());
      //--- Additionally print out descriptions on all uncertainties
      // for ( int iUnco = 0 ; iUnco<BlockB_Data->Nuncorrel ; iUnco++ ){
      //   printf("%30s   ",BlockB_Data->UncDescr[iUnco].c_str());
      // }
      // for ( int iCorr = 0 ; iCorr<BlockB_Data->Ncorrel ; iCorr++ ){
      //   printf("%30s   ",BlockB_Data->CorDescr[iCorr].c_str());
      // }
      // cout<<endl;

      cout << SSEP << endl;
      for (unsigned int i=0; i<xs.size(); i++) {
         for (int j=0; j<NDim; j++) {
            if (i==0) {
               NDimBins[j] = 1;
            } else if (LoBin[i-1][j] < LoBin[i][j]) {
               NDimBins[j]++;
            } else if (LoBin[i][j] < LoBin[i-1][j]) {
               NDimBins[j] = 1;
            }
         }
         printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E   % -#10.3E",
                i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                NDimBins[1],LoBin[i][1],UpBin[i][1],xs[i]);

         //--- For now add sources quadratically
         double DXsUncorLo = 0.;
         double DXsUncorHi = 0.;
         for (int iUnco = 0 ; iUnco<BlockB_Data->Nuncorrel ; iUnco++) {
            DXsUncorLo += BlockB_Data->UncorLo[i][iUnco] * BlockB_Data->UncorLo[i][iUnco];
            DXsUncorHi += BlockB_Data->UncorHi[i][iUnco] * BlockB_Data->UncorHi[i][iUnco];
         }
         DXsUncorLo = -sqrt(DXsUncorLo);
         DXsUncorHi = +sqrt(DXsUncorHi);
         printf(" %+8.2E %+8.2E",DXsUncorHi,DXsUncorLo);

         double DXsCorrLo = 0.;
         double DXsCorrHi = 0.;
         for (int iCorr = 0 ; iCorr<BlockB_Data->Ncorrel ; iCorr++) {
            DXsCorrLo += BlockB_Data->CorrLo[i][iCorr] * BlockB_Data->CorrLo[i][iCorr];
            DXsCorrHi += BlockB_Data->CorrHi[i][iCorr] * BlockB_Data->CorrHi[i][iCorr];
         }
         DXsCorrLo = -sqrt(DXsCorrLo);
         DXsCorrHi = +sqrt(DXsCorrHi);
         printf(" %+8.2E %+8.2E",DXsCorrHi,DXsCorrLo);
         cout << endl;
      }
   } else {
      printf("WARNING! Sorry, currently no output for more than 2 dimensions!\n");
   }
}


//______________________________________________________________________________


void FastNLOReader::PrintCrossSectionsWithReference() {
   //
   //  Print Cross sections in NLO, k-factors and Reference table cross sections
   //
   //  Please mention, that the reference cross section can be easily deviating
   //  more than 20% (scales, pdfs, alpha_s, etc...). This does not mean that
   //  the table is wrong!
   //


   if (XSection.empty()) {
      CalcCrossSection();
   }
   if (XSectionRef.empty() && XSectionRef_s1.empty()) {
      CalcReferenceCrossSection();
   }

   vector < double > xs = XSection;
   vector < double > xsref;

   if (GetIsFlexibleScaleTable()) {
      if (fMuFFunc == kScale1 && fMuRFunc == kScale1)   {
         printf(" *  FastNLOReader::PrintCrossSectionsWithReference. Info. Taking reference cross sections 's1'\n");
         xsref = XSectionRef_s1;
      } else if (fMuFFunc == kScale2 && fMuRFunc == kScale2) {
         printf(" *  FastNLOReader::PrintCrossSectionsWithReference. Info. Taking reference cross sections 's2'\n");
         xsref = XSectionRef_s2;
      } else if (fMuFFunc == kQuadraticMean && fMuRFunc == kQuadraticMean) {
         printf(" *  FastNLOReader::PrintCrossSectionsWithReference. Info. Taking reference cross sections 'mixed'\n");
         xsref = XSectionRefMixed;
      } else {
         xsref = XSectionRefMixed;
         printf(" *  FastNLOReader::PrintCrossSectionsWithReference. Info. Taking reference cross sections 'mixed'\n");
      }
   } else xsref = XSectionRef;


   printf(" *  \n");
   printf(" *  FastNLO Cross sections for\n");
   for (unsigned int i = 0 ; i < ScDescript.size() ; i++) {
      printf(" *     %s\n",ScDescript[i].c_str());
   }
   printf(" *  at sqrt(s) = %8.2f GeV\n", Ecms);
   printf(" *  \n");
   printf(" *  This is a %s-differential table in %s", ((NDim==1)?"single":"double"),DimLabel[0].c_str());
   if (NDim==2) printf(" and %s",DimLabel[1].c_str());
   printf(" *  \n");
   printf(" *  Please mention, that the reference cross section can easily deviating up to more\n *  than 20%% due to different scale choices, alhpa_s value/evolution, PDFs, etc.");
   printf(" *  This does not mean, that this FastNLO table is wrong!\n\n");
   printf(" *  There are three reference cross sections stored for different scale choices.\n");
   printf(" *  If you have choosen mu_r=mu_f=%s, or mu_r=mu_f=%s or mu_r=mu_f=sqrt((%s^2+%s^2)/2), then you access automatically the corresponding reference cross section.\n",
          BBlocksSMCalc[0][1]->ScaleDescript[0][0].c_str(),BBlocksSMCalc[0][1]->ScaleDescript[0][1].c_str(),BBlocksSMCalc[0][1]->ScaleDescript[0][0].c_str(),BBlocksSMCalc[0][1]->ScaleDescript[0][1].c_str());
   printf(" *  In any other case your reference cross section is calculated using mu_r=mu_f=sqrt((%s^2+%s^2)/2).\n",BBlocksSMCalc[0][1]->ScaleDescript[0][0].c_str(),BBlocksSMCalc[0][1]->ScaleDescript[0][1].c_str());
   printf(" *  To be fully consistent with the nlojet++ reference cross section, you also have to adjust alpha_s and the alpha_s evolution accordingly.\n\n");

   printf("\n");
   printf(" *\n");

   string Aunits[16]  = { "[b] --   ","","","[mb] --  ","","","[mu b] --","","","[nb] --  ","","","[pb] --  ","","","[fb] --  "};
   string Nounits[16] = { " --      ","",""," --      ","",""," --      ","",""," --      ","",""," --      ","",""," --      "};
   string* unit = fUnits==kAbsoluteUnits ? Aunits : Nounits;

   if (NDim == 2) {
      double lobindim2 = -321312;
      printf(" *  - Bin - |   ---  %5s  ---        -- XS-FNLO %s -- k-factor -- |  -- XS-ref (NLOJET++) --    Diff [%%]\n",DimLabel[0].c_str(),unit[Ipublunits].c_str());
      printf(" *  -----------------------------------------------------------------------------------------------------------\n");
      for (unsigned int i=0; i<xs.size(); i++) {
         if (LoBin[i][1] != lobindim2) {
            printf(" *                    ---->  from %9.3f to %9.3f in %s  <----\n",LoBin[i][1],UpBin[i][1],DimLabel[1].c_str());
            lobindim2 = LoBin[i][1];
         }
         printf(" *   %4.0f   | %9.3f - %9.3f      % 9.4e           % 5.3f      |     % 9.4e            % 5.4f\n",i*1.,LoBin[i][0],UpBin[i][0],xs[i],kFactor[i],xsref[i],(xs[i]-xsref[i])/xsref[i]*100.);
      }
   }

   else {
      printf("FastNLOReader::PrintCrossSections( ). Info. Single differential printing of cross sections not yet nicely implemented.\n");
      printf("   ---  %s  ---        - Bin -    -- XS-FNLO  --       -- XS-ref (NLOJET++) --    Diff [%%]\n",DimLabel[0].c_str());
      for (unsigned int i=0; i<xs.size(); i++) {
         printf("  %9.3f - %9.3f   %3.0f         % 9.4e           % 9.4e          % 5.4f\n",LoBin[i][0],UpBin[i][0],i*1.,xs[i],xsref[i],(xs[i]-xsref[i])/xsref[i]*100.);
      }
   }
   printf(" *  ------------------------------------------------------------------------------------------------------------\n");
}


//______________________________________________________________________________


int FastNLOReader::GetNScaleVariations() const {
   if (GetIsFlexibleScaleTable()) {
      info["GetNScaleVariations"]<<"This is a 'flexible-scale' table, therefore you can choose all desired scale variations."<<endl;
      return 0;
   }

   // Check for maximal scale variation of all rel. and active SM calcs
   // Assume a maximum of 10!
   unsigned int scalevarmax = 10;
   for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
      if (!BBlocksSMCalc.empty()) {
         for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
            // Do not check pQCD LO or mult. corrections
            if (bUseSMCalc[j][i] && !BBlocksSMCalc[j][i]->IAddMultFlag &&
                  !(j==kFixedOrder && i==kLeading)) {
               if (BBlocksSMCalc[j][i]->Nscalevar[0] < (int)scalevarmax) {
                  scalevarmax = BBlocksSMCalc[j][i]->Nscalevar[0];
               }
            }
         }
      }
   }
   debug["GetNScaleVariations"]<<"Found "<<scalevarmax<<" scale variations."<<endl;
   return scalevarmax;
}


//______________________________________________________________________________


vector < double > FastNLOReader::GetScaleFactors() const {
   if (GetIsFlexibleScaleTable()) {
      info["GetScaleFactors"]<<"This is a 'flexible scale table', therefore you can choose all desired scale variations."<<endl;
      return vector<double>();
   }
   return BBlocksSMCalc[0][1]->ScaleFac[0];
}


//______________________________________________________________________________


vector < double > FastNLOReader::GetCrossSection() {
   // Get fast calculated NLO cross section
   if (XSection.empty()) CalcCrossSection();
   return XSection;
}


//______________________________________________________________________________


vector < double > FastNLOReader::GetKFactors() {
   // Get ratio of fast calculated NLO to LO cross section
   if (XSection.empty()) CalcCrossSection();
   return kFactor;
}

//______________________________________________________________________________

vector < double > FastNLOReader::GetQScales(int irelord) {
   // Get XSection weighted Q scale in bin
   if (XSection.empty()) CalcCrossSection();
   if (irelord == 0) {
      return QScale_LO;
   } else {
      return QScale;
   }
}


//______________________________________________________________________________


vector < double > FastNLOReader::GetReferenceCrossSection() {
   // Get reference cross section from direct nlojet++ calculation
   if (XSectionRef.empty() && XSectionRef_s1.empty()) {
      CalcReferenceCrossSection();
   }
   if (GetIsFlexibleScaleTable()) {
      if (fMuFFunc == kScale1 && fMuRFunc == kScale1)                   return XSectionRef_s1;
      else if (fMuFFunc == kScale2 && fMuRFunc == kScale2)              return XSectionRef_s2;
      else if (fMuFFunc == kQuadraticMean && fMuRFunc == kQuadraticMean)return XSectionRefMixed;
      else return XSectionRefMixed;
   } else return XSectionRef; // XSectionRef from BlockB-Ref
}


//______________________________________________________________________________


void FastNLOReader::CalcReferenceCrossSection() {
   //
   //  Initialize the internal arrays for the reference cross
   //  sections with the information from the FastNLO file
   //

   XSectionRef.clear();
   XSectionRef.resize(NObsBin);

   XSectionRefMixed.clear();
   XSectionRef_s1.clear();
   XSectionRef_s2.clear();
   XSectionRefMixed.resize(NObsBin);
   XSectionRef_s1.resize(NObsBin);
   XSectionRef_s2.resize(NObsBin);

   if (BlockB_LO_Ref && BlockB_NLO_Ref) {

      for (int i=0; i<NObsBin; i++) {
         double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
         for (int l=0; l<BlockB_LO_Ref->NSubproc; l++) {
            XSectionRef[i] +=  BlockB_LO_Ref->SigmaTilde[i][0][0][0][l] * unit; // no scalevariations in LO tables
         }
         for (int l=0; l<BlockB_NLO_Ref->NSubproc; l++) {
            XSectionRef[i] +=  BlockB_NLO_Ref->SigmaTilde[i][fScalevar][0][0][l] * unit;
         }
      }

   }

   if (GetIsFlexibleScaleTable()) {
      for (int i=0; i<NObsBin; i++) {
         double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
         for (int n=0; n<BBlocksSMCalc[0][1]->NSubproc; n++) {
            XSectionRefMixed[i]             += BBlocksSMCalc[0][0] ->SigmaRefMixed[i][n] * unit;
            XSectionRef_s1[i]               += BBlocksSMCalc[0][0] ->SigmaRef_s1[i][n] * unit;
            XSectionRef_s2[i]               += BBlocksSMCalc[0][0] ->SigmaRef_s2[i][n] * unit;
         }
         for (int n=0; n<BBlocksSMCalc[0][1]->NSubproc; n++) {
            XSectionRefMixed[i]             += BBlocksSMCalc[0][1]->SigmaRefMixed[i][n] * unit;
            XSectionRef_s1[i]               += BBlocksSMCalc[0][1]->SigmaRef_s1[i][n] * unit;
            XSectionRef_s2[i]               += BBlocksSMCalc[0][1]->SigmaRef_s2[i][n] * unit;
         }
      }
   }

   if (!GetIsFlexibleScaleTable() && (BlockB_NLO_Ref==NULL))
      warn["CalcReferenceCrossSection"]<<"No reference cross sections available.\n";

}


//______________________________________________________________________________
bool FastNLOReader::PrepareCache() {
   // check pdf cache
   const double PDFcks = CalcNewPDFChecksum();
   if (fPDFCached==0. || (fPDFCached!=0. && fabs(PDFcks/fPDFCached -1.) > 1.e-7)) {
      debug["PrepareCache"]<<"Need to refill PDFCache, since PDFCecksum="<<PDFcks<<" and fPDFCached="<<fPDFCached<<endl;
      FillPDFCache(PDFcks);
   } else {
      debug["PrepareCache"]<<"No need to refill PDFCache."<<endl;
   }
   // check pdf cache
   if (!fPDFSuccess) {
      error["PrepareCache"]<<"Cannot calculate cross sections. PDF has not been initalized successfully."<<endl;
      return false;
   }

   // check alpha_s cache
   const double asref = CalcReferenceAlphas();
   if (fAlphasCached == 0. || fAlphasCached != asref) {
      debug["PrepareCache"]<<"Need to refill AlphasCache, since fAlphasCached="<<fAlphasCached<<endl;
      FillAlphasCache();
   }
   // do we now have an alphas?
   if (fAlphasCached==0. || fAlphasCached != asref) {
      error["PrepareCache"]<<"Filling of alpha_s cache failed. fAlphasCached="<<fAlphasCached<<"\tasref="<<asref<<endl;
      return false;
   }
   return true;
}


//______________________________________________________________________________


void FastNLOReader::CalcCrossSection() {
   debug["CalcCrossSection"]<<endl;
   //
   //  Initialize the internal arrays with the NLO cross sections
   //  with the information from the FastNLO file, the pdf and
   //  the defined alpha_s
   //

   XSection_LO.clear();
   XSection.clear();
   XSection_LO.resize(NObsBin);
   XSection.resize(NObsBin);
   kFactor.clear();
   kFactor.resize(NObsBin);
   QScale_LO.clear();
   QScale.clear();
   QScale_LO.resize(NObsBin);
   QScale.resize(NObsBin);

   // handle alpha_s and PDF Cache
   bool CacheOK = PrepareCache();
   if (!CacheOK) {
      error["CalcCrossSection"]<<"Caching failed. Cannot calculate cross sections."<<endl;
      return;
   }

   // perturbative (additive) contributions
   for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
      if (!BBlocksSMCalc[j].empty()) {
         for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
            if (bUseSMCalc[j][i] && !BBlocksSMCalc[j][i]->IAddMultFlag) {
               if (!GetIsFlexibleScaleTable())
                  CalcCrossSectionv20(BBlocksSMCalc[j][i]);
               else
                  CalcCrossSectionv21(BBlocksSMCalc[j][i]);
            }
         }
      }
   }


   // contributions from the a-posteriori scale variation
   if (!GetIsFlexibleScaleTable()) {
      if (abs(fScaleFacMuR-B_NLO()->ScaleFac[0][fScalevar]) > DBL_MIN) {
         CalcAposterioriScaleVariation();
      }
   }

   // calculate LO cross sections
   if (!GetIsFlexibleScaleTable())
      CalcCrossSectionv20(B_LO(),true);
   else
      CalcCrossSectionv21(B_LO(),true);


   // non-perturbative corrections (multiplicative corrections)
   for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
      if (!BBlocksSMCalc[j].empty()) {
         for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
            if (bUseSMCalc[j][i] && BBlocksSMCalc[j][i]->IAddMultFlag && BBlocksSMCalc[j][i]->IContrFlag1 == 4) {
               if (BBlocksSMCalc[j][i]->IContrFlag2 == 1) {
                  debug["CalcCrossSection"]<<"Adding multiplicative non-perturbative correction."<<endl;
                  for (int iB=0; iB<NObsBin; iB++) {
                     XSection[iB] *= BBlocksSMCalc[j][i]->fact[iB];
                     //            XSection_LO[iB]     *= BBlocksSMCalc[j][i]->fact[iB];
                  }
               }
            }
         }
      }
   }

   // ---- k-factor calculation ---- //
   debug["CalcCrossSection"]<<"Calculate k-factors: xs/xs_LO"<<endl;
   for (int i=0; i<NObsBin; i++) {
      kFactor[i] = XSection[i] / XSection_LO[i];
   }

   // ---- Q-scale calculation ---- //
   debug["CalcCrossSection"]<<"Calculate Q-scales: xsQ/xs"<<endl;
   for (int i=0; i<NObsBin; i++) {
      QScale_LO[i] = QScale_LO[i]/XSection_LO[i];
      QScale[i]    = QScale[i]/XSection[i];
   }
}

//______________________________________________________________________________


void FastNLOReader::CalcAposterioriScaleVariation() {
   double scalefac       = fScaleFacMuR/fScaleFacMuF;
   debug["CalcAposterioriScaleVariation"]<<"scalefac="<<scalefac<<endl;
   vector<double>* XS    = &XSection;
   vector<double>* QS    = &QScale;
   const double n     = B_LO()->Npow;
   const double L     = std::log(scalefac);
   const double beta0 = (11.*3.-2.*5)/3.;
   for (int i=0; i<NObsBin; i++) {
      int nxmax = B_LO()->GetNxmax(i);
      double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
      for (int j=0; j<B_LO()->GetTotalScalenodes(); j++) {
         double asnp1 = pow(B_LO()->AlphasTwoPi_v20[i][j],(n+1)/n);//as^n+1
         for (int k=0; k<nxmax; k++) {
            for (int l=0; l<B_LO()->NSubproc; l++) {
               double clo  = B_LO()->SigmaTilde[i][0][j][k][l] *  B_LO()->PdfLc[i][j][k][l] * unit;
               double xsci = asnp1 * clo * n * L * beta0;
               double mur  = fScaleFacMuR * B_LO()->ScaleNode[i][0][0][j];
               XS->at(i) +=  xsci;
               QS->at(i) +=  xsci*mur;
            }
         }
      }
   }
}


//______________________________________________________________________________


void FastNLOReader::CalcCrossSectionv21(FastNLOBlockB* B , bool IsLO) {
   debug["CalcCrossSectionv21"]<<"B->fname="<<B->fname<<"\tNpow="<<B->Npow<<"\tIsLO="<<IsLO<<endl;
   //
   //  Cross section calculation for DIS and HHC tables in v2.1 format
   //

   vector<double>* XS = IsLO ? &XSection_LO : &XSection;
   vector<double>* QS = IsLO ? &QScale_LO : &QScale;
   B->fact.resize(NObsBin);
   for (int i=0; i<NObsBin; i++) {
      B->fact[i]=0;
      int nxmax = B->GetNxmax(i);
      double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
      for (unsigned int jS1=0; jS1<B->ScaleNodeScale1[i].size(); jS1++) {
         for (unsigned int kS2=0; kS2<B->ScaleNodeScale2[i].size(); kS2++) {
            double Q2   = B->ScaleNodeScale1[i][jS1]*B->ScaleNodeScale1[i][jS1];
            double mur      = CalcMu(kMuR , B->ScaleNodeScale1[i][jS1] ,  B->ScaleNodeScale2[i][kS2] , fScaleFacMuR);
            double muf      = CalcMu(kMuF , B->ScaleNodeScale1[i][jS1] ,  B->ScaleNodeScale2[i][kS2] , fScaleFacMuF);
            double mur2 = pow(mur,2);
            double muf2 = pow(muf,2);
            for (int x=0; x<nxmax; x++) {
               for (int n=0; n<B->NSubproc; n++) {
                  double as   = B->AlphasTwoPi[i][jS1][kS2];
                  double pdflc        = B->PdfLcMuVar[i][x][jS1][kS2][n];
                  if (pdflc == 0.) continue;
                  double fac  = as * pdflc * unit;
                  double xsci =  B->SigmaTildeMuIndep[i][x][jS1][kS2][n] *                  fac;
                  if (B->Npow!=ILOord) {
                     xsci             += B->SigmaTildeMuFDep [i][x][jS1][kS2][n] * std::log(muf2) * fac;
                     xsci             += B->SigmaTildeMuRDep [i][x][jS1][kS2][n] * std::log(mur2) * fac;
                     if (BBlocksSMCalc[0][0]->IPDFdef1 == 2) {   // DIS tables use log(mu/Q2) instead of log(mu)
                        xsci -= B->SigmaTildeMuFDep [i][x][jS1][kS2][n] * std::log(Q2) * fac;
                        xsci -= B->SigmaTildeMuRDep [i][x][jS1][kS2][n] * std::log(Q2) * fac;
                     }
                  }
                  XS->at(i)   += xsci;
                  B->fact[i]  += xsci;
                  QS->at(i)   += xsci*mur;
               }
            }
         }
      }
   }
}

//______________________________________________________________________________


void FastNLOReader::CalcCrossSectionv20(FastNLOBlockB* B , bool IsLO) {
   debug["CalcCrossSectionv20"]<<"B->fname="<<B->fname<<"\tNpow="<<B->Npow<<"\tIsLO="<<IsLO<<endl;
   //
   //  Cross section calculation in v2.0 format
   //

   int scaleVar          = B->Npow == ILOord ? 0 : fScalevar;
   vector<double>* XS    = IsLO ? &XSection_LO : &XSection;
   vector<double>* QS    = IsLO ? &QScale_LO : &QScale;
   B->fact.resize(NObsBin);
   for (int i=0; i<NObsBin; i++) {
      B->fact[i] = 0;
      int nxmax = B->GetNxmax(i);
      double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
      for (int j=0; j<B->GetTotalScalenodes(); j++) {
         for (int k=0; k<nxmax; k++) {
            for (int l=0; l<B->NSubproc; l++) {
               double xsci     = B->SigmaTilde[i][scaleVar][j][k][l] *  B->AlphasTwoPi_v20[i][j]  *  B->PdfLc[i][j][k][l] * unit;
               double scalefac = fScaleFacMuR/B->ScaleFac[0][scaleVar];
               double mur      = scalefac * B->ScaleNode[i][0][scaleVar][j];
               XS->at(i)      +=  xsci;
               B->fact[i]     +=  xsci;
               QS->at(i)      +=  xsci*mur;
            }
         }
      }
   }
}


//______________________________________________________________________________


void FastNLOReader::SetUnits(EUnits Unit) {
   if (fUnits != Unit) {
      fUnits  = Unit;
      //CalcCrossSection();
   } else {
      // nothing todo
   }
}

//______________________________________________________________________________


void FastNLOReader::FillAlphasCache() {
   debug["FillAlphasCache"]<<endl;
   //
   //  Fill the internal alpha_s cache.
   //  This is usally called automatically. Only if you
   //  make use of ReFillCache==false options, you have
   //  to take care of this filling by yourself.
   //

   // check if the alpha_s value is somehow reasonable
   debug["FillAlphasCache"]<<"Sanity check!"<<endl;
   TestAlphas();

   // is there a need for a recalclation?
   const double asNew = CalcReferenceAlphas();
   if (asNew == fAlphasCached) {
      debug["FillAlphasCache"]<<"No need for a refilling of AlphasCache. asNew==fAlphasCached="<<asNew<<endl;
   } else {
      fAlphasCached = asNew;
      for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
         if (!BBlocksSMCalc.empty()) {
            for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
               // Check that this contribution type j and no. i should actually be used
               // Otherwise deactivation of e.g. threshold corr. is not respected here
               if (bUseSMCalc[j][i] && !BBlocksSMCalc[j][i]->IAddMultFlag) {
                  if (!GetIsFlexibleScaleTable()) {
                     FillAlphasCacheInBlockBv20(BBlocksSMCalc[j][i]);
                  } else if (GetIsFlexibleScaleTable()) {
                     FillAlphasCacheInBlockBv21(BBlocksSMCalc[j][i]);
                  }
               }
            }
         }
      }
   }
   return ;
}


//______________________________________________________________________________


void FastNLOReader::FillAlphasCacheInBlockBv20(FastNLOBlockB* B) {
   //
   //  Internal method for filling alpha_s cache
   //

   int scaleVar          = B->Npow == ILOord ? 0 : fScalevar;
   // Sanity check that scaleVar is in allowed range
   // For thresh. corr. can otherwise lead to inf and then segfault!
   if (scaleVar >= GetNScaleVariations()) {
      error<<"Trying to refresh  cache for non-existing scale variation no. "<<scaleVar<<" while only "<<GetNScaleVariations()<<" exist in total. Aborted."<<endl;
      exit(1);
   }
   double scalefac       = fScaleFacMuR/B->ScaleFac[0][scaleVar];
   debug["FillAlphasCacheInBlockBv20"]<<"scalefac="<<scalefac<<"\tscaleVar="<<scaleVar<<endl;

   for (int i=0; i<NObsBin; i++) {
      for (int j=0; j<B->GetTotalScalenodes(); j++) {
         double mur        = scalefac * B->ScaleNode[i][0][scaleVar][j];
         double as         = CalcAlphas(mur);
         B->AlphasTwoPi_v20[i][j] = pow(as/TWOPI , B->Npow);
      }
   }
}


//______________________________________________________________________________


void FastNLOReader::FillAlphasCacheInBlockBv21(FastNLOBlockB* B) {
   //
   //  Internal method for filling alpha_s cache
   //

   for (int i=0; i<NObsBin; i++) {
      for (unsigned int jS1=0; jS1<B->ScaleNodeScale1[i].size(); jS1++) {
         for (unsigned int kS2=0; kS2<B->ScaleNodeScale2[i].size(); kS2++) {
            double mur              = CalcMu(kMuR , BBlocksSMCalc[0][0]->ScaleNodeScale1[i][jS1] ,  BBlocksSMCalc[0][0]->ScaleNodeScale2[i][kS2] , fScaleFacMuR);
            double as               = CalcAlphas(mur);
            double alphastwopi      = pow(as/TWOPI, B->Npow);
            B->AlphasTwoPi[i][jS1][kS2] = alphastwopi;
         }
      }
   }
}


//______________________________________________________________________________


double FastNLOReader::CalcAlphas(double Q) {
   //
   //  Internal method for calculating the alpha_s(mu)
   //
   return EvolveAlphas(Q);
}


//______________________________________________________________________________


double FastNLOReader::CalcReferenceAlphas() {
   double mu = 0;
   if (GetIsFlexibleScaleTable()) {
      if (fMuRFunc==kExtern) mu = (*Fct_MuR)(91.,1.)*(fScaleFacMuR+0.1);
      else mu = 91.1876111111+(fMuRFunc*0.1)+(fScaleFacMuR);
   } else mu = 91.187611111115*(fScaleFacMuR+0.1)+fScalevar*0.1;
   double as = CalcAlphas(mu);
   if (isnan(as)) {
      error["CalcReferenceAlphas"]<<"Reference alphas is a 'nan' for scale mu="<<mu<<endl;
      //exit(1);
   }
   return as;
}


//______________________________________________________________________________


double FastNLOReader::CalcNewPDFChecksum() {
   // calculate a PDF checksum to
   // decide, whether PDF cache has to be refilled

   // init PDF and check success
   debug["CalcNewPDFChecksum"]<<"Call InitPDF() in user module."<<endl;
   fPDFSuccess = InitPDF();
   debug["CalcNewPDFChecksum"]<<"Return value InitPDF() = "<<fPDFSuccess<<endl;
   if (!fPDFSuccess) {
      warn["CalcPDFChecksum"]<<"PDF initialization failed. Please check PDF interface in your FastNLO user module."<<endl;
      return 0.;
   }

   // calculate checksum for some scales and flavors
   double muf = 0;
   if (GetIsFlexibleScaleTable()) {
      if (fMuFFunc==kExtern) muf = (*Fct_MuF)(91.,1.)/91.*(fScaleFacMuF+0.1) ;
      else muf = 91.1+0.1*fMuFFunc+fScaleFacMuF;
   } else muf=(fScaleFacMuF+0.1)+fScalevar*0.1;
   double cks = CalcChecksum(muf);
   return cks;
}


//______________________________________________________________________________

double FastNLOReader::CalcChecksum(double mufac) {
   debug["CalcChecksum"]<<"Calculate checksum of 13 flavors, 3 mu_f values, and 3 x-values, for scalefac="<<mufac<<endl;
   double cks = 0;
   vector<double> xfx(13);
   const double mf[3] = { 3,10,91.18};
   const double x[3] = {1.e-1,1.e-2,1.e-3};
   for (int jf = 0 ; jf<3 ; jf++) {
      double mu = mf[jf]* mufac;//(fScaleFacMuF+0.1)+fScalevar*0.1;
      for (int ix = 0 ; ix<3 ; ix++) {
         xfx = GetXFX(x[ix],mu);
         for (unsigned int fl = 0 ; fl<xfx.size() ; fl++) {
            cks+=xfx[fl];
         }
      }
   }
   debug["CalcChecksum"]<<"Calculated checksum = "<<cks<<endl;
   return cks;
}


//______________________________________________________________________________


bool FastNLOReader::TestAlphas() {
   const double as = CalcAlphas(91.18);
   if (as < 0.01 || as > 0.5) {
      warn["TestAlphas"]<<"The alphas value, returned by the user class seems to be unreasonably small/large."<<endl;
      warn["TestAlphas"]<<"The evolution code calculated alphas(Mz~91.18GeV) = "<<as<<endl;
      return false;
   }
   debug["TestAlphas"]<<"Sanity check of alpha_s(MZ=91.18) = "<<as<<endl;
   return true;
}


//______________________________________________________________________________


bool FastNLOReader::TestXFX() {
   vector<double> pdftest = GetXFX(1.e-2,10);
   if (pdftest.size() != 13) {
      error["TestXFX"]<<"The pdf array must have the size of 13 flavors.\n";
      return false;
   }
   // if ( pdftest[6] == 0. )printf("FastNLOReader. Warning. There seems to be no gluon in the pdf.\n");
   // double sum = 0;
   // for ( int i = 0 ; i<13 ; i++ ) sum+=fabs(pdftest[i]);
   // if ( sum== 0. ) printf("FastNLOReader. Error. All 13 pdf probabilities are 0. There might be sth. wrong in the pdf interface. Please check FastNLOUser::GetXFX().\n");
   for (int i = 0 ; i<13 ; i++) {
      if (pdftest[i] > 1.e10 || (pdftest[i] < 1.e-10 && pdftest[i] > 1.e-15)) {
         warn["TestXFX"]<<"The pdf probability of the "<<i<<"'s flavor seeems to be unreasonably large/small (pdf="<<pdftest[i]<<").\n";
      }
   }
   return true;
}

//______________________________________________________________________________



void FastNLOReader::FillPDFCache(double chksum) {
   debug["FillPDFCache"]<<"Passed chksum="<<chksum<<". Do not recalculate checksum (which calls InitPDF()) if chksum!=0."<<endl;
   //
   //  Fill the internal pdf cache.
   //  This function has to be called by the user, since the
   //  pdf parameters and evolutions are calculated externally.
   //

   // reset checknum
   // check if the alpha_s value is somehow reasonable
   double PDFnew = chksum;
   if (chksum == 0.) {
      debug["FillPDFCache"]<<"Calculate Checksum!"<<endl;
      PDFnew = CalcNewPDFChecksum();
      if (PDFnew==0.) {
         warn["FillPDFCache"]<<"PDF Checksum is zero."<<endl;
      }
      debug["FillPDFCache"]<<"PDF Checksum = "<<PDFnew<<endl;
   }

   // is there a need for a recalculation?
   if (fPDFCached != 0. && fabs(PDFnew/fPDFCached - 1.) < 1.e-7) {
      debug["FillPDFCache"]<<"No need for a refilling of PDFCache. fPDFCached=RefreshPDFChecksum()"<<PDFnew<<endl;
   } else {
      debug["FillPDFCache"]<<"Refilling PDF cache"<<endl;
      fPDFCached = PDFnew;

      // check (or not) if the pdf is somehow reasonable
      TestXFX();

      for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
         if (!BBlocksSMCalc.empty()) {
            for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
               // Check that this contribution type j and no. i should actually be used
               // Otherwise deactivation of e.g. threshold corr. is not respected here
               if (bUseSMCalc[j][i] && !BBlocksSMCalc[j][i]->IAddMultFlag) {
                  if (BBlocksSMCalc[j][i]->NScaleDim>1) {
                     error<<"WOW! NScaleDim>1! This is usually not the case!\n";
                     //scaleindex2 = 1; // If we use multiple scales, then mu_f is by convention the second scale -> index=1
                     //fScalevar2 = fScalevar % NfScalevar[1];
                  }

                  // linear: DIS-case
                  // ---- DIS ---- //
                  if (BBlocksSMCalc[j][i]->IPDFdef1 == 2) {
                     if (BBlocksSMCalc[j][i]->NPDFDim == 0) {
                        if (!GetIsFlexibleScaleTable()) FillBlockBPDFLCsDISv20(BBlocksSMCalc[j][i]);
                        else if (GetIsFlexibleScaleTable())   FillBlockBPDFLCsDISv21(BBlocksSMCalc[j][i]);
                     }
                  }
                  // ---- pp ---- //
                  else if (BBlocksSMCalc[j][i]->IPDFdef1 == 3) {
                     if (BBlocksSMCalc[j][i]->NPDFDim == 1) {
                        if (!GetIsFlexibleScaleTable()) FillBlockBPDFLCsHHCv20(BBlocksSMCalc[j][i]);
                        else                                            FillBlockBPDFLCsHHCv21(BBlocksSMCalc[j][i]);
                     } else {
                        error<<"Only half matrices for hh is implemented.\n";
                        exit(1);
                     }
                  } else {
                     error<<"Tables not yet implemented.\n";
                  }
               }
            }
         }
      }
   }
}


//______________________________________________________________________________


void FastNLOReader::FillBlockBPDFLCsDISv20(FastNLOBlockB* B) {
   debug["FillBlockBPDFLCsDISv20"]<<"BlockB = "<<B->fname<<endl;
   int scaleVar          = B->Npow == ILOord ? 0 : fScalevar;
   double scalefac       = B->ScaleFac[0][scaleVar] == fScaleFacMuF ? 1. : fScaleFacMuF;
   vector<double> xfx(13); // PDFs of all partons
   if (!GetIsFlexibleScaleTable()) {
      for (int i=0; i<NObsBin; i++) {
         int nxmax = B->GetNxmax(i);
         for (int j=0; j<B->Nscalenode[0]; j++) {
            for (int k=0; k<nxmax; k++) {
               double xp     = B->XNode1[i][k];
               double muf    = scalefac * B->ScaleNode[i][0][scaleVar][j];
               xfx = GetXFX(xp,muf);
               vector < double > buffer = CalcPDFLinearCombDIS(xfx , B->NSubproc);
               for (int l=0; l<B->NSubproc; l++) {
                  B->PdfLc[i][j][k][l] = buffer[l];
               }
            }
         }
      }
   }
}

//______________________________________________________________________________


void FastNLOReader::FillBlockBPDFLCsDISv21(FastNLOBlockB* B) {
   debug["FillBlockBPDFLCsDISv21"]<<"BlockB = "<<B->fname<<endl;

   if (B->PdfLcMuVar.empty()) {
      error<< "empty."<<endl;
      exit(1);
   }

   for (int i=0; i<NObsBin; i++) {
      // speed up! if mu_f is only dependent on one variable, we can safe the loop over the other one
      for (int x=0; x<B->GetNxmax(i); x++) {
         double xp = B->XNode1[i][x];
         if (fMuFFunc != kScale1 &&  fMuFFunc != kScale2) {   // that't the standard case!
            for (unsigned int jS1=0; jS1<B->ScaleNodeScale1[i].size(); jS1++) {
               for (unsigned int kS2=0; kS2<B->ScaleNodeScale2[i].size(); kS2++) {
                  double muf = CalcMu(kMuF , BBlocksSMCalc[0][0]->ScaleNodeScale1[i][jS1] ,  BBlocksSMCalc[0][0]->ScaleNodeScale2[i][kS2] , fScaleFacMuF);
                  B->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombDIS(GetXFX(xp,muf) , B->NSubproc);
               }
            }
         } else if (fMuFFunc == kScale2) { // speed up
            for (unsigned int kS2=0; kS2<B->ScaleNodeScale2[i].size(); kS2++) {
               double muf = CalcMu(kMuF , 0 ,  BBlocksSMCalc[0][0]->ScaleNodeScale2[i][kS2] , fScaleFacMuF);
               vector < double > buffer = CalcPDFLinearCombDIS(GetXFX(xp,muf) , B->NSubproc);
               for (unsigned int jS1=0; jS1<B->ScaleNodeScale1[i].size(); jS1++) {
                  B->PdfLcMuVar[i][x][jS1][kS2] = buffer;
               }
            }
         } else if (fMuFFunc == kScale1) { // speed up
            for (unsigned int jS1=0; jS1<B->ScaleNodeScale1[i].size(); jS1++) {
               double muf = CalcMu(kMuF , BBlocksSMCalc[0][0]->ScaleNodeScale1[i][jS1] , 0 , fScaleFacMuF);
               vector < double > buffer = CalcPDFLinearCombDIS(GetXFX(xp,muf) , B->NSubproc);
               for (unsigned int kS2=0; kS2<B->ScaleNodeScale2[i].size(); kS2++) {
                  B->PdfLcMuVar[i][x][jS1][kS2] = buffer;
               }
            }
         }
      }
   }
   debug["FillBlockBPDFLCsDISv21"]<<"BlockB = "<<B->fname<<" done." <<endl;
}


//______________________________________________________________________________


void FastNLOReader::FillBlockBPDFLCsHHCv20(FastNLOBlockB* B) {
   int scaleVar          = B->Npow == ILOord ? 0 : fScalevar;
   double scalefac       = fScaleFacMuF/B->ScaleFac[0][scaleVar];
   debug["FillBlockBPDFLCsHHCv20"]<<"scalefac="<<scalefac<<"\tBlockB="<<B<<endl;

   vector < vector < double > > xfx; // PDFs of all partons
   if (!GetIsFlexibleScaleTable()) {
      for (int i=0; i<NObsBin; i++) {
         int nxmax = B->GetNxmax(i);
         int nxbins1 = B->Nxtot1[i]; // number of columns in half matrix
         xfx.resize(nxbins1);
         for (int j=0; j<B->Nscalenode[0]; j++) {
            // determine all pdfs of hadron1
            for (int k=0; k<nxbins1; k++) {
               double xp     = B->XNode1[i][k];
               double muf    = scalefac * B->ScaleNode[i][0][scaleVar][j];
               xfx[k]        = GetXFX(xp,muf);
            }
            int x1bin = 0;
            int x2bin = 0;
            for (int k=0; k<nxmax; k++) {
               // ----- if pp ---- //
               if (B->NPDFPDG[0] == B->NPDFPDG[1]) {
                  B->PdfLc[i][j][k] = CalcPDFLinearCombHHC(xfx[x2bin], xfx[x1bin], B->NSubproc);
               }
               // ----- if ppbar ---- //
               else if (B->NPDFPDG[0] == -B->NPDFPDG[1]) {
                  vector < double > xfxbar(13);
                  for (unsigned int p = 0 ; p<13 ; p++) {
                     xfxbar[p] = xfx[x1bin][12-p];
                  }
                  B->PdfLc[i][j][k] = CalcPDFLinearCombHHC(xfx[x2bin], xfxbar, B->NSubproc);
               } else {
                  printf("FastNLOReader::FillBlockBPDFLCsHHCv20(). This is not pp, nor ppbar, nor pbarpbar!\n");
                  exit(1);
               }
               x1bin++;
               if (x1bin>x2bin) {
                  x1bin = 0;
                  x2bin++;
               }
            }
         }
      }
   }
}


//______________________________________________________________________________



void FastNLOReader::FillBlockBPDFLCsHHCv21(FastNLOBlockB* B) {
   debug["FillBlockBPDFLCsHHCv210"]<<"BlockB="<<B<<endl;
   if (B->PdfLcMuVar.empty()) {
      cout<< "empty."<<endl;
      exit(1);
   }
   vector < vector < double > > xfx; // PDFs of all partons
   for (int i=0; i<NObsBin; i++) {
      int nxmax = B->GetNxmax(i);
      int nxbins1 = B->Nxtot1[i]; // number of columns in half matrix
      xfx.resize(nxbins1);

      if (fMuFFunc != kScale1 &&  fMuFFunc != kScale2)  {   // that't the standard case!
         for (unsigned int jS1=0; jS1<B->ScaleNodeScale1[i].size(); jS1++) {
            for (unsigned int kS2=0; kS2<B->ScaleNodeScale2[i].size(); kS2++) {
               // determine all pdfs of hadron1
               for (int k=0; k<nxbins1; k++) {
                  double muf = CalcMu(kMuF , BBlocksSMCalc[0][0]->ScaleNodeScale1[i][jS1] ,  BBlocksSMCalc[0][0]->ScaleNodeScale2[i][kS2] , fScaleFacMuF);
                  double xp   = B->XNode1[i][k];
                  xfx[k] = GetXFX(xp,muf);
               }
               int x1bin = 0;
               int x2bin = 0;

               for (int x=0; x<nxmax; x++) {
                  // ----- if pp ---- //
                  if (B->NPDFPDG[0] == B->NPDFPDG[1]) {
                     B->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombHHC(xfx[x2bin], xfx[x1bin], B->NSubproc);
                  }
                  // ----- if ppbar ---- //
                  else if (B->NPDFPDG[0] == -B->NPDFPDG[1]) {
                     vector < double > xfxbar(13);
                     for (unsigned int p = 0 ; p<13 ; p++) {
                        xfxbar[p] = xfx[x1bin][12-p];
                     }
                     B->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombHHC(xfx[x2bin], xfxbar, B->NSubproc);
                  } else {
                     error["FillBlockBPDFLCsHHCv21"]<<"This is not pp, nor ppbar, nor pbarpbar!"<<endl;
                     exit(1);
                  }
                  x1bin++;
                  if (x1bin>x2bin) {
                     x1bin = 0;
                     x2bin++;
                  }
               }
            }
         }
      } else if (fMuFFunc == kScale2) {   // speed up
         for (unsigned int kS2=0; kS2<B->ScaleNodeScale2[i].size(); kS2++) {
            // determine all pdfs of hadron1
            for (int k=0; k<nxbins1; k++) {
               double muf = CalcMu(kMuF , 0 ,  BBlocksSMCalc[0][0]->ScaleNodeScale2[i][kS2] , fScaleFacMuF);
               double xp     = B->XNode1[i][k];
               xfx[k] = GetXFX(xp,muf);
            }
            for (unsigned int jS1=0; jS1<B->ScaleNodeScale1[i].size(); jS1++) {
               int x1bin = 0;
               int x2bin = 0;
               for (int x=0; x<nxmax; x++) {
                  // ----- if pp ---- //
                  if (B->NPDFPDG[0] == B->NPDFPDG[1]) {
                     B->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombHHC(xfx[x2bin], xfx[x1bin], B->NSubproc);
                  }
                  // ----- if ppbar ---- //
                  else if (B->NPDFPDG[0] == -B->NPDFPDG[1]) {
                     vector < double > xfxbar(13);
                     for (unsigned int p = 0 ; p<13 ; p++) {
                        xfxbar[p] = xfx[x1bin][12-p];
                     }
                     B->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombHHC(xfx[x2bin], xfxbar, B->NSubproc);
                  } else {
                     error["FillBlockBPDFLCsHHCv21"]<<"This is not pp, nor ppbar, nor pbarpbar!"<<endl;
                     exit(1);
                  }
                  x1bin++;
                  if (x1bin>x2bin) {
                     x1bin = 0;
                     x2bin++;
                  }
               }
            }
         }
      } else if (fMuFFunc == kScale1) {   // speed up
         for (unsigned int jS1=0; jS1<B->ScaleNodeScale1[i].size(); jS1++) {
            // determine all pdfs of hadron1
            for (int k=0; k<nxbins1; k++) {
               double muf = CalcMu(kMuF , BBlocksSMCalc[0][0]->ScaleNodeScale1[i][jS1] , 0 , fScaleFacMuF);
               double xp     = B->XNode1[i][k];
               xfx[k] = GetXFX(xp,muf);
            }
            for (unsigned int kS2=0; kS2<B->ScaleNodeScale2[i].size(); kS2++) {
               int x1bin = 0;
               int x2bin = 0;
               for (int x=0; x<nxmax; x++) {
                  // ----- if pp ---- //
                  if (B->NPDFPDG[0] == B->NPDFPDG[1]) {
                     B->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombHHC(xfx[x2bin], xfx[x1bin], B->NSubproc);
                  }
                  // ----- if ppbar ---- //
                  else if (B->NPDFPDG[0] == -B->NPDFPDG[1]) {
                     vector < double > xfxbar(13);
                     for (unsigned int p = 0 ; p<13 ; p++) {
                        xfxbar[p] = xfx[x1bin][12-p];
                     }
                     B->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombHHC(xfx[x2bin], xfxbar, B->NSubproc);
                  } else {
                     error["FillBlockBPDFLCsHHCv21"]<<"This is not pp, nor ppbar, nor pbarpbar!"<<endl;
                     exit(1);
                  }
                  x1bin++;
                  if (x1bin>x2bin) {
                     x1bin = 0;
                     x2bin++;
                  }
               }
            }
         }
      }
   }
}


//______________________________________________________________________________


vector<double> FastNLOReader::CalcPDFLinearCombDIS(vector<double> pdfx1 , int NSubproc) {
   //
   //  Internal method.
   //  CalcPDFLinearComb is used to calculate the
   //  linear combinations of different parton flavors
   //  according to the used subprocesses of your
   //  calculation/table.
   //

   vector < double > pdflc;
   pdflc.resize(3);
   pdflc[1] = pdfx1[6]; //gluon
   for (int l=0; l<13; l++) {
      double temp = (l==6 ? 0.0 : pdfx1[l]);
      if (!(l&1)) temp *= 4.;
      pdflc[0] += temp; // delta
   }
   pdflc[0] /= 9.;
   if (NSubproc>2) { // only from NLO
      for (int l=0; l<6; l++) {
         pdflc[2] += pdfx1[5-l] + pdfx1[l+7]; // sigma
      }
   }
   return pdflc;
}



//______________________________________________________________________________


vector<double> FastNLOReader::CalcPDFLinearCombHHC(vector<double> pdfx1 , vector<double> pdfx2 , int NSubproc) {
   //
   //  Internal method.
   //  CalcPDFLinearComb is used to calculate the
   //  linear combinations of different parton flavors
   //  according to the used subprocesses of your
   //  calculation/table.
   //

   double SumQ1  = 0;
   double SumQB1 = 0;
   double SumQ2  = 0;
   double SumQB2 = 0;
   vector <double> Q1(6);
   vector <double> QB1(6);
   vector <double> Q2(6);
   vector <double> QB2(6);
   for (int k = 0 ; k<6 ; k++) {
      Q1[k]  = pdfx1[k+7];  //! read 1st PDF at x1
      QB1[k] = pdfx1[5-k];
      SumQ1  += Q1[k];
      SumQB1 += QB1[k];
      Q2[k]  = pdfx2[k+7];//  ! read 2nd PDF at x2
      QB2[k] = pdfx2[5-k];
      SumQ2  += Q2[k];
      SumQB2 += QB2[k];
   }
   double G1     = pdfx1[6];
   double G2     = pdfx2[6];

   //   - compute S,A
   double S = 0;
   double A = 0;
   for (int k = 0 ; k<6 ; k++) {
      S += (Q1[k]*Q2[k]) + (QB1[k]*QB2[k]);
      A += (Q1[k]*QB2[k]) + (QB1[k]*Q2[k]);
   }

   //c   - compute seven combinations
   vector <double> H(7);
   H[0]  = G1*G2;
   H[1] = SumQ1*SumQ2 + SumQB1*SumQB2 - S;
   H[2] = S;
   H[3] = A;
   H[4] = SumQ1*SumQB2 + SumQB1*SumQ2 - A;
   H[5] = (SumQ1+SumQB1)*G2;
   H[6] = G1*(SumQ2+SumQB2);

   if (NSubproc == 6) {
      H[5] += H[6];
      H.resize(6);
   }
   return H;

}


//______________________________________________________________________________


void FastNLOReader::SetExternalFuncForMuR(double(*Func)(double,double)) {
   if (!GetIsFlexibleScaleTable()) {
      warn["SetExternalFuncForMuR"]<<"This is not a flexible-scale table and SetExternalFuncForMuR has no impact.\n";
      man<<"Please use a flexible-scale table, if you want to change your scale definition.\n";
      return;
   }

   Fct_MuR = Func;
   SetFunctionalForm(kExtern , kMuR);
   info["SetExternalFuncForMuR"]<<"Testing external function:"<<endl;
   info<<"Scale1 = 1 ,      Scale2 = 1        ->  mu = func(1,1)             = "<<(*Fct_MuR)(1,1)<<endl;
   info<<"Scale1 = 91.1876, Scale2 = 91.1876  ->  mu = func(91.1876,91.1876) = "<<(*Fct_MuR)(91.1876,91.1876)<<endl;
   info<<"Scale1 = 1,       Scale2 = 91.1876  ->  mu = func(1,91.1876)       = "<<(*Fct_MuR)(1,91.1876)<<endl;
   info<<"Scale1 = 91.1876, Scale2 = 1        ->  mu = func(91.1876,1)       = "<<(*Fct_MuR)(91.1876,1)<<endl;
}


//______________________________________________________________________________


void FastNLOReader::SetExternalFuncForMuF(double(*Func)(double,double)) {
   if (!GetIsFlexibleScaleTable()) {
      warn["SetExternalFuncForMuF"]<<"This is not a flexible-scale table and SetExternalFuncForMuF has no impact.\n";
      man<<"Please use a flexible-scale table, if you want to change your scale definition.\n";
      return;
   }

   Fct_MuF = Func;
   SetFunctionalForm(kExtern , kMuF);
   info["SetExternalFuncForMuF"]<<"Testing external function:"<<endl;
   info<<"Scale1 = 1 ,      Scale2 = 1        ->  mu = func(1,1)             = "<<(*Fct_MuF)(1,1)<<endl;
   info<<"Scale1 = 91.1876, Scale2 = 91.1876  ->  mu = func(91.1876,91.1876) = "<<(*Fct_MuF)(91.1876,91.1876)<<endl;
   info<<"Scale1 = 1,       Scale2 = 91.1876  ->  mu = func(1,91.1876)       = "<<(*Fct_MuF)(1,91.1876)<<endl;
   info<<"Scale1 = 91.1876, Scale2 = 1        ->  mu = func(91.1876,1)       = "<<(*Fct_MuF)(91.1876,1)<<endl;
}



//______________________________________________________________________________


void FastNLOReader::StripWhitespace(string* s) {
   string fastlast = &(*s)[s->size()-1];
   while (!fastlast.compare(" ")) {
      string::iterator it = s->end();
      s->erase(it-1);
      fastlast = &(*s)[s->size()-1];
   }
}

//______________________________________________________________________________


int FastNLOReader::ContrId(const ESMCalculation eCalc, const ESMOrder eOrder) const {
   int Id = -1;
   if (BBlocksSMCalc.empty() || bUseSMCalc[eCalc].empty()) {
      return Id;
   }

   // Requested order
   string requested = fOrdName[eCalc][eOrder];
   // Loop over all available orders of contribution type eCalc
   for (unsigned int i=0; i<BBlocksSMCalc[eCalc].size(); i++) {
      // Found order
      string available = fOrdName[BBlocksSMCalc[eCalc][i]->IContrFlag1-1][BBlocksSMCalc[eCalc][i]->IContrFlag2-1];
      if (available == requested) {
         Id = i;
      }
   }
   return Id;
}

//______________________________________________________________________________

