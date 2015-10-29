#include <cstdlib>
#include <unistd.h>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include "fastnlotk/fastNLOConstants.h"
#include "fastnlotk/fastNLOBase.h"
#include "fastnlotk/fastNLOTools.h"

using namespace std;

//______________________________________________________________________________
bool fastNLOBase::fWelcomeOnce = false;


//______________________________________________________________________________
fastNLOBase::fastNLOBase() : fPrecision(8), logger("fastNLOBase") {
   if (!fWelcomeOnce) PrintWelcomeMessage();
}


//______________________________________________________________________________
fastNLOBase::fastNLOBase(string name) : ffilename(name), fPrecision(8), logger("fastNLOBase") {
   if (!fWelcomeOnce) PrintWelcomeMessage();
}


//______________________________________________________________________________
fastNLOBase::fastNLOBase(const fastNLOBase& other) :
   ffilename(other.ffilename), fPrecision(other.fPrecision),
   Itabversion(other.Itabversion), ScenName(other.ScenName),
   Ncontrib(other.Ncontrib), Nmult(other.Nmult),
   Ndata(other.Ndata), NuserString(other.NuserString),
   NuserInt(other.NuserInt), NuserFloat(other.NuserFloat),
   Imachine(other.Imachine),
   logger("fastNLOBase") {
   //! copy constructor
}


//______________________________________________________________________________
fastNLOBase::~fastNLOBase() {
}


//______________________________________________________________________________
ifstream* fastNLOBase::OpenFileRead() {
   //! Open file-stream for reading table
   // does file exist?
   if (access(ffilename.c_str(), R_OK) != 0) {
      logger.error["OpenFileRead"]<<"File does not exist! Was looking for: "<<ffilename<<". Exiting."<<endl;
      exit(1);
   }
   ifstream* strm = new ifstream(ffilename.c_str(),ios::in);
   return strm;
}


//______________________________________________________________________________
void fastNLOBase::CloseFileRead(ifstream& strm) {
   //! Close file-stream
   strm.close();
   delete &strm;
}


//______________________________________________________________________________
void fastNLOBase::ReadTable() {
   // does file exist?
   // open file
   ifstream* strm = OpenFileRead();
   // read header
   ReadHeader(*strm);
   // close stream
   CloseFileRead(*strm);
}


//______________________________________________________________________________
void fastNLOBase::ReadHeader(istream& table) {
   table.peek();
   if (table.eof()) {
      logger.error["ReadHeader"]<<"Cannot read from stream."<<endl;
   }

   if (!fastNLOTools::ReadMagicNo(table)) {
      logger.error["ReadHeader"]<<"Did not find initial magic number, aborting!"<<endl;
      logger.error["ReadHeader"]<<"Please check compatibility of tables and program version!"<<endl;
      exit(1);
   }
   table >> Itabversion;
   table >> ScenName;
   table >> Ncontrib;
   table >> Nmult;
   table >> Ndata;
   table >> NuserString;
   table >> NuserInt;
   for (int i = 0 ; i<NuserInt ; i++) {
      int IUserLines;
      table >> IUserLines;
      // future code if 'user-blocks' are used ...
      // int NUserFlag;
      // string NUserBlockDescr;
      // table >> NUserFlag;
      // table >> NUserBlockDescr;;
      // if ( known-user-block ) { read-known-userblock... }
      // else { // skip meaningful reading
      //    for ( int i = 2 ; i<NuserInt ; i++ ) {
      //       double devnull;
      //       table >> devnull;
      //    }
      // }
      // ...sofar skip reading
      for (int i = 0 ; i<NuserInt ; i++) {
         double devnull;
         table >> devnull;
      }
   }
   table >> NuserFloat;
   table >> Imachine;
   if (!fastNLOTools::ReadMagicNo(table)) {
      logger.error["ReadHeader"]<<"Did not find final magic number, aborting!"<<endl;
      logger.error["ReadHeader"]<<"Please check compatibility of tables and program version!"<<endl;
      exit(1);
   }
   fastNLOTools::PutBackMagicNo(table);
}


//______________________________________________________________________________
void fastNLOBase::WriteTable() {
   //!
   //! WriteTable(). writes the full FastNLO table to
   //! the previously defined ffilename on disk.
   //!
   //! this function is overwritten by
   //! fastNLOTable::WriteTable();
   //!
   ofstream* table = OpenFileWrite();
   WriteHeader(*table);
   CloseFileWrite(*table);
}


//______________________________________________________________________________
void fastNLOBase::WriteHeader(ostream& table) {
   table << fastNLO::tablemagicno << endl;
   table << Itabversion << endl;
   table << ScenName << endl;
   table << Ncontrib << endl;
   table << Nmult << endl;
   table << Ndata << endl;
   table << NuserString << endl;
   table << NuserInt << endl;
   table << NuserFloat << endl;
   table << Imachine << endl;
}


//______________________________________________________________________________
ofstream* fastNLOBase::OpenFileWrite() {
   //! open ofstream for writing tables
   //! do overwrite existing table
   if (access(ffilename.c_str(), F_OK) == 0) {
      logger.info["OpenFileWrite"]<<"Overwriting the already existing table file: " << ffilename << endl;
   }
   ofstream* stream = new ofstream(ffilename.c_str(),ios::out);
   if (!stream->good()) {
      logger.error["OpenFileWrite"]<<"Cannot open file '"<<ffilename<<"' for writing. Aborting."<<endl;
      exit(2);
   }
   stream->precision(fPrecision);
   return stream;
}


//______________________________________________________________________________
void fastNLOBase::CloseFileWrite(ofstream& table) {
   //! close stream and delete object;
   table << fastNLO::tablemagicno << endl;
   table << fastNLO::tablemagicno << endl;
   table.close();
   delete &table;
}


//______________________________________________________________________________
bool fastNLOBase::IsCompatibleHeader(const fastNLOBase& other) const {
   if (Itabversion!= other.GetItabversion()) {
      logger.warn["IsCompatibleHeader"]<<"Differing versions of table format: "<<Itabversion<<" and "<< other.GetItabversion()<<endl;
      return false;
   }
   if (Ndata + other.GetNdata() > 1) {
      logger.warn["IsCompatibleHeader"]<<"Two tables containing both experimental data are incompatible"<<endl;
      return false;
   }
   if (ScenName!= other.GetScenName()) {
      logger.warn["IsCompatibleHeader"]<<"Differing names of scenarios: "<<ScenName.c_str()<<" and "<<other.ScenName.c_str()<<endl;
      // continue...
   }
   return true;
}


//______________________________________________________________________________
void fastNLOBase::SetHeaderDefaults() {
   // TableMagicNo and ITabVersion are defined as constant in fastNLOConstants.h
   SetScenName("tns2000");
   SetContributionHeader();
}


//______________________________________________________________________________
void fastNLOBase::SetContributionHeader() {
   SetNcontrib(1);
   SetNmult(0);
   SetNdata(0);
   SetNuserString(0);
   SetNuserInt(0);
   SetNuserFloat(0);
   SetImachine(0);
}


//______________________________________________________________________________
void fastNLOBase::ResetHeader() {
   logger.debug["ResetHeader"]<<endl;
   SetNcontrib(0);
   SetNmult(0);
   SetNdata(0);
   SetNuserString(0);
   SetNuserInt(0);
   SetNuserFloat(0);
   SetImachine(0);
}


//______________________________________________________________________________
void fastNLOBase::Print() const {
   PrintHeader();
}


//______________________________________________________________________________
void fastNLOBase::PrintHeader() const {
   printf("\n **************** FastNLO Table Header ******************\n\n");
   printf("   tablemagicno                  %d\n",fastNLO::tablemagicno);
   printf("   Itabversion                   %d\n",Itabversion);
   printf("   ScenName                      %s\n",ScenName.data());
   printf("   Ncontrib                      %d\n",Ncontrib);
   printf("   Nmult                         %d\n",Nmult);
   printf("   Ndata                         %d\n",Ndata);
   printf("   NuserString                   %d\n",NuserString);
   printf("   NuserInt                      %d\n",NuserInt);
   printf("   NuserFloat                    %d\n",NuserFloat);
   printf("   Imachine                      %d\n",Imachine);
   printf("\n ********************************************************\n\n");
}


//______________________________________________________________________________
void fastNLOBase::PrintWelcomeMessage() {

   char fnlo[100];
   sprintf(fnlo,"%c[%d;%dmfast%c[%d;%dmNLO\033[0m",27,0,31,27,0,34);
   char subproject[100]      = FNLO_SUBPROJECT;
   char package_version[100] = FNLO_VERSION;
   char svnrev[100]          = FNLO_SVNREV;
   char authors[500]         = FNLO_AUTHORS;
   char webpage[500]         = FNLO_WEBPAGE;
   char authorsv14[200]      = FNLO_AUTHORSv14;
   char quotev14[200]        = FNLO_QUOTEv14;
   char authorsv2[200]       = FNLO_AUTHORSv2;
   char quotev2[200]         = FNLO_QUOTEv2;
   char years[100]           = FNLO_YEARS;

   cout  << endl;
   cout  << fastNLO::_CSEPSC << endl;
   speaker &shout = logger.shout;
   shout << "" << endl;
   shout << fnlo << "_" << subproject << endl;
   shout << "Version " << package_version << "_" << svnrev << endl;
   shout << "" << endl;
   shout << "C++ program and toolkit to read and create fastNLO v2 tables and" << endl;
   shout << "derive QCD cross sections using PDFs, e.g. from LHAPDF" << endl;
   shout << "" << endl;
   cout  << fastNLO::_SSEPSC << endl;
   shout << "" << endl;
   shout << "Copyright © " << years << " " << fnlo << " Collaboration" << endl;
   shout << authors << endl;
   shout << "" << endl;
   shout << "This program is free software: you can redistribute it and/or modify" << endl;
   shout << "it under the terms of the GNU General Public License as published by" << endl;
   shout << "the Free Software Foundation, either version 3 of the License, or" << endl;
   shout << "(at your option) any later version." << endl;
   shout << "" << endl;
   shout << "This program is distributed in the hope that it will be useful," << endl;
   shout << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
   shout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the" << endl;
   shout << "GNU General Public License for more details." << endl;
   shout << "" << endl;
   shout << "You should have received a copy of the GNU General Public License" << endl;
   shout << "along with this program. If not, see <http://www.gnu.org/licenses/>." << endl;
   shout << "" << endl;
   cout  << fastNLO::_SSEPSC << endl;
   shout << "" << endl;
   shout << "The projects web page can be found at:" << endl;
   shout << "  " << webpage << endl;
   shout << "" << endl;
   shout << "If you use this code, please cite:" << endl;
   shout << "  " << authorsv14 << ", " << quotev14 << endl;
   shout << "  " << authorsv2 << ", " << quotev2 << endl;
   shout << "" << endl;
   cout  << fastNLO::_CSEPSC << endl;
   fWelcomeOnce = true;
}
