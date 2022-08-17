#include <cfloat>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#include "fastnlotk/fastNLOTools.h"

using namespace std;
using namespace say;
using namespace fastNLO;

namespace fastNLOTools {

   //________________________________________________________________________________________________________________ //
   void PrintFastnloVersion() {
      char fnlo[100];
      sprintf(fnlo,"%c[%d;%dmfast%c[%d;%dmNLO\033[0m",27,0,31,27,0,34);
      char subproject[100]      = FNLO_SUBPROJECT;
      char package_version[100] = FNLO_VERSION;
      char gitrev[100]          = FNLO_GITREV;
      cout << fnlo << "_" << subproject << " Version " << package_version << "_" << gitrev << endl;
      return;
   }

   //________________________________________________________________________________________________________________ //
   bool CheckVersion(int version ){
      if ( fastNLO::CompatibleVersions.count(version) == 0 ) {
         error["fastNLOTools::CheckVersion"]<<"This table version ("<<version<<") is incompatible with this fastNLO code."<<endl;
         error["fastNLOTools::CheckVersion"]<<"Supported table versions are:";
         for ( auto i : fastNLO::CompatibleVersions ) error>>" "<<i;
         error>>""<<endl;
         error["fastNLOTools::CheckVersion"]<<"Exiting."<<endl;
         exit(1);
         return false;
      }
      return true;
   }


   //________________________________________________________________________________________________________________ //
   int ReadVector(vector<double >& v, istream& table , double nevts ){
      //! Read values according to the size() of the given vector
      //! from table (v2.0 format).
      const bool ReadBinary = false;
      if ( !ReadBinary ) {
         for( unsigned int i=0 ; i<v.size() ; i++){
            table >> v[i];
            v[i] *= nevts;
            if ( !isfinite(v[i]) ) {
               error["ReadVector"]<<"Non-finite number read from table, aborted! value = " << v[i] << endl;
               error["ReadVector"]<<"Please check the table content." << endl;
               exit(1);
            }
         }
      }
      else {
         table.get();
         float f;
         for ( unsigned int k = 0; k < v.size(); ++k ) {
            table.read(reinterpret_cast<char *>(&f), sizeof(f));
            v[k] = f*nevts;
         }
      }
      return v.size();
   }

   //________________________________________________________________________________________________________________ //
   int ReadUnused(istream& table ){
      //! Read values, which are not known to the current code.
      int nLines = 0;
      table >> nLines;
      if ( nLines==fastNLO::tablemagicno ) {
         error["ReadUnused"]<<"Number of lines identical to magic number. Exiting."<<endl; exit(3);
      }
      string sUnused;
      if ( nLines > 0 ) std::getline(table,sUnused); // discard empty space due to precendent >>
      for( int i=0 ; i<nLines ; i++)
         std::getline(table,sUnused);
      return nLines;
   }


   //________________________________________________________________________________________________________________ //
   int ReadFlexibleVector(vector<double >& v, istream& table , int nProcLast , double nevts ){
      int nn = 0;
      if ( nProcLast == 0 ) {
         table >> nProcLast;
         nn++;
      }
      v.resize(nProcLast);
      for(unsigned int i0=0;i0<v.size();i0++){
         table >> v[i0];
         v[i0] *= nevts;
         nn++;
         if ( !isfinite(v[i0]) ) {
            error["ReadFlexibleVector"]<<"Non-finite number read from table, aborted! value = " << v[i0] << endl;
            error["ReadFlexibleVector"]<<"Please check the table content." << endl;
            exit(1);
         }
      }
      return nn;
   }

   //________________________________________________________________________________________________________________ //
   int ReadFlexibleVector(vector<unsigned long long >& v, istream& table , int nProcLast , double nevts ){
      int nn = 0;
      if ( nProcLast == 0 ) {
         table >> nProcLast;
         nn++;
      }
      v.resize(nProcLast);
      for(unsigned int i0=0;i0<v.size();i0++){
         char buffer[256];
         table >> buffer;
         double value = atof(buffer);
         v[i0]  = value;
         v[i0] *= nevts;
         nn++;
      }
      return nn;
   }


   //________________________________________________________________________________________________________________ //
   int ReadFlexibleVector(vector<std::string >& v, istream& table , int size , double nevts ){
      if ( size == 0 ) table >> size;
      v.resize(size);
      if ( size > 0 ) std::getline(table,v[0]); // discard empty space due to precendent >>
      for( auto& i : v) {
         std::getline(table,i);
      }
      return v.size() + 1;
   }


   //________________________________________________________________________________________________________________ //
   int ReadFlexibleVector(vector<int >& v, istream& table , int size , double nevts ){
      if ( size == 0 ) table >> size;
      v.resize(size);
      for( auto& i : v) {
         table >> i;
      }
      return v.size() + 1;
   }


   //________________________________________________________________________________________________________________ //
   void ResizeVector( v7d& v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5, int dim6 ){
      if ( dim0 > 0 ){
         v.resize(dim0);
         for ( int i= 0 ; i<dim0 ; i++)
            ResizeVector( v[i] , dim1, dim2, dim3, dim4, dim5, dim6 );
      } else {
         error["fastNLOTools::ResizeVector"]<<"Cannot resize table, because dimension is <= zero (dim0="<<dim0<<"). Exiting"<<endl;
         exit(1);
      }
   }


   //________________________________________________________________________________________________________________ //
   void ResizeVector( v6d& v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5 ){
      if ( dim0 > 0 ){
         v.resize(dim0);
         for ( int i= 0 ; i<dim0 ; i++)
            ResizeVector( v[i] , dim1, dim2, dim3, dim4, dim5 );
      } else {
         error["fastNLOTools::ResizeVector"]<<"Cannot resize table, because dimension is <= zero (dim0="<<dim0<<"). Exiting"<<endl;
         exit(1);
      }
   }




   //________________________________________________________________________________________________________________ //
   void ResizeVector( v5d& v, int dim0 , int dim1, int dim2, int dim3, int dim4 ){
      if ( dim0 > 0 ){
         v.resize(dim0);
         for ( int i= 0 ; i<dim0 ; i++)
            ResizeVector( v[i] , dim1, dim2, dim3, dim4 );
      } else {
         error["fastNLOTools::ResizeVector"]<<"Cannot resize table, because dimension is <= zero (dim0="<<dim0<<"). Exiting"<<endl;
         exit(1);
      }
   }


   //________________________________________________________________________________________________________________ //
   void ResizeVector( v4d& v, int dim0 , int dim1, int dim2, int dim3 ){
      if ( dim0 > 0 ){
         v.resize(dim0);
         for ( int i= 0 ; i<dim0 ; i++)
            ResizeVector( v[i] , dim1, dim2, dim3 );
      } else {
         error["fastNLOTools::ResizeVector"]<<"Cannot resize table, because dimension is <= zero (dim0="<<dim0<<"). Exiting"<<endl;
         exit(1);
      }
   }


   //________________________________________________________________________________________________________________ //
   void ResizeVector( v3d& v, int dim0 , int dim1, int dim2 ){
      if ( dim0 > 0 ){
         v.resize(dim0);
         for ( int i= 0 ; i<dim0 ; i++)
            ResizeVector( v[i] , dim1, dim2 );
      } else {
         error["fastNLOTools::ResizeVector"]<<"Cannot resize table, because dimension is <= zero (dim0="<<dim0<<"). Exiting"<<endl;
         exit(1);
      }
   }


   //________________________________________________________________________________________________________________ //
   void ResizeVector( v2d&  v, int dim0 , int dim1 ){
      if ( dim0 > 0 ){
         v.resize(dim0);
         for ( int i= 0 ; i<dim0 ; i++)
            ResizeVector( v[i] , dim1 );
      } else {
         error["fastNLOTools::ResizeVector"]<<"Cannot resize table, because dimension is <= zero (dim0="<<dim0<<"). Exiting"<<endl;
         exit(1);
      }
   }


   //________________________________________________________________________________________________________________ //
   void ResizeVector( v1d& v, int dim0 ){
      if ( dim0 > 0 )
         v.resize(dim0);
      else{
         error["fastNLOTools::ResizeVector"]<<"Cannot resize table, because dimension is <= zero (dim0="<<dim0<<"). Exiting"<<endl;
         exit(1);
      }
   }



   //______________________________________________________________________________
   void ResizeFlexibleVector(vector<double >& v, const vector<double >& nom) {
      v.resize(nom.size());
   }


   //______________________________________________________________________________
   int WriteVector( const vector<string >& v, ostream& table , double nevts ) {
      if ( nevts != 0 ) {
         error["fastNLOTools::WriteVector"]<<"Cannot scale a string table by nevts (nevts="<<nevts<<")."<<endl;
         return -1000;
      }
      else return _Write1DVector(v,table);
   }
   //______________________________________________________________________________
   int WriteVector( const vector<double >& v, ostream& table , double nevts ) {
      return _Write1DVectorByN(v,table,nevts);
   }
   //______________________________________________________________________________
   int WriteVector( const vector<int >& v, ostream& table , double nevts ) {
      return _Write1DVectorByN(v,table,nevts);
   }
   //______________________________________________________________________________
   int WriteVector( const vector<unsigned long long >& v, ostream& table , double nevts ) {
      return _Write1DVectorByN(v,table,nevts);
   }


   //______________________________________________________________________________
   void AddVectors( vector<double >& vSum, const vector<double >& vAdd, double w1, double w2 ) {
      _DoAddVectors(vSum,vAdd,w1,w2);
   }
   //______________________________________________________________________________
   void AddVectors( vector<int >& vSum, const vector<int >& vAdd, double w1, double w2 ) {
      _DoAddVectors(vSum,vAdd,w1,w2);
   }
   //______________________________________________________________________________
   void AddVectors( vector<unsigned long long >& vSum, const vector<unsigned long long >& vAdd, double w1, double w2 ) {
      _DoAddVectors(vSum,vAdd,w1,w2);
   }


   //______________________________________________________________________________
   int WriteFlexibleVector(const vector<double >& v, ostream& table, int nProcLast, double nevts) {
      //! Write 1-dimensional flexible table to disk
      //! nevts: Divide all values by nevts
      //! nProcLast: Specify, if the size of the vector should be written in the first line.
      //! if nProcLast != 0, skip the first line
      if ( nevts == 0 ) {
         error["fastNLOTools::WriteFlexibleVector"]<<"Cannot divide by zero. nProcLast ="<<nProcLast<<endl;
         return -1000;
      }
      if ( nProcLast == 0 )
         table << v.size() << sep;
      if ( nProcLast != 0 && nProcLast != (int)v.size() )
         warn["fastNLOTools::WriteFlexibleVector(double)"]
            <<"Dimension of this vector is not compatible with its size (i.e. nProclast ="<<nProcLast<<", v.size()="<<v.size()<<endl;
      int n = _Write1DVectorByN(v,table,nevts);
      return ( nProcLast == 0 ) ? n+1 : n;
   }

   //______________________________________________________________________________
   int WriteFlexibleVector(const vector<string >& v, ostream& table, int nProcLast, double nevts) {
      //! Write 1-dimensional flexible table to disk
      //! nevts: ignoring nevts !!
      //! nProcLast: Specify, if the size of the vector should be written in the first line.
      //! if nProcLast != 0, skip the fist line
      if ( nevts!= 1 ) warn["fastNLOTools::WriteFlexibleVector(string)"]
         <<"String variable cannot be divided by integer number! Ignoring nevts="<<nevts<<endl;
      if ( nProcLast == 0 )
         table << v.size() << sep;
      if ( nProcLast != 0 && nProcLast != (int)v.size() )
         warn["fastNLOTools::WriteFlexibleVector(string)"]
            <<"Dimension of this vector is not compatible with its size (i.e. nProclast ="<<nProcLast<<", v.size()="<<v.size()<<endl;
      int n = _Write1DVector(v,table);
      return ( nProcLast == 0 ) ? n+1 : n;
   }

   //______________________________________________________________________________
   int WriteFlexibleVector(const vector<unsigned long long >& v, ostream& table, int nProcLast, double nevts) {
      //! Write 1-dimensional flexible table to disk
      //! nevts: ignoring nevts !!
      //! nProcLast: Specify, if the size of the vector should be written in the first line.
      //! if nProcLast != 0, skip the fist line
      if ( nevts!= 1 ) warn["fastNLOTools::WriteFlexibleVector(unsigned long long)"]
         <<"String variable cannot be divided by integer number! Ignoring nevts="<<nevts<<endl;
      if ( nProcLast == 0 )
         table << v.size() << sep;
      if ( nProcLast != 0 && nProcLast != (int)v.size() )
         warn["fastNLOTools::WriteFlexibleVector(string)"]
            <<"Dimension of this vector is not compatible with its size (i.e. nProclast ="<<nProcLast<<", v.size()="<<v.size()<<endl;
      int n = _Write1DVector(v,table);
      return ( nProcLast == 0 ) ? n+1 : n;
   }

   //______________________________________________________________________________
   int WriteFlexibleVector(const vector<int >& v, ostream& table, int nProcLast, double nevts) {
      //! Write 1-dimensional flexible table to disk
      //! nevts: ignoring nevts !!
      //! nProcLast: Specify, if the size of the vector should be written in the first line.
      //! if nProcLast != 0, skip the fist line
      if ( nevts!= 1 ) warn["fastNLOTools::WriteFlexibleVector(int)"]
         <<"Refusing dividing integer numbers by each other! Ignoring nevts="<<nevts<<endl;
      if ( nProcLast == 0 )
         table << v.size() << sep;
      if ( nProcLast != 0 && nProcLast != (int)v.size() )
         warn["fastNLOTools::WriteFlexibleVector(int)"]
            <<"Dimension of this vector is not compatible with its size (i.e. nProclast ="<<nProcLast<<", v.size()="<<v.size()<<endl;
      int n = _Write1DVector(v,table);
      return ( nProcLast == 0 ) ? n+1 : n;
   }

   //________________________________________________________________________________________________________________ //
   void StripWhitespace(string& str) {
      //! remove white spaces from string
      for(string::iterator achar = str.end(); achar>str.begin();achar--) {
         if (*achar==0x20 || *achar==0x00){
            str.erase(achar);
         } else {
            break;
         }
      }
   }

   //________________________________________________________________________________________________________________ //
   void PutBackMagicNo(istream& table){
   //! Put magic number back
      for(int i=0;i<(int)(log10((double)tablemagicno)+1);i++){
         table.unget();
      }
      table.unget();
   }

   //______________________________________________________________________________
   bool ReadMagicNo(istream& table) {
      //! read and crosscheck magic number
      if (table.eof()){
         error["ReadMagicNo"]<<"Cannot read from file. Exiting"<<endl;
         exit(3);
      }
      string line;
      std::getline(table,line);
      if ( line=="" ) std::getline(table,line);  // last one was '<<'
      if( line != std::to_string(tablemagicno)){
         error["ReadMagicNo"]<<"Found '"<<line<<"' instead of "<<tablemagicno<<"."<<endl;
         error["ReadMagicNo"]<<"Did not find magic number, aborting!"<<endl;
         error["ReadMagicNo"]<<"Please check compatibility of tables and program version. Exiting."<<endl;
         exit(2);
         return false;
      };
      return true;
   }

   //______________________________________________________________________________
   std::vector <double> ReadUncertaintyFromFile(std::string filename, unsigned int icola, unsigned int icolb) {
      std::string extension = "";
      std::ifstream infile;
      std::string line;
      std::vector <double> Uncertainty;

      //! Determine extension to differentiate for parsing
      //! - fnlo-tk-statunc:  'log' file extension; column numbers not needed, rel. stat. uncertainty = col #4
      //! - NNLOJET dat file: 'dat' file extension; column numbers not needed, rel. stat. uncertainty = (col #5 / col #4)
      //! - Generic txt file: 'txt' file extension; only icola --> rel. stat. uncertainty = col #icola
      //! -                                         icol a & b --> rel. stat. uncertainty = col #icolb / #icola
      if ( filename.find_last_of(".") != std::string::npos ) {
         extension = filename.substr(filename.find_last_of(".")+1);
      }
      if ( extension != "dat" && extension != "log" && extension != "txt" ) {
         error["ReadUncertaintyFromFile"]<<"Unknown filename extension, aborted! filename = " << filename <<endl;
         exit(34);
      } else if ( extension == "txt" && icola == 0) {
         error["ReadUncertaintyFromFile"]<<"'txt' file found, but column specification is missing, aborted! icola " << icola <<endl;
         exit(35);
      } else if ( extension == "txt" && (icola > 10 || icolb > 10) ) {
         error["ReadUncertaintyFromFile"]<<"'txt' file found, but column specification is too large, aborted! icola, icolb = " << icola << ", " << icolb <<endl;
         exit(35);
      } else {
         info["ReadUncertaintyFromFile"]<<"Reading additional uncertainty content from file: " << filename <<endl;
      }

      infile.open(filename);
      if (infile.is_open()) {
         int  iline = 0;
         bool lline = false;
         // Read line-by-line
         while(std::getline(infile, line)) {
            // Put line into stringstream and read word-by-word
            std::istringstream iss(line);
            std::string word, word1, word2;
            iss >> word;
            // For 'dat' extension assume NNLOJET dat file format:
            // - Skip all lines starting with comment symbol '#'
            // - Read cross section and absolute statistical uncertainty from 4th and 5th columns
            if ( extension == "dat" ) {
               if ( word.at(0) != '#' ) {
                  // Skip first three words of each line
                  iss >> word;
                  iss >> word;
                  double xs, dxs;
                  iss >> xs;
                  iss >> dxs;
                  if ( fabs(xs) > DBL_MIN ) {
                     // Is negative, if NLO_only or NNLO_only x section at production was < 0; keep this as additional information.
                     Uncertainty.push_back(dxs/xs);
                     // Only allow positive numbers with maximum value of 1, i.e. = 100% uncertainty maximum
                     //                     Uncertainty.push_back(std::min(fabs(dxs/xs),1.0));
                  } else {
                     Uncertainty.push_back(0.);
                  }
                  iline += 1;
               }
            }
            // For 'log' extension assume fnlo-tk-stat v2.5 log file format:
            // (New v2.5 separator lines starting with #- - - - - - -)
            // - Start at first line with "#-" as 1st word and
            // - stop again at next line starting with "#-" in 1st word
            else if ( extension == "log" ) {
               if ( word == "#-" ) {
                  lline = ! lline;
               } else if ( lline ) {
                  // Skip second & third word of each uncertainty line (x section; lower uncertainty)
                  iss >> word;
                  iss >> word;
                  double dxsrel;
                  iss >> dxsrel;
                  Uncertainty.push_back(dxsrel);
                  iline += 1;
               }
            }
            // For 'txt' extension either read column #icola or divide column #icolb / #icola; max col = 10
            else if ( extension == "txt" ) {
               double a = 0;
               double b = 0;
               for ( unsigned int ic = 1; ic<11; ic++ ) {
                  if ( ic == icola ) a = std::stod(word);
                  if ( ic == icolb ) b = std::stod(word);
                  iss >> word;
               }
               if ( icolb == 0 ) {
                  Uncertainty.push_back(a);
               } else {
                  if ( fabs(a) > DBL_MIN ) {
                     Uncertainty.push_back(b/a);
                  } else {
                     Uncertainty.push_back(0);
                  }
               }
            } else {
               error["ReadUncertaintyFromFile"]<<"Unknown filename extension, aborted! filename = " << filename <<endl;
               exit(34);
            }
         }
      } else {
         error["ReadUncertaintyFromFile"]<<"Cannot read from file, aborted! filename is: " << filename <<endl;
         exit(33);
      }
      return Uncertainty;
   }

} // end namespace fastNLO
