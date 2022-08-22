#ifndef __fnlotools__
#define __fnlotools__

#include <string>
#include <vector>
#include "fastnlotk/fastNLOConstants.h"
#include "speaker.h"

namespace fastNLOTools {

   const bool binary = false;

   //! - Reading std::vectors from disk
   template<typename T> int ReadVector( std::vector<T>& v, std::istream& table , double nevts = 1);
   int ReadVector( std::vector<double>& v, std::istream& table , double nevts = 1);

   template<typename T>  int ReadFlexibleVector(std::vector<T>& v, std::istream& table, int nProcLast=0 , double nevts = 1);
   int ReadFlexibleVector( std::vector<std::string >& v, std::istream& table , int nProcLast = 0 , double nevts = 1 );
   int ReadFlexibleVector( std::vector<double >& v, std::istream& table , int nProcLast = 0 , double nevts = 1 );
   int ReadFlexibleVector( std::vector<int >& v, std::istream& table , int nProcLast = 0 , double nevts = 1 );
   int ReadFlexibleVector( std::vector<unsigned long long >& v, std::istream& table , int nProcLast = 0 , double nevts = 1 );
   int ReadUnused( std::istream& table );

   //! - Resizing tools
   template<typename T> void ResizeFlexibleVector(std::vector<T>& v, const std::vector<T>& nom);
   void ResizeFlexibleVector( std::vector<double >& v, const std::vector<double >& nom);
   void ResizeFlexibleVector( std::vector<unsigned long long >& v, const std::vector<double >& nom);

   //! - Clearing tools
   template<typename T> void ClearVector(std::vector<std::vector<T > >& v);
   template<typename T> void ClearVector(std::vector<T>& v);

   // there are nicer  options in c++11
   void ResizeVector( fastNLO::v1d& v, int dim0 );
   void ResizeVector( fastNLO::v2d& v, int dim0 , int dim1 );
   void ResizeVector( fastNLO::v3d& v, int dim0 , int dim1, int dim2 );
   void ResizeVector( fastNLO::v4d& v, int dim0 , int dim1, int dim2, int dim3 );
   void ResizeVector( fastNLO::v5d& v, int dim0 , int dim1, int dim2, int dim3, int dim4 );
   void ResizeVector( fastNLO::v6d& v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5 );
   void ResizeVector( fastNLO::v7d& v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5, int dim6 );

   //! - Check if std::vector is empty
   template<typename T> bool IsEmptyVector(const std::vector<std::vector<T > >& v);
   template<typename T> bool IsEmptyVector(const std::vector<T>& v);

   //! - Writing tables to disk
   //! use 'fastNLO::WriteVector(std::vector..., *table, nevts=1) to write fastNLO table in v2.0 format to disk
   //! use 'fastNLO::WriteFlexibleVector(std::vector..., *table, int nProcLast=0, nevts=1) to write 'flexible' table
   template<typename T> int WriteVector( const std::vector<T>& v, std::ostream& table , double nevts=1 );
   template<typename T> int _Write1DVectorByN( const std::vector<T>& v, std::ostream& table , double nevts );
   template<typename T> int _Write1DVector( const std::vector<T>& v, std::ostream& table);
   int WriteVector( const std::vector<double >& v, std::ostream& table , double nevts=1 );
   int WriteVector( const std::vector<std::string >& v, std::ostream& table , double nevts=1 );
   int WriteVector( const std::vector<int >& v, std::ostream& table , double nevts=1 ) ;
   int WriteVector( const std::vector<unsigned long long >& v, std::ostream& table , double nevts=1 );

   template<typename T> int WriteFlexibleVector( const std::vector<T>& v, std::ostream& table, int nProcLast = 0, double nevts=1 );
   int WriteFlexibleVector( const std::vector<double >& v, std::ostream& table, int nProcLast = 0 , double nevts=1 );
   int WriteFlexibleVector( const std::vector<std::string >& v, std::ostream& table, int nProcLast = 0 , double nevts=1 );
   int WriteFlexibleVector( const std::vector<int >& v, std::ostream& table, int nProcLast = 0 , double nevts=1 );
   int WriteFlexibleVector( const std::vector<unsigned long long >& v, std::ostream& table, int nProcLast = 0 , double nevts=1 );

   //! - adding std::vectors
   template<typename T> void AddVectors( std::vector<T>& vSum, const std::vector<T>& vAdd, double w1 = 1, double w2 = 1 );
   template<typename T> void _DoAddVectors( std::vector<T>& vSum, const std::vector<T>& vAdd, double w1 = 1, double w2 = 1 );
   void AddVectors( std::vector<double >& vSum, const std::vector<double >& vAdd, double w1 = 1, double w2 = 1 ) ;
   void AddVectors( std::vector<int >& vSum, const std::vector<int >& vAdd, double w1 = 1, double w2 = 1 ) ;
   void AddVectors( std::vector<unsigned long long >& vSum, const std::vector<unsigned long long >& vAdd, double w1=1, double w2=1  ) ;

   //! - std::string modifications
   void StripWhitespace(std::string& s);

   //! - Printout of std::vectors
   template<typename T> void PrintVector( const std::vector<T>& v, std::string name, std::string prefix="");

   //! - useful i/o
   void PrintFastnloVersion(); //!< Print out fastNLO version
   bool CheckVersion(int version); //!< check version and exit if failed.
   //bool CheckVersion(const std::string& version) {return CheckVersion((int)std::stoi(version));} ; //!< check version and exit if failed.
   bool ReadMagicNo(std::istream& table);                                       //!< Read and check magic number from table.
   void PutBackMagicNo(std::istream& table);                                    //!< Reset magic number, such that it can be recognized by other reading routines

   //! Parse filename for uncertainties
   //! - fnlo-tk-statunc:  'log' file extension; column numbers not needed, rel. stat. uncertainty = col #4
   //! - NNLOJET dat file: 'dat' file extension; column numbers not needed, rel. stat. uncertainty = (col #5 / col #4)
   //! - Generic txt file: 'txt' file extension; only icola --> rel. stat. uncertainty = col #icola
   //! -                                         icol a & b --> rel. stat. uncertainty = col #icolb / #icola
   std::vector <double> ReadUncertaintyFromFile(std::string filename, unsigned int icola = 0, unsigned int icolb = 0);

};



//________________________________________________________________________________________________________________
// Reading functions
template<typename T>
int fastNLOTools::ReadVector( std::vector<T>& v, std::istream& table , double nevts){
   //! Read values according to the size() of the given std::vector
   //! from table (v2.0 format).
   int nn = 0;
   for( unsigned int i=0 ; i<v.size() ; i++ ){
      nn += ReadVector(v[i],table, nevts);
   }
   return nn;
};


template<typename T>
int fastNLOTools::ReadFlexibleVector(std::vector<T>& v, std::istream& table, int nProcLast, double nevts ){
   int nn = 0;
   int size = 0;
   table >> size; nn++;
   v.resize(size);
   for(unsigned int i0=0;i0<v.size();i0++){
      nn += ReadFlexibleVector(v[i0],table,nProcLast,nevts);
   }
   return nn;
};

//________________________________________________________________________________________________________________
// Resizing functions
template<typename T>
void fastNLOTools::ResizeFlexibleVector(std::vector<T>& v, const std::vector<T>& nom) {
   v.resize(nom.size());
   for (unsigned int i = 0 ; i<v.size() ; i++) {
      ResizeFlexibleVector(v[i],nom[i]);
   }
};


//________________________________________________________________________________________________________________
// Clearing
template<typename T>
void fastNLOTools::ClearVector(std::vector<std::vector<T > >& v) {
   for (unsigned int i = 0 ; i<v.size() ; i++) {
      ClearVector(v[i]);
   }
};

template<typename T>
void fastNLOTools::ClearVector(std::vector<T >& v) {
   for (unsigned int i = 0 ; i<v.size() ; i++) {
      v[i]=0;
   }
};


//________________________________________________________________________________________________________________
// Check if std::vector is empty
template<typename T>
bool fastNLOTools::IsEmptyVector(const std::vector<std::vector<T > >& v){
   //! check if std::vector is 'empty', or if sum of all elements is 0.
   if ( v.empty() ) return true;
   for (unsigned int i = 0 ; i<v.size() ; i++) {
      if ( !IsEmptyVector(v[i]) ) return false;
   }
   return true;
}

template<typename T>
bool fastNLOTools::IsEmptyVector(const std::vector<T>& v){
   //! check if std::vector is 'empty', or if sum of all elements is 0.
   if ( v.empty() ) return true;
   for (unsigned int i = 0 ; i<v.size() ; i++) {
      if ( v[i] != 0 ) return false;
   }
   return true;
}


//________________________________________________________________________________________________________________
// Writing functions
template<typename T>
int fastNLOTools::WriteVector( const std::vector<T>& v, std::ostream& table , double nevts) {
   //! Write values of std::vector v to table (v2.0 format) .
   int nn = 0;
   for(unsigned int i=0;i<v.size();i++)
      nn += WriteVector( v[i] , table , nevts );
   return nn;
}

template<typename T>
int fastNLOTools::_Write1DVectorByN( const std::vector<T>& v, std::ostream& table , double nevts) {
   if( nevts == 0) return -1000;
   // --- ascii
   if ( !fastNLOTools::binary ) {
      for(unsigned int i0=0;i0<v.size();i0++)
         table << v[i0] / nevts << fastNLO::sep;
   }
   else {
      // --- binary as float
      std::vector<float> ff;
      ff.reserve(v.size());
      for ( auto val : v ) ff.push_back(val/nevts);
      table << 'b';
      //table.flush();
      table.write(reinterpret_cast<const char *>(&ff[0]), ff.size()*sizeof(float));
      //table << std::endl;
   }
   /* static int bb = 0; */
   /* if ( bb++ > 3 )  exit(1); */
   /* std::cout<<"size: " <<v.size() <<std::endl; */
   /* for ( auto val : v ) std::cout<<"\t"<<val; */
   /* std::cout<<std::endl; */

   return v.size();
}

template<typename T>
int fastNLOTools::_Write1DVector( const std::vector<T>& v, std::ostream& table ) {
   //if ( !fastNLOTools::binary ) {
      for(unsigned int i0=0;i0<v.size();i0++)
         table << v[i0] << fastNLO::sep;
//   else
//      table.write(reinterpret_cast<const char *>(&v[0]), v.size()*sizeof(T));
   return v.size();
}


template<typename T>
int fastNLOTools::WriteFlexibleVector( const std::vector<T>& v, std::ostream& table, int nProcLast , double nevts ) {
   if ( nevts == 0 ) {
      say::error["fastNLOTools::WriteFlexibleVector"]<<"Cannot divide by zero."<<std::endl;
      return -1000;
   }
   int nn = 1;
   table << v.size() << fastNLO::sep;
   for(unsigned int i0=0;i0<v.size();i0++){
      nn += WriteFlexibleVector( v[i0] , table , nProcLast , nevts );
   }
   return nn;
};


//________________________________________________________________________________________________________________
// Adding functions
template<typename T>
void fastNLOTools::AddVectors( std::vector<T>& vSum, const std::vector<T>& vAdd, double w1, double w2 ) {
   //! Add the values of the std::vector vAdd to the std::vector vSum
   //! if weights w1 and w1 are specified, the values are weighted accordingly
   //! i.e.: vSum[i] = w1*vSum[i] + w2*vAdd[i];
   if ( vSum.size() != vAdd.size() ) {
      say::error["fastNLOTools::AddVectors"]
      <<"Cannot add tables with different size. s1="
      <<vSum.size()<<", s2="<<vAdd.size()<<std::endl;
      return;
   }
   for ( unsigned int i = 0 ; i<vSum.size() ; i++ )
      AddVectors( vSum[i], vAdd[i], w1 , w2  );
}

template<typename T>
void fastNLOTools::_DoAddVectors( std::vector<T>& vSum, const std::vector<T>& vAdd, double w1, double w2 ) {
   //! This function infact does the addition
   if ( vSum.size() != vAdd.size() ) {
      say::error["fastNLOTools::_DoAddVectors"]
      <<"Cannot add tables with different size. s1="
      <<vSum.size()<<", s2="<<vAdd.size()<<std::endl;
      return;
   }
   if ( w1==1. && w2==1. )
      for ( unsigned int i = 0 ; i<vSum.size() ; i++ )
         vSum[i] += vAdd[i];
   else
      for ( unsigned int i = 0 ; i<vSum.size() ; i++ )
         vSum[i] =  w1*vSum[i] + w2*vAdd[i];
}

template<typename T>
void fastNLOTools::PrintVector( const std::vector<T>& v, std::string name, std::string prefix){
   std::cout<<" "<<prefix<<" "<<name<<std::endl;
   for(unsigned int i=0;i<v.size();i++){
      std::cout<<" "<<prefix<<"   "<<i<<"\t"<<v[i]<<std::endl;
   }
}


#endif
