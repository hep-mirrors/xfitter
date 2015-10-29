#ifndef __fnlotools__
#define __fnlotools__

#include <string>
#include <vector>
#include "fastnlotk/fastNLOConstants.h"
#include "speaker.h"

namespace fastNLOTools {

   using namespace std;
   using namespace say;
   using namespace fastNLO;

   //! - Reading vectors from disk
   template<typename T> int ReadVector( vector<T>& v, istream& table , double nevts = 1);
   int ReadVector( vector<double>& v, istream& table , double nevts = 1);

   template<typename T>  int ReadFlexibleVector(vector<T>& v, istream& table, int nProcLast=0 , double nevts = 1);
   int ReadFlexibleVector( vector<double >& v, istream& table , int nProcLast = 0 , double nevts = 1 );

   //! - Resizing tools
   template<typename T> void ResizeFlexibleVector(vector<T>& v, const vector<T>& nom);
   void ResizeFlexibleVector( vector<double >& v, const vector<double >& nom);

   //! - Clearing tools
   template<typename T> void ClearVector(vector<vector<T > >& v);
   template<typename T> void ClearVector(vector<T>& v);

   // there are nicer  options in c++11
   void ResizeVector( v1d& v, int dim0 );
   void ResizeVector( v2d& v, int dim0 , int dim1 );
   void ResizeVector( v3d& v, int dim0 , int dim1, int dim2 );
   void ResizeVector( v4d& v, int dim0 , int dim1, int dim2, int dim3 );
   void ResizeVector( v5d& v, int dim0 , int dim1, int dim2, int dim3, int dim4 );
   void ResizeVector( v6d& v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5 );
   void ResizeVector( v7d& v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5, int dim6 );

   //! - Check if vector is empty
   template<typename T> bool IsEmptyVector(const vector<vector<T > >& v);
   template<typename T> bool IsEmptyVector(const vector<T>& v);

   //! - Writing tables to disk
   //! use 'fastNLO::WriteVector(vector..., *table, nevts=1) to write fastNLO table in v2.0 format to disk
   //! use 'fastNLO::WriteFlexibleVector(vector..., *table, int nProcLast=0, nevts=1) to write 'flexible' table
   template<typename T> int WriteVector( const vector<T>& v, ostream& table , double nevts=1 );
   template<typename T> int _Write1DVectorByN( const vector<T>& v, ostream& table , double nevts );
   template<typename T> int _Write1DVector( const vector<T>& v, ostream& table);
   int WriteVector( const vector<double >& v, ostream& table , double nevts=1 );
   int WriteVector( const vector<string >& v, ostream& table , double nevts=1 );
   int WriteVector( const vector<int >& v, ostream& table , double nevts=1 ) ;
   int WriteVector( const vector<unsigned long long >& v, ostream& table , double nevts=1 );
   
   template<typename T> int WriteFlexibleVector( const vector<T>& v, ostream& table, int nProcLast = 0, double nevts=1 );
   int WriteFlexibleVector( const vector<double >& v, ostream& table, int nProcLast = 0 , double nevts=1 );
   int WriteFlexibleVector( const vector<string >& v, ostream& table, int nProcLast = 0 , double nevts=1 );
   int WriteFlexibleVector( const vector<int >& v, ostream& table, int nProcLast = 0 , double nevts=1 );
   int WriteFlexibleVector( const vector<unsigned long long >& v, ostream& table, int nProcLast = 0 , double nevts=1 );
   
   //! - adding vectors
   template<typename T> void AddVectors( vector<T>& vSum, const vector<T>& vAdd, double w1 = 1, double w2 = 1 );
   template<typename T> void _DoAddVectors( vector<T>& vSum, const vector<T>& vAdd, double w1 = 1, double w2 = 1 );
   void AddVectors( vector<double >& vSum, const vector<double >& vAdd, double w1 = 1, double w2 = 1 ) ;
   void AddVectors( vector<int >& vSum, const vector<int >& vAdd, double w1 = 1, double w2 = 1 ) ;
   void AddVectors( vector<unsigned long long >& vSum, const vector<unsigned long long >& vAdd, double w1=1, double w2=1  ) ;

   //! - String modifications
   void StripWhitespace(string& s);

   //! - Printout of vectors
   template<typename T> void PrintVector( const vector<T>& v, string name, string prefix="");

   //! - useful i/o
   bool ReadMagicNo(istream& table);					//!< Read and check magic number from table.	
   void PutBackMagicNo(istream& table);					//!< Reset magic number, such that it can be recognized by other reading routines

};



//________________________________________________________________________________________________________________
// Reading functions
template<typename T>
int fastNLOTools::ReadVector( vector<T>& v, istream& table , double nevts){
   //! Read values according to the size() of the given vector
   //! from table (v2.0 format).
   int nn = 0;
   for( unsigned int i=0 ; i<v.size() ; i++ ){
      nn += ReadVector(v[i],table, nevts);
   }
   return nn;
};


template<typename T>
int fastNLOTools::ReadFlexibleVector(vector<T>& v, istream& table, int nProcLast, double nevts ){
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
void fastNLOTools::ResizeFlexibleVector(vector<T>& v, const vector<T>& nom) {
   v.resize(nom.size());
   for (unsigned int i = 0 ; i<v.size() ; i++) {
      ResizeFlexibleVector(v[i],nom[i]);
   }
};


//________________________________________________________________________________________________________________
// Clearing 
template<typename T>
void fastNLOTools::ClearVector(vector<vector<T > >& v) {
   for (unsigned int i = 0 ; i<v.size() ; i++) {
      ClearVector(v[i]);
   }
};

template<typename T>
void fastNLOTools::ClearVector(vector<T >& v) {
   for (unsigned int i = 0 ; i<v.size() ; i++) {
      v[i]=0;
   }
};


//________________________________________________________________________________________________________________
// Check if vector is empty
template<typename T> 
bool fastNLOTools::IsEmptyVector(const vector<vector<T > >& v){
   //! check if vector is 'empty', or if sum of all elements is 0.
   if ( v.empty() ) return true;
   for (unsigned int i = 0 ; i<v.size() ; i++) {
      if ( !IsEmptyVector(v[i]) ) return false;
   }
   return true;
}

template<typename T> 
bool fastNLOTools::IsEmptyVector(const vector<T>& v){
   //! check if vector is 'empty', or if sum of all elements is 0.
   if ( v.empty() ) return true;
   for (unsigned int i = 0 ; i<v.size() ; i++) {
      if ( v[i] != 0 ) return false;
   }
   return true;
}


//________________________________________________________________________________________________________________
// Writing functions
template<typename T>
int fastNLOTools::WriteVector( const vector<T>& v, ostream& table , double nevts) {
   //! Write values of vector v to table (v2.0 format) .
   int nn = 0;
   for(unsigned int i=0;i<v.size();i++)
      nn += WriteVector( v[i] , table , nevts );
   return nn;
}

template<typename T>
int fastNLOTools::_Write1DVectorByN( const vector<T>& v, ostream& table , double nevts) {
   if( nevts == 0) return -1000;
   for(unsigned int i0=0;i0<v.size();i0++)
      table << v[i0] / nevts << endl;
   return v.size();
}

template<typename T>
int fastNLOTools::_Write1DVector( const vector<T>& v, ostream& table ) {
   for(unsigned int i0=0;i0<v.size();i0++)
      table << v[i0] << endl;
   return v.size();
}


template<typename T>
int fastNLOTools::WriteFlexibleVector( const vector<T>& v, ostream& table, int nProcLast , double nevts ) {
   if ( nevts == 0 ) {
      error["fastNLOTools::WriteFlexibleVector"]<<"Cannot divide by zero."<<endl;
      return -1000;
   }
   int nn = 1;
   table << v.size() << endl;
   for(unsigned int i0=0;i0<v.size();i0++){
      nn += WriteFlexibleVector( v[i0] , table , nProcLast , nevts );
   }
   return nn;
};


//________________________________________________________________________________________________________________
// Adding functions
template<typename T>
void fastNLOTools::AddVectors( vector<T>& vSum, const vector<T>& vAdd, double w1, double w2 ) {
   //! Add the values of the vector vAdd to the vector vSum
   //! if weights w1 and w1 are specified, the values are weighted accordingly
   //! i.e.: vSum[i] = w1*vSum[i] + w2*vAdd[i];
   if ( vSum.size() != vAdd.size() ) {
      error["fastNLOTools::AddVectors"]
	 <<"Cannot add tables with different size. s1="
	 <<vSum.size()<<", s2="<<vAdd.size()<<endl;
      return;
   }
   for ( unsigned int i = 0 ; i<vSum.size() ; i++ )
      AddVectors( vSum[i], vAdd[i], w1 , w2  );
}

template<typename T>
void fastNLOTools::_DoAddVectors( vector<T>& vSum, const vector<T>& vAdd, double w1, double w2 ) {
   //! This function infact does the addition
   if ( vSum.size() != vAdd.size() ) {
      error["fastNLOTools::_DoAddVectors"]
	 <<"Cannot add tables with different size. s1="
	 <<vSum.size()<<", s2="<<vAdd.size()<<endl;
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
void fastNLOTools::PrintVector( const vector<T>& v, string name, string prefix){
   cout<<" "<<prefix<<"   "<<name<<endl;
   for(unsigned int i=0;i<v.size();i++){
      cout<<" "<<prefix<<"     "<<i<<"\t"<<v[i]<<endl;
   }
}


#endif
