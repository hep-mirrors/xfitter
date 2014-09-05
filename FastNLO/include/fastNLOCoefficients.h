#ifndef __fastNLOCoefficients__
#define __fastNLOCoefficients__

#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <iostream>

#include "fastNLOConstants.h"

using namespace std;

class fastNLOCoefficients {

   friend class fastNLOTable;
   friend class fastNLOCreate;

public:
   fastNLOCoefficients();
   fastNLOCoefficients(int NObsBin, int iLOord);
   int Read(istream *table);
   int Write(ostream *table, int option = 0);
   int Copy(fastNLOCoefficients* other);

   void StripWhitespace(string& s) const;

   //    template<typename T> void ResizeFlexibleVector(vector<T>* v, vector<T>* nom);
   //    void ResizeFlexibleVector(vector<double >* v, vector<double >*nom ){ v->resize(nom->size());};

   void ResizeTable( vector<double >* v, int dim0 );
   void ResizeTable( vector<vector<double > >*  v, int dim0 , int dim1 );
   void ResizeTable( vector<vector<double > >*  v, int dim0 , int* dim1GetNxmaxFromDimI );
   void ResizeTable( vector<vector<vector<double > > >* v, int dim0 , int* dim1GetNxmaxFromDimI, int dim2 );
   void ResizeTable( vector<vector<vector<double > > >* v, int dim0 , int dim1, int dim2 );
   void ResizeTable( vector<vector<vector<vector<double > > > >* v, int dim0 , int dim1, int dim2, int dim3 );
   void ResizeTable( vector<vector<vector<vector<double > > > >* v, int dim0 , int dim1, int* dim2GetNxmaxFromDimI, int dim3 );
   void ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4 );
   void ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int dim1, int dim2, int* dim3GetNxmaxFromDimI, int dim4 );
   void ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int* dim1GetNxmaxFromDimI, int dim2, int dim3, int dim4 );
   void ResizeTable( vector<vector<vector<vector<vector<vector<double > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5 );
   void ResizeTable( vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5, int dim6 );
   void ResizeTable( vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int* dim5GetNxmaxFromDimI , int dim6 );


   template<typename T> int ReadTable( vector<T>* v, istream *table );
   int ReadTable( vector<double>* v, istream *table );

   template<typename T>  int ReadFlexibleVector(vector<T>* v, istream* table, bool nProcLast=false);
   int ReadFlexibleVector( vector<double >* v, istream *table , bool nProcLast = false );

   int WriteTable( vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, ostream *table , bool DivByNevt=false , int Nevt=1 );
   int WriteTable( vector<vector<vector<vector<vector<vector<double > > > > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   int WriteTable( vector<vector<vector<vector<vector<double > > > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   int WriteTable( vector<vector<vector<vector<double > > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   int WriteTable( vector<vector<vector<double > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   int WriteTable( vector<vector<double > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   int WriteTable( vector<double >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );

   int WriteFlexibleTable( vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, ostream *table , bool DivByNevt=false , int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( vector<vector<vector<vector<vector<vector<double > > > > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( vector<vector<vector<vector<vector<double > > > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( vector<vector<vector<vector<double > > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( vector<vector<vector<double > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( vector<vector<double > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( vector<double >* v, ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );

   void AddTableToAnotherTable( vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* vSum, vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* vAdd, double w1 = 0, double w2 = 0 );
   void AddTableToAnotherTable( vector<vector<vector<vector<vector<vector<double > > > > > >* vSum, vector<vector<vector<vector<vector<vector<double > > > > > >* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( vector<vector<vector<vector<vector<double > > > > >* vSum, vector<vector<vector<vector<vector<double > > > > >* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( vector<vector<vector<vector<double > > > >* vSum, vector<vector<vector<vector<double > > > >* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( vector<vector<vector<double > > >* vSum, vector<vector<vector<double > > >* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( vector<vector<double > >* vSum, vector<vector<double > >* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( vector<double >* vSum, vector<double >* vAdd, double w1 = 1, double w2 = 1 );

   void Print() const;

   int GetIRef() const {return IRef;}
   int GetIDataFlag() const {return IDataFlag;}
   int GetIAddMultFlag() const {return IAddMultFlag;}
   int GetIContrFlag1() const {return IContrFlag1;}
   int GetIContrFlag2() const {return IContrFlag2;}
   int GetNScaleDep() const {return NScaleDep;}
   int GetNpow() const {return Npow;}
   long long int GetNevt() const {return Nevt;}
   int GetNxmax(int Obsbin) const ;
   int GetXIndex(int Obsbin,int x1bin,int x2bin =0) const ;

   void SetNlojetDefaults();
   void SetIXsectUnits(int n){IXsectUnits = n;}
   void SetIDataFlag(int n){IDataFlag = n;}
   void SetIAddMultFlag(int n){IAddMultFlag = n;}
   void SetIContrFlag1(int n){IContrFlag1 = n;} 
   void SetIContrFlag2(int n){IContrFlag2 = n;} 
   void SetNScaleDep(int n){NScaleDep = n;}
   void SetNlojetDescr(){
      CodeDescript.push_back("NLOJet++_4.1.3");
      CodeDescript.push_back("Z. Nagy, Phys. Rev. Lett. 88, 122003 (2002),");
      CodeDescript.push_back("Z. Nagy, Phys. Rev. D68, 094002 (2003).");
   }
   
   void Add(fastNLOCoefficients* other);
   bool IsLO() const {return IContrFlag1==1 && IContrFlag2==1;}
   bool IsNLO() const {return IContrFlag1==1 && IContrFlag2==2;}
   bool IsReference() const {return IRef>0;};
   int GetTotalScalevars() const ;
   int GetTotalScalenodes() const ;

   static const int DividebyNevt = 1; // shitty definition of a global constant

private:
   int GetScaledimfromvar(int scalevar) const;


protected:

   //fastNLOScenario *fScen;
   int fNObsBins; // obtained from Scenario
   int fILOord;   // obtained from Scenario

   int IXsectUnits;
   int IDataFlag;
   int IAddMultFlag;
   int IContrFlag1;
   int IContrFlag2;
   int NScaleDep;
   // obsolete int NContrDescr;
   vector < string > CtrbDescript;
   // obsolete   int NCodeDescr;
   vector < string > CodeDescript;
   int Nuncorrel;
   vector < string > UncDescr;
   int Ncorrel;
   vector < string > CorDescr;
   vector < double > Xcenter;
   vector < double > Value;
   vector < vector < double > > UncorLo;
   vector < vector < double > > UncorHi;
   vector < vector < double > > CorrLo;
   vector < vector < double > > CorrHi;
   int NErrMatrix;
   vector < vector < double > > matrixelement;
   vector < double > fact;
   int IRef;
   int IScaleDep;
   unsigned long long int Nevt;
   int Npow;
   int NPDF;
   vector < int > NPDFPDG;
   int NPDFDim;
   int NFragFunc;
   vector < int > NFFPDG;
   int NFFDim;
   int NSubproc;
   int IPDFdef1;
   int IPDFdef2;
   int IPDFdef3;
   // Missing: linear PDF combinations for IPDFdef1=0
   //vector < int > Nxtot1;
   vector < double > Hxlim1;
   vector < vector < double > > XNode1;
   //vector < int >  Nxtot2;
   vector < double > Hxlim2;
   vector < vector < double > > XNode2;
   vector < int > Nztot;
   vector < double > Hzlim;
   vector < vector < double > > ZNode;
   int NScales;
   int NScaleDim;
   vector < int > Iscale;
   // obsolete vector < int > NscaleDescript;
   vector < vector < string > > ScaleDescript;


   // ------------------------ not flexible scale --------------------------- //
   vector < int > Nscalevar;
   vector < int > Nscalenode;
   vector < vector < double > > ScaleFac;
   vector < vector < vector < vector < double > > > > ScaleNode;
   //vector < vector < vector < vector < double > > > > HScaleNode;
   vector < vector < vector < vector < vector < double > > > > > SigmaTilde;

   // --------------------------- flexible scale --------------------------- //
   // ---- members to write to disc ---- //
   // SigmaTilde [NObsBins] ['n' x-nodes] [n s1-Nodes] [n s2-Nodes] [nsubproc]
   vector < vector < vector < vector < vector < double > > > > > SigmaTildeMuIndep; 
   vector < vector < vector < vector < vector < double > > > > > SigmaTildeMuFDep; 
   vector < vector < vector < vector < vector < double > > > > > SigmaTildeMuRDep; 
   // SigmaRef [NObsBins] [nsubproc]
   vector < vector < double > > SigmaRefMixed; 
   vector < vector < double > > SigmaRef_s1; 
   vector < vector < double > > SigmaRef_s2; 
   int NscalenodeScale1;
   int NscalenodeScale2;
   // ScaleNodeXY [ObsBin] [NscalenodeScaleX]  
   vector < vector < double > > ScaleNode1;
   vector < vector < double > > ScaleNode2;
   //    vector < vector < double > > HScaleNode1;
   //    vector < vector < double > > HScaleNode2;

};


template<typename T>
int fastNLOCoefficients::ReadFlexibleVector(vector<T>* v, istream* table, bool nProcLast){
   int nn = 0;
   int size = 0;
   *table >> size; nn++;
   v->resize(size);
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += ReadFlexibleVector(&(v->at(i0)),table,nProcLast);
   }
   return nn;
};


// template<typename T>
// void fastNLOCoefficients::ResizeFlexibleVector(vector<T>* v, vector<T>* nom){
//    v->resize(nom->size());
//    for ( unsigned int i = 0 ; i<v->size() ; i++ ){
//       ResizeFlexibleVector(&((*v)[i]),&((*nom)[i]));
//    }
// };

template<typename T> int fastNLOCoefficients::ReadTable( vector<T>* v, istream *table ){
   int nn = 0;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn+= ReadTable(&(*v)[i0],table);
   }
   return nn;
};

#endif
