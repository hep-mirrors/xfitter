#ifndef __fastNLOCoefficients__
#define __fastNLOCoefficients__

#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <iostream>

#include "fastNLOConstants.h"


class fastNLOCoefficients {

   friend class fastNLOTable;
   friend class fastNLOCreate;

public:
   fastNLOCoefficients();
   fastNLOCoefficients(int NObsBin, int iLOord);
   int Read(std::istream *table);
   int Write(std::ostream *table, int option = 0);
   int Copy(fastNLOCoefficients* other);

   void StripWhitespace(std::string& s) const;

   //    template<typename T> void ResizeFlexibleVector(std::vector<T>* v, std::vector<T>* nom);
   //    void ResizeFlexibleVector(std::vector<double >* v, std::vector<double >*nom ){ v->resize(nom->size());};

   void ResizeTable( std::vector<double >* v, int dim0 );
   void ResizeTable( std::vector<std::vector<double > >*  v, int dim0 , int dim1 );
   void ResizeTable( std::vector<std::vector<double > >*  v, int dim0 , int* dim1GetNxmaxFromDimI );
   void ResizeTable( std::vector<std::vector<std::vector<double > > >* v, int dim0 , int* dim1GetNxmaxFromDimI, int dim2 );
   void ResizeTable( std::vector<std::vector<std::vector<double > > >* v, int dim0 , int dim1, int dim2 );
   void ResizeTable( std::vector<std::vector<std::vector<std::vector<double > > > >* v, int dim0 , int dim1, int dim2, int dim3 );
   void ResizeTable( std::vector<std::vector<std::vector<std::vector<double > > > >* v, int dim0 , int dim1, int* dim2GetNxmaxFromDimI, int dim3 );
   void ResizeTable( std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4 );
   void ResizeTable( std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > >* v, int dim0 , int dim1, int dim2, int* dim3GetNxmaxFromDimI, int dim4 );
   void ResizeTable( std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > >* v, int dim0 , int* dim1GetNxmaxFromDimI, int dim2, int dim3, int dim4 );
   void ResizeTable( std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5 );
   void ResizeTable( std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5, int dim6 );
   void ResizeTable( std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int* dim5GetNxmaxFromDimI , int dim6 );


   template<typename T> int ReadTable( std::vector<T>* v, std::istream *table );
   int ReadTable( std::vector<double>* v, std::istream *table );

   template<typename T>  int ReadFlexibleVector(std::vector<T>* v, std::istream* table, bool nProcLast=false);
   int ReadFlexibleVector( std::vector<double >* v, std::istream *table , bool nProcLast = false );

   int WriteTable( std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > > >* v, std::ostream *table , bool DivByNevt=false , int Nevt=1 );
   int WriteTable( std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > >* v, std::ostream *table , bool DivByNevt=false, int Nevt=1 );
   int WriteTable( std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > >* v, std::ostream *table , bool DivByNevt=false, int Nevt=1 );
   int WriteTable( std::vector<std::vector<std::vector<std::vector<double > > > >* v, std::ostream *table , bool DivByNevt=false, int Nevt=1 );
   int WriteTable( std::vector<std::vector<std::vector<double > > >* v, std::ostream *table , bool DivByNevt=false, int Nevt=1 );
   int WriteTable( std::vector<std::vector<double > >* v, std::ostream *table , bool DivByNevt=false, int Nevt=1 );
   int WriteTable( std::vector<double >* v, std::ostream *table , bool DivByNevt=false, int Nevt=1 );

   int WriteFlexibleTable( std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > > >* v, std::ostream *table , bool DivByNevt=false , int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > >* v, std::ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > >* v, std::ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( std::vector<std::vector<std::vector<std::vector<double > > > >* v, std::ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( std::vector<std::vector<std::vector<double > > >* v, std::ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( std::vector<std::vector<double > >* v, std::ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( std::vector<double >* v, std::ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );

   void AddTableToAnotherTable( std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > > >* vSum, std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > > >* vAdd, double w1 = 0, double w2 = 0 );
   void AddTableToAnotherTable( std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > >* vSum, std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > >* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > >* vSum, std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > >* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( std::vector<std::vector<std::vector<std::vector<double > > > >* vSum, std::vector<std::vector<std::vector<std::vector<double > > > >* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( std::vector<std::vector<std::vector<double > > >* vSum, std::vector<std::vector<std::vector<double > > >* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( std::vector<std::vector<double > >* vSum, std::vector<std::vector<double > >* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( std::vector<double >* vSum, std::vector<double >* vAdd, double w1 = 1, double w2 = 1 );

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
   std::vector < std::string > CtrbDescript;
   // obsolete   int NCodeDescr;
   std::vector < std::string > CodeDescript;
   int Nuncorrel;
   std::vector < std::string > UncDescr;
   int Ncorrel;
   std::vector < std::string > CorDescr;
   std::vector < double > Xcenter;
   std::vector < double > Value;
   std::vector < std::vector < double > > UncorLo;
   std::vector < std::vector < double > > UncorHi;
   std::vector < std::vector < double > > CorrLo;
   std::vector < std::vector < double > > CorrHi;
   int NErrMatrix;
   std::vector < std::vector < double > > matrixelement;
   std::vector < double > fact;
   int IRef;
   int IScaleDep;
   unsigned long long int Nevt;
   int Npow;
   int NPDF;
   std::vector < int > NPDFPDG;
   int NPDFDim;
   int NFragFunc;
   std::vector < int > NFFPDG;
   int NFFDim;
   int NSubproc;
   int IPDFdef1;
   int IPDFdef2;
   int IPDFdef3;
   // Missing: linear PDF combinations for IPDFdef1=0
   //std::vector < int > Nxtot1;
   std::vector < double > Hxlim1;
   std::vector < std::vector < double > > XNode1;
   //std::vector < int >  Nxtot2;
   std::vector < double > Hxlim2;
   std::vector < std::vector < double > > XNode2;
   std::vector < int > Nztot;
   std::vector < double > Hzlim;
   std::vector < std::vector < double > > ZNode;
   int NScales;
   int NScaleDim;
   std::vector < int > Iscale;
   // obsolete std::vector < int > NscaleDescript;
   std::vector < std::vector < std::string > > ScaleDescript;


   // ------------------------ not flexible scale --------------------------- //
   std::vector < int > Nscalevar;
   std::vector < int > Nscalenode;
   std::vector < std::vector < double > > ScaleFac;
   std::vector < std::vector < std::vector < std::vector < double > > > > ScaleNode;
   //std::vector < std::vector < std::vector < std::vector < double > > > > HScaleNode;
   std::vector < std::vector < std::vector < std::vector < std::vector < double > > > > > SigmaTilde;

   // --------------------------- flexible scale --------------------------- //
   // ---- members to write to disc ---- //
   // SigmaTilde [NObsBins] ['n' x-nodes] [n s1-Nodes] [n s2-Nodes] [nsubproc]
   std::vector < std::vector < std::vector < std::vector < std::vector < double > > > > > SigmaTildeMuIndep; 
   std::vector < std::vector < std::vector < std::vector < std::vector < double > > > > > SigmaTildeMuFDep; 
   std::vector < std::vector < std::vector < std::vector < std::vector < double > > > > > SigmaTildeMuRDep; 
   // SigmaRef [NObsBins] [nsubproc]
   std::vector < std::vector < double > > SigmaRefMixed; 
   std::vector < std::vector < double > > SigmaRef_s1; 
   std::vector < std::vector < double > > SigmaRef_s2; 
   int NscalenodeScale1;
   int NscalenodeScale2;
   // ScaleNodeXY [ObsBin] [NscalenodeScaleX]  
   std::vector < std::vector < double > > ScaleNode1;
   std::vector < std::vector < double > > ScaleNode2;
   //    std::vector < std::vector < double > > HScaleNode1;
   //    std::vector < std::vector < double > > HScaleNode2;

};


template<typename T>
int fastNLOCoefficients::ReadFlexibleVector(std::vector<T>* v, std::istream* table, bool nProcLast){
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
// void fastNLOCoefficients::ResizeFlexibleVector(std::vector<T>* v, std::vector<T>* nom){
//    v->resize(nom->size());
//    for ( unsigned int i = 0 ; i<v->size() ; i++ ){
//       ResizeFlexibleVector(&((*v)[i]),&((*nom)[i]));
//    }
// };

template<typename T> int fastNLOCoefficients::ReadTable( std::vector<T>* v, std::istream *table ){
   int nn = 0;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn+= ReadTable(&(*v)[i0],table);
   }
   return nn;
};

#endif
