#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffAddBase.h"
#include "fastnlotk/fastNLOTools.h"

using namespace std;
using namespace fastNLO;

//________________________________________________________________________________________________________________ //
bool fastNLOCoeffAddBase::CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet) {
   if ( c->GetIAddMultFlag()==0 && c->GetIDataFlag()==0 ) {
      // Additive contribution
      return true;
   } else if ( c->GetIAddMultFlag()==1 && c->GetIDataFlag()==0 ) {
      // Multiplicative contribution
      return false;
   } else if ( c->GetIAddMultFlag()==0 && c->GetIDataFlag()==1 ) {
      // Data contribution
      return false;
   } else {
      // Unknown contribution
      say::error["fastNLOCoeffAddBase::CheckCoeffConstants"]
         << "Unknown contribution type, aborting! "
         << "IAddMultFlag = " << c->GetIAddMultFlag()
         << ", IDataFlag ="   << c->GetIDataFlag() <<endl;
      exit(1);
   }
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddBase::fastNLOCoeffAddBase(){
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddBase::fastNLOCoeffAddBase(int NObsBin) : fastNLOCoeffBase(NObsBin) {
   NScaleDep = 0;
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddBase::fastNLOCoeffAddBase(const fastNLOCoeffBase& base ) : fastNLOCoeffBase(base) {
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffBase* fastNLOCoeffAddBase::Clone() const {
   //! Use has to take care to delete this object later
   return static_cast<fastNLOCoeffBase*>(new fastNLOCoeffAddBase(*this));
}


///________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::Read(istream& table){
   fastNLOCoeffBase::ReadBase(table);
   CheckCoeffConstants(this);
   ReadCoeffAddBase(table);
   EndReadCoeff(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::ReadCoeffAddBase(istream& table){
   CheckCoeffConstants(this);
   char buffer[5257];
   table >> IRef;
   table >> IScaleDep;
   table >> Nevt;
   table >> Npow;
   int NPDF;
   table >> NPDF;
   if(NPDF>0){
      NPDFPDG.resize(NPDF);
      for(int i=0;i<NPDF;i++){
         table >>  NPDFPDG[i];
      }
   }
   table >> NPDFDim;
   int NFragFunc;
   table >> NFragFunc;
   if(NFragFunc>0){
      NFFPDG.resize(NFragFunc);
      for(int i=0;i<NFragFunc;i++){
         table >>  NFFPDG[i];
      }
   }
   table >> NFFDim;
   table >> NSubproc;
   table >> IPDFdef1;
   table >> IPDFdef2;
   table >> IPDFdef3;
   //printf("  *  fastNLOCoeffAddBase::Read(). IRef : %d, IScaleDep: %d, Nevt: %d, Npow: %d, NPDF: %d, NPDFDim: %d\n", IRef ,IScaleDep  ,Nevt  , Npow ,NPDF , NPDFDim  );

   if(IPDFdef2==0){ // PDF linear combinations are stored herewith
      if ( IPDFdef3 != NSubproc ){
         error["ReadCoeffAddBase"]<<"IPDFdef3 must be equal to NSubproc. (IPDFdef3="<<IPDFdef3<<", NSubproc="<<NSubproc<<"). Exiting."<<endl;
         exit(1);
      }
      int IPDFCoeffFormat = -1;
      table >> IPDFCoeffFormat;
      if ( IPDFCoeffFormat == 0 ) {
         fPDFCoeff.resize(NSubproc);
         for(int k=0;k<NSubproc;k++){
            int NPartonPairs = -1;
            table >> NPartonPairs;
            for(int n=0;n<NPartonPairs;n++){
               int PDF1Flavor=-100, PDF2Flavor=-100;
               if ( IPDFdef1>=3 ) {
                  table >> PDF1Flavor;
                  table >> PDF2Flavor;
               }
               else if ( IPDFdef1>=3 ) {
                  table >> PDF1Flavor;
                  PDF2Flavor = PDF1Flavor;
               }
               fPDFCoeff[k].push_back(make_pair(PDF1Flavor,PDF2Flavor));
            }
         }
      }
      else {
         error["ReadCoeffAddBase"]<<"Only IPDFCoeffFormat==0 is implemented, but IPDFCoeffFormat="<<IPDFCoeffFormat<<". Exiting."<<endl;
         exit(1);
      }
   }
   if(IPDFdef1==0){
      for(int i=0;i<NSubproc;i++){
         // Missing: linear PDF combinations for IPDFdef1=0
         if(NPDF==1){
         }else{
            if(NPDF==2){
            }
         }
      }
   }
   //Nxtot1.resize(fNObsBins);
   XNode1.resize(fNObsBins);
   for(int i=0;i<fNObsBins;i++){
      int xtot;
      table >> xtot;
      //table >> Nxtot1[i];
      //XNode1[i].resize(Nxtot1[i]);
      XNode1[i].resize(xtot);
      for(int j=0;j<xtot;j++){
         table >> XNode1[i][j];
      }
   }
   if(NPDFDim==2){
      //Nxtot2.resize(fNObsBins);
      XNode2.resize(fNObsBins);
      for(int i=0;i<fNObsBins;i++){
         int xtot;
         table >> xtot;
         XNode2[i].resize(xtot);
         //table >> Nxtot2[i];
         //XNode2[i].resize(Nxtot2[i]);
         for(int j=0;j<xtot;j++){
            table >> XNode2[i][j];
         }
      }
   }
   if(NFragFunc>0){
      Nztot.resize(fNObsBins);
      ZNode.resize(fNObsBins);
      for(int i=0;i<fNObsBins;i++){
         table >> Nztot[i];
         ZNode[i].resize(Nztot[i]);
         for(int j=0;j<Nztot[i];j++){
            table >> ZNode[i][j];
         }
      }
   }

   table >> NScales;
   table >> NScaleDim;
   Iscale.resize(NScales);
   for(int i=0;i<NScales;i++){
      table >> Iscale[i];
   }
   int NscaleDescript;
   ScaleDescript.resize(NScaleDim);
   for(int i=0;i<NScaleDim;i++){
      table >> NscaleDescript;
      ScaleDescript[i].resize(NscaleDescript);
      table.getline(buffer,256);
      for(int j=0;j<NscaleDescript;j++){
         table.getline(buffer,256);
         ScaleDescript[i][j] = buffer;
         //            StripWhitespace(ScaleDescript[i][j]);
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::Write(ostream& table) {
   debug["Write"]<<"Calling fastNLOCoeffBase::Write()"<<endl;
   fastNLOCoeffBase::Write(table);
   CheckCoeffConstants(this);
   table << IRef << endl;
   table << IScaleDep << endl;
   table << Nevt << endl;
   table << Npow << endl;
   table << NPDFPDG.size() << endl;
   for(unsigned int i=0;i<NPDFPDG.size();i++){
      table <<  NPDFPDG[i] << endl;
   }
   table << NPDFDim << endl;
   int NFragFunc = NFFPDG.size();
   table << NFragFunc << endl;
   if(NFragFunc>0){
      for(int i=0;i<NFragFunc;i++){
         table <<  NFFPDG[i] << endl;
      }
   }
   table << NFFDim << endl;
   table << NSubproc << endl;
   table << IPDFdef1 << endl;
   table << IPDFdef2 << endl;
   table << IPDFdef3 << endl;

   if(IPDFdef2==0){ // PDF linear combinations are stored herewith
      if ( IPDFdef3 != NSubproc ){
         error["Write"]<<"IPDFdef3 must be equal to NSubproc. (IPDFdef3="<<IPDFdef3<<", NSubproc="<<NSubproc<<"). Exiting."<<endl;
         exit(1);
      }
      int IPDFCoeffFormat = 0 ; // this is format style 0
      table <<  IPDFCoeffFormat << endl;
      for(int k=0;k<NSubproc;k++){
         table << fPDFCoeff[k].size() <<endl; // NPartonParis
         for( unsigned int n=0;n<fPDFCoeff[k].size();n++){
            if ( IPDFdef1>=3 ) {
               table << fPDFCoeff[k][n].first << endl;
               table << fPDFCoeff[k][n].second << endl;
            }
            else if ( IPDFdef1>=3 ) {
               table << fPDFCoeff[k][n].first << endl;
            }
         }
      }
   }

   if(IPDFdef1==0){
      for(int i=0;i<NSubproc;i++){
         // Missing: linear PDF combinations for IPDFdef1=0
         if(NPDFPDG.size()==1){
         }else{
            if(NPDFPDG.size()==2){
            }
         }
      }
   }
   for(int i=0;i<fNObsBins;i++){
      table << XNode1[i].size() << endl;
      for(unsigned int j=0;j<XNode1[i].size();j++){
         table << XNode1[i][j] << endl;
      }
   }
   if(NPDFDim==2){
      for(int i=0;i<fNObsBins;i++){
         table << XNode2[i].size() << endl;
         for(unsigned int j=0;j<XNode2[i].size();j++){
            table << XNode2[i][j] << endl;
         }
      }
   }
   if(NFragFunc>0){
      for(int i=0;i<fNObsBins;i++){
         table << Nztot[i] << endl;
         for(int j=0;j<Nztot[i];j++){
            table << ZNode[i][j] << endl;
         }
      }
   }
   int NScales = Iscale.size();
   table << NScales << endl;
   table << NScaleDim << endl;
   for(int i=0;i<NScales;i++){
      table << Iscale[i] << endl;
   }
   for(int i=0;i<NScaleDim;i++){
      table << ScaleDescript[i].size() << endl;
      for(unsigned int j=0;j<ScaleDescript[i].size();j++){
         table << ScaleDescript[i][j] << endl;
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::Add(const fastNLOCoeffAddBase& other){
   //    double w1 = (double)Nevt / (Nevt+other.Nevt);
   //    double w2 = (double)other.Nevt / (Nevt+other.Nevt);
   Nevt += other.Nevt;
}


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffAddBase::IsCompatible(const fastNLOCoeffAddBase& other) const {
   // chek CoeffBase variables
   if ( ! ((fastNLOCoeffBase*)this)->IsCompatible(other)) {
      say::debug["fastNLOCoeffAddBase::IsCompatible"]<<"fastNLOCoeffBase not compatible."<<endl;
      return false;
   }
   if ( IRef != other.GetIRef() ) {
      //warn["IsCompatible"]<<""<<endl;
      say::warn["fastNLOCoeffAddBase::IsCompatible"]<<"Different number of IRef detected."<<endl;
      return false;
   }
   if ( IScaleDep != other.GetIScaleDep() ) {
      say::warn["fastNLOCoeffAddBase::IsCompatible"]<<"Different number of IScaleDep detected."<<endl;
      return false;
   }
   if ( Npow != other.GetNpow() ) {
      say::warn["fastNLOCoeffAddBase::IsCompatible"]<<"Different number of NPow detected."<<endl;
      //warn["IsCompatible"]<<""<<endl;
      return false;
   }
   if ( GetNPDF() != other.GetNPDF() ) {
      say::warn["fastNLOCoeffAddBase::IsCompatible"]<<"Different number of NPDF detected."<<endl;
      //warn["IsCompatible"]<<""<<endl;
      return false;
   }
   if ( NSubproc != other.GetNSubproc() ) {
      say::warn["fastNLOCoeffAddBase::IsCompatible"]<<"Different numbers for NSubproc detected."<<endl;
      //warn["IsCompatible"]<<""<<endl;
      return false;
   }
   // check x-nodes briefly
   if ( fNObsBins != other.GetNObsBin() ){
      say::warn["IsCompatible"]<<"Different number of bins detected."<<endl;
      return false;
   }
   // check x-nodes briefly
   for ( int i = 0 ; i< fNObsBins ;i++ ){
      if ( GetNxmax(i) != other.GetNxmax(i) ){
	 say::warn["fastNLOCoeffAddBase::IsCompatible"]<<"Different number of x-nodes detected."<<endl;
	 return false;
      }
      if ( GetNxtot1(i) != other.GetNxtot1(i) ){
	 say::warn["fastNLOCoeffAddBase::IsCompatible"]<<"Different number of x-nodes detected."<<endl;
	 return false;
      }
      if ( GetXNode1(i,0) != other.GetXNode1(i,0) ){
	 say::warn["fastNLOCoeffAddBase::IsCompatible"]<<"Different values for x-nodes detected."<<endl;
	 return false;
      }
      if ( GetXNode1(i,1) != other.GetXNode1(i,1) ){
	 say::warn["fastNLOCoeffAddBase::IsCompatible"]<<"Different values for x-nodes detected."<<endl;
	 return false;
      }
   }
   // succesful!
   return true;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::Clear() {
   //! Clear all coefficients and event counts
   Nevt = 0;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::NormalizeCoefficients() {
   Nevt = 1;
}


//________________________________________________________________________________________________________________ //
int fastNLOCoeffAddBase::GetNxmax(int i) const {
   int nxmax = 0;
   switch (NPDFDim) {
   case 0: nxmax = (int)XNode1[i].size();
      break;
      //   case 1: nxmax = ((int)pow((double)Nxtot1[i],2)+Nxtot1[i])/2;
   case 1: nxmax = ((int)pow((double)XNode1[i].size(),2)+XNode1[i].size())/2;
      break;
   case 2: nxmax = XNode1[i].size()*XNode2[i].size();
      break;
   default: ;
   }
   return nxmax;
};


//________________________________________________________________________________________________________________ //
int fastNLOCoeffAddBase::GetXIndex(int Obsbin,int x1bin,int x2bin) const {
   int xbin = 0;
   switch (NPDFDim) {
   case 0: xbin = x1bin; // linear
      break;
   case 1: xbin = x1bin + (x2bin*(x2bin+1)/2);    // half matrix
      break;
   case 2: xbin = x1bin + x2bin * XNode1[Obsbin].size(); // full matrix
      break;
   default: ;
   }
   return xbin;
};


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::Print() const {
   fastNLOCoeffBase::Print();
   printf(" **************** FastNLO Table: fastNLOCoeffAddBase ****************\n");
   printf(" B   IRef                          %d\n",IRef);
   printf(" B   NScaleDep                     %d\n",NScaleDep);
   printf(" B   IScaleDep                     %d\n",IScaleDep);
   printf(" B   Nevt                          %f\n",Nevt);
   printf(" B   Npow                          %d\n",Npow);
   printf(" B   NPDF                          %lu\n",NPDFPDG.size());
   for(unsigned int i=0;i<NPDFPDG.size();i++){
      printf(" B    - NPDFPDG[%d]                 %d\n",i,NPDFPDG[i]);
   }
   printf(" B   NPDFDim                       %d\n",NPDFDim);
   printf(" B   NFragFunc                     %lu\n",NFFPDG.size());
   for(unsigned int i=0;i<NFFPDG.size();i++){
      printf(" B    - NFFPDG[%d]               %d\n",i,NFFPDG[i]);
   }
   printf(" B   NFFDim                        %d\n",NFFDim);
   printf(" B   NSubproc                      %d\n",NSubproc);
   printf(" B   IPDFdef1                      %d\n",IPDFdef1);
   printf(" B   IPDFdef2                      %d\n",IPDFdef2);
   printf(" B   IPDFdef3                      %d\n",IPDFdef3);
   printf(" B   Nxtot1[0-%d]             ",fNObsBins);
   for(int i=0;i<fNObsBins;i++){
      printf("%lu ,",XNode1[i].size());
   }
   printf("\n");
   //     for(int i=0;i<fNObsBins;i++){
   //       printf(" B    XNode1[%d]             ",i);
   //       for(int j=0;j<Nxtot1[i];j++){
   //   printf(" B   %8.4f ,",XNode1[i][j]);
   //       }
   //       printf(" B   \n");
   //     }
   printf(" B   if (NPDFDim==2), you could print xnodes2 here. (NPDFDim = %d)\n",NPDFDim);
   printf(" B   if (NFragFunc>0), you could print xnodes2 here. (NFragFunc = %lu)\n",NFFPDG.size());
   printf(" B   NScales                       %lu\n",Iscale.size());
   for(unsigned int i=0;i<Iscale.size();i++){
      printf(" B    - Iscale[%d]                  %d\n",i,Iscale[i]);
   }
   printf(" B   NScaleDim                     %d\n",NScaleDim);
   for(int i=0;i<NScaleDim;i++){
      for(unsigned int j=0;j<ScaleDescript[i].size();j++){
         printf(" B    -  - ScaleDescript[%d][%d]     %s\n",i,j,ScaleDescript[i][j].data());
      }
   }
   printf(" *******************************************************\n");
}


//________________________________________________________________________________________________________________ //
