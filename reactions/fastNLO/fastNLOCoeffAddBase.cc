#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffAddBase.h"
#include "fastnlotk/fastNLOTools.h"

using namespace std;
using namespace fastNLO;

//________________________________________________________________________________________________________________ //
fastNLOCoeffAddBase::fastNLOCoeffAddBase(int NObsBin)
   : fastNLOCoeffBase(NObsBin), IRef(), IScaleDep(), Nevt(), Npow(), NPDFPDG(),
     NPDFDim(), NFFPDG(), NFFDim(), NSubproc(), IPDFdef1(), IPDFdef2(), IPDFdef3(),
     fPDFCoeff(), Hxlim1(), XNode1(), Hxlim2(), XNode2(), Nztot(), Hzlim(), ZNode(),
     NScales(), NScaleDim(), Iscale(), ScaleDescript() {

}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddBase::fastNLOCoeffAddBase(const fastNLOCoeffBase& base ) : fastNLOCoeffBase(base) {
}


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
      say::error["fastNLOCoeffAddBase::CheckCoeffConstants"] << "Unknown contribution type, aborting! "
                                                             << "IAddMultFlag = " << c->GetIAddMultFlag()
                                                             << ", IDataFlag ="   << c->GetIDataFlag() <<endl;
      exit(1);
   }
}

//________________________________________________________________________________________________________________ //
fastNLOCoeffAddBase* fastNLOCoeffAddBase::Clone() const {
   //! Use has to take care to delete this object later
   return new fastNLOCoeffAddBase(*this);
}


///________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::Read(istream& table, int ITabVersionRead){
   debug["Read"]<<"Start reading table ..."<<endl;
   fastNLOCoeffBase::ReadBase(table, ITabVersionRead);
   CheckCoeffConstants(this);
   ReadCoeffAddBase(table, ITabVersionRead);
   fastNLOCoeffBase::ReadCoeffInfoBlocks(table, ITabVersionRead);
   EndReadCoeff(table, ITabVersionRead);
   debug["Read"]<<"Finished reading table ..."<<endl;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::ReadCoeffAddBase(istream& table, int ITabVersionRead){
   CheckCoeffConstants(this);
   char buffer[5257];
   //   string stest;
   // if ( fVersionRead>=24000 ) table >> stest; //"fastNLO_CoeffAddBase"
   // if ( fVersionRead>=24000 ) fastNLOTools::ReadUnused(table);
   table >> IRef;
   table >> IScaleDep;
   // if ( fVersionRead >= 24000 ) {
   //    table >> Nevt;
   //    table >> fWgt.WgtNevt;
   //    table >> fWgt.NumTable;
   //    table >> fWgt.WgtNumEv;
   //    table >> fWgt.WgtSumW2;
   //    table >> fWgt.SigSumW2;
   //    table >> fWgt.SigSum;
   //    fastNLOTools::ReadFlexibleVector ( fWgt.WgtObsSumW2, table );
   //    fastNLOTools::ReadFlexibleVector ( fWgt.SigObsSumW2, table );
   //    fastNLOTools::ReadFlexibleVector ( fWgt.SigObsSum, table );
   //    fastNLOTools::ReadFlexibleVector ( fWgt.WgtObsNumEv, table );
   // }
   // else {
      table >> Nevt;
      double readNevt = Nevt;
      if ( Nevt <= 0 ) { // v2300
         table >> Nevt;
         table >> fWgt.WgtNevt;
         if ( readNevt<=-2 ) table >> fWgt.NumTable;
         table >> fWgt.WgtNumEv;
         table >> fWgt.WgtSumW2;
         table >> fWgt.SigSumW2;
         table >> fWgt.SigSum;
         fastNLOTools::ReadFlexibleVector ( fWgt.WgtObsSumW2, table );
         fastNLOTools::ReadFlexibleVector ( fWgt.SigObsSumW2, table );
         fastNLOTools::ReadFlexibleVector ( fWgt.SigObsSum, table );
         fastNLOTools::ReadFlexibleVector ( fWgt.WgtObsNumEv, table );
      }
   // }
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

   sub_enabled.clear();
   sub_enabled.resize(NSubproc, true); // enable all subprocesses by default

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
               else if ( IPDFdef1==2 ) {
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

   // if ( fVersionRead>=24000 ) fastNLOTools::ReadUnused(table);
   // if ( fVersionRead>=24000 ) fastNLOTools::ReadUnused(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::Write(ostream& table, int itabversion) {
   debug["Write"]<<"Calling fastNLOCoeffBase::Write()"<<endl;
   fastNLOCoeffBase::Write(table,itabversion);
   CheckCoeffConstants(this);
   // if ( itabversion >= 24000 ) table << "fastNLO_CoeffAddBase" << sep;
   // if ( itabversion >= 24000 ) table << 0 << sep; // v2.4, but yet unused
   table << IRef << sep;
   table << IScaleDep << sep;
   if ( itabversion==23000 || itabversion==23500 || itabversion==23600 || itabversion==25000 ) { // detailed storage of weights
      if ( itabversion==23000 || itabversion==23500 ) table << -1 << sep; // -1: read the values below
      else table << -2 << sep; // -1: read the values below
      table << Nevt << sep;
      table << fWgt.WgtNevt << sep;
      if ( itabversion >= 23600 ) table << fWgt.NumTable << sep;
      table << fWgt.WgtNumEv << sep;
      table << fWgt.WgtSumW2 << sep;
      table << fWgt.SigSumW2 << sep;
      table << fWgt.SigSum << sep;
      fastNLOTools::WriteFlexibleVector ( fWgt.WgtObsSumW2, table );
      fastNLOTools::WriteFlexibleVector ( fWgt.SigObsSumW2, table );
      fastNLOTools::WriteFlexibleVector ( fWgt.SigObsSum, table );
      fastNLOTools::WriteFlexibleVector ( fWgt.WgtObsNumEv, table );
   }
   // else if ( itabversion>=24000 ) { // detailed storage of weights
   //    table << Nevt << sep;
   //    table << fWgt.WgtNevt << sep;
   //    table << fWgt.NumTable << sep;
   //    table << fWgt.WgtNumEv << sep;
   //    table << fWgt.WgtSumW2 << sep;
   //    table << fWgt.SigSumW2 << sep;
   //    table << fWgt.SigSum << sep;
   //    fastNLOTools::WriteFlexibleVector ( fWgt.WgtObsSumW2, table );
   //    fastNLOTools::WriteFlexibleVector ( fWgt.SigObsSumW2, table );
   //    fastNLOTools::WriteFlexibleVector ( fWgt.SigObsSum, table );
   //    fastNLOTools::WriteFlexibleVector ( fWgt.WgtObsNumEv, table );
   // }
   else {
      table << Nevt << sep;
   }
   table << Npow << sep;
   table << NPDFPDG.size() << sep;
   for(unsigned int i=0;i<NPDFPDG.size();i++){
      table <<  NPDFPDG[i] << sep;
   }
   table << NPDFDim << sep;
   int NFragFunc = NFFPDG.size();
   table << NFragFunc << sep;
   if(NFragFunc>0){
      for(int i=0;i<NFragFunc;i++){
         table <<  NFFPDG[i] << sep;
      }
   }
   table << NFFDim << sep;
   table << NSubproc << sep;
   table << IPDFdef1 << sep;
   table << IPDFdef2 << sep;
   table << IPDFdef3 << sep;

   if(IPDFdef2==0){ // PDF linear combinations are stored herewith
      info["Write"]<<"Writing PDF coefficients into table."<<endl;
      if ( IPDFdef3 != NSubproc ){
         error["Write"]<<"IPDFdef3 must be equal to NSubproc. (IPDFdef3="<<IPDFdef3<<", NSubproc="<<NSubproc<<"). Exiting."<<endl;
         exit(1);
      }
      int IPDFCoeffFormat = 0 ; // this is format style 0
      table <<  IPDFCoeffFormat << sep;
      for(int k=0;k<NSubproc;k++){
         table << fPDFCoeff[k].size() <<sep; // NPartonParis
         for( unsigned int n=0;n<fPDFCoeff[k].size();n++){
            if ( IPDFdef1>=3 ) {
               table << fPDFCoeff[k][n].first << sep;
               table << fPDFCoeff[k][n].second << sep;
            }
            else if ( IPDFdef1==2 ) {
               table << fPDFCoeff[k][n].first << sep;
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
      table << XNode1[i].size() << sep;
      for(unsigned int j=0;j<XNode1[i].size();j++){
         table << XNode1[i][j] << sep;
      }
   }
   if(NPDFDim==2){
      for(int i=0;i<fNObsBins;i++){
         table << XNode2[i].size() << sep;
         for(unsigned int j=0;j<XNode2[i].size();j++){
            table << XNode2[i][j] << sep;
         }
      }
   }
   if(NFragFunc>0){
      for(int i=0;i<fNObsBins;i++){
         table << Nztot[i] << sep;
         for(int j=0;j<Nztot[i];j++){
            table << ZNode[i][j] << sep;
         }
      }
   }
   int NScales = Iscale.size();
   table << NScales << sep;
   table << NScaleDim << sep;
   for(int i=0;i<NScales;i++){
      table << Iscale[i] << sep;
   }
   for(int i=0;i<NScaleDim;i++){
      table << ScaleDescript[i].size() << sep;
      for(unsigned int j=0;j<ScaleDescript[i].size();j++){
         table << ScaleDescript[i][j] << sep;
      }
   }
   // if ( itabversion>=24000 ) table << 0 << sep; // v2.4, but yet unused
   // if ( itabversion>=24000 ) table << 0 << sep; // v2.4, but yet unused

}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::Add(const fastNLOCoeffAddBase& other, fastNLO::EMerge moption){
   //    double w1 = (double)Nevt / (Nevt+other.Nevt);
   //    double w2 = (double)other.Nevt / (Nevt+other.Nevt);
   if ( Nevt==1 || other.GetNevt()==1 ) {
      if ( moption != fastNLO::kAdd && moption != fastNLO::kUnweighted && moption != fastNLO::kAttach ) {
         error["Add"]<<"Table has weight 1, which is invalid for mergeing purposes."<<endl;
         error["Add"]<<"Possibly, the table is a result from a previous 'append' or 'unweighted' mergeing."<<endl;
         exit(4);
      }
   }


   if ( moption==fastNLO::kAttach ) {
      //Nevt = Nevt;// stays!
      fWgt.WgtNevt  = 1;
      fWgt.NumTable += other.fWgt.NumTable;
      fWgt.WgtNumEv += other.fWgt.WgtNumEv;
      fWgt.WgtSumW2 += other.fWgt.WgtSumW2;
      fWgt.SigSumW2 += other.fWgt.SigSumW2;
      fWgt.SigSum   += other.fWgt.SigSum;
      for ( unsigned int iAddProc = 0 ; iAddProc<other.fWgt.WgtObsSumW2.size() ; iAddProc++ ) {
         fWgt.WgtObsSumW2.push_back(other.fWgt.WgtObsSumW2[iAddProc]);
         fWgt.SigObsSumW2.push_back(other.fWgt.SigObsSumW2[iAddProc]);
         fWgt.SigObsSum.  push_back(other.fWgt.SigObsSum[iAddProc]);
         fWgt.WgtObsNumEv.push_back(other.fWgt.WgtObsNumEv[iAddProc]);
      }
      if ( other.GetPDFCoeff().size() ==0  || this->GetPDFCoeff().size()==0 ) {
         error["Add"]<<"Mergeing option 'kAttach' requires that PDF coefficients are stored in table.!"<<endl;
      }
      for ( unsigned int iAddProc = 0 ; iAddProc<other.GetPDFCoeff().size() ;iAddProc++ ) {
         fPDFCoeff.push_back(other.GetPDFCoeff()[iAddProc]);
      }
      NSubproc += other.GetNSubproc();
      IPDFdef3 += other.GetIPDFdef3();
   }
   else {
      Nevt += other.Nevt;
      fWgt.WgtNevt  += other.fWgt.WgtNevt;
      fWgt.NumTable += other.fWgt.NumTable;
      fWgt.WgtNumEv += other.fWgt.WgtNumEv;
      fWgt.WgtSumW2 += other.fWgt.WgtSumW2;
      fWgt.SigSumW2 += other.fWgt.SigSumW2;
      fWgt.SigSum   += other.fWgt.SigSum;
      fastNLOTools::AddVectors( fWgt.WgtObsSumW2, other.fWgt.WgtObsSumW2 );
      fastNLOTools::AddVectors( fWgt.SigObsSumW2, other.fWgt.SigObsSumW2 );
      fastNLOTools::AddVectors( fWgt.SigObsSum,   other.fWgt.SigObsSum );
      fastNLOTools::AddVectors( fWgt.WgtObsNumEv, other.fWgt.WgtObsNumEv );
   }

}


//________________________________________________________________________________________________________________ //
double fastNLOCoeffAddBase::GetMergeWeight(fastNLO::EMerge moption, int proc, int bin) const {

   //!< Get a bin and subprocess dependent weight for merging puprposes.
   if      ( moption == kMerge    )   return fWgt.WgtNevt; // Nevt
   else if ( moption == kUnweighted ) return fWgt.NumTable;
   else if ( moption == kAdd )     return 1;
   else if ( moption == kNumEvent )   return double(fWgt.WgtNumEv);
   else if ( moption == kSumW2    )   return fWgt.WgtSumW2;
   else if ( moption == kSumSig2  )   return fWgt.SigSumW2;
   else if ( moption == kSumUser  )   return fWgt.SigSum;
   else if ( moption == kNumEventBinProc ) return double(fWgt.WgtObsNumEv[proc][bin]);
   else if ( moption == kSumW2BinProc    ) return fWgt.WgtObsSumW2[proc][bin];
   else if ( moption == kSumSig2BinProc  ) return fWgt.SigObsSumW2[proc][bin];
   else if ( moption == kSumUserBinProc  ) return fWgt.SigObsSum[proc][bin];
   error["GetMergeWeight"]<<"Weighting option not recognized: "<<moption<<endl;
   exit(4);
   return 0;
}


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffAddBase::IsCompatible(const fastNLOCoeffAddBase& other) const {
   // chek CoeffBase variables
   if ( ! ((fastNLOCoeffBase*)this)->IsCompatible(other)) {
      debug["IsCompatible"]<<"fastNLOCoeffBase not compatible."<<endl;
      return false;
   }
   if ( IRef != other.GetIRef() ) {
      //warn["IsCompatible"]<<""<<endl;
      warn["IsCompatible"]<<"Different number of IRef detected."<<endl;
      return false;
   }
   if ( Nevt * other.Nevt < 0 ) {
      // skip, if the two tables store the event weights in different formats
      // If this is needed, simple solutions are thinkable
      warn["IsCompatible"]<<"Tables use different format for normalisation."<<endl;
      return false;
   }
   if ( IScaleDep != other.GetIScaleDep() ) {
      warn["IsCompatible"]<<"Different number of IScaleDep detected."<<endl;
      return false;
   }
   if ( Npow != other.GetNpow() ) {
      warn["IsCompatible"]<<"Different number of NPow detected."<<endl;
      return false;
   }
   if ( GetNPDF() != other.GetNPDF() ) {
      warn["IsCompatible"]<<"Different number of NPDF detected."<<endl;
      return false;
   }
   if ( NSubproc != other.GetNSubproc() ) {
      warn["IsCompatible"]<<"Different numbers for NSubproc detected."<<endl;
      //return false;
      warn["IsCompatible"]<<"Continuing! (experimental: This is needed for kAttach, but may causes bugs otherwise. Please report!)"<<endl;
   }
   // check x-nodes briefly
   if ( fNObsBins != other.GetNObsBin() ){
      warn["IsCompatible"]<<"Different number of bins detected."<<endl;
      return false;
   }
   // check x-nodes briefly
   for ( int i = 0 ; i< fNObsBins ;i++ ){
      if ( GetNxmax(i) != other.GetNxmax(i) ){
         error["IsCompatible"]<<"Different number of x-nodes detected: "<<GetNxmax(i)<<" <-> "<<other.GetNxmax(i)<<endl;
         return false;
      }
      if ( GetNxtot1(i) != other.GetNxtot1(i) ){
         error["IsCompatible"]<<"Different number of x-nodes detected: "<<GetNxtot1(i)<<" <-> "<<other.GetNxtot1(i)<<endl;
         return false;
      }
      if ( GetXNode1(i,0) != other.GetXNode1(i,0) ){
         warn["IsCompatible"]<<"Different values for x-nodes detected. Lowest x-node: "<<GetXNode1(i,0)<<" <-> "<<other.GetXNode1(i,0)<<endl;
         return false;
      }
      // if ( GetXNode1(i,1) != other.GetXNode1(i,1) ){
      //    warn["IsCompatible"]<<"Different values for x-nodes detected."<<endl;
      //    return false;
      // }
   }
   // succesful!
   return true;
}


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffAddBase::IsCatenable(const fastNLOCoeffAddBase& other) const {
   // check CoeffBase variables
   if ( ! ((fastNLOCoeffBase*)this)->IsCatenable(other)) {
      debug["IsCatenable"]<<"fastNLOCoeffBase not compatible. Skipped."<<endl;
      return false;
   }
   if ( Nevt * other.Nevt < 0 ) {
      // skip, if the two tables store the event weights in different formats
      // If this is needed, simple solutions are thinkable
      debug["IsCatenable"]<<"Tables use different format for table normalisation. Skipped."<<endl;
      return false;
   }
   if ( IRef != other.GetIRef() ) {
      debug["IsCatenable"]<<"Different number of IRef detected. Skipped."<<endl;
      return false;
   }
   if ( IScaleDep != other.GetIScaleDep() ) {
      debug["IsCatenable"]<<"Different number of IScaleDep detected. Skipped."<<endl;
      return false;
   }
   if ( Npow != other.GetNpow() ) {
      debug["IsCatenable"]<<"Different number of NPow detected. Skipped."<<endl;
      return false;
   }
   if ( GetNPDF() != other.GetNPDF() ) {
      debug["IsCatenable"]<<"Different number of NPDF detected. Skipped."<<endl;
      return false;
   }
   if ( NSubproc != other.GetNSubproc() ) {
      debug["IsCatenable"]<<"Different numbers for NSubproc detected. Skipped."<<endl;
      return false;
   }
   // check x-nodes briefly
   // for ( int i = 0 ; i< fNObsBins ;i++ ){
   //    if ( GetNxmax(i) != other.GetNxmax(i) ){
   //       debug["IsCatenable"]<<"Different number of x-nodes detected."<<endl;
   //       return false;
   //    }
   //    if ( GetNxtot1(i) != other.GetNxtot1(i) ){
   //       debug["IsCatenable"]<<"Different number of x-nodes detected."<<endl;
   //       return false;
   //    }
   //    if ( GetXNode1(i,0) != other.GetXNode1(i,0) ){
   //       debug["IsCatenable"]<<"Different values for x-nodes detected."<<endl;
   //       return false;
   //    }
   //    if ( GetXNode1(i,1) != other.GetXNode1(i,1) ){
   //       debug["IsCatenable"]<<"Different values for x-nodes detected."<<endl;
   //       return false;
   //    }
   // }
   info["IsCatenable"]<<"Base parameters of additive contribution allow catenation"<<endl;
   return true;
}

//________________________________________________________________________________________________________________
bool fastNLOCoeffAddBase::SubSelect( vector< pair<int,int> > processes, bool on ) {
   vector<int> s;
   s.clear();
   for ( unsigned int k=0; k<processes.size(); k++ ) {
      pair<int,int> p = processes[k];
      // search for selected process in fPDFCoeff
      for ( unsigned int i = 0; i<fPDFCoeff.size(); i++ ) {
         for ( unsigned int j = 0; j<fPDFCoeff[i].size(); j++ ) {
            if ( p == fPDFCoeff[i][j] ) {
               // found! now check if the other proesses in this subcontribution are also to be selected
               vector< pair<int,int> > p_list = fPDFCoeff[i];
               bool f = true;
               for ( unsigned int n = 0; n<p_list.size(); n++ ) {
                  bool ff = false;
                  for ( unsigned int m = 0; m<processes.size(); m++ )
                     if ( p_list[n] == processes[m] )
                        ff = true;
                  f &= ff;
               }
               if (!f)
                  return false;
               s.push_back(i);
            }
         }
      }
      // uncomment to throw error on non-existing pairs
      //if (!fff)
      //         return false;
   }
   // now activate the selected subcontributions and return succes (true)
   for ( unsigned int k = 0; k<s.size(); k++ )
      SubEnable(s[k], on);

   return true;
}



//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::Clear() {
   //! Clear all coefficients and event counts
   Nevt = 0;
   fWgt.WgtNevt = 0;
   fWgt.NumTable = 1;
   fWgt.WgtNumEv = 0;
   fWgt.WgtSumW2 = 0;
   fWgt.SigSumW2 = 0;
   fWgt.SigSum   = 0;
   fastNLOTools::ClearVector(fWgt.WgtObsSumW2);
   fastNLOTools::ClearVector(fWgt.SigObsSumW2);
   fastNLOTools::ClearVector(fWgt.SigObsSum);
   fastNLOTools::ClearVector(fWgt.WgtObsNumEv);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::NormalizeCoefficients(double wgt) {
   Nevt = wgt;
   // Don't touch other weights.
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
double fastNLOCoeffAddBase::GetX1(int iObsBin, int iXnode) const {
   // return x-value of PDF1 at node iXnode
   switch (NPDFDim) {
   case 0:
      return GetXNode1(iObsBin,iXnode);
   case 1:
      //
      // cout<<"GetX1 not implemented for half-matrix notation!"<<endl;
      //
      return 1;
   case 2:
      return GetXNode1(iObsBin, iXnode % GetNxtot1(iObsBin) );
   default: return 1;
   }
   return 1;
}

//________________________________________________________________________________________________________________ //
double fastNLOCoeffAddBase::GetX2(int iObsBin, int iXnode) const {
   // return x-value of PDF1 at node iXnode
   switch (NPDFDim) {
   case 0:
      return 1;
   case 1:
      //
      //cout<<"GetX2 not implemented for half-matrix notation!"<<endl;
      //
      return 1;
   case 2:
      return GetXNode2(iObsBin, iXnode / GetNxtot1(iObsBin) );
   default: return 1;
   }
   return 1;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::Print(int iprint) const {
   if ( !(iprint < 0) ) {
      fastNLOCoeffBase::Print(iprint);
      cout << fastNLO::_DSEP20C << " fastNLO Table: CoeffAddBase " << fastNLO::_DSEP20 << endl;
   } else {
      cout << endl << fastNLO::_CSEP20C << " fastNLO Table: CoeffAddBase " << fastNLO::_CSEP20 << endl;
   }
   printf(" # No. of events (Nevt)                %f\n",Nevt);
   if ( fWgt.WgtNevt!= 0 || fWgt.WgtSumW2!= 0 ) {
      printf(" # Weight of table [=Nevt] (fWgtNevt)  %f\n",fWgt.WgtNevt);
      printf(" # Number of tables merged together    %d\n",fWgt.NumTable);
      printf(" # No. of filled events (WgtNumEv)     %llu\n",fWgt.WgtNumEv);
      printf(" # Sum of weights squared (WgtSumW2)   %f\n",fWgt.WgtSumW2);
      printf(" # Sum of sigma squared (SigSumW2)     %f\n",fWgt.SigSumW2);
      printf(" # Sum of sigma (SigSum)               %f\n",fWgt.SigSum);
      printf(" # Sigma / Nevt (SigSum/WgtNevt)       %f\n",fWgt.SigSum/fWgt.WgtNevt);
   }

   printf(" # Abs. order in a_s (Npow)            %d\n",Npow);
   printf(" # No. of hadrons involved (NPDF)      %lu\n",NPDFPDG.size());
   fastNLOTools::PrintVector(NPDFPDG,"Type(s) of hadrons (NPDFPDG)","#");

   printf(" # PDF storage format (NPDFDim)        %d\n",NPDFDim);
   if ( NPDFDim == 0 ) {
      printf(" #   --> x-interpolation storage format: Linear\n");
   } else if ( NPDFDim == 1 ) {
      printf(" #   --> x-interpolation storage format: Half-Matrix\n");
   } else if ( NPDFDim == 2 ) {
      printf(" #   --> x-interpolation storage format: Full-Matrix\n");
   } else {
      error["Print"] << "Unknown interpolation storage structure, aborting! "
                                          << " NPDFDim = " << NPDFDim << endl;
   }

   for (int i=0; i<fNObsBins; i++) {
      // Print only for first and last observable bin
      if (i==0 || i==fNObsBins-1) {
         printf(" # Observable bin no. %d\n",i+1);
         printf(" #   No. of X1 nodes (XNode1[i].size())          %d\n",(int)GetXNodes1(i).size());
      }
   }
   printf(" # No. of scales involved (NScales)    %lu\n",Iscale.size());
   for(int i=0;i<NScaleDim;i++){
      fastNLOTools::PrintVector(ScaleDescript[i],"Scale descriptions (ScaleDescript)","#");
   }
   if ( std::abs(iprint) > 0 ) {
      cout << fastNLO::_SSEP20C << " Extended information (iprint > 0) " << fastNLO::_SSEP20 << endl;
      if ( NScales > 0 ) {fastNLOTools::PrintVector(Iscale,"Iscale (Unused, always 0) (Iscale)","#  ");}
      printf(" #   IRef                              %d\n",IRef);
      printf(" #   IScaleDep (Unused, always 0)      %d\n",IScaleDep);
      printf(" #   NFragFunc                         %lu\n",NFFPDG.size());
      if ( NFFPDG.size() > 0 ) {fastNLOTools::PrintVector(NFFPDG,"Type(s) of hadrons (NFFPDG)","#");}
      printf(" #   NFFDim                            %d\n",NFFDim);
      printf(" #   NScaleDim                         %d\n",NScaleDim);
      printf(" #   NSubproc                          %d\n",NSubproc);
      printf(" #   IPDFdef1                          %d\n",IPDFdef1);
      printf(" #   IPDFdef2                          %d\n",IPDFdef2);
      printf(" #   IPDFdef3                          %d\n",IPDFdef3);
      char buffer[1024];
      for (int i=0; i<fNObsBins; i++) {
         // Print only for first and last observable bin
         if (i==0 || i==fNObsBins-1) {
            printf(" #   Observable bin no. %d\n",i+1);
            snprintf(buffer, sizeof(buffer), "X1 nodes (XNode1[%d][])",i);
            fastNLOTools::PrintVector(GetXNodes1(i),buffer,"#    ");
            if ( NPDFDim == 2 ) {
               snprintf(buffer, sizeof(buffer), "X2 nodes (XNode2[%d][])",i);
               fastNLOTools::PrintVector(GetXNodes2(i),buffer,"#    ");
            }
         }
      }
   }
   cout << " # No. of subprocess contributions     " << fPDFCoeff.size() << endl;
   cout << " #ISub  " << "Processes" << endl;
   for ( unsigned int i = 0; i<fPDFCoeff.size(); i++ ) {
      cout << " #";
      int w = cout.width();
      cout.width(4);
      cout << i;
      cout.width(w);
      cout << ": ";
      for ( unsigned int j = 0; j<fPDFCoeff[i].size(); j++ ) {
         cout << "(" << fPDFCoeff[i][j].first << "," << fPDFCoeff[i][j].second << ")" << ", ";
      }
      cout << endl;
   }
   if ( iprint < 0 ) {
      cout << fastNLO::_CSEPSC << endl;
   } else {
      //      cout << fastNLO::_DSEPSC << endl;
   }
}


//________________________________________________________________________________________________________________ //

// Erase observable bin
void fastNLOCoeffAddBase::EraseBin(unsigned int iObsIdx, int ITabVersionRead) {
   debug["EraseBin"]<<"Erasing observable bin in CoeffAddBase with bin index " << iObsIdx << endl;
   if ( XNode1.size() == 0 ) {
      error["EraseBin"]<<"All additive contribution bins deleted already. Aborted!" << endl;
      exit(1);
   }
   if ( XNode1.size() != 0 ) XNode1.erase(XNode1.begin()+iObsIdx);
   if ( NPDFDim==2 && XNode2.size() != 0 ) XNode2.erase(XNode2.begin()+iObsIdx);
   for ( unsigned int i = 0 ; i<fWgt.WgtObsSumW2.size() ; i++ ) {
      fWgt.WgtObsSumW2[i].erase(fWgt.WgtObsSumW2[i].begin()+iObsIdx);
      fWgt.SigObsSumW2[i].erase(fWgt.SigObsSumW2[i].begin()+iObsIdx);
      fWgt.SigObsSum[i].  erase(fWgt.SigObsSum  [i].begin()+iObsIdx);
      fWgt.WgtObsNumEv[i].erase(fWgt.WgtObsNumEv[i].begin()+iObsIdx);
   }
   if ( ! ( ITabVersionRead < 25000 ) ) {
      if ( NCoeffInfoBlocks > 0 ) {
         debug["EraseBin"]<<"Found " << NCoeffInfoBlocks << " InfoBlocks with bins to be erased, too." << endl;
         for ( int i=0; i < NCoeffInfoBlocks; i++ ) {
            if ( ICoeffInfoBlockFlag1[i] == 0 && ICoeffInfoBlockFlag2[i] == 0 ) {
               CoeffInfoBlockContent[i].erase(CoeffInfoBlockContent[i].begin()+iObsIdx);
               NCoeffInfoBlockCont[i] = NCoeffInfoBlockCont[i] - 1;
            } else {
               error["EraseBin"]<<"Erase bin not yet implemented for InfoBlocks other than with flags 1,2 = 0, 0:" <<
                  ICoeffInfoBlockFlag1[i] << ", " << ICoeffInfoBlockFlag2[i] << ", aborted!" << endl;
               exit(567);
            }
         }
      }
   }
   fastNLOCoeffBase::EraseBin(iObsIdx);
}

// Catenate observable bin
void fastNLOCoeffAddBase::CatBin(const fastNLOCoeffAddBase& other, unsigned int iObsIdx, int ITabVersionRead) {
   debug["CatBin"]<<"Catenating observable bin to CoeffAddBase with bin index " << iObsIdx << endl;
   if ( XNode1.size() == 0 ) {
      error["CatBin"]<<"Initial additive table is empty. Aborted!" << endl;
      exit(1);
   }
   if ( XNode1.size() != 0 ) XNode1.push_back(other.XNode1[iObsIdx]);
   if ( NPDFDim==2 &&  XNode2.size() != 0 ) XNode2.push_back(other.XNode2[iObsIdx]);
   for ( unsigned int i = 0 ; i<fWgt.WgtObsSumW2.size() ; i++ ) {
      fWgt.WgtObsSumW2[i].push_back(other.fWgt.WgtObsSumW2[i][iObsIdx]);
      fWgt.SigObsSumW2[i].push_back(other.fWgt.SigObsSumW2[i][iObsIdx]);
      fWgt.SigObsSum[i].  push_back(other.fWgt.SigObsSum  [i][iObsIdx]);
      fWgt.WgtObsNumEv[i].push_back(other.fWgt.WgtObsNumEv[i][iObsIdx]);
   }
   if ( ! ( ITabVersionRead < 25000 ) ) {
      if ( NCoeffInfoBlocks > 0 ) {
         debug["CatBin"]<<"Found " << NCoeffInfoBlocks << " InfoBlocks with bins to be catenated, too." << endl;
         for ( int i=0; i < NCoeffInfoBlocks; i++ ) {
            if ( ICoeffInfoBlockFlag1[i] == 0 && ICoeffInfoBlockFlag2[i] == 0 ) {
               CoeffInfoBlockContent[i].push_back(other.CoeffInfoBlockContent[i][iObsIdx]);
               NCoeffInfoBlockCont[i] = NCoeffInfoBlockCont[i] + 1;
            } else {
               error["CatBin"]<<"Catenate bins not yet implemented for InfoBlocks other than with flags 1,2 = 0, 0:" <<
                  ICoeffInfoBlockFlag1[i] << ", " << ICoeffInfoBlockFlag2[i] << ", aborted!" << endl;
               exit(678);
            }
         }
      }
   }
   fastNLOCoeffBase::CatBin(other, iObsIdx);
}
