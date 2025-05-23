#include <cstdlib>
#include <fastnlotk/fastNLOEvent.h>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "fastnlotk/fastNLOCoeffAddFix.h"
#include "fastnlotk/fastNLOTools.h"
#include "fastnlotk/speaker.h"

using namespace std;
using namespace fastNLO;

//________________________________________________________________________________________________________________ //
fastNLOCoeffAddFix::fastNLOCoeffAddFix(int NObsBin)
   : fastNLOCoeffAddBase(NObsBin), Nscalevar(), ScaleFac(), ScaleNode(), SigmaTilde(),
     AlphasTwoPi_v20(), PdfLc(), PdfSplLc1(), PdfSplLc2() {
   SetClassName("fastNLOCoeffAddFix");
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddFix::fastNLOCoeffAddFix(const fastNLOCoeffBase& base) : fastNLOCoeffAddBase(base) {
   SetClassName("fastNLOCoeffAddFix");
}



//________________________________________________________________________________________________________________ //
bool fastNLOCoeffAddFix::CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet) {
   bool ret = fastNLOCoeffAddBase::CheckCoeffConstants(c,quiet);
   if ( ret && c->GetNScaleDep() == 0 ) return true;
   else if ( c->GetNScaleDep() >= 3 ) {
      if ( !quiet)
         say::error["fastNLOCoeffAddFix::CheckCoeffConstants"]
            <<"This is not a fixed order v2.0  table. NScaleDep must be equal 0 but is NScaleDep="
            <<c->GetNScaleDep()<<endl;
      return false;
   }
   else return false;
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddFix* fastNLOCoeffAddFix::Clone() const {
   //! Use has to take care to delete this object later
   return new fastNLOCoeffAddFix(*this);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Read(istream& table, int ITabVersionRead){
   fastNLOCoeffBase::ReadBase(table, ITabVersionRead);
   ReadRest(table, ITabVersionRead);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ReadRest(istream& table, int ITabVersionRead){
   CheckCoeffConstants(this);
   fastNLOCoeffAddBase::ReadCoeffAddBase(table, ITabVersionRead);
   ReadCoeffAddFix(table, ITabVersionRead);
   fastNLOCoeffBase::ReadCoeffInfoBlocks(table, ITabVersionRead);
   EndReadCoeff(table, ITabVersionRead);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ReadCoeffAddFix(istream& table, int ITabVersionRead){
   debug["ReadCoeffAddFix::ReadCoeffAddFix"]<<"Start reading coefficients for table version "<<ITabVersionRead<<endl;
   CheckCoeffConstants(this);
   Nscalevar.resize(NScaleDim);
   vector<int> Nscalenode(NScaleDim);
   for(int i=0;i<NScaleDim;i++){
      table >> Nscalevar[i];
      table >> Nscalenode[i];
   }
   // printf("  *  fastNLOCoeffAddFix::Read().bins %d, NScalevar[0] %d, Nscalenode[0] %d,  NScaleDim %d  \n",
   // fNObsBins, Nscalevar[0] , Nscalenode[0] , NScaleDim );
   // pre-binary
   // ScaleFac.resize(NScaleDim);
   // for(int i=0;i<NScaleDim;i++){
   //    ScaleFac[i].resize(Nscalevar[i]);
   //    for(int j=0;j<Nscalevar[i];j++){
   //       table >> ScaleFac[i][j];
   //    }
   // }

   ScaleFac.resize(NScaleDim);
   for(int i=0;i<NScaleDim;i++){
      ScaleFac[i].resize(Nscalevar[i]);
   }
   fastNLOTools::ReadVector( ScaleFac , table , 1 );

   fastNLOTools::ResizeVector( ScaleNode , fNObsBins, 1 , Nscalevar[0] , Nscalenode[0] );
   int nsn = fastNLOTools::ReadVector( ScaleNode , table );
   debug["fastNLOCoeffAddFix::Read()"]<<"Read "<<nsn<<" lines of ScaleNode."<<endl;

   ResizeSigmaTilde();
   ResizePdfLC();
   ResizePdfSplLC();
   int nst = fastNLOTools::ReadVector( SigmaTilde , table , Nevt);
   debug["fastNLOCoeffAddFix::Read()"]<<"Read "<<nsn+nst<<" lines of fastNLO v2 tables."<<endl;

   // prepare members for evaluation
   fastNLOTools::ResizeVector(AlphasTwoPi_v20 , fNObsBins, GetTotalScalenodes());
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ResizeSigmaTilde(){
   //! resize SigmaTilde and ensure that all entries are empty
   SigmaTilde.resize(fNObsBins);
   for( int i=0 ; i<fNObsBins ; i++ ){
      int nxmax = GetNxmax(i);
      SigmaTilde[i].resize(GetTotalScalevars());
      for( int k=0 ; k<GetTotalScalevars() ; k++ ){
         SigmaTilde[i][k].resize(GetTotalScalenodes());
         for( int l=0 ; l<GetTotalScalenodes() ; l++ ){
            //ResizeVector(SigmaTilde[i][k][l],nxmax,NSubproc);
            SigmaTilde[i][k][l].resize(nxmax);
            for( int m=0 ; m<nxmax ; m++ ){
               SigmaTilde[i][k][l][m].resize(NSubproc);
               for( int n=0 ; n<NSubproc ; n++ ){
                  SigmaTilde[i][k][l][m][n] = 0.;
               }
            }
         }
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ResizePdfLC(){
   //! resize PdfLC
   PdfLc.resize(fNObsBins);
   for( int i=0 ; i<fNObsBins ; i++ ){
      int nxmax = GetNxmax(i);
      int totalscalenodes = GetTotalScalenodes();
      PdfLc[i].resize(totalscalenodes);
      for( int l=0 ; l<totalscalenodes ; l++ ){
         fastNLOTools::ResizeVector(PdfLc[i][l],nxmax,NSubproc);
      }
   }
}

//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ResizePdfSplLC(){
   //! resize PdfSplLC
   PdfSplLc1.resize(fNObsBins);
   PdfSplLc2.resize(fNObsBins);
   for( int i=0 ; i<fNObsBins ; i++ ){
      int nxmax = GetNxmax(i);
      int totalscalenodes = GetTotalScalenodes();
      PdfSplLc1[i].resize(totalscalenodes);
      PdfSplLc2[i].resize(totalscalenodes);
      for( int l=0 ; l<totalscalenodes ; l++ ){
         fastNLOTools::ResizeVector(PdfSplLc1[i][l],nxmax,NSubproc);
         fastNLOTools::ResizeVector(PdfSplLc2[i][l],nxmax,NSubproc);
      }
   }
}



//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Write(ostream& table, int itabversion){
   //! Write coefficient table to disk (ostream)
   CheckCoeffConstants(this);
   fastNLOCoeffAddBase::Write(table,itabversion);

   for(int i=0;i<NScaleDim;i++){
      table << Nscalevar[i] << sep;
      table << GetNScaleNode() << sep;
   }
   fastNLOTools::WriteVector( ScaleFac , table );
   int nsn = fastNLOTools::WriteVector( ScaleNode , table );
   int nst = fastNLOTools::WriteVector( SigmaTilde , table , Nevt);
   info["Write"]<<"Wrote "<<nst+nsn<<" lines into fastNLO table."<<endl;

   if ( itabversion >= 25000 ) {
      fastNLOCoeffBase::WriteCoeffInfoBlocks(table, itabversion);
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Add(const fastNLOCoeffAddBase& other, fastNLO::EMerge moption){
   //! Add another coefficient table to this table
   bool ok = CheckCoeffConstants(this);
   if ( !ok ) {
      error["Add"]<<"Cannot add tables."<<endl;
      return;
   }
   const fastNLOCoeffAddFix& othfix = (const fastNLOCoeffAddFix&)other;

   fastNLO::v5d otherSigmaTilde = othfix.SigmaTilde;
   for(int obsBin=0 ; obsBin<fNObsBins ; obsBin++){
      for(int scalevar=0 ; scalevar<GetTotalScalevars() ; scalevar++){
         for(int scalenode=0 ; scalenode<GetTotalScalenodes() ; scalenode++){
            fastNLO::v2d& thisX = SigmaTilde[obsBin][scalevar][scalenode];
            fastNLO::v2d& otherX = otherSigmaTilde[obsBin][scalevar][scalenode];
            int thisNxtot1 = GetNxtot1(obsBin);
            int thisNxtot2 = GetNxtot2(obsBin);
            int otherNxtot1 = othfix.GetNxtot1(obsBin);
            int otherNxtot2 = othfix.GetNxtot2(obsBin);
            if (thisNxtot1 < otherNxtot1) {
               fastNLOTools::ExtendSigmaTildeX(thisX, thisNxtot1, otherNxtot1,
                  thisNxtot2, thisNxtot2, NPDFDim, fastNLO::v1d(NSubproc));
            } else if (thisNxtot1 > otherNxtot1) {
               fastNLOTools::ExtendSigmaTildeX(otherX, otherNxtot1, thisNxtot1,
                  otherNxtot2, otherNxtot2, NPDFDim, fastNLO::v1d(NSubproc));
            }
            int maxNxtot1 = max(thisNxtot1, otherNxtot1);
            if (thisNxtot2 < otherNxtot2) {
               fastNLOTools::ExtendSigmaTildeX(thisX, maxNxtot1, maxNxtot1,
                  thisNxtot2, otherNxtot2, NPDFDim, fastNLO::v1d(NSubproc));
            } else if (thisNxtot2 > otherNxtot2) {
               fastNLOTools::ExtendSigmaTildeX(otherX, maxNxtot1, maxNxtot1,
                  otherNxtot2, thisNxtot2, NPDFDim, fastNLO::v1d(NSubproc));
            }
         }
      }
   }

   if ( moption==fastNLO::kMerge ) {
      fastNLOTools::AddVectors(SigmaTilde, otherSigmaTilde);
   } else if ( moption==fastNLO::kAttach ) {
      for( int i=0 ; i<fNObsBins ; i++ ){
         int nxmax = GetNxmax(i);
         for( int k=0 ; k<GetTotalScalevars() ; k++ ){
            for( int l=0 ; l<GetTotalScalenodes() ; l++ ){
               for( int m=0 ; m<nxmax ; m++ ){
                  for( int n=0 ; n<other.GetNSubproc() ; n++ ){ // attach all other subprocesses
                     double s2  = otherSigmaTilde[i][k][l][m][n];
                     s2 *= this->Nevt/other.GetNevt();
                     this->SigmaTilde[i][k][l][m].push_back(s2);
                  }
               }
            }
         }
      }
   }
   else {
      for( int i=0 ; i<fNObsBins ; i++ ){
         int nxmax = GetNxmax(i);
         for( int k=0 ; k<GetTotalScalevars() ; k++ ){
            for( int l=0 ; l<GetTotalScalenodes() ; l++ ){
               for( int m=0 ; m<nxmax ; m++ ){
                  for( int n=0 ; n<NSubproc ; n++ ){
                     double w1  = this->GetMergeWeight(moption,n,i);
                     double w2  = other.GetMergeWeight(moption,n,i);
                     double& s1 = this->SigmaTilde[i][k][l][m][n];
                     double s2  = otherSigmaTilde[i][k][l][m][n];
                     if ( s1!=0 || s2!=0 ) {
                        if ( w1==0 || w2==0 ) {
                           error["fastNLOCoeffAddFix"]<<"Mergeing weight is 0, but sigma tilde is non-zero. Cannot proceed!"<<endl;
                           exit(3);
                        }
                        s1 = ( w1*s1/Nevt + w2*s2/other.GetNevt() ) / (w1 + w2 ) * ( Nevt + other.GetNevt() ) ;
                     }
                  }
               }
            }
         }
      }
   }
   //Nevt += othfix.Nevt;
   fastNLOCoeffAddBase::Add(other,moption);
   if ( moption==fastNLO::kAdd ) {
      NormalizeCoefficients(2);
      Nevt = 1;
      fWgt.WgtNevt = 1;
   }
   else if ( moption==fastNLO::kUnweighted ) {
      NormalizeCoefficients(1);
   }
   else if ( moption==fastNLO::kAttach ) {
      NormalizeCoefficients(1);
   }
}


//________________________________________________________________________________________________________________ //
int fastNLOCoeffAddFix::GetTotalScalevars() const {
   //! Get nuber of scale-variations
   int totalscalevars=1;
   for(int scaledim=0;scaledim<NScaleDim;scaledim++){
      totalscalevars *= Nscalevar[scaledim];
   }
   return totalscalevars;
}


//________________________________________________________________________________________________________________ //
int fastNLOCoeffAddFix::GetTotalScalenodes() const {
   //! Get number of scale nodes
   if ( !ScaleNode.empty() ) return ScaleNode[0][0][0].size();
   else return 0;
   //    int totalscalenodes=1;
   //    for(int scaledim=0;scaledim<NScaleDim;scaledim++){
   //       if ( !ScaleNode.empty() )
   //  totalscalenodes *= ScaleNode[0][0][0].size();
   //       //Nscalenode[scaledim];
   //    }
   //    return totalscalenodes;
}


//________________________________________________________________________________________________________________ //
bool  fastNLOCoeffAddFix::IsCompatible(const fastNLOCoeffAddFix& other) const {
   //! Check for compatibility for merging/adding of two contributions
   if ( ! ((fastNLOCoeffAddBase*)this)->IsCompatible(other)) return false;
   if ( GetNScaleNode() != other.GetNScaleNode() ) {
      say::warn["IsCompatible"]<<"Incompatible number of scale nodes found."<<endl;
      return false;
   }
   if ( GetNScalevar() != other.GetNScalevar() ) {
      say::warn["IsCompatible"]<<"Incompatible number of scale variations found."<<endl;
      return false;
   }
   if ( GetAvailableScaleFactors()[GetNScalevar()-1] != other.GetAvailableScaleFactors()[GetNScalevar()-1] ) {
      say::warn["IsCompatible"]<<"Incompatible scale variations found."<<endl;
      return false;
   }
   for ( int i=0 ; i<fNObsBins ; i++ ){
      for ( int is=0 ; is<GetNScaleNode() ; is++ ){
         if ( GetScaleNode(i,0,is) != other.GetScaleNode(i,0,is) ) {
            say::warn["IsCompatible"]<<"Incompatible scale node found."<<endl;
            return false;
         }
      }
   }
   return true;
}


//________________________________________________________________________________________________________________ //
bool  fastNLOCoeffAddFix::IsCatenable(const fastNLOCoeffAddFix& other) const {
   //! Check for compatibility of catenating observable bins
   if ( ! ((fastNLOCoeffAddBase*)this)->IsCatenable(other)) return false;
   if ( GetNScaleNode() != other.GetNScaleNode() ) {
      debug["IsCatenable"]<<"Incompatible number of scale nodes found. Skipped."<<endl;
      return false;
   }
   if ( GetNScalevar() != other.GetNScalevar() ) {
      debug["IsCatenable"]<<"Incompatible number of scale variations found. Skipped."<<endl;
      return false;
   }
   if ( GetAvailableScaleFactors()[GetNScalevar()-1] != other.GetAvailableScaleFactors()[GetNScalevar()-1] ) {
      debug["IsCatenable"]<<"Incompatible scale variations found. Skipped."<<endl;
      return false;
   }
   info["IsCatenable"]<<"Fix-scale contributions are catenable"<<endl;
   return true;
}


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffAddFix::IsEquivalent(const fastNLOCoeffBase& other, double rtol) const {
   const fastNLOCoeffAddFix* op = dynamic_cast<const fastNLOCoeffAddFix*>(&other);
   if (op == nullptr) {
      debug["IsEquivalent"] << "other is not of type fastNLOCoeffAddFix." << endl;
      return false;
   }

   if (NPDFDim != op->GetNPDFDim()) {
      debug["IsEquivalent"] << "NPDFDim is different: " << NPDFDim << " <-> " << op->GetNPDFDim() << endl;
      return false;
   }
   if (!fastNLOTools::SameTails(GetAllXNodes1(), op->GetAllXNodes1(), rtol)) {
      debug["IsEquivalent"] << "XNode1 not equivalent, see above." << endl;
      return false;
   }
   if (!fastNLOTools::SameTails(GetAllXNodes2(), op->GetAllXNodes2(), rtol)) {
      debug["IsEquivalent"] << "XNode2 not equivalent, see above." << endl;
      return false;
   }
   fastNLO::v5d ost5 = op->SigmaTilde;
   if (SigmaTilde.size() != ost5.size()) {
      debug["IsEquivalent"] << "Number of observable bins is different." << endl;
      return false;
   }
   for (unsigned int obsBin = 0; obsBin < SigmaTilde.size(); obsBin++) {
      unsigned int tOffsetXN1, oOffsetXN1, tOffsetXN2, oOffsetXN2;
      int tnb1 = GetNxtot1(obsBin);
      int onb1 = op->GetNxtot1(obsBin);
      int tnb2 = GetNxtot2(obsBin);
      int onb2 = op->GetNxtot2(obsBin);
      int x1Max = std::min(tnb1, onb1);
      int x2Max = std::min(tnb2, onb2);
      if (tnb1 > onb1) {
         tOffsetXN1 = tnb1 - onb1;
         oOffsetXN1 = 0;
      } else {
         tOffsetXN1 = 0;
         oOffsetXN1 = onb1 - tnb1;
      }
      if (NPDFDim > 1) {
         if (tnb2 > onb2) {
            tOffsetXN2 = tnb2 - onb2;
            oOffsetXN2 = 0;
         } else {
            tOffsetXN2 = 0;
            oOffsetXN2 = onb2 - tnb2;
         }
      } else {
         tOffsetXN2 = tOffsetXN1;
         oOffsetXN2 = oOffsetXN1;
      }

      fastNLO::v4d tst4 = SigmaTilde[obsBin];
      fastNLO::v4d ost4 = ost5[obsBin];
      if (tst4.size() != ost4.size()) {
         return false;
      }
      for (unsigned int scaleVar = 0; scaleVar < tst4.size(); scaleVar++) {
         fastNLO::v3d tst3 = tst4[scaleVar];
         fastNLO::v3d ost3 = ost4[scaleVar];
         if (tst3.size() != ost3.size()) {
            return false;
         }
         for (unsigned int scaleNode = 0; scaleNode < tst3.size(); scaleNode++) {
            fastNLO::v2d tst2 = tst3[scaleNode];
            fastNLO::v2d ost2 = ost3[scaleNode];
            int xIt2Max;
            switch (NPDFDim) {
               case 0:
                  xIt2Max = 1;
                  break;
               case 1:
                  xIt2Max = x1Max;
                  break;
               case 2:
                  xIt2Max = x2Max;
                  break;
               default:
                  error["IsEquivalent"] << "Unsupported NPDFDim: " << NPDFDim << endl;
                  exit(1);
            }
            for (int xIt2 = 0; xIt2 < xIt2Max; xIt2++) {
               int xIt1Max;
               switch (NPDFDim) {
                  case 0:
                  case 2:
                     xIt1Max = x1Max;
                     break;
                  case 1:
                     xIt1Max = xIt2 + 1;
                     break;
                  default:
                     error["IsEquivalent"] << "Unsupported NPDFDim: " << NPDFDim << endl;
                     exit(1);
               }
               for (int xIt1 = 0; xIt1 < xIt1Max; xIt1++) {
                  fastNLO::v1d tst1 = tst2[GetXIndex(obsBin, xIt1 + tOffsetXN1, xIt2 + tOffsetXN2)];
                  fastNLO::v1d ost1 = ost2[op->GetXIndex(obsBin, xIt1 + oOffsetXN1, xIt2 + oOffsetXN2)];
                  if (tst1.size() != ost1.size()) {
                     return false;
                  }
                  for (unsigned int subProcess = 0; subProcess < tst1.size(); subProcess++) {
                     double rdiff = std::abs((tst1[subProcess] - ost1[subProcess]) / tst1[subProcess]);
                     if (rdiff > rtol) {
                        debug["IsEquivalent"] << "SigmaTilde rdiff too high: obsBin=" << obsBin << " scaleVar=" << scaleVar
                           << " scaleNode=" << scaleNode << " x1=" << (xIt1 - x1Max) << " x2=" << (xIt2 - x2Max)
                           << " subProcess=" << subProcess << " t=" << tst1[subProcess] << " o=" << ost1[subProcess]
                           << " rdiff=" << rdiff << endl;
                        return false;
                     }
                  }
               }
            }
         }
      }
   }
   return true;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Clear() {
   //! Set all elelments of SigmaTilde to zero.
   fastNLOCoeffAddBase::Clear();
   fastNLOTools::ClearVector(SigmaTilde);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::NormalizeCoefficients(double wgt){
   //!< Set number of events to wgt (default=1) and normalize coefficients accordingly.
   //! This means, that the information about the
   //! number of events is essentially lost
   if ( wgt == Nevt ) return;
   MultiplyCoefficientsByConstant(wgt/Nevt);
   fastNLOCoeffAddBase::NormalizeCoefficients(wgt); //Nevt = wgt;
}

//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::NormalizeCoefficients(const std::vector<std::vector<double> >& wgtProcBin){
   //!< Set number of events to wgtProcBin for each subprocess and bin
   //!< and normalize coefficients accordingly.
   if ( int(wgtProcBin.size()) != GetNSubproc() ) {//NObs
      error["NormalizeCoefficients"]<<"Dimension of weights (iObs) incompatible with table (wgtProcBin must have dimension [iProc][iBin])."<<endl;
      exit(4);
   }

   for ( int iProc = 0 ; iProc<GetNSubproc(); iProc++ ) {
      if ( int(wgtProcBin[iProc].size()) != GetNObsBin() ) {
         error["NormalizeCoefficients"]<<"Dimension of weights (iProc) incompatible with table (wgtProcBin must have dimension [iProc][iBin])."<<endl;
         exit(4);
      }
      for ( int iObs = 0 ; iObs<GetNObsBin(); iObs++ ) {
         MultiplyBinProc(iObs, iProc, wgtProcBin[iProc][iObs]/Nevt);
      }
   }
   fastNLOCoeffAddBase::NormalizeCoefficients(wgtProcBin);
   //Nevt = 1;
}

//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::MultiplyCoefficientsByConstant(double fact) {
   for (unsigned int i=0; i<SigmaTilde.size(); i++) {
      MultiplyBin(i,fact);
   }
}

//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::MultiplyBin(unsigned int iObsIdx, double fact) {
   //! Multiply observable bin
   for (int m=0 ; m<GetNSubproc() ; m++)
      MultiplyBinProc(iObsIdx,m,fact);
}

//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::MultiplyBinProc(unsigned int iObsIdx, unsigned int iProc, double fact) {
   //! Multiply observable bin for a single subprocess
   debug["fastNLOCoeffAddFix::MultiplyBin"]<<"Multiplying table entries in CoeffAddFix for bin index " << iObsIdx << " by factor " << fact << endl;
   for (unsigned int s=0 ; s<SigmaTilde[iObsIdx].size() ; s++) {
      for (unsigned int x=0 ; x<SigmaTilde[iObsIdx][s].size() ; x++) {
         for (unsigned int l=0 ; l<SigmaTilde[iObsIdx][s][x].size() ; l++) {
            SigmaTilde[iObsIdx][s][x][l][iProc] *= fact;
         }
      }
   }
   fastNLOCoeffAddBase::MultiplyBinProc(iObsIdx, iProc, fact);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Print(int iprint) const {
   if ( !(iprint < 0) ) {
      fastNLOCoeffAddBase::Print(iprint);
      cout << fastNLO::_DSEP20C << " fastNLO Table: CoeffAddFix " << fastNLO::_DSEP20 << endl;
   } else {
      cout << endl << fastNLO::_CSEP20C << " fastNLO Table: CoeffAddFix " << fastNLO::_CSEP20 << endl;
   }
   for (int i=0; i<NScaleDim; i++) {
      printf(" # No. of scale variations (Nscalevar) %d\n",GetNScalevar());
      fastNLOTools::PrintVector(GetAvailableScaleFactors(),"Scale factors (ScaleFac[0][])","#");
      printf(" # No. of scale nodes (Nscalenode)     %d\n",GetNScaleNode());
   }
   if ( std::abs(iprint) > 0 ) {
      cout << fastNLO::_SSEP20C << " Extended information (iprint > 0) " << fastNLO::_SSEP20 << endl;
      char buffer[1024];
      for (int i=0; i<fNObsBins; i++) {
         // Print only for first and last observable bin
         if (i==0 || i==fNObsBins-1) {
            printf(" #   Observable bin no. %d\n",i+1);
            snprintf(buffer, sizeof(buffer), "Scale nodes (ScaleNode[%d][0][0][])",i);
            fastNLOTools::PrintVector(GetScaleNodes(i,0),buffer,"#  ");
         }
      }
   }
   if ( std::abs(iprint) > 1 ) {
      cout << fastNLO::_SSEP20C << " Extended information (iprint > 1) " << fastNLO::_SSEP20 << endl;
      printf(" #   Printing of SigmaTilde not yet implemented.\n");
   }
   cout << fastNLO::_CSEPSC << endl;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::EraseBin(unsigned int iObsIdx, int ITabVersionRead) {
   //! Erase observable bin
   debug["fastNLOCoeffAddFix::EraseBin"]<<"Erasing table entries in CoeffAddFix for bin index " << iObsIdx << endl;
   if ( ScaleNode.size() == 0 ) {
      say::error["EraseBin"]<<"All fix-scale bins deleted already. Aborted!" << endl;
      exit(1);
   }
   if ( ScaleNode.size() != 0 ) ScaleNode.erase(ScaleNode.begin()+iObsIdx);
   if ( SigmaTilde.size() != 0 ) SigmaTilde.erase(SigmaTilde.begin()+iObsIdx);
   fastNLOCoeffAddBase::EraseBin(iObsIdx,ITabVersionRead);
}

// Catenate observable bin
void fastNLOCoeffAddFix::CatBin(const fastNLOCoeffAddFix& other, unsigned int iObsIdx, int ITabVersionRead) {
   debug["fastNLOCoeffAddFix::CatBin"]<<"Catenating observable bin in CoeffAddFix corresponding to bin index " << iObsIdx << endl;
   if ( ScaleNode.size() == 0 ) {
      say::error["CatBin"]<<"Initial fix-scale table is empty. Aborted!" << endl;
      exit(1);
   }
   unsigned int nold = ScaleNode.size();
   if ( ScaleNode.size() != 0 ) {
      ScaleNode.resize(nold+1);
      ScaleNode[nold] = other.ScaleNode[iObsIdx];
   }
   if ( SigmaTilde.size() != 0 ) {
      SigmaTilde.resize(nold+1);
      SigmaTilde[nold] = other.SigmaTilde[iObsIdx];
   }
   fastNLOCoeffAddBase::CatBin(other, iObsIdx, ITabVersionRead);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ExtendSigmaTildeX(int ObsBin, unsigned int OldXSize1, unsigned int OldXSize2) {
   int nScalevar = GetNScalevar();
   int nScaleNode = GetNScaleNode();
   int newXSize1 = GetNxtot1(ObsBin);
   int newXSize2 = GetNxtot2(ObsBin);
   for (int scaleVar = 0; scaleVar < nScalevar; scaleVar++) {
      for (int scaleNode = 0; scaleNode < nScaleNode; scaleNode++) {
         fastNLOTools::ExtendSigmaTildeX(
            SigmaTilde[ObsBin][scaleVar][scaleNode], OldXSize1, newXSize1, OldXSize2, newXSize2,
            NPDFDim, fastNLO::v1d(NSubproc));
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Fill(fnloEvent& Event, int ObsBin, int X, int scalevar, const vector<pair<int, double>>& nmu1,
      const vector<pair<int, double>>& nmu2, int SubProcess, double w) {
   for (unsigned int m1 = 0 ; m1<nmu1.size() ; m1++) {
      SigmaTilde[ObsBin][scalevar][nmu1[m1].first][X][SubProcess] += Event._w * nmu1[m1].second * w;
   }
}
