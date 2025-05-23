#ifndef __fastNLOCoeffAddBase__
#define __fastNLOCoeffAddBase__


#include "fastNLOCoeffBase.h"
#include "fastNLOConstants.h"
#include "fastNLOInterpolBase.h"
#include "fastNLOEvent.h"

namespace fastNLO {
   struct WgtStat {
      double WgtNevt = 0; //!< 'number of events', i.e. normalisation as suggested by generator (identical to previously use 'Nevt')
      int NumTable = 1; //!< Number of tables merged into this table
      unsigned long long WgtNumEv = 0; //!< number of entries
      double WgtSumW2 = 0; //!< Sum of all weight**2
      double SigSumW2 = 0; //!< Sum of all sigma**2 (i.e. (wgt*alpha*pdf)**2 )
      double SigSum = 0; //!< Sum of all sigma (i.e. (wgt*alpha*pdf)**2 )
      fastNLO::v2d WgtObsSumW2; //!< sumw2[proc][obs]
      fastNLO::v2d SigObsSumW2; //!< sumw2[proc][obs]
      fastNLO::v2d SigObsSum;   //!< sum[proc][obs]
      std::vector < std::vector < unsigned long long > > WgtObsNumEv; //!< Nentries[proc][obs]
      //fastNLO::v2d SigObsUser;   //!< sum[proc][obs]

      void Erase() {
         WgtNevt=0;
         NumTable=1;
         WgtNumEv=0;
         WgtSumW2=0;
         SigSumW2=0;
         SigSum=0;
         for ( auto& i : WgtObsSumW2 ) for ( auto& j : i ) j=0;
         for ( auto& i : SigObsSumW2 ) for ( auto& j : i ) j=0;
         for ( auto& i : SigObsSum   ) for ( auto& j : i ) j=0;
         for ( auto& i : WgtObsNumEv ) for ( auto& j : i ) j=0;
         //for ( auto& i : SigObsUser  ) for ( auto& j : i ) j=0;
      };

      void Add(const WgtStat& other ) {
         this->WgtNevt  += other.WgtNevt;
         this->NumTable += other.NumTable;
         this->WgtNumEv += other.WgtNumEv;
         this->WgtSumW2 += other.WgtSumW2;
         this->SigSumW2 += other.SigSumW2;
         this->SigSum   += other.SigSum;
         if ( this->WgtObsNumEv.size() != other.WgtObsNumEv.size() ) exit(8);
         for ( unsigned int i = 0 ; i<WgtObsNumEv.size() ; i++ ) {
            if ( this->WgtObsNumEv[i].size() != other.WgtObsNumEv[i].size() ) exit(8);
            for ( unsigned int j = 0 ; j<WgtObsNumEv[i].size() ; j++ ) {
               this->WgtObsSumW2[i][j] += other.WgtObsSumW2[i][j];
               this->SigObsSumW2[i][j] += other.SigObsSumW2[i][j];
               this->SigObsSum[i][j]   += other.SigObsSum[i][j];
               this->WgtObsNumEv[i][j] += other.WgtObsNumEv[i][j];
               //this->SigObsUser[i][j]  += other.SigObsUser[i][j];
            }
         }
      };
      WgtStat& operator+=(const WgtStat& other) { this->Add(other); return *this;}

      void SetWgtUser( const fastNLO::v2d& wgtUser ) {
         //! Set user weights
         //! wgtUser[proc][obs]
         if ( wgtUser.size()!=WgtObsNumEv.size() ) {
            std::cerr<<"Error [fastNLO::WgtStat::SetWgtUser()] Array with wrong size (wrong number of subprocesses. Must be "
                     <<WgtObsNumEv.size()<<", but is "<< wgtUser.size()<<std::endl;
            exit(8);
         }
         if ( wgtUser[0].size()!=WgtObsNumEv[0].size() ) {
            std::cerr<<"Error [fastNLO::WgtStat::SetWgtUser()] Array with wrong size (wrong number of obs bins. Must be "
                     <<WgtObsNumEv[0].size()<<", but is "<< wgtUser[0].size()<<std::endl;
            exit(8);
         }
         SigObsSum = wgtUser;
         //SigObsUser = wgtUser;
      };

      void SetWgtUser( double wgtUser ) {
         //! Set user weight valid for all entries
         SigSum = wgtUser;
         for ( auto& proc : SigObsSum ) {
            for ( auto& w : proc ) w=wgtUser;
         }
      };

      void SetWgtUser( const std::vector<double>& wgtUserObs ) {
         //! Set user weights
         //! wgtUserObs[obs]
         if ( wgtUserObs.size()!=WgtObsNumEv[0].size() ) {
            std::cerr<<"Error [fastNLO::WgtStat::SetWgtUser()] Array with wrong size (wrong number of obs bins. Must be "
                     <<WgtObsNumEv[0].size()<<", but is "<< wgtUserObs.size()<<std::endl;
            exit(8);
         }
         for ( std::vector<double>& proc : SigObsSum ) {
            proc = wgtUserObs;
         }
      };
   };

}

class fastNLOCoeffAddBase : public fastNLOCoeffBase {

   friend class fastNLOTable;
   friend class fastNLOCreate;

public:
   fastNLOCoeffAddBase() = delete;
   fastNLOCoeffAddBase(int NObsBin);
   explicit fastNLOCoeffAddBase(const fastNLOCoeffBase& base);
   virtual ~fastNLOCoeffAddBase() {}
   virtual fastNLOCoeffAddBase* Clone() const;                                     //!< returns 'new' copy of this instance.
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false);
   void Read(std::istream& table, int ITabVersionRead);
   virtual void Write(std::ostream& table, int ITabVersionWrite);
   virtual void Add(const fastNLOCoeffAddBase& other, fastNLO::EMerge moption=fastNLO::kMerge);
   virtual void Print(int iprint) const;

   // Manipulate coefficient bins
   virtual void Clear();//!< Clear all coefficients and event counters
   virtual void NormalizeCoefficients(double wgt=1); //!< Set number of events to unity and normalize coefficients accordingly
   virtual void NormalizeCoefficients(const std::vector<std::vector<double> >& wgtProcBin) {};
   virtual void MultiplyCoefficientsByConstant(double fact) {};//!< Multiply all coefficients of all bins by a constant factor
   virtual void MultiplyBin(unsigned int iObsIdx, double fact) {};  //!< Multiply coefficients of one observable bin a factor
   virtual void MultiplyBinProc(unsigned int iObsIdx, unsigned int iProc, double fact) {}; //!< Multiply coefficients of one observable bin a factor (idx starting from 0)
   virtual void EraseBin(unsigned int iObsIdx, int ITabVersionRead);//!< Erase observable bin from table
   virtual void CatBin(const fastNLOCoeffAddBase& other, unsigned int iObsIdx, int ITabVersionRead); //!< Catenate observable to table

   int GetIRef() const {return IRef;}
   void SetIRef(int iref=1) {IRef=iref;}
   double GetNevt() const { return Nevt; }
   double GetNevt(int NObsBin, int NSubproc) const {
      //return fWgt.WgtObsSumW2[NSubproc][NObsBin];
      if (Nevt > 0) return Nevt;
      else {std::cout<<"Todo. Preparation for v2.3."<< std::endl; return Nevt;}
   }
   void SetNevt(double nevt) { Nevt = nevt;}                                    //!< Set number of events
   int GetNxmax(int Obsbin) const ;
   int GetXIndex(int Obsbin,int x1bin,int x2bin =0) const ;
   int GetNSubproc() const { return NSubproc;}
   int GetIScaleDep() const { return IScaleDep;}
   int GetNPDF() const {return NPDFPDG.size();}
   int GetPDFPDG(int iPDF) const {return NPDFPDG[iPDF];}
   int GetNPDFDim() const {return NPDFDim;}
   int GetIPDFdef1() const { return IPDFdef1; }
   int GetIPDFdef2() const { return IPDFdef2; }
   int GetIPDFdef3() const { return IPDFdef3; }
   int GetNpow() const {return Npow;}
   int GetNScales() const {return NScales;}
   int GetNScaleDim() const {return NScaleDim;}
   //std::vector<std::string > GetScaleDescript(int iScale=0) const { return ScaleDescript[iScale]; };
   std::string GetScaleDescription(int iScale=0) const { return ScaleDescript[0][iScale]; };         // getter for scale description of scale iScale
   const std::vector<std::vector<std::string > >& GetScaleDescr() const { return ScaleDescript; }
   int GetNxtot1(int iBin) const { return XNode1[iBin].size(); }
   int GetNxtot2(int iBin) const { return XNode2.size() > 0 ? XNode2[iBin].size() : -1; }

   double GetXNode1(int iObsBin, int iNode) const { return XNode1[iObsBin][iNode]; }
   double GetXNode2(int iObsBin, int iNode) const { return XNode2[iObsBin][iNode]; }
   double GetX1(int iObsBin, int iXnode) const; //! return x value of pdf1 for x-node 1
   double GetX2(int iObsBin, int iXnode) const; //! return x value of pdf1 for x-node 1

   const std::vector < double >& GetXNodes1(int iObsBin) const { return XNode1[iObsBin]; }
   const std::vector < double >& GetXNodes2(int iObsBin) const { return XNode2[iObsBin]; }

   fastNLO::v2d GetAllXNodes1() const { return XNode1; }
   fastNLO::v2d GetAllXNodes2() const { return XNode2; }
   fastNLO::v2d* AccessAllXNodes1() { return &XNode1; }
   fastNLO::v2d* AccessAllXNodes2() { return &XNode2; }
   bool IsReference() const {return IRef>0;};
   bool IsCompatible(const fastNLOCoeffAddBase& other) const;
   bool IsCatenable(const fastNLOCoeffAddBase& other) const;

   const std::vector<std::vector<std::pair<int,int> > >& GetPDFCoeff() const { return fPDFCoeff;}

   bool SubIsEnabled(int index) const { return index < NSubproc && sub_enabled[index]; }  //!< Returns true if a single subcontribution is enabled in this contribution
   void SubEnable(int index, bool on = true) { if ( index < NSubproc ) sub_enabled[index] = on; } //!< Enables/Disables a single subcontribution
   void SubEnableAll( bool on = true ) { sub_enabled.assign(sub_enabled.size(), on); } //!< Enables/Disables all subcontributions
   bool SubSelect( std::vector< std::pair<int,int> > processes, bool on = true); //!< checks if the given subprocess list is compatible with the subcontribution list and enables/disables the corresponding subcontributions. Returns true on sucess and false on failure.

   const fastNLO::WgtStat& GetWgtStat() const { return fWgt;} //!< Get weight and event counts
   fastNLO::WgtStat& AccessWgtStat() { return fWgt;} //!< Get weight and event counts
   double GetMergeWeight(fastNLO::EMerge moption, int proc, int bin) const ; //!< Get merge weight for a given bin and subprocess

   void ExtendX(int ObsBin, std::vector<fastNLOInterpolBase*>& KernX1, std::vector<fastNLOInterpolBase*>& KernX2);
   virtual void ExtendSigmaTildeX(int ObsBin, unsigned int OldXSize1, unsigned int OldXSize2);
   virtual void Fill(fnloEvent& Event, int ObsBin, int X, int scalevar, const std::vector<std::pair<int, double>>& nmu1,
      const std::vector<std::pair<int, double>>& nmu2, int SubProcess, double w);

protected:
   void ReadCoeffAddBase(std::istream& table, int ITabVersionRead);
   int GetScaledimfromvar(int scalevar) const;

   int IRef = 0;
   int IScaleDep = 0;
   double Nevt = 0;
   int Npow = 0;
   std::vector < int > NPDFPDG;
   int NPDFDim = 0;
   std::vector < int > NFFPDG;
   int NFFDim = 0;
   int NSubproc = 0;
   int IPDFdef1 = 0;
   int IPDFdef2 = 0;
   int IPDFdef3 = 0;
   std::vector<std::vector<std::pair<int,int> > > fPDFCoeff;                                                   //! fPDFCoeff[iSubProc][iPartonPair][pair]
   // Missing: linear PDF combinations for IPDFdef1=0
   std::vector < double > Hxlim1;
   fastNLO::v2d XNode1;
   std::vector < double > Hxlim2;
   fastNLO::v2d XNode2;
   std::vector < int > Nztot;
   std::vector < double > Hzlim;
   fastNLO::v2d ZNode;
   int NScales = 0;
   int NScaleDim = 0;
   std::vector < int > Iscale;                                                                       // not used
   std::vector < std::vector < std::string > > ScaleDescript;

   std::vector < bool > sub_enabled;

   fastNLO::WgtStat fWgt; //!< event and weight counts
   // double fWgtNevt = 0; //!< 'number of events', i.e. normalisation as suggested by generator (identical to previously use 'Nevt')
   // unsigned long long fWgtNumEv = 0; //!< number of entries
   // double fWgtSumW2 = 0; //!< Sum of all weight**2
   // double fSigSumW2 = 0; //!< Sum of all sigma**2 (i.e. (wgt*alpha*pdf)**2 )
   // double fSigSum = 0; //!< Sum of all sigma (i.e. (wgt*alpha*pdf)**2 )
   // fastNLO::v2d fWgtObsSumW2; //!< sumw2[proc][obs]
   // fastNLO::v2d fSigObsSumW2; //!< sumw2[proc][obs]
   // fastNLO::v2d fSigObsSum;   //!< sum[proc][obs]
   // std::vector < std::vector < unsigned long long > > fWgtObsNumEv; //!< Nentries[proc][obs]


};

#endif
