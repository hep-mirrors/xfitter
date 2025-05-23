#ifndef __fastNLOCoeffAddFix__
#define __fastNLOCoeffAddFix__

#include "fastNLOCoeffAddBase.h"
#include "fastNLOConstants.h"
#include "fastNLOEvent.h"


class fastNLOCoeffAddFix : public fastNLOCoeffAddBase {

   friend class fastNLOTable;
   friend class fastNLOReader;
   friend class fastNLOCreate;

public:
   fastNLOCoeffAddFix() = delete;
   fastNLOCoeffAddFix(int NObsBin);
   explicit fastNLOCoeffAddFix(const fastNLOCoeffBase&);
   virtual ~fastNLOCoeffAddFix(){;}
   virtual fastNLOCoeffAddFix* Clone() const;                                     //!< returns 'new' copy of this instance.
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false);
   virtual void Read(std::istream&table, int ITabVersionRead);
   void ReadRest(std::istream& table, int ITabVersionRead);
   virtual void Write(std::ostream& table, int ITabVersionWrite);
   virtual void Add(const fastNLOCoeffAddBase& other, fastNLO::EMerge moption = fastNLO::kMerge);
   virtual void Print(int iprint) const;

   // Manipulate coefficient bins
   virtual void Clear();//!< Clear all coefficients and event counters
   virtual void NormalizeCoefficients(double wgt=1);//!<a Set number of events to wgt and re-normalize coefficients accordingly
   virtual void NormalizeCoefficients(const std::vector<std::vector<double> >& wgtProcBin);
   virtual void MultiplyCoefficientsByConstant(double fact);//!< Multiply all coefficients of all bins by a constant factor
   virtual void MultiplyBin(unsigned int iObsIdx, double fact); //!< Multiply coefficients of one bin a factor
   virtual void MultiplyBinProc(unsigned int iObsIdx, unsigned int iProc, double fact); //!< Multiply coefficients of one bin and subprocess by a factor
   // Erase observable bin from table
   virtual void EraseBin(unsigned int iObsIdx, int ITabVersionRead);
   // Catenate observable to table
   virtual void CatBin(const fastNLOCoeffAddFix& other, unsigned int iObsIdx, int ITabVersionRead);

   int GetTotalScalevars() const ;
   int GetTotalScalenodes() const ;
   int GetNScaleNode() const { return GetTotalScalenodes(); }
   int GetNScalevar() const { return Nscalevar[0];}
   fastNLO::v1d GetAvailableScaleFactors() const { return ScaleFac[0]; }
   double GetScaleFactor(int iVar) const {
      if ( iVar >= (int)ScaleFac[0].size() )
         this->error["GetScaleFactor"]<<"Scalevariation no. "<<iVar<<" not available. There are only "<<GetNScalevar()<<" available in this table."<< std::endl;
      return ScaleFac[0][iVar];
   }

   double GetSigmaTilde(int iObs, int iSvar, int ix, int is, int iN ) const { return SigmaTilde[iObs][iSvar][ix][is][iN];}
   double GetScaleNode(int iObs, int iSvar, int iNode ) const { return ScaleNode[iObs][0][iSvar][iNode]; }
   std::vector < double > GetScaleNodes(int iObs, int iSvar) const { return ScaleNode[iObs][0][iSvar]; }

   void ResizePdfLC();
   void ResizePdfSplLC();
   void ResizeSigmaTilde();
   bool IsCompatible(const fastNLOCoeffAddFix& other) const;                   //!< Check for compatibility of two contributions for merging/adding
   bool IsCatenable(const fastNLOCoeffAddFix& other) const;        //!< Check for compatibility of two contributions for merging/adding
   bool IsEquivalent(const fastNLOCoeffBase& other, double rtol) const;

   void ExtendSigmaTildeX(int ObsBin, unsigned int OldXSize1, unsigned int OldXSize2);
   void Fill(fnloEvent& Event, int ObsBin, int X, int scalevar, const std::vector<std::pair<int, double>>& nmu1,
      const std::vector<std::pair<int, double>>& nmu2, int SubProcess, double w);

protected:
   void ReadCoeffAddFix(std::istream& table, int ITabVersionRead);

   std::vector < int > Nscalevar;
   //std::vector < int > Nscalenode;
   fastNLO::v2d ScaleFac;
   fastNLO::v4d ScaleNode;
   fastNLO::v5d SigmaTilde; // units are (p)barn * Nevt / BinSize

public:
   fastNLO::v2d AlphasTwoPi_v20;
   fastNLO::v4d PdfLc;
   fastNLO::v4d PdfSplLc1;
   fastNLO::v4d PdfSplLc2;
};

#endif
