#ifndef __fastNLOCoeffAddFlex__
#define __fastNLOCoeffAddFlex__

#include "fastNLOCoeffAddBase.h"
#include "fastNLOConstants.h"


class fastNLOCoeffAddFlex : public fastNLOCoeffAddBase {

   friend class fastNLOTable;
   friend class fastNLOReader;
   friend class fastNLOCreate;

public:
   fastNLOCoeffAddFlex() = delete;
   fastNLOCoeffAddFlex(int NObsBin, int iLOord);
   explicit fastNLOCoeffAddFlex(const fastNLOCoeffBase& base , int iLOord);
   virtual ~fastNLOCoeffAddFlex(){;}
   virtual fastNLOCoeffAddFlex* Clone() const;                                     //!< returns 'new' copy of this instance.
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false) ;
   virtual void Read(std::istream& table, int ITabVersionRead);
   void ReadRest(std::istream& table, int ITabVersionRead);
   virtual void Write(std::ostream& table, int ITabVersionWrite);
   virtual void Print(int iprint) const;
   virtual void Add(const fastNLOCoeffAddBase& other, fastNLO::EMerge moption = fastNLO::kMerge);

   // Manipulate coefficient bins
   virtual void Clear(); //!< Clear all coefficients and event counters
   virtual void NormalizeCoefficients(double wgt=1); //!< Set number of events to given wgt and re-normalize coefficients accordingly
   virtual void NormalizeCoefficients(const std::vector<std::vector<double> >& wgtProcBin);
   virtual void MultiplyCoefficientsByConstant(double fact); //!< Multiply all coefficients of all bins by a constant factor
   virtual void MultiplyBin(unsigned int iObsIdx, double fact); //!< Multiply coefficients of one bin a factor (iObsIdx starting with index 0)
   virtual void MultiplyBinProc(unsigned int iObsIdx, unsigned int iProc, double fact); //!< Multiply coefficients of one bin and subprocess a factor
   virtual void EraseBin(unsigned int iObsIdx, int ITabVersionRead); //!< Erase observable bin from table
   virtual void CatBin(const fastNLOCoeffAddFlex& other, unsigned int iObsIdx, int ITabVersionRead); //!< Catenate observable to table

   unsigned int GetNScaleNode1(int iObsBin) const { return ScaleNode1[iObsBin].size(); };
   unsigned int GetNScaleNode2(int iObsBin) const { return ScaleNode2[iObsBin].size(); };
   double GetScaleNode1(int iObsBin, int iNode) const { return ScaleNode1[iObsBin][iNode]; };
   double GetScaleNode2(int iObsBin, int iNode) const { return ScaleNode2[iObsBin][iNode]; };
   std::vector < double > GetScaleNodes1(int iObsBin) const { return ScaleNode1[iObsBin]; };
   std::vector < double > GetScaleNodes2(int iObsBin) const { return ScaleNode2[iObsBin]; };
   bool IsCompatible(const fastNLOCoeffAddFlex& other) const;                   //!< check for compatibilty for adding/merging of two tables
   bool IsCatenable(const fastNLOCoeffAddFlex& other) const;        //!< Check for compatibility of two contributions for merging/adding
   std::vector<fastNLO::v5d*>  AccessSigmaTildes() {
      return {&SigmaTildeMuIndep,&SigmaTildeMuRDep,&SigmaTildeMuFDep,&SigmaTildeMuRRDep,&SigmaTildeMuFFDep,&SigmaTildeMuRFDep};
   };//!< Get access to sigma tilde
   std::vector<const fastNLO::v5d*> GetSigmaTildes() const {
      return {&SigmaTildeMuIndep,&SigmaTildeMuRDep,&SigmaTildeMuFDep,&SigmaTildeMuRRDep,&SigmaTildeMuFFDep,&SigmaTildeMuRFDep};
   };//!< Get access to sigma tilde

protected:

   void ReadCoeffAddFlex(std::istream& table, int ITabVersionRead);

   int fILOord;   // obtained from Scenario
   int fSTildeDISFormat = 1; // format of sigma-tilde coefficients (0: log(mu2/q2), 1: log(mu2))

   // SigmaTilde [NObsBins] ['n' x-nodes] [n s1-Nodes] [n s2-Nodes] [nsubproc]
   fastNLO::v5d SigmaTildeMuIndep; // units are (p)barn * Nevt / BinSize
   fastNLO::v5d SigmaTildeMuFDep;
   fastNLO::v5d SigmaTildeMuRDep;
   fastNLO::v5d SigmaTildeMuRRDep;
   fastNLO::v5d SigmaTildeMuFFDep;
   fastNLO::v5d SigmaTildeMuRFDep;
   // SigmaRef [NObsBins] [nsubproc]
   fastNLO::v2d SigmaRefMixed;  // units are (p)barn * Nevt / BinSize
   fastNLO::v2d SigmaRef_s1;
   fastNLO::v2d SigmaRef_s2;
   //int NscalenodeScale1;
   //int NscalenodeScale2;
   // ScaleNodeXY [ObsBin] [NscalenodeScaleX]
   fastNLO::v2d ScaleNode1;
   fastNLO::v2d ScaleNode2;

public:
   fastNLO::v3d AlphasTwoPi;
   fastNLO::v5d PdfLcMuVar;
   fastNLO::v5d PdfXfx;

};

#endif
