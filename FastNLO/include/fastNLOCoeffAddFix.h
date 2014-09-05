#ifndef __fastNLOCoeffAddFix__
#define __fastNLOCoeffAddFix__

#include "fastNLOCoeffAddBase.h"
#include "fastNLOConstants.h"

using namespace std;

class fastNLOCoeffAddFix : public fastNLOCoeffAddBase {

   friend class fastNLOReader;
   friend class fastNLOCreate;

public:
   fastNLOCoeffAddFix(int NObsBin);
   fastNLOCoeffAddFix(const fastNLOCoeffBase&);
   virtual ~fastNLOCoeffAddFix(){;}
   virtual fastNLOCoeffBase* Clone() const;                                     //!< returns 'new' copy of this instance.
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false);
   virtual void Read(istream&table);
   void ReadRest(istream& table);
   virtual void Write(ostream& table);
   virtual void Add(const fastNLOCoeffAddBase& other);
   virtual void Print() const;
   virtual void Clear();                                                        //!< Clear all coefficients and event counters
   virtual void NormalizeCoefficients();                                        //!< Set number of events to 1 and normalize coefficients accordingly.
   virtual void MultiplyCoefficientsByConstant(double coef);                    //!< Multiply all coefficients by constant coef

   int GetTotalScalevars() const ;
   int GetTotalScalenodes() const ;
   int GetNScaleNode() const { return GetTotalScalenodes(); }
   int GetNScalevar() const { return Nscalevar[0];}
   v1d GetAvailableScaleFactors() const { return ScaleFac[0]; }
   double GetScaleFactor(int iVar) const {
      if ( iVar >= (int)ScaleFac[0].size() )
         this->error["GetScaleFactor"]<<"Scalevariation no. "<<iVar<<" not available. There are only "<<GetNScalevar()<<" available in this table."<<endl;
      return ScaleFac[0][iVar];
   }

   double GetSigmaTilde(int iObs, int iSvar, int ix, int is, int iN ) const { return SigmaTilde[iObs][iSvar][ix][is][iN];}
   double GetScaleNode(int iObs, int iSvar, int iNode ) const { return ScaleNode[iObs][0][iSvar][iNode];}

   void ResizePdfLC();
   void ResizePdfSplLC();
   void ResizeSigmaTilde();


protected:
   fastNLOCoeffAddFix();
   void ReadCoeffAddFix(istream& table);

   vector < int > Nscalevar;
   //vector < int > Nscalenode;
   v2d ScaleFac;
   v4d ScaleNode;
   v5d SigmaTilde; // units are (p)barn * Nevt / BinSize

public:
   v2d AlphasTwoPi_v20;
   v4d PdfLc;
   v4d PdfSplLc1;
   v4d PdfSplLc2;
};

#endif
