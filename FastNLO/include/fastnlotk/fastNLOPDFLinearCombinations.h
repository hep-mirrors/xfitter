#ifndef __fastNLOLinearCombinations__
#define __fastNLOLinearCombinations__

#include "speaker.h"
#include "fastNLOCoeffAddBase.h"


class fastNLOPDFLinearCombinations {

public:
   fastNLOPDFLinearCombinations();
   ~fastNLOPDFLinearCombinations();

   std::vector<double > CalcPDFLinearCombination(const fastNLOCoeffAddBase* c, const std::vector<double>& pdfx1 = std::vector<double>(), const std::vector<double>& pdfx2 = std::vector<double>() , bool pdf2IsAntiParticle = false) const;

protected:
   std::vector<double > MakeAntiHadron(const std::vector<double >& hadron) const;

private: 
   std::vector<double > CalcPDFLCTwoHadrons(const fastNLOCoeffAddBase* c, const std::vector<double>& pdfx1, const std::vector<double>& pdfx2 ) const ;
   std::vector<double > CalcPDFLCOneHadron(const fastNLOCoeffAddBase* c, const std::vector<double>& pdfx1 ) const;

   std::vector<double> CalcPDFDIS(const fastNLOCoeffAddBase* c, const std::vector<double>& pdfx1) const;
   std::vector<double> CalcPDFDISFromTable(const fastNLOCoeffAddBase* c, const std::vector<double>& pdfx1) const ; // DIS. PDFLiCos are stored in table
   std::vector<double> CalcPDFHHCFromTable(const fastNLOCoeffAddBase* c, const std::vector<double>& pdfx1 , const std::vector<double>& pdfx2) const ; // hh collisions. PDFLiCos are stored in table
   std::vector<double> CalcPDFHHC(const fastNLOCoeffAddBase* c, const std::vector<double>& pdfx1 , const std::vector<double>& pdfx2) const ; // jets in hh
   std::vector<double> CalcDefaultPDFLiCos(const fastNLOCoeffAddBase* c, const std::vector<double>& pdfx1 , const std::vector<double>& pdfx2) const ; // jets in hh
   std::vector<double> CalcPDFttbar(const fastNLOCoeffAddBase* c, const std::vector<double>& pdfx1 , const std::vector<double>& pdfx2) const ; // ttbar
   std::vector<double> CalcPDFThreshold(const fastNLOCoeffAddBase* c, const std::vector<double>& pdfx1 , const std::vector<double>& pdfx2) const ; // pp->2jets 
};
#endif
