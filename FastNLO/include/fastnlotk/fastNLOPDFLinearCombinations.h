#ifndef __fastNLOLinearCombinations__
#define __fastNLOLinearCombinations__

#include <string>
#include "speaker.h"
#include "fastNLOCoeffAddBase.h"

using namespace std;

class fastNLOPDFLinearCombinations {

public:
   fastNLOPDFLinearCombinations();
   ~fastNLOPDFLinearCombinations();

   vector<double > CalcPDFLinearCombination(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 = vector<double>(), const vector<double>& pdfx2 = vector<double>() , bool pdf2IsAntiParticle = false) const;

protected:
   vector<double > MakeAntiHadron(const vector<double >& hadron) const;

private: 
   vector<double > CalcPDFLCTwoHadrons(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1, const vector<double>& pdfx2 ) const ;
   vector<double > CalcPDFLCOneHadron(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 ) const;

   vector<double> CalcPDFDIS(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1) const;
   vector<double> CalcPDFHHCFromTable(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 , const vector<double>& pdfx2) const ; // hh collisions. PDFLiCos are stored in table
   vector<double> CalcPDFHHC(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 , const vector<double>& pdfx2) const ; // jets in hh
   vector<double> CalcDefaultPDFLiCos(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 , const vector<double>& pdfx2) const ; // jets in hh
   vector<double> CalcPDFttbar(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 , const vector<double>& pdfx2) const ; // ttbar
   vector<double> CalcPDFThreshold(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 , const vector<double>& pdfx2) const ; // pp->2jets 
};
#endif
