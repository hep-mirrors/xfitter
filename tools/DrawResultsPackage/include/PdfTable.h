#ifndef PdfTable_h
#define PdfTable_h

#include <TString.h>
#include <vector>
#include <string>
#include <TGraphAsymmErrors.h>

using std::vector;
using std::string;

//
// @brief Interface class 
// 


//
// @brief Reader of the PDF file
//
class PdfTable  {
 public:
 PdfTable():fQ2value(0),fXmin(0),fXmax(0),fNxValues(0),fNPdfs(0),fTable(NULL)
{}

// @brief Constructor, reads the file
  PdfTable(TString fName);
  virtual ~PdfTable();

  PdfTable* CreatePdfTable(const Char_t* filename);

  // @brief Get number of X points
  virtual const int GetNx(){return fNxValues;}
  // @brief Get number of PDfs
  const int GetNPdfs(){return fNPdfs;}
  // @brief Get Q2 value
  virtual const double GetQ2(){return fQ2value;}

  // @brief Get Xmin value
  const double GetXmin(){return fXmin;}
  // @brief Get Xmax value
  const double GetXmax(){return fXmax;}
  
  // @brief Get PDF iPdf  at a given ix grid point 
  const double GetPDF(int iX, int iPdf);

  // @brief Get PDF iPdf  at a given ix grid point 
  const double GetPDF(int iX, string name);

  // @brief Get PDF as a vector 
  const vector<double> GetPDF(string name);

  // @brief Get TGraph pdf for PDF with the name vs x
  virtual TGraphAsymmErrors* GetPDFGraph(string name);

  // @brief get index of the array
  const int GetIndex(string name);

  const double *GetTable(){return fTable;}

  const string GetColumnName(int i){return fColumnNames[i];}

 protected:
  double fQ2value;
  double fXmin;
  double fXmax;
  int    fNxValues;
  int    fNPdfs;
  double *fTable;
  vector<string> fColumnNames;
};


//
// @brief Helper class to read PDF error sets
//

class PdfErrorTables : public PdfTable{
 public:
 PdfErrorTables():PdfTable(),fErrorTables(NULL)
    {}
  // @brief Constructor to read the files from directory base for given Q2 set
  PdfErrorTables(string base, int iQ2=0, Bool_t SymErrors = kFALSE, TString option = TString("b"));

  // @brief Return TGraph 
  virtual TGraphAsymmErrors* GetPDFGraph(string name);

 private:
  // @brief List of PDF tables to store error sets
  vector<PdfTable*> fErrorTables;
  Bool_t fSymmetric;
  void GetPDFError(int ix, int iPDF, double* eminus, double* eplus);
  
  TString sErrOpt;
};


#endif
