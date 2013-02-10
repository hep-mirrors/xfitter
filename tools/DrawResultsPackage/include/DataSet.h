#ifndef  DataSet_h
#define  DataSet_h

#include <TString.h>
#include <TObject.h>
#include <iostream>
#include <TGraphErrors.h>
#include <TObjArray.h>
#include <vector>
#include <TObjString.h>
#include <TOrdCollection.h>
#include <map>
#include <TH1F.h>

using std::cout;
using std::cerr;
using std::endl;
using std::map;

struct SubPlot {
  TString* fTitle;
  TString* fPlotOption;
  TGraphErrors* fDUnc;
  TGraphErrors* fDTot;
  TGraphErrors* fTheo;
  TGraphErrors* fTMod;
  TH1F* fHistogram;
  bool fXlog;
  bool fYlog;
  SubPlot(const char* descriptor);
  virtual ~SubPlot();
  void AddPoint(double x, double data, double uncorrerr, double toterr, double theory, double theory_mod);
  void PrepareHistogram(bool RatioToData);
};


class  DataSet {
 private:    
  TString* fName;
  Int_t fSetId;
  map<int, SubPlot*> fSubPlots;

 public: 
  DataSet();
  DataSet(Int_t SetId, const Char_t* name);
  virtual ~DataSet();
  Int_t inline GetSetId() {return fSetId;};
  inline const Char_t* GetName() {return fName->Data();}
  
  void AddNewPlot(const char* descriptor);
  void AddPoint(const char* s, double data, double uncorrerr, double toterr, double theory, double theory_mod);
  inline int GetNSubPlots() {return fSubPlots.size();}
  
  TH1F* GetHistogram(int i, bool RatioToData);
  TGraphErrors* GetDataUnc(int i);
  TGraphErrors* GetDataTot(int i);
  TGraphErrors* GetTheo(int i);
  TGraphErrors* GetTMod(int i);
  bool          GetXlog(int i);
  bool          GetYlog(int i);
  
 private:
  SubPlot* GetSubPlot(int idx);
  
      };

#endif
