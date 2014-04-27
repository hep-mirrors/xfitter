#include "DataSet.h"
#include <TH1F.h>
#include <stdlib.h>
#include <TMath.h>

SubPlot::SubPlot(const char* descriptor) {
  fPlotOption = new TString(descriptor);
  fDUnc = new TGraphErrors;
  fDTot = new TGraphErrors;
  fTheo = new TGraphErrors;
  fTMod = new TGraphErrors;
  fHistogram = NULL;
  fXlog = kFALSE;
  fYlog = kFALSE;
  fXTitle = "";
  fYTitle = "";
  fTitle = "";
  fLabel = "";
  fExp = "";
  fXmin = 0;
  fXmax = 0;
  fYminR = 0;
  fYmaxR = 0;
}
SubPlot::~SubPlot() {
  delete fPlotOption;
  delete fDUnc;
  delete fDTot;
  delete fTheo;
  delete fTMod;
  if(fHistogram) delete fHistogram;
}


void SubPlot::AddPoint(double x, double bin1, double bin2, double data, double uncorrerr, double toterr, double theory, double theory_mod, double therr_up, double therr_down, double pull) {
  Int_t N = fDUnc->GetN();
  Bins1.push_back(bin1);
  Bins2.push_back(bin2);
  Data.push_back(data);
  Uncor.push_back(uncorrerr);
  Toterr.push_back(toterr);
  Th.push_back(theory);
  Thshift.push_back(theory_mod);
  ThErrUp.push_back(therr_up);
  ThErrDown.push_back(therr_down);
  Pull.push_back(pull);
  fDUnc->Set(N+1);
  fDTot->Set(N+1);
  fTheo->Set(N+1);
  fTMod->Set(N+1);
  fDUnc->SetPoint(N, x, data);
  fDUnc->SetPointError(N, 0., uncorrerr);
  fDTot->SetPoint(N, x, data);
  fDTot->SetPointError(N, 0., toterr);
  fTheo->SetPoint(N, x, theory);
  fTheo->SetPointError(N, 0., 0.);
  fTMod->SetPoint(N, x, theory_mod);
  fTMod->SetPointError(N, 0., 0.);
  //  fDUnc->Sort();
  //  fDTot->Sort();
  //  fTheo->Sort();
  //  fTMod->Sort();
};

void SubPlot::PrepareHistogram(bool RatioToData) {
  
  double Xmin = -9999.;
  double Xmax = 9999.;
  fXlog = kFALSE;
  fYlog = kFALSE;
  
  TObjArray* array = fPlotOption->Tokenize("@");
  
  for(int i=0; i<array->GetEntries(); i++) { // first loop to detect x axis
    TString str( ((TObjString*)array->At(i))->GetString().Data());
    if(str.BeginsWith("Xmin:")) {
      str.ReplaceAll("Xmin:","");
      Xmin = str.Atof();
      fXmin = str.Atof();
    }
    else if(str.BeginsWith("Xmax:")) {
      str.ReplaceAll("Xmax:","");
      Xmax = str.Atof();
      fXmax = str.Atof();
    }
  }
  if(Xmin == -9999.) {
    double min = 99999.;
    double max = -99999.;
    for(int i=0; i<fDUnc->GetN(); i++) {
      if(min > fDUnc->GetX()[i]) min = fDUnc->GetX()[i];
      if(max < fDUnc->GetX()[i]) max = fDUnc->GetX()[i];
    }
    Xmin = min - (max-min)*0.1;
  }

  if(Xmax == 9999.) {
    double min = 99999.;
    double max = -99999.;
    for(int i=0; i<fDUnc->GetN(); i++) {
      if(min > fDUnc->GetX()[i]) min = fDUnc->GetX()[i];
      if(max < fDUnc->GetX()[i]) max = fDUnc->GetX()[i];
    }
    Xmax = max + (max-min)*0.1;
  }
  if(Xmin == Xmax) {  // if there is only one point on the x axis
    Xmin *= 0.9;
    Xmax *= 1.1;
  }


  if(fHistogram) delete fHistogram;
  fHistogram = new TH1F("","",2, Xmin, Xmax);
  fHistogram->SetStats(kFALSE);
  fHistogram->SetTitle("unknown");


  // set default 
  if(!RatioToData) {
    fHistogram->SetMaximum(fDTot->GetHistogram()->GetMaximum());
    fHistogram->SetMaximum( TMath::Max(fTheo->GetHistogram()->GetMaximum(), fHistogram->GetMaximum()));
    fHistogram->SetMaximum( TMath::Max(fTMod->GetHistogram()->GetMaximum(), fHistogram->GetMaximum()));
    fHistogram->SetMaximum( fHistogram->GetMaximum() * 1.1);
    
    fHistogram->SetMinimum(fDTot->GetHistogram()->GetMinimum());
    fHistogram->SetMinimum( TMath::Max(fTheo->GetHistogram()->GetMinimum(), fHistogram->GetMinimum()));
    fHistogram->SetMinimum( TMath::Max(fTMod->GetHistogram()->GetMinimum(), fHistogram->GetMinimum()));
    fHistogram->SetMinimum( fHistogram->GetMinimum() * 0.9);
  }   
  else {
    fHistogram->SetMaximum(0.5);
    fHistogram->SetMinimum(-0.5);
  }

  for(int i=0; i<array->GetEntries(); i++) {
    TString str( ((TObjString*)array->At(i))->GetString().Data());

    if(str.BeginsWith("ExtraLabel:")) {
      str.ReplaceAll("ExtraLabel:","");
      fLabel = str.Data();
    }

    if(str.BeginsWith("Experiment:")) {
      str.ReplaceAll("Experiment:","");
      fExp = str.Data();
    }

    if(str.BeginsWith("Title:")) {
      str.ReplaceAll("Title:","");
      fHistogram->SetTitle(str.Data());
      fTitle = str.Data();
    }
    else if(str.BeginsWith("XTitle:")) {
      str.ReplaceAll("XTitle:","");
      fHistogram->GetXaxis()->SetTitle(str.Data());
      fXTitle = str.Data();
    }
    else if(str.BeginsWith("YTitle:")) {
      str.ReplaceAll("YTitle:","");
      fHistogram->GetYaxis()->SetTitle(str.Data());
      fYTitle = str.Data();
    }
    else if(str.BeginsWith("YminR:")) {
      str.ReplaceAll("YminR:","");
      fHistogram->SetMinimum(str.Atof());
      fYminR = str.Atof();
    }
    else if(str.BeginsWith("YmaxR:")) {
      str.ReplaceAll("YmaxR:","");
      fHistogram->SetMaximum(str.Atof());
      fYmaxR = str.Atof();
    }
    else if(str.BeginsWith("Ymin:")) {
      if(!RatioToData) {
	str.ReplaceAll("Ymin:","");
	fHistogram->SetMinimum(str.Atof());
      }
    }
    else if(str.BeginsWith("Ymax:")) {
      if(!RatioToData) {
	str.ReplaceAll("Ymax:","");
	fHistogram->SetMaximum(str.Atof());
      }
    }
    else if(str.BeginsWith("Xlog")) {
      fXlog = kTRUE;
    }
    else if(str.BeginsWith("Ylog")) {
      if(!RatioToData) {
	fYlog = kTRUE;
      }
    }

  }
  delete array;
}



DataSet::DataSet(Int_t SetId, const Char_t* name) {
  fName = new TString(name);
  fSetId = SetId;
}
DataSet::DataSet() {
  fName = new TString;
}

DataSet::~DataSet(){
  delete fName;
}
 
 SubPlot* DataSet::GetSubPlot(int idx) {
  int i=0;
  if(idx >= fSubPlots.size()) return NULL;
  map<int, SubPlot*>::iterator it = fSubPlots.begin();
  while(i!=idx) {it++; i++;}
  return (*it).second;
}
 
TH1F* DataSet::GetHistogram(int idx, bool RatioToData) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) return NULL;
  s->PrepareHistogram(RatioToData);
  return s->fHistogram;
}

TGraphErrors*  DataSet::GetDataUnc(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) return NULL;
  return s->fDUnc;
}
TGraphErrors*  DataSet::GetDataTot(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) return NULL;
  return s->fDTot;
}
TGraphErrors*  DataSet::GetTheo(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) return NULL;
  return s->fTheo;
}
TGraphErrors*  DataSet::GetTMod(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) return NULL;
  return s->fTMod;
}
bool  DataSet::GetXlog(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) return kFALSE;
  return s->fXlog;
}
bool  DataSet::GetYlog(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) return kFALSE;
  return s->fYlog;
}
float  DataSet::GetXmin(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) return 0;
  return s->fXmin;
}
float  DataSet::GetXmax(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) return 0;
  return s->fXmax;
}
float  DataSet::GetYminR(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) return 0;
  return s->fYminR;
}
float  DataSet::GetYmaxR(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) return 0;
  return s->fYmaxR;
}
string  DataSet::GetXTitle(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) return "";
  return s->fXTitle;
}
string  DataSet::GetYTitle(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) return "";
  return s->fYTitle;
}
string  DataSet::GetTitle(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) return "";
  return s->fTitle;
}
string  DataSet::GetLabel(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) return "";
  return s->fLabel;
}
string  DataSet::GetExperiment(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) return "";
  return s->fExp;
}
vector <float>  DataSet::getbins1(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) {vector<float> v; return v;}
  return s->Bins1;
}
vector <float>  DataSet::getbins2(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) {vector<float> v; return v;}
  return s->Bins2;
}
vector <float>  DataSet::getdata(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) {vector<float> v; return v;}
  return s->Data;
}
vector <float>  DataSet::getuncor(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) {vector<float> v; return v;}
  return s->Uncor;
}
vector <float>  DataSet::gettoterr(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) {vector<float> v; return v;}
  return s->Toterr;
}
vector <float>  DataSet::gettheory(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) {vector<float> v; return v;}
  return s->Th;
}
vector <float>  DataSet::gettheoryshifted(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) {vector<float> v; return v;}
  return s->Thshift;
}
vector <float>  DataSet::gettherrup(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) {vector<float> v; return v;}
  return s->ThErrUp;
}
vector <float>  DataSet::gettherrdown(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) {vector<float> v; return v;}
  return s->ThErrDown;
}
vector <float>  DataSet::getpulls(int idx) {
  SubPlot* s = GetSubPlot(idx);
  if(!s) {vector<float> v; return v;}
  return s->Pull;
}

void DataSet::AddNewPlot(const char* descriptor) {
  TString temp(descriptor);
  TObjArray* array = temp.Tokenize("@");
  temp.Form(((TObjString*)  array->At(0))->GetString().Data());
  delete array;
  temp.ReplaceAll("Plot","");
  int iplot = temp.Atoi();
  fSubPlots.erase(iplot);
  fSubPlots[iplot] = new SubPlot(descriptor);  
}

void DataSet::AddPoint(const char* s, double bin1, double bin2, double data, double uncorrerr, double toterr, double theory, double theory_mod, double therr_up, double therr_down, double pull) {
  TString str(s);
  TObjArray* array = str.Tokenize("/");
  int iplot = ((TObjString*) array->At(0))->GetString().Atoi();
  double x = ((TObjString*) array->At(1))->GetString().Atof();
  map<int,SubPlot*>::iterator it = fSubPlots.find(iplot);
  if(it == fSubPlots.end()) fSubPlots[iplot] = new SubPlot("");
  fSubPlots[iplot]->AddPoint(x, bin1, bin2, data, uncorrerr, toterr, theory, theory_mod, therr_up, therr_down, pull);
}
