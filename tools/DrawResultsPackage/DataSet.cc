#include "DataSet.h"
#include <TH1F.h>
#include <stdlib.h>
#include <TMath.h>

SubPlot::SubPlot(const char* descriptor) {
  fTitle = new TString;
  fPlotOption = new TString(descriptor);
  fDUnc = new TGraphErrors;
  fDTot = new TGraphErrors;
  fTheo = new TGraphErrors;
  fTMod = new TGraphErrors;
  fHistogram = NULL;
  fXlog = kFALSE;
  fYlog = kFALSE;
}
SubPlot::~SubPlot() {
  delete fTitle;
  delete fPlotOption;
  delete fDUnc;
  delete fDTot;
  delete fTheo;
  delete fTMod;
  if(fHistogram) delete fHistogram;
}


void SubPlot::AddPoint(double x, double data, double uncorrerr, double toterr, double theory, double theory_mod) {
  Int_t N = fDUnc->GetN();
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
  fDUnc->Sort();
  fDTot->Sort();
  fTheo->Sort();
  fTMod->Sort();
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
    }
    else if(str.BeginsWith("Xmax:")) {
      str.ReplaceAll("Xmax:","");
      Xmax = str.Atof();
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

    if(str.BeginsWith("Title:")) {
      str.ReplaceAll("Title:","");
      fHistogram->SetTitle(str.Data());
    }
    else if(str.BeginsWith("XTitle:")) {
      str.ReplaceAll("XTitle:","");
      fHistogram->GetXaxis()->SetTitle(str.Data());
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
    else if(str.BeginsWith("YminR:")) {
      if(RatioToData) {
	str.ReplaceAll("YminR:","");
	fHistogram->SetMinimum(str.Atof());
      }
    }
    else if(str.BeginsWith("YmaxR:")) {
      if(RatioToData) {
	str.ReplaceAll("YmaxR:","");
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

void DataSet::AddPoint(const char* s, double data, double uncorrerr, double toterr, double theory, double theory_mod) {
  TString str(s);
  TObjArray* array = str.Tokenize("/");
  int iplot = ((TObjString*) array->At(0))->GetString().Atoi();
  double x = ((TObjString*) array->At(1))->GetString().Atof();
  map<int,SubPlot*>::iterator it = fSubPlots.find(iplot);
  if(it == fSubPlots.end()) fSubPlots[iplot] = new SubPlot("");
  fSubPlots[iplot]->AddPoint(x, data, uncorrerr, toterr, theory, theory_mod);
}
