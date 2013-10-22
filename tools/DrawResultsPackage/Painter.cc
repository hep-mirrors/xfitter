#include <Painter.h>
#include <TPaveLabel.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TGraphAsymmErrors.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TMath.h>


Painter::Painter(bool DrawBands, bool DrawBase, bool DrawExp, bool DrawModel, bool DrawParam){
  fPath = new TString("../../output/");
  fPathRef = new TString("../../output/");
  fOutput = NULL;
  fOutputRef = NULL;
  fPsFileName = new TString("DrawResults.ps");
  gROOT->SetStyle("Plain");
  cout << endl;
  //cout << "TO DO: in fittedresults.txt q2 and x for sets 61-64 are switched"<<endl;
  fColor = kRed;
  fColorRef = kBlue;
  fFillStyle = 1001;
  fFillStyleRef = 0; //3010
  fBands = DrawBands;
  fHighlight = kRed;
  
  fDrawBase = DrawBase;
  fDrawExp = DrawExp;
  fDrawModel = DrawModel;
  fDrawParam = DrawParam;
  fPathBase = new TString("../../output/");
  fPathExp = new TString("../../output/");
  fPathModel = new TString("../../output/");
  fPathParam = new TString("../../output/");
  fOutputBase = NULL;
  fOutputExp = NULL;
  fOutputModel = NULL;
  fOutputParam = NULL;
  fColorBase = kBlue;
  fColorExp = kRed;
  fColorModel = kYellow;
  fColorParam = kGreen;
  fFillStyleBase = 0;
  fFillStyleExp = 1001;
  fFillStyleModel = 1001;
  fFillStyleParam = 1001;
}

Painter::~Painter(){ 
  delete fPath;
  delete fPathRef;
  if(fOutput) delete fOutput;
  if(fOutputRef) delete fOutputRef;

  TCanvas* can = new TCanvas;
  cout << "Output stored in " << fPsFileName->Data() << " file"<<endl;
  fPsFileName->Append(")");
  can->Print(fPsFileName->Data());
  delete can;

  delete fPsFileName;
  
  delete fPathBase;
  delete fPathExp;
  delete fPathModel;
  delete fPathParam;
  if(fOutputBase) delete fOutputBase;
  if(fOutputExp) delete fOutputExp;
  if(fOutputModel) delete fOutputModel;
  if(fOutputParam) delete fOutputParam;
}

Int_t Painter::Prepare() {
  
 if ( ! (fDrawExp||fDrawModel||fDrawParam) ) {
   
  fOutput = new Output(fPath->Data());
  if(fPath->CompareTo(fPathRef->Data()))   //if fPath and fPathRef are different
    fOutputRef = new Output(fPathRef->Data());
  else 
    fOutputRef = NULL;

  if(fOutputRef == NULL) fPsFileName->Form("%s/DrawResults.ps", fPathRef->Data());
  fPsFileName->Append("(");

  fOutput->Prepare(fBands);
  if(fOutputRef) fOutputRef->Prepare(fBands);
  
 } else {
   
  fPsFileName->Append("(");
   
  if(fDrawBase) { 
    fOutputBase = new Output(fPathBase->Data());
    fOutputBase->Prepare(false);
  }
  if(fDrawExp) {
    fOutputExp = new Output(fPathExp->Data());
    fOutputExp->Prepare(true, TString("e"));
  }
  if(fDrawModel) {
    fOutputModel = new Output(fPathModel->Data());
    fOutputModel->Prepare(true, TString("m"));
  }
  if(fDrawParam) {
    fOutputParam = new Output(fPathParam->Data());
    fOutputParam->Prepare(true, TString("p"));
  }
   
 } 
 return 0;
}

Int_t Painter::Draw() {
  
  this->Prepare();
  
 if ( ! (fDrawExp||fDrawModel||fDrawParam) ) {

  // Get number of Q2 files:
  Int_t NQ2Files = fOutput->GetNQ2Files();

  for(Int_t i=0; i<NQ2Files; i++) {
    this->DrawPDF(i);
  }
  
//  this->DrawPull();
  
  Int_t NRefDataSets = 0;
  if(fOutputRef) NRefDataSets = fOutputRef->GetNsets();
  Bool_t RefSetsDrawn[NRefDataSets];
  for(Int_t i=0; i<NRefDataSets; i++) RefSetsDrawn[i] = kFALSE;

  for (Int_t i=0; i<fOutput->GetNsets(); i++) {
    DataSet* dataset = fOutput->GetSet(i);
    DataSet* datasetref = NULL;
    for (Int_t iref=0; iref<NRefDataSets; iref++) {
      if(fOutputRef->GetSet(iref)->GetSetId() == dataset->GetSetId()) {
	datasetref = fOutputRef->GetSet(iref);
	RefSetsDrawn[iref] = kTRUE;
	break;
      }
    } 
    DrawDataSet(dataset, datasetref, kFALSE, fColor);
    DrawDataSet(dataset, datasetref, kTRUE, fColor);
  }
  for (Int_t iref=0; iref<NRefDataSets; iref++) {
    if(!RefSetsDrawn[iref])
      DrawDataSet(fOutputRef->GetSet(iref), NULL, kFALSE, fColorRef);
  }

  this->DrawFitResults();

  this->DrawCorrelations(fOutput);
  if(fOutputRef) 
    this->DrawCorrelations(fOutputRef);

  this->DrawMessages(fOutput);
  if(fOutputRef) 
    this->DrawMessages(fOutputRef);

 } else {
   
   Int_t NQ2Files = 0;
   
   if (fDrawBase) {
     NQ2Files = fOutputBase->GetNQ2Files();
   } else if (fDrawExp) {
     NQ2Files = fOutputExp->GetNQ2Files();
   } else if (fDrawModel) {
     NQ2Files = fOutputModel->GetNQ2Files();
   } else if (fDrawParam) {
     NQ2Files = fOutputParam->GetNQ2Files();
   }
   
   for(Int_t i=0; i<NQ2Files; i++) {
    this->DrawPDF(i);
   }
   
   
   Int_t NMainDataSets = 0;
   
   if (fDrawBase) NMainDataSets = fOutputBase->GetNsets();
    else if (fDrawExp) NMainDataSets = fOutputExp->GetNsets();
     else if (fDrawModel) NMainDataSets = fOutputModel->GetNsets();
      else if (fDrawParam) NMainDataSets = fOutputParam->GetNsets();
   
   Int_t MainDataSetsIds[NMainDataSets];
   
   if (fDrawBase) for (int i=0; i<NMainDataSets; i++) MainDataSetsIds[i] = fOutputBase->GetSet(i)->GetSetId();
    else if (fDrawExp) for (int i=0; i<NMainDataSets; i++) MainDataSetsIds[i] = fOutputExp->GetSet(i)->GetSetId();
     else if (fDrawModel) for (int i=0; i<NMainDataSets; i++) MainDataSetsIds[i] = fOutputModel->GetSet(i)->GetSetId();
      else if (fDrawParam) for (int i=0; i<NMainDataSets; i++) MainDataSetsIds[i] = fOutputParam->GetSet(i)->GetSetId();
   
   Int_t NDataSetsBase = 0;
   if (fDrawBase) NDataSetsBase = fOutputBase->GetNsets();
   Bool_t SetsDrawnBase[NDataSetsBase];
   for(Int_t i=0; i<NDataSetsBase; i++) SetsDrawnBase[i] = kFALSE;
   
   Int_t NDataSetsExp = 0;
   if (fDrawExp) NDataSetsExp = fOutputExp->GetNsets();
   Bool_t SetsDrawnExp[NDataSetsExp];
   for(Int_t i=0; i<NDataSetsExp; i++) SetsDrawnExp[i] = kFALSE;
   
   Int_t NDataSetsModel = 0;
   if (fDrawModel) NDataSetsModel = fOutputModel->GetNsets();
   Bool_t SetsDrawnModel[NDataSetsModel];
   for(Int_t i=0; i<NDataSetsModel; i++) SetsDrawnModel[i] = kFALSE;
   
   Int_t NDataSetsParam = 0;
   if (fDrawParam) NDataSetsParam = fOutputParam->GetNsets();
   Bool_t SetsDrawnParam[NDataSetsParam];
   for(Int_t i=0; i<NDataSetsParam; i++) SetsDrawnParam[i] = kFALSE;
   
   DataSet* datasetBase = NULL;
   DataSet* datasetExp = NULL;
   DataSet* datasetModel = NULL;
   DataSet* datasetParam = NULL;
   
   for (Int_t i=0; i<NMainDataSets; i++) {
    for (Int_t iref=0; iref<NDataSetsBase; iref++) {
      if(fOutputBase->GetSet(iref)->GetSetId() == MainDataSetsIds[i]) {
	datasetBase = fOutputBase->GetSet(iref);
	SetsDrawnBase[iref] = kTRUE;
	break;
      }
    }
    for (Int_t iref=0; iref<NDataSetsExp; iref++) {
      if(fOutputExp->GetSet(iref)->GetSetId() == MainDataSetsIds[i]) {
	datasetExp = fOutputExp->GetSet(iref);
	SetsDrawnExp[iref] = kTRUE;
	break;
      }
    }
    for (Int_t iref=0; iref<NDataSetsModel; iref++) {
      if(fOutputModel->GetSet(iref)->GetSetId() == MainDataSetsIds[i]) {
	datasetModel = fOutputModel->GetSet(iref);
	SetsDrawnModel[iref] = kTRUE;
	break;
      }
    }
    for (Int_t iref=0; iref<NDataSetsParam; iref++) {
      if(fOutputParam->GetSet(iref)->GetSetId() == MainDataSetsIds[i]) {
	datasetParam = fOutputParam->GetSet(iref);
	SetsDrawnParam[iref] = kTRUE;
	break;
      }
    }
    DrawDataSetEMP(datasetBase, datasetExp, datasetModel, datasetParam, kFALSE);
    DrawDataSetEMP(datasetBase, datasetExp, datasetModel, datasetParam, kTRUE);
  }
  for (Int_t iref=0; iref<NDataSetsBase; iref++) {
    if(!SetsDrawnBase[iref])
      DrawDataSetEMP(datasetBase, NULL, NULL, NULL, kFALSE);
  }
  for (Int_t iref=0; iref<NDataSetsExp; iref++) {
    if(!SetsDrawnExp[iref])
      DrawDataSetEMP(NULL, datasetExp, NULL, NULL, kFALSE);
  }
  for (Int_t iref=0; iref<NDataSetsModel; iref++) {
    if(!SetsDrawnModel[iref])
      DrawDataSetEMP(NULL, NULL, datasetModel, NULL, kFALSE);
  }
  for (Int_t iref=0; iref<NDataSetsParam; iref++) {
    if(!SetsDrawnParam[iref])
      DrawDataSetEMP(NULL, NULL, NULL, datasetParam, kFALSE);
  }
  
  this->DrawFitResultsEMP();

  //   cout << "I should be able to draw more fancy plots soon." << endl;
  
  if (fDrawBase) this->DrawCorrelations(fOutputBase);
  if (fDrawExp) this->DrawCorrelations(fOutputExp);
  if (fDrawModel) this->DrawCorrelations(fOutputModel);
  if (fDrawParam) this->DrawCorrelations(fOutputParam);

  if (fDrawBase) this->DrawMessages(fOutputBase);
  if (fDrawExp) this->DrawMessages(fOutputExp);
  if (fDrawModel) this->DrawMessages(fOutputModel);
  if (fDrawParam) this->DrawMessages(fOutputParam);
   
 }
 
 return 0;
}

TText* Painter::AddLineToPave(TObjArray* paves, float& yposition, const char* text, const char* option) {
  static int NLines = 30;
  //  static float yposition = 1.;

  TString Option(option);
  TPaveText* pave = (TPaveText*) paves->At(paves->GetEntries()-1);
  if(!pave) return NULL;

  if(pave->GetListOfLines()->GetEntries() >= NLines-1) {
    TPaveText* pavenew = new TPaveText(pave->GetX1(), pave->GetY1(), pave->GetX2(), pave->GetY2());
    paves->AddLast(pavenew);
    pave = pavenew;
    yposition = 1.;
  }

  yposition -= 1./(float)NLines;
  TText* T = pave->AddText(0.1, yposition, text);
  T->SetTextAlign(12);
  
  if (Option.Contains("B"))   T->SetTextFont(102);
  else                        T->SetTextFont(82);
  if      (Option.Contains("N"))   T->SetTextColor(fColor);
  else if (Option.Contains("R"))   T->SetTextColor(fColorRef);
  else if (Option.Contains("A"))   T->SetTextColor(fColorBase);
  else if (Option.Contains("E"))   T->SetTextColor(fColorExp);
  else if (Option.Contains("M"))   T->SetTextColor(fColorModel);
  else if (Option.Contains("P"))   T->SetTextColor(fColorParam);

  T->SetTextSize(0.025);
  return T;
}

void Painter::FillPavesWithFitResults(TObjArray* paves, Output* output) {

  //  TObjArray* names = output->GetFittedParametersNames();
  //TObjArray* namesNuisance= output->GetNuisanceParNames();
  Int_t NNuisance = output->GetNNuisanceParameters();
  TString str;
  
  float ypos = 1.;
  int NChar = 40;
  
  if(!output->GetParametersCheck()) return;
  Int_t NFittedParameters = output->GetNFittedParameters();
  
  AddLineToPave(paves, ypos, "Results for:","B");
  
  for(int i=0; i<=output->GetName()->Length()/NChar; i++) {
    if(output==fOutput)
      AddLineToPave(paves, ypos, TString((*(output->GetName()))(i*NChar, NChar)).Data() ,"BN");
    else if(output==fOutputBase)
      AddLineToPave(paves, ypos, TString((*(output->GetName()))(i*NChar, NChar)).Data() ,"BA");
    else if(output==fOutputExp)
      AddLineToPave(paves, ypos, TString((*(output->GetName()))(i*NChar, NChar)).Data() ,"BE");
    else if(output==fOutputModel)
      AddLineToPave(paves, ypos, TString((*(output->GetName()))(i*NChar, NChar)).Data() ,"BM");
    else if(output==fOutputParam)
      AddLineToPave(paves, ypos, TString((*(output->GetName()))(i*NChar, NChar)).Data() ,"BP");
    else
      AddLineToPave(paves, ypos, TString((*(output->GetName()))(i*NChar, NChar)).Data() ,"BR");
  }
  
  AddLineToPave(paves, ypos, "","");
  
  str.Form("Fitted %d parameters:", NFittedParameters);; 
  AddLineToPave(paves, ypos, str.Data(),"");
  str.Form("(most reliable available method: %s" , output->GetErrorCalculationMethod()->Data()); 
  AddLineToPave(paves, ypos, str.Data(),"");
  str.Form("giving confidence in errors: %s)", output->GetErrorTrustLevel()->Data()); 
  AddLineToPave(paves, ypos, str.Data(),"");
  
  
  for(int i=0; i<NFittedParameters; i++) {
    str.Form("%4d:\t%6s = %6.3f  #pm%6.3f", i+1, output->GetFittedParametersName(i)->Data(),
	     output->GetFittedParameter(i, false), output->GetFittedParameter(i, true));
    AddLineToPave(paves, ypos, str.Data(),"");
  }
  if(!output->GetNuisanceCheck()) return;
  
  AddLineToPave(paves, ypos, "","");
  AddLineToPave(paves, ypos, "Nuisance Parameters:","");
  for(int i=0; i<NNuisance; i++) {
    if(!output->GetNuisanceParNames(i)) continue;
    str.Form("%4d:\t%17s = %5.2f  #pm%5.2f", i+1, output->GetNuisanceParNames(i)->Data(),
	     output->GetNuisancePar(i, false), output->GetNuisancePar(i, true));
    TText* T = AddLineToPave(paves, ypos, str.Data(),""); 
    if(output->GetNuisancePar(i, false) > 2.5 || output->GetNuisancePar(i, false) < -2.5 ) {
      T->SetTextColor(fHighlight);
      T->SetTextFont(102);
    }
  }
}

Int_t Painter::DrawFitResults() {

  TCanvas* can = new TCanvas;
  TObjArray* pavesL = new TObjArray; pavesL->SetOwner();
  TObjArray* pavesR = new TObjArray; pavesR->SetOwner();

  pavesL->AddLast(new TPaveText(0.05, 0.05, 0.5, 0.95, "br"));
  FillPavesWithFitResults(pavesL, fOutput);
  if(fOutputRef) {
    pavesR->AddLast(new TPaveText(0.5, 0.05, 0.95, 0.95, "br"));
    FillPavesWithFitResults(pavesR, fOutputRef);
  }

  for(int j=0; j<TMath::Max(pavesL->GetEntries(), pavesR->GetEntries()); j++) {
    TPaveText* paveL = (TPaveText*) pavesL->At(j);
    TPaveText* paveR = (TPaveText*) pavesR->At(j);
    TCanvas* can = new TCanvas;
    if(paveL) paveL->Draw();
    if(paveR) paveR->Draw();
    PrintCanvas(can);
    delete can;
  }
  delete pavesL; 
  delete pavesR; 
  delete can;
}

Int_t Painter::DrawFitResultsEMP() {

  bool isSamplDrawn[4]; //0 - base, 1 - exp, 2 - model, 3 - param.
  bool isPavesInUse[2]; //0 - pavesL, 1 - pavesR.
  
  for (int i=0; i<4; i++) isSamplDrawn[i] = false;
  
      TCanvas* can = new TCanvas;
    TObjArray* pavesL = new TObjArray; pavesL->SetOwner();
    TObjArray* pavesR = new TObjArray; pavesR->SetOwner();
  
  for (int jc=0; jc<2; jc++) {
    
    for (int i=0; i<2; i++) isPavesInUse[i] = false;
    
    if (fDrawBase && !isSamplDrawn[0]) {
      if (!isPavesInUse[0]) {
	pavesL->AddLast(new TPaveText(0.05, 0.05, 0.5, 0.95, "br"));
	FillPavesWithFitResults(pavesL, fOutputBase);
	isSamplDrawn[0] = true;
	isPavesInUse[0] = true;
      } else if (!isPavesInUse[1]) {
	pavesR->AddLast(new TPaveText(0.5, 0.05, 0.95, 0.95, "br"));
	FillPavesWithFitResults(pavesR, fOutputBase);
	isSamplDrawn[0] = true;
	isPavesInUse[1] = true;
      }
    }
    if (fDrawExp && !isSamplDrawn[1]) {
      if (!isPavesInUse[0]) {
	pavesL->AddLast(new TPaveText(0.05, 0.05, 0.5, 0.95, "br"));
	FillPavesWithFitResults(pavesL, fOutputExp);
	isSamplDrawn[1] = true;
	isPavesInUse[0] = true;
      } else if (!isPavesInUse[1]) {
	pavesR->AddLast(new TPaveText(0.5, 0.05, 0.95, 0.95, "br"));
	FillPavesWithFitResults(pavesR, fOutputExp);
	isSamplDrawn[1] = true;
	isPavesInUse[1] = true;
      }
    }
    if (fDrawModel && !isSamplDrawn[2]) {
      if (!isPavesInUse[0]) {
	pavesL->AddLast(new TPaveText(0.05, 0.05, 0.5, 0.95, "br"));
	FillPavesWithFitResults(pavesL, fOutputModel);
	isSamplDrawn[2] = true;
	isPavesInUse[0] = true;
      } else if (!isPavesInUse[1]) {
	pavesR->AddLast(new TPaveText(0.5, 0.05, 0.95, 0.95, "br"));
	FillPavesWithFitResults(pavesR, fOutputModel);
	isSamplDrawn[2] = true;
	isPavesInUse[1] = true;
      }
    }
    
    if (fDrawParam && !isSamplDrawn[3]) {
      if (!isPavesInUse[0]) {
	pavesL->AddLast(new TPaveText(0.05, 0.05, 0.5, 0.95, "br"));
	FillPavesWithFitResults(pavesL, fOutputParam);
	isSamplDrawn[3] = true;
	isPavesInUse[0] = true;
      } else if (!isPavesInUse[1]) {
	pavesR->AddLast(new TPaveText(0.5, 0.05, 0.95, 0.95, "br"));
	FillPavesWithFitResults(pavesR, fOutputParam);
	isSamplDrawn[3] = true;
	isPavesInUse[1] = true;
      }
    }
  }
  
    for(int j=0; j<TMath::Max(pavesL->GetEntries(), pavesR->GetEntries()); j++) {
      TPaveText* paveL = (TPaveText*) pavesL->At(j);
      TPaveText* paveR = (TPaveText*) pavesR->At(j);
      TCanvas* can2 = new TCanvas;
      if(paveL) paveL->Draw();
      if(paveR) paveR->Draw();
      PrintCanvas(can2);
      delete can2;
    }
    delete pavesL; 
    delete pavesR; 
    delete can;


  return 0;
}

void Painter::DrawCorrelations(Output* output) {

  if(!output->GetCovarianceCheck()) return;
  TCanvas* can = new TCanvas;
  TPaveText* pave = new TPaveText(0.05, 0.05, 0.95, 0.95);
  TString* str = new TString;
    
  //TObjArray* names = output->GetFittedParametersNames();
  Int_t NFittedParameters = output->GetNFittedParameters();
  float xpos = 0.05;
  float ypos = 0.95;
  float ystep = 0.9 / (NFittedParameters+4);
  float xstep = 0.8 / (NFittedParameters+3);
  TText* T;

  T = pave->AddText(xpos, ypos, "Estimated correlation factors for"); ypos -= ystep;
  T->SetTextFont(102); 
  
  T = pave->AddText(xpos, ypos, output->GetName()->Data()); T->SetTextFont(102); 
  if(output==fOutput) T->SetTextColor(fColor);
  else                        T->SetTextColor(fColorRef);
  ypos -= ystep;

  str->Form("(most reliable available method: %s, giving confidence in correlation estimates: %s)", 
	    output->GetCorrelationCalculationMethod()->Data(),output->GetErrorTrustLevel()->Data());
  T = pave->AddText(xpos, ypos, str->Data());
  
  ypos -= ystep;
  ypos -= ystep;
  
  for(int i=0; i<NFittedParameters; i++) {
    str->Form("%d. %s", i+1, output->GetFittedParametersName(i)->Data());
    pave->AddText(xpos, ypos, str->Data());
    for(int j=0; j<=i; j++) {
      xpos = 0.2+j*xstep;
      
      if(i==NFittedParameters-1) {
	str->Form("   %d.", j+1);
	T = pave->AddText(xpos, ypos-ystep, str->Data());
      }
      
      Double_t cor = output->GetCorPar(i, j);
      str->Form("%6.2f", cor);
      T = pave->AddText(xpos, ypos, str->Data());
      if(i!=j && (cor>0.9 || cor<-0.9)) {
	T->SetTextColor(fHighlight);
	T->SetTextFont(102);
      }
    }
    xpos = 0.05;
    ypos -= ystep;
  }
  
  for(int i=0; i<pave->GetListOfLines()->GetEntries(); i++) {
    T = (TText*) pave->GetListOfLines()->At(i);
    T->SetTextAlign(11);  
    T->SetTextSize(0.02);
    if(T->GetTextFont()==102)     T->SetTextFont(102); 
    else                          T->SetTextFont(82); 
  }

  pave->Draw();
  PrintCanvas(can);
  delete can;
  delete pave;
  delete str;
}

void Painter::DrawMessages(Output* output) {
  if(!output->GetMessagesCheck()) return;
  TObjArray* paves = new TObjArray; paves->SetOwner();
  paves->AddLast(new TPaveText(0.05, 0.05, 0.95, 0.95, "br"));
  float ypos = 1.;
  int NChar = 80;
  
  AddLineToPave(paves, ypos, "Fit messages for:","B");
  for(int i=0; i<=output->GetName()->Length()/NChar; i++) {
    if(output==fOutput)
      AddLineToPave(paves, ypos, TString((*(output->GetName()))(i*NChar, NChar)).Data() ,"BN");
    else if(output==fOutputBase)
      AddLineToPave(paves, ypos, TString((*(output->GetName()))(i*NChar, NChar)).Data() ,"BA");
    else if(output==fOutputExp)
      AddLineToPave(paves, ypos, TString((*(output->GetName()))(i*NChar, NChar)).Data() ,"BE");
    else if(output==fOutputModel)
      AddLineToPave(paves, ypos, TString((*(output->GetName()))(i*NChar, NChar)).Data() ,"BM");
    else if(output==fOutputParam)
      AddLineToPave(paves, ypos, TString((*(output->GetName()))(i*NChar, NChar)).Data() ,"BP");
    else
      AddLineToPave(paves, ypos, TString((*(output->GetName()))(i*NChar, NChar)).Data() ,"BR");
  }

  AddLineToPave(paves, ypos, "","");


  TObjArray* messages = output->GetMessages();
  for(int i=0; i<messages->GetEntries(); i++) {
    AddLineToPave(paves, ypos, ((TObjString*)messages->At(i))->GetString().Data(),"");
  }
  
  for(int j=0; j<paves->GetEntries(); j++) {
    TPaveText* pave = (TPaveText*) paves->At(j);
    TCanvas* can = new TCanvas;
    pave->Draw();
    PrintCanvas(can);
    delete can;
  }
  delete paves;
}


Int_t Painter::DrawPull() {
//  TCanvas* can = new TCanvas;
//
//  TH1F* h = fOutput->GetPull();
//  TH1F* hRef = NULL;
//  if(fOutputRef) hRef = fOutputRef->GetPull();
//
//  h->SetTitle("Pull histogram");
//  h->SetLineColor(kRed);
//
//  h->Draw();
//
//  if(hRef) {
//    hRef->SetLineColor(kBlue);
//    hRef->Draw("same");
//  }
//
//  PrintCanvas(can);
//  delete can;
//
}

Int_t Painter::DrawPDF(Int_t ival) {
  TObjArray* TrashBin = new TObjArray; TrashBin->SetOwner();

  TPaveLabel* Q2Label = new TPaveLabel(0.7, 0.08, 1.0, 0.16, "","NDC"); TrashBin->AddLast(Q2Label);
  Q2Label->SetBorderSize(1); Q2Label->SetFillColor(kWhite); Q2Label->SetBorderSize(0);

  
  // Q2 value of the bin:
  Double_t fQ2val;
  
  if ( ! (fDrawExp||fDrawModel||fDrawParam) ) {
    fQ2val = fOutput->GetQ2Value(ival);
  } else {  
    if (fDrawBase) {
     fQ2val = fOutputBase->GetQ2Value(ival);
   } else if (fDrawExp) {
     fQ2val = fOutputExp->GetQ2Value(ival);
   } else if (fDrawModel) {
     fQ2val = fOutputModel->GetQ2Value(ival);
   } else if (fDrawParam) {
     fQ2val = fOutputParam->GetQ2Value(ival);
   }  
  }
  
  // Label:
  char label[32]; 
  sprintf (label,"Q^{2} = %12.2f GeV^{2}",fQ2val);
  Q2Label->SetLabel(label);
  
  TCanvas* can = new TCanvas("can","can",600, 400);
  can->Divide(3,2);
  
  if ( ! (fDrawExp||fDrawModel||fDrawParam) ) {
    
    PlotPdfSub(can->cd(1), fOutput, fOutputRef, ival, "xU, xu_{V}",
	     Output::kU,     Output::kUv, can->cd(6), TrashBin);
    PlotPdfSub(can->cd(2), fOutput, fOutputRef, ival, "xD, xd_{V}",
	     Output::kD,     Output::kDv, NULL,  TrashBin);
    PlotPdfSub(can->cd(3), fOutput, fOutputRef, ival, "xg",
	     Output::kGluon, Output::kNULL,NULL,  TrashBin);
    PlotPdfSub(can->cd(4), fOutput, fOutputRef, ival, "x#bar{U}",
	     Output::kUb,   Output::kNULL, NULL, TrashBin);
    PlotPdfSub(can->cd(5), fOutput, fOutputRef, ival, "x#bar{D}",
	     Output::kDb,    Output::kNULL, NULL, TrashBin);
	     
  } else {
    
    PlotPdfSubEMP(can->cd(1), ival, "xU, xu_{V}", Output::kU, Output::kUv, can->cd(6), TrashBin);
    PlotPdfSubEMP(can->cd(2), ival, "xD, xd_{V}", Output::kD, Output::kDv, NULL, TrashBin);
    PlotPdfSubEMP(can->cd(3), ival, "xg", Output::kGluon, Output::kNULL, NULL, TrashBin);
    PlotPdfSubEMP(can->cd(4), ival, "x#bar{U}", Output::kUb, Output::kNULL, NULL, TrashBin);
    PlotPdfSubEMP(can->cd(5), ival, "x#bar{D}", Output::kDb, Output::kNULL, NULL, TrashBin);
    
  }
	     
  can->cd();
  Q2Label->Draw();
  PrintCanvas(can);
  delete can;
//  delete TrashBin;           // I am sorry, but this delete makes "invalid pointer: 0x089a8850" error for DrawResults --bands <1 option>
  return 0;
}


Int_t Painter::PlotPdfSub(TVirtualPad* pad, Output* FitterOut, Output* FitterRef, 
				  Int_t Q2Bin, const Char_t* Title, Output::pdf pdf1, Output::pdf pdf2,
				  TVirtualPad* legend, TObjArray* TrashBin) {
  TGraphAsymmErrors* graph =   FitterOut->GetPdf(pdf1, Q2Bin);
  TGraphAsymmErrors* graphR = NULL; 
  TGraphAsymmErrors* graph2 = NULL;
  TGraphAsymmErrors* graphR2 = NULL;
  Double_t RatioSize = 0.;

  if(FitterRef) {
    graphR = FitterRef->GetPdf(pdf1, Q2Bin);
    RatioSize = 0.3;
  }
  if(pdf2 != Output::kNULL) {
    graph2 =   FitterOut->GetPdf(pdf2, Q2Bin);
    if(FitterRef) graphR2 =  FitterRef->GetPdf(pdf2, Q2Bin);
  }

  graph->SetLineColor(fColor);
  graph->SetFillColor(fColor);
  graph->SetFillStyle(fFillStyle);

  if(graph2) {
    graph2->SetLineColor(fColor);
    graph2->SetFillColor(fColor);
    graph2->SetFillStyle(fFillStyle);
  }
  if(graphR) {
    graphR->SetLineColor(fColorRef);
    graphR->SetFillColor(fColorRef);
    graphR->SetFillStyle(fFillStyleRef);
  }
  if(graphR2) {
    graphR2->SetLineColor(fColorRef);
    graphR2->SetFillColor(fColorRef);
    graphR2->SetFillStyle(fFillStyleRef);
  }

  pad->cd();
  TPad* pad1 = new TPad("pad1","pad1", 0., RatioSize, 1., 1.);   TrashBin->AddLast(pad1);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogx();

  graph->SetTitle("");
  graph->GetYaxis()->SetTitle("xP(x)");
  graph->GetYaxis()->SetTitleOffset(1.5);  
  graph->GetYaxis()->SetLabelSize(0.05); 
  graph->GetYaxis()->SetNdivisions(506); 
  graph->GetYaxis()->SetLabelOffset(0.02);
  graph->SetLineColor(fColor);
  graph->GetXaxis()->SetTitle("x");
  graph->GetXaxis()->SetTitleOffset(0.5);  
  graph->GetXaxis()->SetTitleSize(.06);  
  graph->GetXaxis()->SetLabelSize(0.04); 
  graph->GetXaxis()->SetLabelOffset(0.015);


  TGraphAsymmErrors* graphRL = NULL;
  TGraphAsymmErrors* graphR2L = NULL;

  graph->Draw("AC3");
  if(graphR) {
    graphR->DrawClone("3 same");
    } 
  if(graph2)  graph2->Draw("L3 same");
  if(graphR2) {
    graphR2->Draw("3 same");
  }

  TPaveLabel* box = new TPaveLabel(0.0, 0.75, 0.09, 0.91, "xP(x)", "NDC"); TrashBin->AddLast(box);
  box->SetFillColor(kWhite); box->SetBorderSize(0); box->SetTextAngle(90.); box->SetTextSize(0.35);
  box->Draw();

  // Plotting the ratio pad if needed
  if(graphR) {
    Int_t Color_Ratio = kBlack;
    Int_t Style_Ratio = 2;

    graph->GetYaxis()->SetTitle("xP(x)");
    graph->GetYaxis()->SetTitleOffset(1.);  
    graph->GetYaxis()->SetLabelSize(0.04); 
    graph->GetYaxis()->SetNdivisions(506); 
    graph->GetYaxis()->SetLabelOffset(0.02);
 
    TGraph* ref_ratio = (TGraph*) (graph->Clone()); TrashBin->AddLast(ref_ratio);
    TGraphAsymmErrors* graph_ratio = (TGraphAsymmErrors*) (graph->Clone()); TrashBin->AddLast(graph_ratio);
    TGraphAsymmErrors* graphR_ratio = (TGraphAsymmErrors*) (graphR->Clone()); TrashBin->AddLast(graphR_ratio);
    for(Int_t i=0; i<graph->GetN(); i++) {
      ref_ratio->SetPoint(i, graph->GetX()[i], graph->GetY()[i] / graphR->GetY()[i]);

      graph_ratio->SetPoint(i, graph->GetX()[i], 1.);
      if(graph->GetY()[i] > 0.) 
        graph_ratio->SetPointError(i, 0., 0., graph->GetEYlow()[i] / graph->GetY()[i], graph->GetEYhigh()[i] / graph->GetY()[i]);

      graphR_ratio->SetPoint(i, graphR->GetX()[i], 1.);
      if(graphR->GetY()[i] > 0.)
        graphR_ratio->SetPointError(i, 0., 0., graphR->GetEYlow()[i] / graphR->GetY()[i], graphR->GetEYhigh()[i] / graphR->GetY()[i]);
    }

    pad->cd();
    TPad* pad2 = new TPad("pad2","pad2", 0., 0., 1., RatioSize);   TrashBin->AddLast(pad2);
    pad2->Draw();

    pad1->SetBottomMargin(0.);
    pad2->SetTopMargin(0.);
    pad2->SetBottomMargin(0.2);

    pad2->cd();
    pad2->SetLogx();

    ref_ratio->SetLineColor(Color_Ratio);  
    ref_ratio->SetLineStyle(Style_Ratio);  
    ref_ratio->SetFillColor(Color_Ratio);
    ref_ratio->GetYaxis()->SetTitle("ratio                ");
    ref_ratio->GetYaxis()->SetTitleSize(0.13);
    ref_ratio->GetYaxis()->SetTitleOffset(0.35);
    ref_ratio->GetYaxis()->SetNdivisions(505);
    ref_ratio->GetYaxis()->SetLabelSize(0.11);
    ref_ratio->GetXaxis()->SetLabelSize(0.11);
    ref_ratio->GetXaxis()->SetLabelOffset(0.03);
    ref_ratio->GetXaxis()->SetTitle("x  ");
    ref_ratio->GetXaxis()->SetTitleOffset(0.5);  
    ref_ratio->GetXaxis()->SetTitleSize(.15);  
    
    ref_ratio->SetMaximum(1.19);
    ref_ratio->SetMinimum(0.81);
    ref_ratio->Draw("ALX");

    TGraphAsymmErrors* graphRL_ratio = (TGraphAsymmErrors*) graphR_ratio->Clone(); TrashBin->AddLast(graphRL_ratio);
    graphRL_ratio->SetFillStyle(0);

    graph_ratio->Draw("3");
    graphR_ratio->Draw("3");
    graphRL_ratio->Draw("3");

    ref_ratio->Draw("LX");
  }


  TPaveLabel* label = new TPaveLabel(0.49, 0.85, 0.51, 0.87, Title,"NDC"); TrashBin->AddLast(label);
  label->SetBorderSize(0); label->SetFillColor(kWhite); label->SetTextSize(3.0);
  pad->cd();
  label->Draw();

  if(legend) {
    legend->cd();
    TPaveLabel* lab1 = new TPaveLabel(0., 0.4, 1.0, 0.5, FitterOut->GetName()->Data(), "NDC");
    TrashBin->AddLast(lab1); lab1->SetFillColor(kWhite); lab1->SetBorderSize(0); lab1->SetTextColor(fColor);
    lab1->Draw();
    if(FitterRef) {
      TPaveLabel* lab2 = new TPaveLabel(0., 0.6, 1.0, 0.7, FitterRef->GetName()->Data(), "NDC");
      TrashBin->AddLast(lab2); lab2->SetFillColor(kWhite); lab2->SetBorderSize(0); lab2->SetTextColor(fColorRef);
      lab2->Draw();
    }
  }
  return 0;
}

Int_t Painter::PlotPdfSubEMP(TVirtualPad* pad, Int_t Q2Bin, const Char_t* Title, Output::pdf pdf1, Output::pdf pdf2,
				  TVirtualPad* legend, TObjArray* TrashBin) {
  
  TGraphAsymmErrors* graphBase = NULL;
  TGraphAsymmErrors* graphExp = NULL;
  TGraphAsymmErrors* graphModel = NULL;
  TGraphAsymmErrors* graphParam = NULL;
  TGraphAsymmErrors* graphBase2 = NULL;
  TGraphAsymmErrors* graphExp2 = NULL;
  TGraphAsymmErrors* graphModel2 = NULL;
  TGraphAsymmErrors* graphParam2 = NULL;
  
  TGraphAsymmErrors* graphModelsumm = NULL;
  TGraphAsymmErrors* graphParamsumm = NULL;
  TGraphAsymmErrors* graphModelsumm2 = NULL;
  TGraphAsymmErrors* graphParamsumm2 = NULL;
  
  TGraph* graph_lineUp = NULL;
  TGraph* graph_lineDown = NULL;
  TGraph* graphRatio_lineUp = NULL;
  TGraph* graphRatio_lineDown = NULL;
  TGraph* graph2_lineUp = NULL;
  TGraph* graph2_lineDown = NULL;
  
  
  Double_t RatioSize = 0.;

  bool isSmthDrawn = false;
  
  if (fDrawBase) {
    graphBase = fOutputBase->GetPdf(pdf1, Q2Bin);
    graphBase->SetLineColor(fColorBase);
    graphBase->SetFillColor(fColorBase);
    graphBase->SetFillStyle(fFillStyleBase);
  }
  if (fDrawExp) {
    graphExp   = fOutputExp->GetPdf(pdf1, Q2Bin);
    RatioSize = 0.3;
    graphExp->SetLineColor(fColorExp);
    graphExp->SetFillColor(fColorExp);
    graphExp->SetFillStyle(fFillStyleExp);
  }
  if (fDrawModel) {
    graphModel = fOutputModel->GetPdf(pdf1, Q2Bin);
    RatioSize = 0.3;
    graphModelsumm = fOutputModel->GetPdf(pdf1, Q2Bin);
    graphModelsumm->SetLineColor(fColorModel);
    graphModelsumm->SetFillColor(fColorModel);
    graphModelsumm->SetFillStyle(fFillStyleModel);
  }
  if (fDrawParam) {
    graphParam = fOutputParam->GetPdf(pdf1, Q2Bin);
    RatioSize = 0.3;
    graphParamsumm = fOutputParam->GetPdf(pdf1, Q2Bin);
    graphParamsumm->SetLineColor(fColorParam);
    graphParamsumm->SetFillColor(fColorParam);
    graphParamsumm->SetFillStyle(fFillStyleParam);
  }
  
  if (fDrawParam) {
    graph_lineUp = new TGraph(graphParamsumm->GetN());
    graph_lineDown = new TGraph(graphParamsumm->GetN());
  } else if (fDrawModel) {
    graph_lineUp = new TGraph(graphModelsumm->GetN());
    graph_lineDown = new TGraph(graphModelsumm->GetN());
  } else if (fDrawExp) {
    graph_lineUp = new TGraph(graphExp->GetN());
    graph_lineDown = new TGraph(graphExp->GetN());
  }
  
  if(pdf2 != Output::kNULL) {
    if (fDrawBase) {
      graphBase2 = fOutputBase->GetPdf(pdf2, Q2Bin);
      graphBase2->SetLineColor(fColorBase);
      graphBase2->SetFillColor(fColorBase);
      graphBase2->SetFillStyle(fFillStyleBase);
    }
    if (fDrawExp) {
      graphExp2 = fOutputExp->GetPdf(pdf2, Q2Bin);
      graphExp2->SetLineColor(fColorExp);
      graphExp2->SetFillColor(fColorExp);
      graphExp2->SetFillStyle(fFillStyleExp);
    }
    if (fDrawModel) {
      graphModel2 = fOutputModel->GetPdf(pdf2, Q2Bin);
      graphModelsumm2 = fOutputModel->GetPdf(pdf2, Q2Bin);
      graphModelsumm2->SetLineColor(fColorModel);
      graphModelsumm2->SetFillColor(fColorModel);
      graphModelsumm2->SetFillStyle(fFillStyleModel);
    }
    if (fDrawParam) {
      graphParam2 = fOutputParam->GetPdf(pdf2, Q2Bin);
      graphParamsumm2 = fOutputParam->GetPdf(pdf2, Q2Bin);
      graphParamsumm2->SetLineColor(fColorParam);
      graphParamsumm2->SetFillColor(fColorParam);
      graphParamsumm2->SetFillStyle(fFillStyleParam);
    }
    
    if (fDrawParam) {
      graph2_lineUp = new TGraph(graphParamsumm2->GetN());
      graph2_lineDown = new TGraph(graphParamsumm2->GetN());
    } else if (fDrawModel) {
      graph2_lineUp = new TGraph(graphModelsumm2->GetN());
      graph2_lineDown = new TGraph(graphModelsumm2->GetN());
    } else if (fDrawExp) {
      graph2_lineUp = new TGraph(graphExp2->GetN());
      graph2_lineDown = new TGraph(graphExp2->GetN());
    }
  
  }
  
  if (graph_lineUp && graph_lineDown) {
    graph_lineUp->SetLineColor(kBlack);
    graph_lineDown->SetLineColor(kBlack);
    graph_lineUp->SetLineWidth(0.1);
    graph_lineDown->SetLineWidth(0.1);
  }
  
  if (graph2_lineUp && graph2_lineDown) {
    graph2_lineUp->SetLineColor(kBlack);
    graph2_lineDown->SetLineColor(kBlack);
    graph2_lineUp->SetLineWidth(0.1);
    graph2_lineDown->SetLineWidth(0.1);
  }

  pad->cd();
  TPad* pad1 = new TPad("pad1","pad1", 0., RatioSize, 1., 1.);   TrashBin->AddLast(pad1);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogx();
  
  
  if (fDrawBase) {
    graphBase->SetTitle("");
    graphBase->GetYaxis()->SetTitle("xP(x)");
    graphBase->GetYaxis()->SetTitleOffset(1.5);  
    graphBase->GetYaxis()->SetLabelSize(0.05); 
    graphBase->GetYaxis()->SetNdivisions(506); 
    graphBase->GetYaxis()->SetLabelOffset(0.02);
    graphBase->GetXaxis()->SetTitle("x");
    graphBase->GetXaxis()->SetTitleOffset(0.5);  
    graphBase->GetXaxis()->SetTitleSize(.06);  
    graphBase->GetXaxis()->SetLabelSize(0.04); 
    graphBase->GetXaxis()->SetLabelOffset(0.015);
  }
  if (fDrawExp) {
    graphExp->SetTitle("");
    graphExp->GetYaxis()->SetTitle("xP(x)");
    graphExp->GetYaxis()->SetTitleOffset(1.5);  
    graphExp->GetYaxis()->SetLabelSize(0.05); 
    graphExp->GetYaxis()->SetNdivisions(506); 
    graphExp->GetYaxis()->SetLabelOffset(0.02);
    graphExp->GetXaxis()->SetTitle("x");
    graphExp->GetXaxis()->SetTitleOffset(0.5);  
    graphExp->GetXaxis()->SetTitleSize(.06);  
    graphExp->GetXaxis()->SetLabelSize(0.04); 
    graphExp->GetXaxis()->SetLabelOffset(0.015);
  }
  if (fDrawModel) {
    graphModelsumm->SetTitle("");
    graphModelsumm->GetYaxis()->SetTitle("xP(x)");
    graphModelsumm->GetYaxis()->SetTitleOffset(1.5);  
    graphModelsumm->GetYaxis()->SetLabelSize(0.05); 
    graphModelsumm->GetYaxis()->SetNdivisions(506); 
    graphModelsumm->GetYaxis()->SetLabelOffset(0.02);
    graphModelsumm->GetXaxis()->SetTitle("x");
    graphModelsumm->GetXaxis()->SetTitleOffset(0.5);  
    graphModelsumm->GetXaxis()->SetTitleSize(.06);  
    graphModelsumm->GetXaxis()->SetLabelSize(0.04); 
    graphModelsumm->GetXaxis()->SetLabelOffset(0.015);
  }
  if (fDrawParam) {
    graphParamsumm->SetTitle("");
    graphParamsumm->GetYaxis()->SetTitle("xP(x)");
    graphParamsumm->GetYaxis()->SetTitleOffset(1.5);  
    graphParamsumm->GetYaxis()->SetLabelSize(0.05); 
    graphParamsumm->GetYaxis()->SetNdivisions(506); 
    graphParamsumm->GetYaxis()->SetLabelOffset(0.02);
    graphParamsumm->GetXaxis()->SetTitle("x");
    graphParamsumm->GetXaxis()->SetTitleOffset(0.5);  
    graphParamsumm->GetXaxis()->SetTitleSize(.06);  
    graphParamsumm->GetXaxis()->SetLabelSize(0.04); 
    graphParamsumm->GetXaxis()->SetLabelOffset(0.015);
  }
  
//====================  summ of exp+model, exp+model+param  ====================
  
  int Npoints = 0;
  Double_t xa, xb, xc, ya, yb, yc;
  Double_t xela, xeha, yela, yeha, xelb, xehb, yelb, yehb, xels, xehs, yels, yehs;
  if (fDrawModel || fDrawParam) {
    if (fDrawModel) {
      Npoints = graphModelsumm->GetN();
    } else {
      Npoints = graphParamsumm->GetN();
    }
    for (int ibin=0; ibin<Npoints; ibin++) {
      if (fDrawModel) {
        if (fDrawExp) {
	  graphModelsumm->GetPoint(ibin, xa, ya);
	  graphExp->GetPoint(ibin, xb, yb);
	  if (xa == xb) {
	    xela = graphModelsumm->GetErrorXlow(ibin);
	    xeha = graphModelsumm->GetErrorXhigh(ibin);
	    yela = graphModelsumm->GetErrorYlow(ibin);
	    yeha = graphModelsumm->GetErrorYhigh(ibin);
	    
	    xelb = graphExp->GetErrorXlow(ibin);
	    xehb = graphExp->GetErrorXhigh(ibin);
	    yelb = graphExp->GetErrorYlow(ibin);
	    yehb = graphExp->GetErrorYhigh(ibin);
	    
	    xels = TMath::Sqrt(xela*xela + xelb*xelb);
	    xehs = TMath::Sqrt(xeha*xeha + xehb*xehb);
	    yels = TMath::Sqrt(yela*yela + yelb*yelb);
	    yehs = TMath::Sqrt(yeha*yeha + yehb*yehb);
	    
	    graphModelsumm->SetPointError(ibin, xels, xehs, yels, yehs);
	  } else {
	    cout<<"WARNING: The x values for the Exp and Model points are different !"<<endl;
	  }
	}
      }
      if (fDrawParam) {
        if (fDrawModel) {
	  graphParamsumm->GetPoint(ibin, xa, ya);
	  graphModelsumm->GetPoint(ibin, xb, yb);
	  if (xa == xb) {
	    xela = graphParamsumm->GetErrorXlow(ibin);
	    xeha = graphParamsumm->GetErrorXhigh(ibin);
	    yela = graphParamsumm->GetErrorYlow(ibin);
	    yeha = graphParamsumm->GetErrorYhigh(ibin);
	    
	    xelb = graphModelsumm->GetErrorXlow(ibin);
	    xehb = graphModelsumm->GetErrorXhigh(ibin);
	    yelb = graphModelsumm->GetErrorYlow(ibin);
	    yehb = graphModelsumm->GetErrorYhigh(ibin);
	    
	    xels = TMath::Sqrt(xela*xela + xelb*xelb);
	    xehs = TMath::Sqrt(xeha*xeha + xehb*xehb);
	    yels = TMath::Sqrt(yela*yela + yelb*yelb);
	    yehs = TMath::Sqrt(yeha*yeha + yehb*yehb);
	    
	    graphParamsumm->SetPointError(ibin, xels, xehs, yels, yehs);
	  } else {
	    cout<<"WARNING: The x values for the Model and Param points are different !"<<endl;
	  }
	} else if (fDrawExp) {
	  graphParamsumm->GetPoint(ibin, xa, ya);
	  graphExp->GetPoint(ibin, xb, yb);
	  if (xa == xb) {
	    xela = graphParamsumm->GetErrorXlow(ibin);
	    xeha = graphParamsumm->GetErrorXhigh(ibin);
	    yela = graphParamsumm->GetErrorYlow(ibin);
	    yeha = graphParamsumm->GetErrorYhigh(ibin);
	    
	    xelb = graphExp->GetErrorXlow(ibin);
	    xehb = graphExp->GetErrorXhigh(ibin);
	    yelb = graphExp->GetErrorYlow(ibin);
	    yehb = graphExp->GetErrorYhigh(ibin);
	    
	    xels = TMath::Sqrt(xela*xela + xelb*xelb);
	    xehs = TMath::Sqrt(xeha*xeha + xehb*xehb);
	    yels = TMath::Sqrt(yela*yela + yelb*yelb);
	    yehs = TMath::Sqrt(yeha*yeha + yehb*yehb);
	    
	    graphParamsumm->SetPointError(ibin, xels, xehs, yels, yehs);
	  } else {
	    cout<<"WARNING: The x values for the Exp and Param points are different !"<<endl;
	  }
	}
      }
      
        if(pdf2 != Output::kNULL) {
          if (fDrawModel) {
          if (fDrawExp) {
	    graphModelsumm2->GetPoint(ibin, xa, ya);
	    graphExp2->GetPoint(ibin, xb, yb);
	    if (xa == xb) {
	      xela = graphModelsumm2->GetErrorXlow(ibin);
	      xeha = graphModelsumm2->GetErrorXhigh(ibin);
	      yela = graphModelsumm2->GetErrorYlow(ibin);
	      yeha = graphModelsumm2->GetErrorYhigh(ibin);
	    
	      xelb = graphExp2->GetErrorXlow(ibin);
	      xehb = graphExp2->GetErrorXhigh(ibin);
	      yelb = graphExp2->GetErrorYlow(ibin);
	      yehb = graphExp2->GetErrorYhigh(ibin);
	    
	      xels = TMath::Sqrt(xela*xela + xelb*xelb);
	      xehs = TMath::Sqrt(xeha*xeha + xehb*xehb);
	      yels = TMath::Sqrt(yela*yela + yelb*yelb);
	      yehs = TMath::Sqrt(yeha*yeha + yehb*yehb);
	    
	      graphModelsumm2->SetPointError(ibin, xels, xehs, yels, yehs);
	    } else {
	      cout<<"WARNING: The x values for the Exp2 and Model2 points are different !"<<endl;
	    }
	  }
        }
        if (fDrawParam) {
          if (fDrawModel) {
	    graphParamsumm2->GetPoint(ibin, xa, ya);
	    graphModelsumm2->GetPoint(ibin, xb, yb);
	    if (xa == xb) {
	      xela = graphParamsumm2->GetErrorXlow(ibin);
	      xeha = graphParamsumm2->GetErrorXhigh(ibin);
	      yela = graphParamsumm2->GetErrorYlow(ibin);
	      yeha = graphParamsumm2->GetErrorYhigh(ibin);
	    
	      xelb = graphModelsumm2->GetErrorXlow(ibin);
	      xehb = graphModelsumm2->GetErrorXhigh(ibin);
	      yelb = graphModelsumm2->GetErrorYlow(ibin);
	      yehb = graphModelsumm2->GetErrorYhigh(ibin);
	    
	      xels = TMath::Sqrt(xela*xela + xelb*xelb);
	      xehs = TMath::Sqrt(xeha*xeha + xehb*xehb);
	      yels = TMath::Sqrt(yela*yela + yelb*yelb);
	      yehs = TMath::Sqrt(yeha*yeha + yehb*yehb);
	    
	      graphParamsumm2->SetPointError(ibin, xels, xehs, yels, yehs);
	    } else {
	      cout<<"WARNING: The x values for the Model2 and Param2 points are different !"<<endl;
	    }
	  } else if (fDrawExp) {
	    graphParamsumm2->GetPoint(ibin, xa, ya);
	    graphExp2->GetPoint(ibin, xb, yb);
	    if (xa == xb) {
	      xela = graphParamsumm2->GetErrorXlow(ibin);
	      xeha = graphParamsumm2->GetErrorXhigh(ibin);
	      yela = graphParamsumm2->GetErrorYlow(ibin);
	      yeha = graphParamsumm2->GetErrorYhigh(ibin);
	    
	      xelb = graphExp2->GetErrorXlow(ibin);
	      xehb = graphExp2->GetErrorXhigh(ibin);
	      yelb = graphExp2->GetErrorYlow(ibin);
	      yehb = graphExp2->GetErrorYhigh(ibin);
	    
	      xels = TMath::Sqrt(xela*xela + xelb*xelb);
	      xehs = TMath::Sqrt(xeha*xeha + xehb*xehb);
	      yels = TMath::Sqrt(yela*yela + yelb*yelb);
	      yehs = TMath::Sqrt(yeha*yeha + yehb*yehb);
	    
	      graphParamsumm2->SetPointError(ibin, xels, xehs, yels, yehs);
	    } else {
	      cout<<"WARNING: The x values for the Exp2 and Param2 points are different !"<<endl;
	    }
	  }
        }
      }
    }
  
  }
//====================  end of exp+model, exp+model+param  ====================

//====================  black lines for bands  ====================

  if (fDrawParam) {
    for (int ibin=0; ibin<graphParamsumm->GetN(); ibin++) {
      graphParamsumm->GetPoint(ibin, xa, ya);
      yela = graphParamsumm->GetErrorYlow(ibin);
      yeha = graphParamsumm->GetErrorYhigh(ibin);
      
      graph_lineUp->SetPoint(ibin, xa, (ya+yeha));
      graph_lineDown->SetPoint(ibin, xa, (ya-yela));
    }
  } else if (fDrawModel) {
    for (int ibin=0; ibin<graphModelsumm->GetN(); ibin++) {
      graphModelsumm->GetPoint(ibin, xa, ya);
      yela = graphModelsumm->GetErrorYlow(ibin);
      yeha = graphModelsumm->GetErrorYhigh(ibin);
      
      graph_lineUp->SetPoint(ibin, xa, (ya+yeha));
      graph_lineDown->SetPoint(ibin, xa, (ya-yela));
    }
  } else if (fDrawExp) {
    for (int ibin=0; ibin<graphExp->GetN(); ibin++) {
      graphExp->GetPoint(ibin, xa, ya);
      yela = graphExp->GetErrorYlow(ibin);
      yeha = graphExp->GetErrorYhigh(ibin);
      
      graph_lineUp->SetPoint(ibin, xa, (ya+yeha));
      graph_lineDown->SetPoint(ibin, xa, (ya-yela));
    }
  }
  
  if(pdf2 != Output::kNULL) {
    if (fDrawParam) {
      for (int ibin=0; ibin<graphParamsumm2->GetN(); ibin++) {
        graphParamsumm2->GetPoint(ibin, xa, ya);
        yela = graphParamsumm2->GetErrorYlow(ibin);
        yeha = graphParamsumm2->GetErrorYhigh(ibin);
      
        graph2_lineUp->SetPoint(ibin, xa, (ya+yeha));
        graph2_lineDown->SetPoint(ibin, xa, (ya-yela));
      }
    } else if (fDrawModel) {
      for (int ibin=0; ibin<graphModelsumm2->GetN(); ibin++) {
        graphModelsumm2->GetPoint(ibin, xa, ya);
        yela = graphModelsumm2->GetErrorYlow(ibin);
        yeha = graphModelsumm2->GetErrorYhigh(ibin);
      
        graph2_lineUp->SetPoint(ibin, xa, (ya+yeha));
        graph2_lineDown->SetPoint(ibin, xa, (ya-yela));
      }
    } else if (fDrawExp) {
      for (int ibin=0; ibin<graphExp2->GetN(); ibin++) {
        graphExp2->GetPoint(ibin, xa, ya);
        yela = graphExp2->GetErrorYlow(ibin);
        yeha = graphExp2->GetErrorYhigh(ibin);
      
        graph2_lineUp->SetPoint(ibin, xa, (ya+yeha));
        graph2_lineDown->SetPoint(ibin, xa, (ya-yela));
      }
    }
  }

//=================  end of black lines for bands  =================
  
  isSmthDrawn = false;
  if (fDrawParam) {
    if (isSmthDrawn) {
      graphParamsumm->Draw("3 same");
    } else {
      graphParamsumm->Draw("AC3");
      isSmthDrawn = true;
    }
  }
  if (fDrawModel) {
    if (isSmthDrawn) {
      graphModelsumm->Draw("3 same");
    } else {
      graphModelsumm->Draw("AC3");
      isSmthDrawn = true;
    }
  }
  if (fDrawExp) {
    if (isSmthDrawn) {
      graphExp->Draw("3 same");
    } else {
      graphExp->Draw("AC3");
      isSmthDrawn = true;
    }
  }
  if (fDrawBase) {
    if (isSmthDrawn) {
      graphBase->Draw("3 same");
    } else {
      graphBase->Draw("3");
      isSmthDrawn = true;
    }
  }
  
  if (graph_lineUp && graph_lineDown) {
    graph_lineUp->Draw("C same");
    graph_lineDown->Draw("C same");
  }
  
  if(pdf2 != Output::kNULL) {
    if (isSmthDrawn) {
      if (fDrawParam) graphParamsumm2->Draw("L3 same");
      if (fDrawModel) graphModelsumm2->Draw("L3 same");
      if (fDrawExp)   graphExp2->Draw("L3 same");
      if (fDrawBase)  graphBase2->Draw("3 same");
      
      if (graph2_lineUp && graph2_lineDown) {
        graph2_lineUp->Draw("C same");
        graph2_lineDown->Draw("C same");
      }
    }
  }

  TPaveLabel* box = new TPaveLabel(0.0, 0.75, 0.09, 0.91, "xP(x)", "NDC"); TrashBin->AddLast(box);
  box->SetFillColor(kWhite); box->SetBorderSize(0); box->SetTextAngle(90.); box->SetTextSize(0.35);
  box->Draw();

  
  // Plotting the ratio pad if needed
  
  bool doDrawRatio = false;
  int nGraphs = 0;
  if (fDrawBase)  nGraphs++;
  if (fDrawExp)   nGraphs++;
  if (fDrawModel) nGraphs++;
  if (fDrawParam) nGraphs++;
  
  if (nGraphs > 1) doDrawRatio = true;
  
  TGraphAsymmErrors* graphMainRatio = NULL;
  
  if (doDrawRatio) {
    if (fDrawBase) {
      graphMainRatio = (TGraphAsymmErrors*) (graphBase->Clone());
    } else if (fDrawExp) {
      graphMainRatio = (TGraphAsymmErrors*) (graphExp->Clone());
    } else if (fDrawModel) {
      graphMainRatio = (TGraphAsymmErrors*) (graphModelsumm->Clone());
    } else {
     cout << "Cannot draw grahs ratio correctly - no main graph." << endl;
     doDrawRatio = false;
    }
    
    if (fDrawParam) {
      graphRatio_lineUp = new TGraph(graphParamsumm->GetN());
      graphRatio_lineDown = new TGraph(graphParamsumm->GetN());
    } else if (fDrawModel) {
      graphRatio_lineUp = new TGraph(graphModelsumm->GetN());
      graphRatio_lineDown = new TGraph(graphModelsumm->GetN());
    } else if (fDrawExp) {
      graphRatio_lineUp = new TGraph(graphExp->GetN());
      graphRatio_lineDown = new TGraph(graphExp->GetN());
    }
    
    if (fDrawParam) {
      graphRatio_lineUp->SetLineColor(kBlack);
      graphRatio_lineDown->SetLineColor(kBlack);
      graphRatio_lineUp->SetLineWidth(0.1);
      graphRatio_lineDown->SetLineWidth(0.1);
    }
  }
  
  if (doDrawRatio) {
    
    TGraph* Base_ratio;
    TGraph* Exp_ratio;
    TGraph* Model_ratio;
    TGraph* Param_ratio;
    TGraphAsymmErrors* graphExp_ratio;
    TGraphAsymmErrors* graphModel_ratio;
    TGraphAsymmErrors* graphParam_ratio;
    
    Int_t Style_Ratio = 2;
    
    if (fDrawBase) {
      graphBase->GetYaxis()->SetTitle("xP(x)");
      graphBase->GetYaxis()->SetTitleOffset(1.);  
      graphBase->GetYaxis()->SetLabelSize(0.04); 
      graphBase->GetYaxis()->SetNdivisions(506); 
      graphBase->GetYaxis()->SetLabelOffset(0.02);
      
      Base_ratio = (TGraph*) (graphBase->Clone());
      
      for(Int_t i=0; i<graphBase->GetN(); i++) {
      Base_ratio->SetPoint(i, graphBase->GetX()[i], graphBase->GetY()[i] / graphMainRatio->GetY()[i]);
      }
    }
    if (fDrawExp) {
      graphExp->GetYaxis()->SetTitle("xP(x)");
      graphExp->GetYaxis()->SetTitleOffset(1.);  
      graphExp->GetYaxis()->SetLabelSize(0.04); 
      graphExp->GetYaxis()->SetNdivisions(506); 
      graphExp->GetYaxis()->SetLabelOffset(0.02);
      
      Exp_ratio = (TGraph*) (graphExp->Clone());
      graphExp_ratio = (TGraphAsymmErrors*) (graphExp->Clone());
      
      for(Int_t i=0; i<graphExp->GetN(); i++) {
      Exp_ratio->SetPoint(i, graphExp->GetX()[i], graphExp->GetY()[i] / graphMainRatio->GetY()[i]);
      
      graphExp_ratio->SetPoint(i, graphExp->GetX()[i], graphExp->GetY()[i] / graphMainRatio->GetY()[i]);
      if(graphExp->GetY()[i] > 0.) 
        graphExp_ratio->SetPointError(i, 0., 0., graphExp->GetEYlow()[i] / graphMainRatio->GetY()[i], graphExp->GetEYhigh()[i] / graphMainRatio->GetY()[i]);
      }
    }
    if (fDrawModel) {
      graphModelsumm->GetYaxis()->SetTitle("xP(x)");
      graphModelsumm->GetYaxis()->SetTitleOffset(1.);  
      graphModelsumm->GetYaxis()->SetLabelSize(0.04); 
      graphModelsumm->GetYaxis()->SetNdivisions(506); 
      graphModelsumm->GetYaxis()->SetLabelOffset(0.02);
      
      Model_ratio = (TGraph*) (graphModelsumm->Clone());
      graphModel_ratio = (TGraphAsymmErrors*) (graphModelsumm->Clone());
      
      for(Int_t i=0; i<graphModelsumm->GetN(); i++) {
      Model_ratio->SetPoint(i, graphModelsumm->GetX()[i], graphModelsumm->GetY()[i] / graphMainRatio->GetY()[i]);
      
      graphModel_ratio->SetPoint(i, graphModelsumm->GetX()[i], graphModelsumm->GetY()[i] / graphMainRatio->GetY()[i]);
      if(graphModelsumm->GetY()[i] > 0.) 
        graphModel_ratio->SetPointError(i, 0., 0., graphModelsumm->GetEYlow()[i] / graphMainRatio->GetY()[i], graphModelsumm->GetEYhigh()[i] / graphMainRatio->GetY()[i]);
      }
    }
    if (fDrawParam) {
      graphParamsumm->GetYaxis()->SetTitle("xP(x)");
      graphParamsumm->GetYaxis()->SetTitleOffset(1.);  
      graphParamsumm->GetYaxis()->SetLabelSize(0.04); 
      graphParamsumm->GetYaxis()->SetNdivisions(506); 
      graphParamsumm->GetYaxis()->SetLabelOffset(0.02);
      
      Param_ratio = (TGraph*) (graphParam->Clone());
      graphParam_ratio = (TGraphAsymmErrors*) (graphParamsumm->Clone());
      
      for(Int_t i=0; i<graphParamsumm->GetN(); i++) {
      Param_ratio->SetPoint(i, graphParamsumm->GetX()[i], graphParamsumm->GetY()[i] / graphMainRatio->GetY()[i]);
      
      graphParam_ratio->SetPoint(i, graphParamsumm->GetX()[i], graphParamsumm->GetY()[i] / graphMainRatio->GetY()[i]);
      if(graphParamsumm->GetY()[i] > 0.) 
        graphParam_ratio->SetPointError(i, 0., 0., graphParamsumm->GetEYlow()[i] / graphMainRatio->GetY()[i], graphParamsumm->GetEYhigh()[i] / graphMainRatio->GetY()[i]);
      }
    }
    
    pad->cd();
    TPad* pad2 = new TPad("pad2","pad2", 0., 0., 1., RatioSize);   TrashBin->AddLast(pad2);
    pad2->Draw();

    pad1->SetBottomMargin(0.);
    pad2->SetTopMargin(0.);
    pad2->SetBottomMargin(0.2);

    pad2->cd();
    pad2->SetLogx();
    
    if (fDrawBase) {
     Base_ratio->SetLineColor(fColorBase);  
//     Base_ratio->SetLineStyle(Style_Ratio);  
     Base_ratio->SetFillColor(fColorBase);
     Base_ratio->GetYaxis()->SetTitle("ratio");
     Base_ratio->GetYaxis()->SetTitleSize(0.13);
     Base_ratio->GetYaxis()->SetTitleOffset(0.35);
     Base_ratio->GetYaxis()->SetNdivisions(505);
     Base_ratio->GetYaxis()->SetLabelSize(0.11);
     Base_ratio->GetXaxis()->SetLabelSize(0.11);
     Base_ratio->GetXaxis()->SetLabelOffset(0.03);
     Base_ratio->GetXaxis()->SetTitle("x  ");
     Base_ratio->GetXaxis()->SetTitleOffset(0.5);  
     Base_ratio->GetXaxis()->SetTitleSize(.15);
     Base_ratio->SetMaximum(1.19);
     Base_ratio->SetMinimum(0.81);
     Base_ratio->Draw("ALX");
    }
    if (fDrawExp) {
     Exp_ratio->SetLineColor(fColorExp + 2);  
     Exp_ratio->SetLineStyle(Style_Ratio);  
     Exp_ratio->SetFillColor(fColorExp);
     Exp_ratio->GetYaxis()->SetTitle("ratio");
     Exp_ratio->GetYaxis()->SetTitleSize(0.13);
     Exp_ratio->GetYaxis()->SetTitleOffset(0.35);
     Exp_ratio->GetYaxis()->SetNdivisions(505);
     Exp_ratio->GetYaxis()->SetLabelSize(0.11);
     Exp_ratio->GetXaxis()->SetLabelSize(0.11);
     Exp_ratio->GetXaxis()->SetLabelOffset(0.03);
     Exp_ratio->GetXaxis()->SetTitle("x  ");
     Exp_ratio->GetXaxis()->SetTitleOffset(0.5);  
     Exp_ratio->GetXaxis()->SetTitleSize(.15);
     Exp_ratio->SetMaximum(1.19);
     Exp_ratio->SetMinimum(0.81);
     Exp_ratio->Draw("ALX");
    }
    if (fDrawModel) {
     Model_ratio->SetLineColor(fColorModel + 2);  
     Model_ratio->SetLineStyle(Style_Ratio);  
     Model_ratio->SetFillColor(fColorModel);
     Model_ratio->GetYaxis()->SetTitle("ratio");
     Model_ratio->GetYaxis()->SetTitleSize(0.13);
     Model_ratio->GetYaxis()->SetTitleOffset(0.35);
     Model_ratio->GetYaxis()->SetNdivisions(505);
     Model_ratio->GetYaxis()->SetLabelSize(0.11);
     Model_ratio->GetXaxis()->SetLabelSize(0.11);
     Model_ratio->GetXaxis()->SetLabelOffset(0.03);
     Model_ratio->GetXaxis()->SetTitle("x  ");
     Model_ratio->GetXaxis()->SetTitleOffset(0.5);  
     Model_ratio->GetXaxis()->SetTitleSize(.15);
     Model_ratio->SetMaximum(1.19);
     Model_ratio->SetMinimum(0.81);
     Model_ratio->Draw("ALX");
    }
    if (fDrawParam) {
     Param_ratio->SetLineColor(fColorParam + 2);  
     Param_ratio->SetLineStyle(Style_Ratio);  
     Param_ratio->SetFillColor(fColorParam);
     Param_ratio->GetYaxis()->SetTitle("ratio");
     Param_ratio->GetYaxis()->SetTitleSize(0.13);
     Param_ratio->GetYaxis()->SetTitleOffset(0.35);
     Param_ratio->GetYaxis()->SetNdivisions(505);
     Param_ratio->GetYaxis()->SetLabelSize(0.11);
     Param_ratio->GetXaxis()->SetLabelSize(0.11);
     Param_ratio->GetXaxis()->SetLabelOffset(0.03);
     Param_ratio->GetXaxis()->SetTitle("x  ");
     Param_ratio->GetXaxis()->SetTitleOffset(0.5);  
     Param_ratio->GetXaxis()->SetTitleSize(.15);
     Param_ratio->SetMaximum(1.19);
     Param_ratio->SetMinimum(0.81);
     Param_ratio->Draw("ALX");
    }
    
  if (fDrawParam) {
    for (int ibin=0; ibin<graphParamsumm->GetN(); ibin++) {
      graphParam_ratio->GetPoint(ibin, xa, ya);
      yela = graphParam_ratio->GetErrorYlow(ibin);
      yeha = graphParam_ratio->GetErrorYhigh(ibin);
      
      graphRatio_lineUp->SetPoint(ibin, xa, (ya+yeha));
      graphRatio_lineDown->SetPoint(ibin, xa, (ya-yela));
    }
  } else if (fDrawModel) {
    for (int ibin=0; ibin<graphModelsumm->GetN(); ibin++) {
      graphModel_ratio->GetPoint(ibin, xa, ya);
      yela = graphModel_ratio->GetErrorYlow(ibin);
      yeha = graphModel_ratio->GetErrorYhigh(ibin);
      
      graphRatio_lineUp->SetPoint(ibin, xa, (ya+yeha));
      graphRatio_lineDown->SetPoint(ibin, xa, (ya-yela));
    }
  } else if (fDrawExp) {
    for (int ibin=0; ibin<graphExp->GetN(); ibin++) {
      graphExp_ratio->GetPoint(ibin, xa, ya);
      yela = graphExp_ratio->GetErrorYlow(ibin);
      yeha = graphExp_ratio->GetErrorYhigh(ibin);
      
      graphRatio_lineUp->SetPoint(ibin, xa, (ya+yeha));
      graphRatio_lineDown->SetPoint(ibin, xa, (ya-yela));
    }
  }
    
    if (fDrawParam) {
      graphParam_ratio->Draw("3");
      Param_ratio->Draw("LX");
    }
    if (fDrawModel) {
      graphModel_ratio->Draw("3");
      Model_ratio->Draw("LX");
    }
    if (fDrawExp) {
      graphExp_ratio->Draw("3");
      Exp_ratio->Draw("LX");
    }
    if (fDrawBase) {
      Base_ratio->Draw("LX");
    } 
    
    if (graphRatio_lineUp && graphRatio_lineDown) {
      graphRatio_lineUp->Draw("C same");
      graphRatio_lineDown->Draw("C same");
    }
       
  }

  TPaveLabel* label = new TPaveLabel(0.49, 0.85, 0.51, 0.87, Title,"NDC"); TrashBin->AddLast(label);
  label->SetBorderSize(0); label->SetFillColor(kWhite); label->SetTextSize(3.0);
  pad->cd();
  label->Draw();

  
  if(legend) {
    double irisey1 = 0.4;
    double irisey2 = 0.5;
    
    legend->cd();
    if (fDrawParam) {
      TPaveLabel* labparam = new TPaveLabel(0., irisey1, 1.0, irisey2, fOutputParam->GetName()->Data(), "NDC");
      labparam->SetTextSize(0.9);
      TrashBin->AddLast(labparam);
      labparam->SetFillColor(kWhite);
      labparam->SetBorderSize(0);
      labparam->SetTextColor(fColorParam);
      labparam->Draw();
      irisey1+=0.12; irisey2+=0.12;
    }
    if (fDrawModel) {
      TPaveLabel* labmodel = new TPaveLabel(0., irisey1, 1.0, irisey2, fOutputModel->GetName()->Data(), "NDC");
      labmodel->SetTextSize(0.9);
      TrashBin->AddLast(labmodel);
      labmodel->SetFillColor(kWhite);
      labmodel->SetBorderSize(0);
      labmodel->SetTextColor(fColorModel);
      labmodel->Draw();
      irisey1+=0.12; irisey2+=0.12;
    }
    if (fDrawExp) {
      TPaveLabel* labexp = new TPaveLabel(0., irisey1, 1.0, irisey2, fOutputExp->GetName()->Data(), "NDC");
      labexp->SetTextSize(0.9);
      TrashBin->AddLast(labexp);
      labexp->SetFillColor(kWhite);
      labexp->SetBorderSize(0);
      labexp->SetTextColor(fColorExp);
      labexp->Draw();
      irisey1+=0.12; irisey2+=0.12;
    }
    if (fDrawBase) {
      TPaveLabel* labbase = new TPaveLabel(0., irisey1, 1.0, irisey2, fOutputBase->GetName()->Data(), "NDC");
      labbase->SetTextSize(0.9);
      TrashBin->AddLast(labbase);
      labbase->SetFillColor(kWhite);
      labbase->SetBorderSize(0);
      labbase->SetTextColor(fColorBase);
      labbase->Draw();
      irisey1+=0.12; irisey2+=0.12;
    }
  }
  
  return 0;
}


void Painter::ScaleGraph2ToGraph1(TGraph* graph1, TGraph* graph2, TLine*& line, TGaxis*& axis, Double_t MeanRatio) {
  Double_t Max1 = -99999.; Double_t Min1 = 99999.;
  Double_t Max2 = -99999.; Double_t Min2 = 99999.;
  
  for(Int_t i=0; i<graph1->GetN(); i++) {
    if(graph1->GetY()[i] > Max1) Max1 = graph1->GetY()[i];
    if(graph2->GetY()[i] > Max2) Max2 = graph2->GetY()[i];
    if(graph1->GetY()[i] < Min1) Min1 = graph1->GetY()[i];
    if(graph2->GetY()[i] < Min2) Min2 = graph2->GetY()[i];
  }
  Max1 = Max1 * 1.2;
  Min1 = 0.;
  Double_t MinMax2 = Max2 - Min2;
  Max2 = Max2 + MinMax2/5.;
  Min2 = Min2 - MinMax2/5.;
  Max2 = MeanRatio + 0.2;
  Min2 = MeanRatio - 0.2;

  graph1->SetMinimum(Min1);
  graph1->SetMaximum(Max1);

  Double_t B = (Max1-Min1) / (Max2-Min2);
  Double_t A = Max1 - B*Max2;
  for(Int_t i=0; i<graph1->GetN(); i++) {
    graph2->SetPoint(i, graph2->GetX()[i], A + B * graph2->GetY()[i]);
  }
  line = new TLine(graph1->GetX()[0], A + B * MeanRatio, graph1->GetX()[graph1->GetN()-1], A + B * MeanRatio);
  axis = new TGaxis(graph1->GetXaxis()->GetXmax(), Min1,graph1->GetXaxis()->GetXmax(), Max1, Min2, Max2, 510,"+L");
  axis->SetLabelSize(0.04); axis->SetNdivisions(506); axis->SetLabelOffset(0.02);

  axis->SetTitle("ratio                                                                                                                        ");
  axis->SetTitleOffset(0.5); 
  axis->SetTextAngle(180.); 
}

Int_t Painter::DrawDataSet(DataSet* dataset, DataSet* datasetref, Bool_t RatioToData, EColor color) {

  int N = dataset->GetNSubPlots();

  TCanvas* can = new TCanvas;
  double MarkerSize, HistTitleY, HistTitleX;
  
  if(N<=0) return 1;
  else if (N==1) {                  MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N==2)  {can->Divide(1,2); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=4)  {can->Divide(2,2); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=6)  {can->Divide(2,3); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=9)  {can->Divide(3,3); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=12) {can->Divide(3,4); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=16) {can->Divide(4,4); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=20) {can->Divide(4,5); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=25) {can->Divide(5,5); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=30) {can->Divide(5,6); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=36) {can->Divide(6,6); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else {
    cout << "DrawDataSet does not support drawing more than 36 plots per dataset" <<endl; 
    delete can;
    return 1;
  }

  for(int i=0; i<N; i++) {
    can->cd(i+1);

    gStyle->SetTitleX(HistTitleX); //title X location
    gStyle->SetTitleY(HistTitleY); //title Y location
    //gStyle->SetTitleW(0.5); //title width
    //gStyle->SetTitleH(0.1); //title height 
    TH1F* h = dataset->GetHistogram(i,RatioToData);
    h->Draw();
    TGraphErrors* gDUnc = dataset->GetDataUnc(i);
    TGraphErrors* gDTot = dataset->GetDataTot(i);
    TGraphErrors* gTheo = dataset->GetTheo(i);
    TGraphErrors* gTMod = dataset->GetTMod(i);

    TGraphErrors* gDUncR = NULL;
    TGraphErrors* gDTotR = NULL;
    TGraphErrors* gTheoR = NULL;
    TGraphErrors* gTModR = NULL;
    
    if(datasetref) {
      gDUncR = datasetref->GetDataUnc(i);
      gDTotR = datasetref->GetDataTot(i);
      gTheoR = datasetref->GetTheo(i);
      gTModR = datasetref->GetTMod(i);
    }    

    if(RatioToData) {
      for(int j=0; j<gDUnc->GetN(); j++) {
	double div = gDUnc->GetY()[j];
	gDTot->GetY()[j] /=  div;
	gDTot->GetEY()[j] /=  div;
	gTheo->GetY()[j] /=  div;
	gTMod->GetY()[j] /=  div;
	gDUnc->GetEY()[j] /= div;
	gDUnc->GetY()[j] /=  div;
	gDTot->GetY()[j] -= 1.;
	gTheo->GetY()[j] -= 1.;
	gTMod->GetY()[j] -= 1.;
	gDUnc->GetY()[j] -= 1.;

	if(datasetref) {
	  div = gDUncR->GetY()[j];
	  gDTotR->GetY()[j] /=  div;
	  gDTotR->GetEY()[j] /=  div;
	  gTheoR->GetY()[j] /=  div;
	  gTModR->GetY()[j] /=  div;
	  gDUncR->GetEY()[j] /= div;
	  gDUncR->GetY()[j] /=  div;
	  gDTotR->GetY()[j] -= 1.;
	  gTheoR->GetY()[j] -= 1.;
	  gTModR->GetY()[j] -= 1.;
	  gDUncR->GetY()[j] -= 1.;
	}
      }
    }

    if(dataset->GetXlog(i)) gPad->SetLogx();
    if(dataset->GetYlog(i)) gPad->SetLogy();

    gDUnc->SetMarkerStyle(20);
    gDUnc->SetMarkerSize(MarkerSize);
    gTheo->SetLineColor(color);
    gTMod->SetLineColor(color);
    gTMod->SetLineStyle(3);

    gDUnc->Draw("same P");
    gDTot->Draw("same P");
    gTheo->Draw("same L");
    gTMod->Draw("same L");

    if(datasetref) {
      gTheoR->SetLineColor(fColorRef);
      gTModR->SetLineColor(fColorRef);
      gTModR->SetLineStyle(3);
      gTheoR->Draw("same L");
      gTModR->Draw("same L");
    }
  }
  
  TPaveLabel label(0.0, 0.985, 0.5, 1.0, dataset->GetName(),"NDC");
  //cout << dataset->GetName() << endl;
  label.SetFillColor(kWhite);
  label.SetBorderSize(0);
  can->cd(0);
  label.Draw();

  TLegend leg(0.5, 0.95, 1.0, 1.0, "","NDC");
  TString* temp = new TString;
  leg.SetBorderSize(0); 
  leg.SetFillColor(kWhite);
  if(datasetref) leg.SetNColumns(2); 

  leg.AddEntry(dataset->GetTheo(0), fOutput->GetName()->Data(),"L");
  if(datasetref) 
    leg.AddEntry(datasetref->GetTheo(0), fOutputRef->GetName()->Data(),"L");

  temp->Form("%s (modified)", fOutput->GetName()->Data());
  leg.AddEntry(dataset->GetTMod(0), temp->Data(),"L");
  if(datasetref) {
    temp->Form("%s (modfied)", fOutputRef->GetName()->Data());
    leg.AddEntry(datasetref->GetTMod(0), temp->Data(),"L");
  }
  //leg.Draw();

  
  PrintCanvas(can);
  delete can; delete temp;
}

Int_t Painter::DrawDataSetEMP(DataSet* datasetB, DataSet* datasetE, DataSet* datasetM, DataSet* datasetP, Bool_t RatioToData) {

  int N;
  
  if (datasetP) N = datasetP->GetNSubPlots();
    else if (datasetM) N = datasetM->GetNSubPlots();
     else if (datasetE) N = datasetE->GetNSubPlots();
      else if (datasetB) N = datasetB->GetNSubPlots();
  
  TCanvas* can = new TCanvas;
  double MarkerSize, HistTitleY, HistTitleX;
  
  if(N<=0) return 1;
  else if (N==1) {                  MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N==2)  {can->Divide(1,2); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=4)  {can->Divide(2,2); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=6)  {can->Divide(2,3); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=9)  {can->Divide(3,3); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=12) {can->Divide(3,4); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=16) {can->Divide(4,4); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=20) {can->Divide(4,5); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=25) {can->Divide(5,5); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=30) {can->Divide(5,6); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else if(N<=36) {can->Divide(6,6); MarkerSize = 0.3; HistTitleY = 0.97; HistTitleX = 0.1;}
  else {
    cout << "DrawDataSetEMP does not support drawing more than 36 plots per dataset" <<endl; 
    delete can;
    return 1;
  }
  
  for(int i=0; i<N; i++) {
    can->cd(i+1);

    gStyle->SetTitleX(HistTitleX); //title X location
    gStyle->SetTitleY(HistTitleY); //title Y location
    
    TH1F* h;
    if (datasetP) h = datasetP->GetHistogram(i,RatioToData);
     else if (datasetM) h = datasetM->GetHistogram(i,RatioToData);
      else if (datasetE) h = datasetE->GetHistogram(i,RatioToData);
       else if (datasetB) h = datasetB->GetHistogram(i,RatioToData);
    h->Draw();
    
    TGraphErrors* gDUncBase = NULL;
    TGraphErrors* gDTotBase = NULL;
    TGraphErrors* gTheoBase = NULL;
    TGraphErrors* gTModBase = NULL;
    
    TGraphErrors* gDUncExp = NULL;
    TGraphErrors* gDTotExp = NULL;
    TGraphErrors* gTheoExp = NULL;
    TGraphErrors* gTModExp = NULL;
    
    TGraphErrors* gDUncModel = NULL;
    TGraphErrors* gDTotModel = NULL;
    TGraphErrors* gTheoModel = NULL;
    TGraphErrors* gTModModel = NULL;
    
    TGraphErrors* gDUncParam = NULL;
    TGraphErrors* gDTotParam = NULL;
    TGraphErrors* gTheoParam = NULL;
    TGraphErrors* gTModParam = NULL;
    
    if(datasetB) {
      gDUncBase = datasetB->GetDataUnc(i);
      gDTotBase = datasetB->GetDataTot(i);
      gTheoBase = datasetB->GetTheo(i);
      gTModBase = datasetB->GetTMod(i);
    }
    if(datasetE) {
      gDUncExp = datasetE->GetDataUnc(i);
      gDTotExp = datasetE->GetDataTot(i);
      gTheoExp = datasetE->GetTheo(i);
      gTModExp = datasetE->GetTMod(i);
    }
    if(datasetM) {
      gDUncModel = datasetM->GetDataUnc(i);
      gDTotModel = datasetM->GetDataTot(i);
      gTheoModel = datasetM->GetTheo(i);
      gTModModel = datasetM->GetTMod(i);
    }
    if(datasetP) {
      gDUncParam = datasetP->GetDataUnc(i);
      gDTotParam = datasetP->GetDataTot(i);
      gTheoParam = datasetP->GetTheo(i);
      gTModParam = datasetP->GetTMod(i);
    }
    
    if(RatioToData) {
    
      int NgDU = 0;
      double div = 0.;
      
      if (datasetP) NgDU = gDUncParam->GetN();
       else if (datasetM) NgDU = gDUncModel->GetN();
        else if (datasetE) NgDU = gDUncExp->GetN();
         else if (datasetB) NgDU = gDUncBase->GetN();
	 
      for(int j=0; j<NgDU; j++) {
	if(datasetB) {
	  div = gDUncBase->GetY()[j];
	  gDTotBase->GetY()[j] /=  div;
	  gDTotBase->GetEY()[j] /=  div;
	  gTheoBase->GetY()[j] /=  div;
	  gTModBase->GetY()[j] /=  div;
	  gDUncBase->GetY()[j] /=  div;
	  gDUncBase->GetEY()[j] /= div;
	  gDTotBase->GetY()[j] -= 1.;
	  gTheoBase->GetY()[j] -= 1.;
	  gTModBase->GetY()[j] -= 1.;
	  gDUncBase->GetY()[j] -= 1.;
	}
	if(datasetE) {
	  div = gDUncExp->GetY()[j];
	  gDTotExp->GetY()[j] /=  div;
	  gDTotExp->GetEY()[j] /=  div;
	  gTheoExp->GetY()[j] /=  div;
	  gTModExp->GetY()[j] /=  div;
	  gDUncExp->GetY()[j] /=  div;
	  gDUncExp->GetEY()[j] /= div;
	  gDTotExp->GetY()[j] -= 1.;
	  gTheoExp->GetY()[j] -= 1.;
	  gTModExp->GetY()[j] -= 1.;
	  gDUncExp->GetY()[j] -= 1.;
	}
	if(datasetM) {
	  div = gDUncModel->GetY()[j];
	  gDTotModel->GetY()[j] /=  div;
	  gDTotModel->GetEY()[j] /=  div;
	  gTheoModel->GetY()[j] /=  div;
	  gTModModel->GetY()[j] /=  div;
	  gDUncModel->GetY()[j] /=  div;
	  gDUncModel->GetEY()[j] /= div;
	  gDTotModel->GetY()[j] -= 1.;
	  gTheoModel->GetY()[j] -= 1.;
	  gTModModel->GetY()[j] -= 1.;
	  gDUncModel->GetY()[j] -= 1.;
	}
	if(datasetP) {
	  div = gDUncParam->GetY()[j];
	  gDTotParam->GetY()[j] /=  div;
	  gDTotParam->GetEY()[j] /=  div;
	  gTheoParam->GetY()[j] /=  div;
	  gTModParam->GetY()[j] /=  div;
	  gDUncParam->GetY()[j] /=  div;
	  gDUncParam->GetEY()[j] /= div;
	  gDTotParam->GetY()[j] -= 1.;
	  gTheoParam->GetY()[j] -= 1.;
	  gTModParam->GetY()[j] -= 1.;
	  gDUncParam->GetY()[j] -= 1.;
	}
      }
    }
    
    if(datasetB) {
      gTheoBase->SetLineColor(fColorBase);
      gTModBase->SetLineColor(fColorBase);
      gTModBase->SetLineStyle(3);
    }
    if(datasetE) {
      gTheoExp->SetLineColor(fColorExp);
      gTModExp->SetLineColor(fColorExp);
      gTModExp->SetLineStyle(3);
    }
    if(datasetM) {
      gTheoModel->SetLineColor(fColorModel);
      gTModModel->SetLineColor(fColorModel);
      gTModModel->SetLineStyle(3);
    }
    if(datasetP) {
      gTheoParam->SetLineColor(fColorParam);
      gTModParam->SetLineColor(fColorParam);
      gTModParam->SetLineStyle(3);
    }
    
    if (datasetP) {
      if(datasetP->GetXlog(i)) gPad->SetLogx();
      if(datasetP->GetYlog(i)) gPad->SetLogy();
      gDUncParam->SetMarkerStyle(20);
      gDUncParam->SetMarkerSize(MarkerSize);
      gDUncParam->Draw("same P");
      gDTotParam->Draw("same P");
    } else if (datasetM) {
      if(datasetM->GetXlog(i)) gPad->SetLogx();
      if(datasetM->GetYlog(i)) gPad->SetLogy();
      gDUncModel->SetMarkerStyle(20);
      gDUncModel->SetMarkerSize(MarkerSize);
      gDUncModel->Draw("same P");
      gDTotModel->Draw("same P");
    } else if (datasetE) {
      if(datasetE->GetXlog(i)) gPad->SetLogx();
      if(datasetE->GetYlog(i)) gPad->SetLogy();
      gDUncExp->SetMarkerStyle(20);
      gDUncExp->SetMarkerSize(MarkerSize);
      gDUncExp->Draw("same P");
      gDTotExp->Draw("same P");
    } else if (datasetB) {
      if(datasetB->GetXlog(i)) gPad->SetLogx();
      if(datasetB->GetYlog(i)) gPad->SetLogy();
      gDUncBase->SetMarkerStyle(20);
      gDUncBase->SetMarkerSize(MarkerSize);
      gDUncBase->Draw("same P");
      gDTotBase->Draw("same P");
    }

    if (datasetP) {
      gTheoParam->Draw("same L");
      gTModParam->Draw("same L");
    }
    if (datasetM) {
      gTheoModel->Draw("same L");
      gTModModel->Draw("same L");
    }
    if (datasetE) {
      gTheoExp->Draw("same L");
      gTModExp->Draw("same L");
    }
    if (datasetB) {
      gTheoBase->Draw("same L");
      gTModBase->Draw("same L");
    }
    
  }
    
  TString labelName;
  
  if (datasetP) labelName = TString(datasetP->GetName());
   else if (datasetM) labelName = TString(datasetM->GetName());
    else if (datasetE) labelName = TString(datasetE->GetName());
     else if (datasetB) labelName = TString(datasetB->GetName());
  
  TPaveLabel label(0.0, 0.985, 0.5, 1.0, labelName,"NDC");
  //cout << labelName << endl;
  label.SetFillColor(kWhite);
  label.SetBorderSize(0);
  can->cd(0);
  label.Draw();

  TLegend leg(0.5, 0.95, 1.0, 1.0, "","NDC");
  TString* temp = new TString;
  leg.SetBorderSize(0); 
  leg.SetFillColor(kWhite);
  
  int legNcol = 0;
  if (datasetP) legNcol++;
  if (datasetM) legNcol++;
  if (datasetE) legNcol++;
  if (datasetB) legNcol++;
  if (legNcol > 1) leg.SetNColumns(legNcol);

  if (datasetP) leg.AddEntry(datasetP->GetTheo(0), fOutputParam->GetName()->Data(),"L");
  if (datasetM) leg.AddEntry(datasetM->GetTheo(0), fOutputModel->GetName()->Data(),"L");
  if (datasetE) leg.AddEntry(datasetE->GetTheo(0), fOutputExp->GetName()->Data(),"L");
  if (datasetB) leg.AddEntry(datasetB->GetTheo(0), fOutputBase->GetName()->Data(),"L");

  if(datasetP) {
    temp->Form("%s (modfied)", fOutputParam->GetName()->Data());
    leg.AddEntry(datasetP->GetTMod(0), temp->Data(),"L");
  }
  if(datasetM) {
    temp->Form("%s (modfied)", fOutputModel->GetName()->Data());
    leg.AddEntry(datasetM->GetTMod(0), temp->Data(),"L");
  }
  if(datasetE) {
    temp->Form("%s (modfied)", fOutputExp->GetName()->Data());
    leg.AddEntry(datasetE->GetTMod(0), temp->Data(),"L");
  }
  if(datasetB) {
    temp->Form("%s (modfied)", fOutputBase->GetName()->Data());
    leg.AddEntry(datasetB->GetTMod(0), temp->Data(),"L");
  }
  //leg.Draw();

  
  PrintCanvas(can);
  delete can; delete temp;
  
  return 0;
}

void Painter::PrintCanvas(TCanvas* can) {
  if(fPsFileName->Contains(".txt"))
    fPsFileName->ReplaceAll(".txt/","_");
  can->Print(fPsFileName->Data());
  fPsFileName->ReplaceAll("(","");
  TString* temp = new TString;
  static Int_t idx = 0;
  idx++;
  temp->Form("DrawResults_%03d.eps",idx);
  can->Print(temp->Data());
  delete temp;
}
