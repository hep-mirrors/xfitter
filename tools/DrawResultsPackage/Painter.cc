#include <Painter.h>
#include <TPaveLabel.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TPaveText.h>
#include <TStyle.h>

Painter::Painter(bool DrawBands){
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
}

Int_t Painter::Prepare() {
  fOutput = new Output(fPath->Data());
  if(fPath->CompareTo(fPathRef->Data())) 
    fOutputRef = new Output(fPathRef->Data());
  else 
    fOutputRef = NULL;

  if(fOutputRef == NULL) fPsFileName->Form("%s/DrawResults.ps", fPathRef->Data());
  fPsFileName->Append("(");

  fOutput->Prepare(fBands);
  if(fOutputRef) fOutputRef->Prepare(fBands); 
}

Int_t Painter::Draw() {
  this->Prepare();

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
  if      (Option.Contains("M"))   T->SetTextColor(fColor);
  else if (Option.Contains("R"))   T->SetTextColor(fColorRef);

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
      AddLineToPave(paves, ypos, TString((*(output->GetName()))(i*NChar, NChar)).Data() ,"BM");
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
      AddLineToPave(paves, ypos, TString((*(output->GetName()))(i*NChar, NChar)).Data() ,"BM");
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
  Double_t fQ2val = fOutput->GetQ2Value(ival);

  // Label:
  char label[32]; 
  sprintf (label,"Q^{2} = %12.2f GeV^{2}",fQ2val);
  Q2Label->SetLabel(label);


  TCanvas* can = new TCanvas("can","can",600, 400);
  can->Divide(3,2);

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

  can->cd();
  Q2Label->Draw();
  PrintCanvas(can);
  delete TrashBin;
  delete can;
}


Int_t Painter::PlotPdfSub(TVirtualPad* pad, Output* FitterOut, Output* FitterRef, 
				  Int_t Q2Bin, const Char_t* Title, Output::pdf pdf1, Output::pdf pdf2,
				  TVirtualPad* legend, TObjArray* TrashBin) {
  TGraph* graph =   FitterOut->GetPdf(pdf1, Q2Bin);
  TGraph* graphR = NULL; 
  TGraph* graph2 = NULL;
  TGraph* graphR2 = NULL;
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
