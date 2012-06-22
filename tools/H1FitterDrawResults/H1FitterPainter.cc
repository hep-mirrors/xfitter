#include <H1FitterPainter.h>
#include <TPaveLabel.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TPaveText.h>

H1FitterPainter::H1FitterPainter(bool DrawBands){
  fPath = new TString("../../output/");
  fPathRef = new TString("../../output/");
  fH1FitterOutput = NULL;
  fH1FitterOutputRef = NULL;
  fPsFileName = new TString("DrawResults.ps");
  gROOT->SetStyle("Plain");
  cout << endl;
  //cout << "TO DO: in fittedresults.txt q2 and x for sets 61-64 are switched"<<endl;
  fColor = kRed;
  fColorRef = kBlue;
  fFillStyle = 1001;
  fFillStyleRef = 0; //3010
  fBands = DrawBands;
}

H1FitterPainter::~H1FitterPainter(){ 
  delete fPath;
  delete fPathRef;
  if(fH1FitterOutput) delete fH1FitterOutput;
  if(fH1FitterOutputRef) delete fH1FitterOutputRef;

  TCanvas* can = new TCanvas;
  cout << "Output stored in " << fPsFileName->Data() << " file"<<endl;
  fPsFileName->Append(")");
  can->Print(fPsFileName->Data());
  delete can;

  delete fPsFileName;
}

Int_t H1FitterPainter::Prepare() {
  fH1FitterOutput = new H1FitterOutput(fPath->Data());
  if(fPath->CompareTo(fPathRef->Data())) 
    fH1FitterOutputRef = new H1FitterOutput(fPathRef->Data());
  else 
    fH1FitterOutputRef = NULL;

  if(fH1FitterOutputRef == NULL) fPsFileName->Form("%s/DrawResults.ps", fPathRef->Data());
  fPsFileName->Append("(");

  fH1FitterOutput->Prepare(fBands);
  if(fH1FitterOutputRef) fH1FitterOutputRef->Prepare(fBands); 
}

Int_t H1FitterPainter::Draw() {
  this->Prepare();

  // Get number of Q2 files:
  Int_t NQ2Files = fH1FitterOutput->GetNQ2Files();

  for(Int_t i=0; i<NQ2Files; i++) {
    this->DrawPDF(i);
  }
  
  this->DrawPull();
  
  Int_t NRefDataSets = 0;
  if(fH1FitterOutputRef) NRefDataSets = fH1FitterOutputRef->GetNsets();
  Bool_t RefSetsDrawn[NRefDataSets];
  for(Int_t i=0; i<NRefDataSets; i++) RefSetsDrawn[i] = kFALSE;

  for (Int_t i=0; i<fH1FitterOutput->GetNsets(); i++) {
    DataSet* dataset = fH1FitterOutput->GetSet(i);
    DataSet* datasetref = NULL;
    for (Int_t iref=0; iref<NRefDataSets; iref++) {
      if(fH1FitterOutputRef->GetSet(iref)->GetSetId() == dataset->GetSetId()) {
	datasetref = fH1FitterOutputRef->GetSet(iref);
	RefSetsDrawn[iref] = kTRUE;
	break;
      }
    } 
    DrawDataSet(dataset, datasetref);
    //DrawDataSetRatio(dataset, datasetref);
  }
  for (Int_t iref=0; iref<NRefDataSets; iref++) {
    if(!RefSetsDrawn[iref])
      DrawDataSet(fH1FitterOutputRef->GetSet(iref), NULL, fColorRef);
  }

  this->DrawFitResults();

  this->DrawMessages(fH1FitterOutput);
  if(fH1FitterOutputRef) 
    this->DrawMessages(fH1FitterOutputRef);
}

void H1FitterPainter::AddLineToPave(TObjArray* paves, float& yposition, const char* text, const char* option) {
  static int NLines = 30;
  //  static float yposition = 1.;

  TString Option(option);
  TPaveText* pave = (TPaveText*) paves->At(paves->GetEntries()-1);
  if(!pave) return;

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
}

void H1FitterPainter::FillPavesWithFitResults(TObjArray* paves, H1FitterOutput* output) {

  TObjArray* names = output->GetFittedParametersNames();
  TObjArray* namesNuisance= output->GetNuisanceParNames();
  TString* str = new TString;

  float ypos = 1.;
  int NChar = 40;

  AddLineToPave(paves, ypos, "Results for:","B");

  for(int i=0; i<=output->GetName()->Length()/NChar; i++) {
    if(output==fH1FitterOutput)
      AddLineToPave(paves, ypos, TString((*(output->GetName()))(i*NChar, NChar)).Data() ,"BM");
    else
      AddLineToPave(paves, ypos, TString((*(output->GetName()))(i*NChar, NChar)).Data() ,"BR");
  }

  AddLineToPave(paves, ypos, "","");

  str->Form("Fitted %d parameters:", names->GetEntries()); 
  AddLineToPave(paves, ypos, str->Data(),"");
  str->Form("(most reliable available method: %s" , output->GetErrorCalculationMethod()->Data()); 
  AddLineToPave(paves, ypos, str->Data(),"");
  str->Form("giving confidence in errors: %s)", output->GetErrorTrustLevel()->Data()); 
  AddLineToPave(paves, ypos, str->Data(),"");


  for(int i=0; i<names->GetEntries(); i++) {
    str->Form("%4d:\t%6s = %6.3f  #pm%6.3f", i+1, ((TObjString*)names->At(i))->GetString().Data(), 
	      output->GetFittedParameter(i, false), output->GetFittedParameter(i, true));
    AddLineToPave(paves, ypos, str->Data(),"");
  }
  
  AddLineToPave(paves, ypos, "","");
  AddLineToPave(paves, ypos, "Nuisance Parameters:","");
  for(int i=0; i<namesNuisance->GetEntries(); i++) {
    str->Form("%4d:\t%17s = %5.2f  #pm%5.2f", i+1, ((TObjString*)namesNuisance->At(i))->GetString().Data(), 
	      output->GetNuisancePar(i, false), output->GetNuisancePar(i, true));
    AddLineToPave(paves, ypos, str->Data(),""); 
  }
  delete str;
}

Int_t H1FitterPainter::DrawFitResults() {

  TCanvas* can = new TCanvas;
  TObjArray* pavesL = new TObjArray; pavesL->SetOwner();
  TObjArray* pavesR = new TObjArray; pavesR->SetOwner();

  pavesL->AddLast(new TPaveText(0.05, 0.05, 0.5, 0.95, "br"));
  FillPavesWithFitResults(pavesL, fH1FitterOutput);
  if(fH1FitterOutputRef) {
    pavesR->AddLast(new TPaveText(0.5, 0.05, 0.95, 0.95, "br"));
    FillPavesWithFitResults(pavesR, fH1FitterOutputRef);
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

void H1FitterPainter::DrawMessages(H1FitterOutput* output) {

  TObjArray* paves = new TObjArray; paves->SetOwner();
  paves->AddLast(new TPaveText(0.05, 0.05, 0.95, 0.95, "br"));
  float ypos = 1.;
  int NChar = 80;
  
  AddLineToPave(paves, ypos, "Fit messages for:","B");
  for(int i=0; i<=output->GetName()->Length()/NChar; i++) {
    if(output==fH1FitterOutput)
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


Int_t H1FitterPainter::DrawPull() {
//  TCanvas* can = new TCanvas;
//
//  TH1F* h = fH1FitterOutput->GetPull();
//  TH1F* hRef = NULL;
//  if(fH1FitterOutputRef) hRef = fH1FitterOutputRef->GetPull();
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

Int_t H1FitterPainter::DrawPDF(Int_t ival) {
  TObjArray* TrashBin = new TObjArray; TrashBin->SetOwner();

  TPaveLabel* Q2Label = new TPaveLabel(0.7, 0.08, 1.0, 0.16, "","NDC"); TrashBin->AddLast(Q2Label);
  Q2Label->SetBorderSize(1); Q2Label->SetFillColor(kWhite); Q2Label->SetBorderSize(0);

  // Q2 value of the bin:
  Double_t fQ2val = fH1FitterOutput->GetQ2Value(ival);

  // Label:
  char label[32]; 
  sprintf (label,"Q^{2} = %12.2f GeV^{2}",fQ2val);
  Q2Label->SetLabel(label);


  TCanvas* can = new TCanvas("can","can",600, 400);
  can->Divide(3,2);

  PlotPdfSub(can->cd(1), fH1FitterOutput, fH1FitterOutputRef, ival, "xU, xu_{V}",
	     H1FitterOutput::kU,     H1FitterOutput::kUv, can->cd(6), TrashBin);
  PlotPdfSub(can->cd(2), fH1FitterOutput, fH1FitterOutputRef, ival, "xD, xd_{V}",
	     H1FitterOutput::kD,     H1FitterOutput::kDv, NULL,  TrashBin);
  PlotPdfSub(can->cd(3), fH1FitterOutput, fH1FitterOutputRef, ival, "xg",
	     H1FitterOutput::kGluon, H1FitterOutput::kNULL,NULL,  TrashBin);
  PlotPdfSub(can->cd(4), fH1FitterOutput, fH1FitterOutputRef, ival, "x#bar{U}",
	     H1FitterOutput::kUb,   H1FitterOutput::kNULL, NULL, TrashBin);
  PlotPdfSub(can->cd(5), fH1FitterOutput, fH1FitterOutputRef, ival, "x#bar{D}",
	     H1FitterOutput::kDb,    H1FitterOutput::kNULL, NULL, TrashBin);

  can->cd();
  Q2Label->Draw();
  PrintCanvas(can);
  delete TrashBin;
  delete can;
}


Int_t H1FitterPainter::PlotPdfSub(TVirtualPad* pad, H1FitterOutput* FitterOut, H1FitterOutput* FitterRef, 
				  Int_t Q2Bin, const Char_t* Title, H1FitterOutput::pdf pdf1, H1FitterOutput::pdf pdf2,
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
  if(pdf2 != H1FitterOutput::kNULL) {
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

void H1FitterPainter::ScaleGraph2ToGraph1(TGraph* graph1, TGraph* graph2, TLine*& line, TGaxis*& axis, Double_t MeanRatio) {
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


Int_t H1FitterPainter::DrawDataSetRatio(DataSet* dataset, DataSet* datasetref) {

  if(datasetref==NULL) return 1;

  TGraphErrors* gTheo;
  TGraphErrors* gTheoRef;
  TObjArray* TrashBin = new TObjArray; TrashBin->SetOwner();
    
  TCanvas* can = new TCanvas;
  if     (dataset->GetNGraphs()<=1)  can->Divide(1,1);
  else if(dataset->GetNGraphs()<=2)  can->Divide(1,2);
  else if(dataset->GetNGraphs()<=4)  can->Divide(2,2);
  else if(dataset->GetNGraphs()<=9)  can->Divide(3,3);
  else if(dataset->GetNGraphs()<=16) can->Divide(4,4);
  else if(dataset->GetNGraphs()<=25) can->Divide(5,5);
  else if(dataset->GetNGraphs()<=36) can->Divide(6,6);
  else {cout << "can not DrawDataSetRatio, too many graphs " <<dataset->GetNGraphs() <<endl; return 1;}

  for(Int_t i=0; i<dataset->GetNGraphs(); i++) {
    gTheo = (TGraphErrors*) dataset->GetTheory(i)->Clone();
    TrashBin->AddLast(gTheo);
    gTheoRef= datasetref->GetTheory(i);

    for(Int_t k=0; k<gTheo->GetN(); k++) {
      gTheo->GetY()[k] = gTheo->GetY()[k] / gTheoRef->GetY()[k];
    }
    gTheo->SetTitle(dataset->GetLabel(i)->GetString().Data());
    gTheo->SetMaximum(1.1);
    gTheo->SetMinimum(0.9);
    gTheo->SetMarkerStyle(20);
    gTheo->SetMarkerSize(0.5);
    can->cd(i+1);
    gTheo->Draw("ALP");
  }
  TPaveLabel* Title = new TPaveLabel(0.0, 0.97, 0.6, 1.0, dataset->GetName(), "NDC");
  Title->SetFillColor(kWhite); Title->SetBorderSize(0);
  can->cd(0);Title->Draw();

  PrintCanvas(can);
  delete can; delete TrashBin; delete Title;
  return 0;
}

Int_t H1FitterPainter::DrawDataSet(DataSet* dataset, DataSet* datasetref, EColor color) {

  if(dataset->GetNGraphs()==0) return 1;

  Int_t SetId = dataset->GetSetId();
  if(SetId==61 || SetId==62 || SetId==63 || SetId==64 || SetId==35) { // Painting for all the supported SetIds
    
    TObjArray* TrashBin = new TObjArray; TrashBin->SetOwner();

    Double_t YMaximum = 0.; // axis borders
    Double_t YMinimum = 0.;
    Double_t XMinimum = 0.;
    Double_t XMaximum = 1.;

    Double_t MarkerSize = 0.005; // size of the data point marker
    Double_t ALabelSize = 0.06;  // size of the axis' label
    Bool_t Logx = kFALSE;
    Bool_t Logy = kFALSE;

    Double_t LabelSize = 1.;  // size of the graph label text
    Double_t Chi2Size = 1.;   // size of the chi2 text

    TCanvas* can = new TCanvas;
    can->SetTopMargin(0.3);
    switch(SetId) { // Canvas definition for supported SetIds
    case 61: 
      can->Divide(4,5,0.); 
      XMinimum = 0.001;   
      XMaximum = 1.; 
      YMaximum = 1.99; 
      YMinimum = 0.0001; 
      MarkerSize = 0.005; 
      ALabelSize = 0.06; 
      Chi2Size = 1.; 
      LabelSize = .6; 
      Logx = kTRUE;
      Logy = kFALSE;
      break;
    case 62: 
      can->Divide(6,6,0.); 
      XMinimum = 0.00003; 
      XMaximum = 1.; 
      YMaximum = 1.59; 
      YMinimum = 0.0001; 
      MarkerSize = 0.005; 
      ALabelSize = 0.06; 
      Chi2Size = 1.; 
      LabelSize = .6; 
      Logx = kTRUE;
      Logy = kFALSE;
      break;
    case 63: 
      can->Divide(3,4,0.); 
      XMinimum = 0.002;   
      XMaximum = 1.; 
      YMaximum = 0.99; 
      YMinimum = 0.0001; 
      MarkerSize = 0.5;   
      ALabelSize = 0.08; 
      Chi2Size = 1.; 
      LabelSize = .6; 
      Logx = kTRUE;
      Logy = kFALSE;
      break;
    case 64: 
      can->Divide(3,3,0.); 
      XMinimum = 0.002;   
      XMaximum = 1.; 
      YMaximum = 1.99; 
      YMinimum = 0.0001; 
      MarkerSize = 0.5;   
      ALabelSize = 0.08; 
      Chi2Size = 1.; 
      LabelSize = .6; 
      Logx = kTRUE;
      Logy = kFALSE;
      break;
    case 35: 
      can->Divide(3,2,0.); 
      XMinimum = 5.;      
      XMaximum = 50.;
      YMaximum = 990.;  
      YMinimum = 0.15;     
      MarkerSize = 0.5;   
      ALabelSize = 0.05; 
      Chi2Size = .6; 
      LabelSize = 0.4; 
      Logx = kTRUE;
      Logy = kTRUE;
      break;
    }
    
    TGraphErrors* gDUnc;
    TGraphErrors* gDTot;
    TGraphErrors* gTheo;
    TGraphErrors* gTheoRef;
    
    for(Int_t i=0; i<dataset->GetNGraphs(); i++) {
      gDUnc = dataset->GetDataUncr(i);
      gDTot = dataset->GetDataTotal(i);
      gTheo = dataset->GetTheory(i);
      gTheoRef = NULL;
      if(datasetref) gTheoRef= datasetref->GetTheory(i);

      Double_t chi2offx = 0.;
      Double_t chi2offy = 0.;
      
      // Additional plotting changes for supported SetIds

      if(SetId==61 && i%4==0) chi2offx = 0.1;
      if(SetId==61 && i>15)   chi2offy = 0.08;


      if(SetId==62) chi2offx = 0.1;
      if(SetId==62 && i%6==0) chi2offx += 0.1;
      if(SetId==62 && i>29)   chi2offy = 0.08;

      if(SetId==63) {chi2offy = -0.05; chi2offx = -0.05;}
      if(SetId==63 && i<3) {chi2offy += 0.4; chi2offx += 0.5;}
      if(SetId==63 && i==0) {chi2offy += 0.0; chi2offx += -0.1;}
      if(SetId==63 && i%3==0) chi2offx += 0.1;
      if(SetId==63 && i>8)   chi2offy += 0.08;

      if(SetId==64) {chi2offy = -0.00; chi2offx = 0.02;}
      if(SetId==64 && i==0) {chi2offy = 0.0; chi2offx = -0.08;}
      if(SetId==64 && i<3) {chi2offy += 0.4; chi2offx += 0.4;}
      if(SetId==64 && i%3==0) chi2offx += 0.08;
      if(SetId==64 && i>5)   {chi2offy += 0.05; chi2offx += 0.05;}


      if(SetId==63 && i>2) YMaximum = 0.299;
      if(SetId==63 && i>5) {YMaximum = 0.0399;; YMinimum = 0.000001;}
      if(SetId==63 && i>8) {YMaximum = 0.000499; YMinimum = 0.000001;}

      if(SetId==64 && i>2) YMaximum = 0.1199;
      if(SetId==64 && i>5) YMaximum = 0.01599;


      if(SetId==35 && i>2)    LabelSize = 0.35;
      if(SetId==35 && i%3==0) chi2offx += 0.08;
       


      TH1F* h = new TH1F("","",1, XMinimum, XMaximum); TrashBin->SetOwner(h);
      h->SetMinimum(YMinimum);
      h->SetMaximum(YMaximum);
      h->SetStats(kFALSE);
      h->GetYaxis()->SetLabelSize(ALabelSize);
      h->GetXaxis()->SetLabelSize(ALabelSize);


      TObjString* lstring = dataset->GetLabel(i);
      TPaveLabel* label = new TPaveLabel(0.4, 0.7, 0.9, 0.9,"null","NDC"); TrashBin->AddLast(label);
      if(lstring)	label->SetLabel(lstring->GetString().Data());
      label->SetFillColor(kWhite);
      label->SetBorderSize(0);
      label->SetTextSize(LabelSize);

      Double_t Chi2 = dataset->GetChi2(i);
      Int_t Npoints = dataset->GetNpts(i);


      TPaveLabel* labelchi2 = new TPaveLabel(0.25+chi2offx, 0.1+chi2offy, 0.3+chi2offx, 0.2+chi2offy,"","NDC"); TrashBin->AddLast(labelchi2);
      TString* temp = new TString; temp->Form("#chi^{2} / npts = %.1f / %d", Chi2, Npoints);
      labelchi2->SetLabel(temp->Data());
      delete temp;
      labelchi2->SetFillColor(kWhite);
      labelchi2->SetBorderSize(0);
      labelchi2->SetTextSize(Chi2Size);
      labelchi2->SetTextColor(fColor);

      TPaveLabel* labelchi2ref = NULL;
      if(datasetref) {
	if(SetId==63 && i < 3) chi2offy += 0.03;
	Chi2 = datasetref->GetChi2(i);
	Npoints = datasetref->GetNpts(i);
	labelchi2ref = new TPaveLabel(0.25+chi2offx, 0.2+chi2offy, 0.3+chi2offx, 0.3+chi2offy,"","NDC"); 
	TrashBin->AddLast(labelchi2ref);
	TString* temp = new TString; temp->Form("#chi^{2} / npts = %.1f / %d", Chi2, Npoints);
	labelchi2ref->SetLabel(temp->Data());
	delete temp;
	labelchi2ref->SetFillColor(kWhite);
	labelchi2ref->SetBorderSize(0);
	labelchi2ref->SetTextSize(Chi2Size);
	labelchi2ref->SetTextColor(fColorRef);
      }


      gDUnc->SetTitle("");
      gDUnc->SetMarkerStyle(20);
      gDUnc->SetMarkerSize(MarkerSize);
      
      gTheo->SetLineColor(color);
      if(gTheoRef) gTheoRef->SetLineColor(fColorRef);

      can->cd(i+1);
      if (Logx) gPad->SetLogx();
      if (Logy) gPad->SetLogy();

      h->Draw();
      gDUnc->Draw("same P");
      gDTot->Draw("same P");
      gTheo->Draw("same L");
      if(gTheoRef)     gTheoRef->Draw("same L");
      label->Draw();
      labelchi2->Draw();
      if(labelchi2ref) labelchi2ref->Draw();
    }
    can->cd(0);
    TPaveLabel* Title = new TPaveLabel(0.0, 0.97, 0.6, 1.0, dataset->GetName(), "NDC");
    Title->SetFillColor(kWhite); Title->SetBorderSize(0);
    Title->Draw();

    TLegend* leg = new TLegend(0.6, 0.97, 0.9, 1.0, "","NDC");
    leg->SetNColumns(2); leg->SetBorderSize(0); leg->SetFillColor(kWhite);
    leg->AddEntry(gTheo, fH1FitterOutput->GetName()->Data(),"L");
    if(fH1FitterOutputRef) leg->AddEntry(gTheoRef, fH1FitterOutputRef->GetName()->Data(),"L");
    leg->Draw();

    PrintCanvas(can);    
    delete can;
    delete TrashBin;
    delete Title; delete leg;
  }
  else { // GENERAL
    TGraphErrors* gDUnc = dataset->GetDataUncr(0);
    TGraphErrors* gDTot = dataset->GetDataTotal(0);
    TGraphErrors* gTheo = dataset->GetTheory(0);
    TGraphErrors* gTheoRef = NULL;
    if(datasetref) gTheoRef= datasetref->GetTheory(0);
    
    gDUnc->SetMarkerStyle(20);
    gDUnc->SetMarkerSize(0.5);
    gDUnc->SetTitle("");
    gDUnc->GetYaxis()->SetLabelSize(0.06);
    
    gDTot->SetMarkerStyle(1);
    gDTot->SetMarkerSize(0.05);
    
    gTheo->SetLineColor(color);
    if(gTheoRef) gTheoRef->SetLineColor(fColorRef);
    
    TH1F* h = new TH1F("","", gDUnc->GetN(), -0.5, gDUnc->GetN()-0.5);
    h->SetMaximum(gDUnc->GetHistogram()->GetMaximum()*1.1);
    h->SetMinimum(dataset->GetMinimum(0));
    h->SetStats(kFALSE);
    h->GetXaxis()->SetLabelSize(0.03);
    dataset->FillLabels(h);
    
    TCanvas* can = new TCanvas;
    can->SetLogy();
    h->Draw();
    gDUnc->Draw("same P");
    gDTot->Draw("same P");
    gTheo->Draw("same L");
    if(gTheoRef)     gTheoRef->Draw("same L");
    
    Double_t Chi2 = dataset->GetChi2(0);
    Int_t Npoints = dataset->GetNpts(0);
    
    TPaveLabel* labelchi2 = new TPaveLabel(0.25, 0.15, 0.3, 0.2,"","NDC"); 
    TString* temp = new TString; temp->Form("#chi^{2} / npts = %.1f / %d", Chi2, Npoints);
    labelchi2->SetLabel(temp->Data());
    delete temp;
    labelchi2->SetFillColor(kWhite);
    labelchi2->SetBorderSize(0);
    labelchi2->SetTextSize(1.0);
    labelchi2->SetTextColor(color);
    
    TPaveLabel* labelchi2ref = NULL;
    if(datasetref) {
      Chi2 = datasetref->GetChi2(0);
      Npoints = Npoints = datasetref->GetNpts(0);
      labelchi2ref = new TPaveLabel(0.25, 0.3, 0.3, 0.4,"","NDC"); 
      TString* temp = new TString; temp->Form("#chi^{2} / npts = %.1f / %d", Chi2, Npoints);
      labelchi2ref->SetLabel(temp->Data());
      delete temp;
      labelchi2ref->SetFillColor(kWhite);
      labelchi2ref->SetBorderSize(0);
      labelchi2ref->SetTextSize(0.5);
      labelchi2ref->SetTextColor(fColorRef);
    }


    TPaveLabel* Title = new TPaveLabel(0.0, 0.91, 0.5, 1.0, dataset->GetName(), "NDC");
    Title->SetFillColor(kWhite); Title->SetBorderSize(0);
    Title->Draw();
    
    TLegend* leg = new TLegend(0.6, 0.97, 0.9, 1.0, "","NDC");
    leg->SetNColumns(2); leg->SetBorderSize(0); leg->SetFillColor(kWhite);
    leg->AddEntry(gTheo, fH1FitterOutput->GetName()->Data(),"L"); 
    if(gTheoRef) leg->AddEntry(gTheoRef, fH1FitterOutputRef->GetName()->Data(),"L");
    leg->Draw();
    labelchi2->Draw();
    if(labelchi2ref) labelchi2ref->Draw();
    PrintCanvas(can);
    
    delete h;
    delete can;
    delete Title; delete leg; delete labelchi2;
    if(labelchi2ref) delete labelchi2ref;
  }
}


void H1FitterPainter::PrintCanvas(TCanvas* can) {
  can->Print(fPsFileName->Data());
  fPsFileName->ReplaceAll("(","");
//  TString* temp = new TString;
//  static Int_t idx = 0;
//  idx++;
//  temp->Form("DrawResults_%03d.eps",idx);
//  can->Print(temp->Data());
//  delete temp;
}
