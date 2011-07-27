#include <H1FitterPainter.h>
#include <TPaveLabel.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TH1F.h>

H1FitterPainter::H1FitterPainter(){
  fPath = new TString("../../output/");
  fPathRef = new TString("../../output/");
  fH1FitterOutput = NULL;
  fH1FitterOutputRef = NULL;
  fPsFileName = new TString("DrawResults.ps");
  gROOT->SetStyle("Plain");
  cout << endl;
  cout << "TO DO: pdfs_q2val_*.txt should provide the information about Q2"<<endl;
  cout << "TO DO: in fittedresults.txt q2 and x for sets 61-64 are switched"<<endl;
}

H1FitterPainter::~H1FitterPainter(){ 
  delete fPath;
  delete fPathRef;
  if(fH1FitterOutput) delete fH1FitterOutput;
  if(fH1FitterOutputRef) delete fH1FitterOutputRef;

  TCanvas* can = new TCanvas;
  fPsFileName->Append(")");
  can->Print(fPsFileName->Data());
  delete can;

  delete fPsFileName;
}

Int_t H1FitterPainter::Prepare() {
  fPsFileName->Append("(");
  fH1FitterOutput = new H1FitterOutput(fPath->Data());
  if(fPath->CompareTo(fPathRef->Data())) 
    fH1FitterOutputRef = new H1FitterOutput(fPathRef->Data());
  else 
    fH1FitterOutputRef = NULL;

  fH1FitterOutput->Prepare();
  if(fH1FitterOutputRef) fH1FitterOutputRef->Prepare(); 
}

Int_t H1FitterPainter::Draw() {
  this->Prepare();
  for(Int_t i=0; i<5; i++)
    this->DrawPDF(i);
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
  }
  for (Int_t iref=0; iref<NRefDataSets; iref++) {
    if(!RefSetsDrawn[iref])
      DrawDataSet(fH1FitterOutputRef->GetSet(iref), NULL, kBlue);
  }
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

  if(ival==0) Q2Label->SetLabel("Q^{2} = 1.9 GeV^{2}");  // THIS SHOULD BE AUTOMATIC!!!
  if(ival==1) Q2Label->SetLabel("Q^{2} = 4.0 GeV^{2}");
  if(ival==2) Q2Label->SetLabel("Q^{2} = 10. GeV^{2}");
  if(ival==3) Q2Label->SetLabel("Q^{2} = 100. GeV^{2}");
  if(ival==4) Q2Label->SetLabel("Q^{2} = 6464 GeV^{2}");
  if(ival==5) Q2Label->SetLabel("Q^{2} = 8317 GeV^{2}");

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

  Int_t Color = kRed;
  Int_t RColor = kBlue;
  Int_t Color_Ratio = kBlack;
  Int_t Style_Ratio = 2;

  pad->cd(); pad->SetLogx();

  TGraph* graph =   FitterOut->GetPdf(pdf1, Q2Bin);
  TGraph* graphR = NULL; if(FitterRef) graphR = FitterRef->GetPdf(pdf1, Q2Bin);

  TGraph* ref_ratio = (TGraph*) (graph->Clone()); TrashBin->AddLast(ref_ratio);
  TGaxis* AxisRatio = NULL;

  if(!graph) {cout << "Graph = NULL, return" << endl; return 1;}
  graph->SetTitle("");
  graph->GetXaxis()->SetTitle("x");   graph->GetYaxis()->SetTitle("xP(x)");
  graph->GetXaxis()->SetTitleOffset(1.2);  graph->GetYaxis()->SetTitleOffset(1.5);  
  graph->GetXaxis()->SetLabelSize(0.05); graph->GetXaxis()->SetLabelOffset(0.02);
  graph->GetYaxis()->SetLabelSize(0.05); graph->GetYaxis()->SetNdivisions(506);  graph->GetYaxis()->SetLabelOffset(0.02);
  graph->SetLineColor(Color);
  if(graphR) {
    graphR->SetLineColor(RColor);
  }

  TGraph* graph2 = NULL;
  TGraph* graphR2 = NULL;
  TLine* line = NULL;
  if(pdf2 != H1FitterOutput::kNULL) {
    graph2 =   FitterOut->GetPdf(pdf2, Q2Bin);
    if(FitterRef) graphR2 =  FitterRef->GetPdf(pdf2, Q2Bin);
    graph2->SetLineColor(Color);
    if(FitterRef)  graphR2->SetLineColor(RColor);
  }
 

  if(graphR) {
    for(Int_t i=0; i<graph->GetN(); i++) {
      ref_ratio->SetPoint(i, graph->GetX()[i], graph->GetY()[i] / graphR->GetY()[i]);
    }    
    ScaleGraph2ToGraph1(graph,    ref_ratio,   line, AxisRatio, 1.0);
    TrashBin->AddLast(line); TrashBin->AddLast(AxisRatio);
    ref_ratio->SetLineColor(Color_Ratio);  ref_ratio->SetLineStyle(Style_Ratio);  ref_ratio->SetFillColor(Color_Ratio);
    line->SetLineColor(Color_Ratio);       line->SetLineStyle(3);
  }

  graph->Draw("ALX");
  if(graphR) graphR->Draw("L3X same");
  if(graph2)  graph2->Draw("L3X same");
  if(graphR2) graphR2->Draw("L3X same");

  if(graphR) {
    ref_ratio->Draw("L same");
    AxisRatio->Draw("C same");
    line->Draw();
  }

  TPaveLabel* label = new TPaveLabel(0.49, 0.85, 0.51, 0.87, Title,"NDC"); TrashBin->AddLast(label);
  label->SetBorderSize(0); label->SetFillColor(kWhite); label->SetTextSize(3.0);
  pad->cd();
  label->Draw();

  TPaveLabel* box = new TPaveLabel(0.0, 0.75, 0.08, 0.91, "xP(x)", "NDC"); TrashBin->AddLast(box);
  box->SetFillColor(kWhite); box->SetBorderSize(0); box->SetTextAngle(90.); box->SetTextSize(0.35);
  box->Draw();
 
  if(legend) {
    legend->cd();
    TPaveLabel* lab1 = new TPaveLabel(0., 0.4, 1.0, 0.5, FitterOut->GetDirectory()->Data(), "NDC");
    TrashBin->AddLast(lab1); lab1->SetFillColor(kWhite); lab1->SetBorderSize(0); lab1->SetTextColor(Color);
    lab1->Draw();
    if(FitterRef) {
      TPaveLabel* lab2 = new TPaveLabel(0., 0.6, 1.0, 0.7, FitterRef->GetDirectory()->Data(), "NDC");
      TrashBin->AddLast(lab2); lab2->SetFillColor(kWhite); lab2->SetBorderSize(0); lab2->SetTextColor(RColor);
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


Int_t H1FitterPainter::DrawDataSet(DataSet* dataset, DataSet* datasetref, EColor color) {

  if(dataset->GetNGraphs()==0) return 1;

  Int_t SetId = dataset->GetSetId();
  if(SetId==61 || SetId==62 || SetId==63 || SetId==64) {

    TObjArray* TrashBin = new TObjArray; TrashBin->SetOwner();
    Double_t Maximum = 0.; Double_t Minimum = 0.;
    Double_t XMinimum = 0.;
    Double_t MarkerSize = 0.005;
    Double_t LabelSize = 0.06;
    TCanvas* can = new TCanvas;
    can->SetTopMargin(0.3);
    switch(SetId) {
    case 61: can->Divide(4,5,0.); XMinimum = 0.001;   Maximum = 1.99; Minimum = 0.0001; MarkerSize = 0.005; LabelSize = 0.06; break;
    case 62: can->Divide(6,6,0.); XMinimum = 0.00003; Maximum = 1.59; Minimum = 0.0001; MarkerSize = 0.005; LabelSize = 0.06; break;
    case 63: can->Divide(3,4,0.); XMinimum = 0.002;   Maximum = 0.99; Minimum = 0.0001; MarkerSize = 0.5;   LabelSize = 0.08; break;
    case 64: can->Divide(3,3,0.); XMinimum = 0.002;   Maximum = 1.99; Minimum = 0.0001; MarkerSize = 0.5;   LabelSize = 0.08; break;
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

      //if(SetId==63 && i>8) {Maximum = 0.000499; Minimum = 0.000001;}

      Double_t chi2offx = 0.;
      Double_t chi2offy = 0.;

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


      if(SetId==63 && i>2) Maximum = 0.299;
      if(SetId==63 && i>5) {Maximum = 0.0399;; Minimum = 0.000001;}
      if(SetId==63 && i>8) {Maximum = 0.000499; Minimum = 0.000001;}

      if(SetId==64 && i>2) Maximum = 0.1199;
      if(SetId==64 && i>5) Maximum = 0.01599;


      TH1F* h = new TH1F("","",1, XMinimum, 1.); TrashBin->SetOwner(h);
      h->SetMinimum(Minimum);
      h->SetMaximum(Maximum);
      h->SetStats(kFALSE);
      h->GetYaxis()->SetLabelSize(LabelSize);
      h->GetXaxis()->SetLabelSize(LabelSize);


      TObjString* lstring = dataset->GetLabel(i);
      TPaveLabel* label = new TPaveLabel(0.4, 0.7, 0.9, 0.9,"null","NDC"); TrashBin->AddLast(label);
      if(lstring)	label->SetLabel(lstring->GetString().Data());
      label->SetFillColor(kWhite);
      label->SetBorderSize(0);
      label->SetTextSize(0.6);

      Double_t Chi2 = dataset->GetChi2(i);
      Int_t Npoints = dataset->GetNpts(i);


      TPaveLabel* labelchi2 = new TPaveLabel(0.25+chi2offx, 0.1+chi2offy, 0.3+chi2offx, 0.2+chi2offy,"","NDC"); TrashBin->AddLast(labelchi2);
      TString* temp = new TString; temp->Form("#chi^{2} / npts = %.1f / %d", Chi2, Npoints);
      labelchi2->SetLabel(temp->Data());
      delete temp;
      labelchi2->SetFillColor(kWhite);
      labelchi2->SetBorderSize(0);
      labelchi2->SetTextSize(1.0);
      labelchi2->SetTextColor(kRed);

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
	labelchi2ref->SetTextSize(1.0);
	labelchi2ref->SetTextColor(kBlue);
      }


      gDUnc->SetTitle("");
      gDUnc->SetMarkerStyle(20);
      gDUnc->SetMarkerSize(MarkerSize);
      
      gTheo->SetLineColor(color);
      if(gTheoRef) gTheoRef->SetLineColor(kBlue);

      can->cd(i+1);
      gPad->SetLogx();

      //gDUnc->GetHistogram()->GetXaxis()->SetRangeUser(0.0001, 1.);
      
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
    leg->SetNColumns(2); leg->SetBorderSize(0);
    leg->AddEntry(gTheo, fH1FitterOutput->GetDirectory()->Data(),"L");
    if(fH1FitterOutputRef) leg->AddEntry(gTheoRef, fH1FitterOutputRef->GetDirectory()->Data(),"L");
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
    if(gTheoRef) gTheoRef->SetLineColor(kBlue);
    
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
      labelchi2ref->SetTextSize(1.5);
      labelchi2ref->SetTextColor(kBlue);
    }


    TPaveLabel* Title = new TPaveLabel(0.0, 0.91, 0.5, 1.0, dataset->GetName(), "NDC");
    Title->SetFillColor(kWhite); Title->SetBorderSize(0);
    Title->Draw();
    
    TLegend* leg = new TLegend(0.6, 0.97, 0.9, 1.0, "","NDC");
    leg->SetNColumns(2); leg->SetBorderSize(0);
    leg->AddEntry(gTheo, fH1FitterOutput->GetDirectory()->Data(),"L");
    if(gTheoRef) leg->AddEntry(gTheoRef, fH1FitterOutputRef->GetDirectory()->Data(),"L");
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
  TString* temp = new TString;
  static Int_t idx = 0;
  idx++;
  temp->Form("EpsFiles/DrawResults_%03d.eps",idx);
  can->Print(temp->Data());
  delete temp;
}
