#include "DataSet.h"
#include <TH1F.h>
#include <stdlib.h>
#include <TMath.h>

DataSet::DataSet(Int_t SetId, const Char_t* name, const Char_t* v1, const Char_t* v2, const Char_t* v3) {
  cout << SetId << " " << name<<endl;
  fName = new TString(name);
  fV1 = new TString(v1);
  fV2 = new TString(v2);
  fV3 = new TString(v3);
  fSetId = SetId;

  fDUnc = new TOrdCollection; fDUnc->SetOwner();
  fDTot = new TOrdCollection; fDTot->SetOwner();
  fTheo = new TOrdCollection; fTheo->SetOwner();
  fLabels = new TOrdCollection; fLabels->SetOwner();
  fDLabels.clear();
}
DataSet::DataSet() {
  fName = new TString;
  fV1 = new TString;
  fV2 = new TString;
  fV3 = new TString;
  fDUnc = new TOrdCollection; fDUnc->SetOwner();
  fDTot = new TOrdCollection; fDTot->SetOwner();
  fTheo = new TOrdCollection; fTheo->SetOwner();
  fLabels = new TOrdCollection; fLabels->SetOwner();
  fDLabels.clear();
}

DataSet::~DataSet(){
  delete fName;
  delete fV1;  delete fV2;  delete fV3;
  delete fDUnc;
  delete fDTot;
  delete fTheo;
  delete fLabels;
}

TGraphErrors* DataSet::GetDataUncr(Int_t i) {
  if(i<0 || i>=fDUnc->GetEntries()) return NULL;
  return (TGraphErrors*) fDUnc->At(i);
}
TGraphErrors* DataSet::GetDataTotal(Int_t i) {
  if(i<0 || i>=fDTot->GetEntries()) return NULL;
  return (TGraphErrors*) fDTot->At(i);
}
TGraphErrors* DataSet::GetTheory(Int_t i) {
  if(i<0 || i>=fTheo->GetEntries()) return NULL;
  return (TGraphErrors*) fTheo->At(i);
}
TObjString* DataSet::GetLabel(Int_t i) {
  if(i<0 || i>=fLabels->GetEntries()) return NULL;
  return (TObjString*) fLabels->At(i);
}

Int_t DataSet::GetNpts(Int_t i) {
  if(i<0 || i>=fTheo->GetEntries()) return 0;
  return ((TGraphErrors*) fDTot->At(i))->GetN();
}
Double_t DataSet::GetChi2(Int_t i) {
  if(i<0 || i>=fTheo->GetEntries()) return 0.;
  TGraphErrors* gDTot = (TGraphErrors*) fDUnc->At(i);
  TGraphErrors* gTheo = (TGraphErrors*) fTheo->At(i);
 
  Double_t* Data = gDTot->GetY();
  Double_t* Error = gDTot->GetEY();
  Double_t* Theo = gTheo->GetY();

  Double_t Chi2 = 0.;
  for(Int_t i = 0; i<gDTot->GetN(); i++) {
    Chi2 += TMath::Power((Data[i] - Theo[i]) / Error[i], 2.);
  }
  return Chi2;
}

Int_t DataSet::FindGraphIndex(Double_t value, const Char_t* label) {

  for(Int_t i=0; i<fDLabels.size(); i++) {
    if(fDLabels[i] < value) continue;
    if(fDLabels[i] > value) {
      fDLabels.insert(fDLabels.begin()+i,value);
      fDUnc->AddBefore(fDUnc->At(i),new TGraphErrors(0));
      fDTot->AddBefore(fDTot->At(i),new TGraphErrors(0));
      fTheo->AddBefore(fTheo->At(i),new TGraphErrors(0));
      
      fLabels->AddBefore(fLabels->At(i),new TObjString(label));
      return i;
    }
    else if(fDLabels[i] == value) {
      return i;
    }
  }
    
  fDLabels.push_back(value);
  fDUnc->AddLast(new TGraphErrors(0));
  fDTot->AddLast(new TGraphErrors(0));
  fTheo->AddLast(new TGraphErrors(0));
    
  fLabels->AddLast(new TObjString(label));
  return fDLabels.size() - 1;
}

void DataSet::AddPoint(Int_t GraphIdx, Double_t x, Double_t data, Double_t uncorrerr, Double_t toterr, Double_t theory) {
  TGraphErrors* gDUnc = (TGraphErrors*) fDUnc->At(GraphIdx);
  TGraphErrors* gDTot = (TGraphErrors*) fDTot->At(GraphIdx);
  TGraphErrors* gTheo = (TGraphErrors*) fTheo->At(GraphIdx);

  Int_t N = gDUnc->GetN();
  gDUnc->Set(N+1);
  gDTot->Set(N+1);
  gTheo->Set(N+1);
  
  gDUnc->SetPoint(N, x, data);
  gDUnc->SetPointError(N, 0., uncorrerr);
    
  gDTot->SetPoint(N, x, data);
  gDTot->SetPointError(N, 0., toterr);
    
  gTheo->SetPoint(N, x, theory);
  gTheo->SetPointError(N, 0., 0.);
  
  gDUnc->Sort();
  gDTot->Sort();
  gTheo->Sort();
}

void DataSet::AddPoint(Double_t v1, Double_t v2, Double_t v3, Double_t data, Double_t uncorrerr, Double_t toterr, Double_t theory, Double_t pull) {

  if(fSetId==61 || fSetId==62 || fSetId==63 || fSetId==64) {
    // different graph for each Q2 bin (column #2, v2)
    // graphs as a function of x (column #1, v1)
    // one needs to create the proper label too

    TString* temp = new TString;
    if(v2>10.) temp->Form("Q^{2} = %.0f GeV^{2}", v2);
    else       temp->Form("Q^{2} = %.1f GeV^{2}", v2);
    Int_t GraphIdx = FindGraphIndex(v2, temp->Data()); // find proper graph depending on the v2 value (or create if necessery)
    delete temp;
    
    AddPoint(GraphIdx, v1, data, uncorrerr, toterr, theory); // fill the graph as a function of v1 value

  }
  else if(fSetId==35) {
    // different graph for each Q2 bin (column #1 (v1) and column #2 (v2))
    // graphs as a function of x (column #3, v3)
    // one needs to create the proper label too
    
    TString* temp = new TString;
    temp->Form("%.0f < Q^{2} < %.0f GeV^{2}", v1, v2);
    Int_t GraphIdx = FindGraphIndex(v1, temp->Data()); // find proper graph depending on the v1 value (or create if necessery)
    delete temp;
    
    AddPoint(GraphIdx, v3, data, uncorrerr, toterr, theory); // fill the graph as a function of v1 value

  }
  else {                                // GENERAL:  

    TGraphErrors* gDUnc = NULL;
    TGraphErrors* gDTot = NULL;
    TGraphErrors* gTheo = NULL;

    if(fDUnc->GetEntries()==0) {
      gDUnc = new TGraphErrors(0);
      gDTot = new TGraphErrors(0);
      gTheo = new TGraphErrors(0);
      fDUnc->AddLast(gDUnc);
      fDTot->AddLast(gDTot);
      fTheo->AddLast(gTheo);
    }
    else {
      gDUnc = (TGraphErrors*) fDUnc->At(0);
      gDTot = (TGraphErrors*) fDTot->At(0);
      gTheo = (TGraphErrors*) fTheo->At(0);
    }
    
    Int_t N = gDUnc->GetN();
    gDUnc->Set(N+1);
    gDTot->Set(N+1);
    gTheo->Set(N+1);
    
    gDUnc->SetPoint(N, (Double_t)N, data);
    gDUnc->SetPointError(N, 0., uncorrerr);
    
    gDTot->SetPoint(N, (Double_t)N, data);
    gDTot->SetPointError(N, 0., toterr);
    
    gTheo->SetPoint(N, (Double_t)N, theory);
    gTheo->SetPointError(N, 0., 0.);
    
    TString* temp = new TString;
    temp->Form("%s:%.2f_%s:%.2f_%s:%.2f",fV1->Data(), v1, fV2->Data(), v2, fV3->Data(), v3);
    fLabels->AddLast(new TObjString(temp->Data()));
    delete temp;
  }
}

Double_t DataSet::GetMinimum(Int_t i) {
  Double_t Min = 99999.;
  TGraphErrors* gDUnc = (TGraphErrors*) fDUnc->At(i);
  TGraphErrors* gDTot = (TGraphErrors*) fDTot->At(i);
  TGraphErrors* gTheo = (TGraphErrors*) fTheo->At(i);
  Double_t x,y;

  for(Int_t i=0; i<gDUnc->GetN(); i++) {
    gDUnc->GetPoint(i, x, y);
    if(y < Min) Min = y;
    gDTot->GetPoint(i, x, y);
    if(y < Min) Min = y;
    gTheo->GetPoint(i, x, y);
    if(y < Min) Min = y;
  }
  return 0.5 * Min;
}

Int_t DataSet::FillLabels(TH1F* h) {
  for(Int_t i=0; i<fLabels->GetEntries(); i++)
    h->GetXaxis()->SetBinLabel(i+1, ((TObjString*)fLabels->At(i))->GetString().Data());
  return 0;
}

//void DataSet::PlotDataSet(TCanvas* can, const Char_t* name, DataSet* SetRef, const Char_t* nameRef) {
//  Float_t Chi2=0.; Int_t NPoints=0;
//  Float_t Chi2Ref=0.;
//  
//  if(SetRef) {
//    if(fSetId != SetRef->GetSetId()) {
//      Error("PlotDataSet","SetId mismatch with SetRef %d, %d", fSetId, SetRef->GetSetId());   return;
//    }
//  }
//  cout << "Plot Set " << fSetId << endl;
//  TH1F* h = NULL;
//  can->cd();
//  TPaveLabel* Label = NULL;
//  TPaveLabel* Label2 = NULL;
//  TLegend* Legend = NULL;
//  Int_t Nq2=0;
//  TString* temp = new TString;
//  Int_t n;
//  switch(fSetId) {
//  case 34:
//  case 35:
//  case 36:
//  case 38:
//  case 41:
//  case 42:
//    can->Divide(3, 2, 0., 0.);
//    for(Int_t iq2=0; iq2<6; iq2++) {
//      for(Int_t i=0; i<fData[iq2]->GetN(); i++) {
//	//cout << fData[iq2]->GetY()[i] << " " << fTheory[iq2]->GetY()[i] << " " << fData[iq2]->GetEYhigh()[i]<<endl;
//	Chi2 += TMath::Power((fData[iq2]->GetY()[i] - fTheory[iq2]->GetY()[i]),2.) / TMath::Power(fData[iq2]->GetEYhigh()[i],2.) ;
//	if(SetRef)
//	  Chi2Ref += TMath::Power((fData[iq2]->GetY()[i] - SetRef->GetTheory(iq2)->GetY()[i]),2.) / TMath::Power(fData[iq2]->GetEYhigh()[i],2.) ;
//	NPoints++;
//      }
//      if(fQ2Bins[iq2] == -1) {cout << "Wrong q2 bin, empty bin " << iq2<<endl; continue;}
//      if(iq2==0)  temp->Form("150 < Q^{2} < 200 GeV^{2}");
//      if(iq2==1)  temp->Form("200 < Q^{2} < 270 GeV^{2}");
//      if(iq2==2)  temp->Form("270 < Q^{2} < 400 GeV^{2}");
//      if(iq2==3)  temp->Form("400 < Q^{2} < 700 GeV^{2}");
//      if(iq2==4)  temp->Form("700 < Q^{2} < 5000 GeV^{2}");
//      if(iq2==5)  temp->Form("5000 < Q^{2} < 15000 GeV^{2}");
//      if(fSetId==41 || fSetId==42) { // ZEUS VALUES
//	if(iq2==0)  temp->Form(" 125 < Q^{2} < 250 GeV^{2} ");
//	if(iq2==1)  temp->Form(" 250 < Q^{2} < 500 GeV^{2} ");
//	if(iq2==2)  temp->Form(" 500 < Q^{2} < 1000 GeV^{2}");
//	if(iq2==3)  temp->Form("1000 < Q^{2} < 2000 GeV^{2}");
//	if(iq2==4)  temp->Form("2000 < Q^{2} < 5000 GeV^{2}");
//	if(iq2==5)  temp->Form("   Q^{2} > 5000 GeV^{2}    ");
//      }
//
//      if(iq2<3)   {Label = new TPaveLabel(0.25, 0.10, 0.65, 0.15, temp->Data(), "NDC");Label->SetTextSize(1.65);}
//      else        {Label = new TPaveLabel(0.4, 0.85, 0.8, 0.9,  temp->Data(), "NDC");Label->SetTextSize(1.3);}
//      Label->SetFillColor(kWhite); Label->SetBorderSize(0);
//      fTrash->AddLast(Label);
//
//      can->cd(iq2+1); 
//      if(!fNormalizeToTheory) {
//	gPad->SetLogy();
//	if(fSetId==34) {
//	  if(iq2<3) {fData[iq2]->SetMaximum(0.5); fData[iq2]->SetMinimum(0.0002);}
//	  else if(iq2<4)      {fData[iq2]->SetMaximum(0.05); fData[iq2]->SetMinimum(0.0002);}
//	  else if(iq2<5)      {fData[iq2]->SetMaximum(0.005); fData[iq2]->SetMinimum(0.00002);}
//	  else                {fData[iq2]->SetMaximum(0.00009); fData[iq2]->SetMinimum(0.000001);}
//	}
//	else if(fSetId==35 || fSetId==36) {
//	  if(iq2<3)      {fData[iq2]->SetMaximum(0.08); fData[iq2]->SetMinimum(0.00005);}
//	  else           {fData[iq2]->SetMaximum(0.2); fData[iq2]->SetMinimum(0.0002);}
//	}
//	else if(fSetId==38) {
//	  if(iq2<3)      {fData[iq2]->SetMaximum(0.02); fData[iq2]->SetMinimum(0.00002);}
//	  else           {fData[iq2]->SetMaximum(0.05); fData[iq2]->SetMinimum(0.00005);}
//	}
//	else if(fSetId==41 || fSetId==42) {
//	  if(iq2<3)      {fData[iq2]->SetMaximum(50.); fData[iq2]->SetMinimum(0.01);}
//	  else           {fData[iq2]->SetMaximum(5.); fData[iq2]->SetMinimum(0.002);}
//	}
//      }
//      else {
//	fData[iq2]->SetMaximum(1.5); 
//	fData[iq2]->SetMinimum(0.5);
//      }
//      fData[iq2]->SetMarkerStyle(20);
//      fData[iq2]->SetMarkerSize(0.5);
//      fData[iq2]->SetTitle("");
//      fData[iq2]->GetYaxis()->SetLabelSize(0.06);
//      fData[iq2]->GetXaxis()->SetLabelSize(0.06);
//      fData[iq2]->GetXaxis()->SetTitle("E_{T} [GeV]");
//      if(fSetId==38) fData[iq2]->GetXaxis()->SetTitle("<P_{T}> [GeV]");
//      fData[iq2]->GetXaxis()->SetTitleSize(0.05);
//      fData[iq2]->Draw("AP");
//      fDataUncorrErr[iq2]->Draw("same P");
//
//
//      TLegend* Legend1 = new TLegend(0.1, 0.65, 0.95, 0.95, "", "NDC");
//      Legend1->AddEntry(fTheoryModel[iq2], "Model/Parametrization","f");
//      Legend1->AddEntry(fTheoryMur[iq2], "Scale","f");
//      Legend1->SetTextSize(0.06);
//      Legend1->SetFillColor(kWhite);
//      Legend1->SetBorderSize(0);
//
//      fTheoryPar[iq2]->SetLineColor(kGreen); fTheoryPar[iq2]->SetFillColor(kGreen);
//      fTheoryModel[iq2]->SetLineColor(kBlue); fTheoryModel[iq2]->SetFillColor(kBlue);
//      fTheoryMur[iq2]->SetLineColor(kRed); fTheoryMur[iq2]->SetFillColor(kRed);
//
//
//      fTheoryMur[iq2]->Draw("3");
//      fTheoryModel[iq2]->Draw("3");
//    //fTheoryPar[iq2]->Draw("3");
//    
//      fData[iq2]->SetLineWidth(2);
//      fData[iq2]->SetMarkerSize(1.);
//
//
//      fData[iq2]->Draw("same P");
//      fDataUncorrErr[iq2]->Draw("same P");
//
//      if(fSetId==36 && iq2==2) Legend1->Draw();
//
////      fTheoryExp[iq2]->SetLineColor(kBlue); fTheoryExp[iq2]->SetFillColor(kBlue);
////      fTheoryExp[iq2]->Draw("3");
//
//
//
//      fTheory[iq2]->SetLineColor(kRed); fTheory[iq2]->SetFillColor(kWhite);
//      //fTheory[iq2]->Draw("L");
//      if(SetRef) {
//	SetRef->GetTheory(iq2)->SetLineColor(kBlue); SetRef->GetTheory(iq2)->SetFillColor(kWhite);
//	SetRef->GetTheory(iq2)->Draw("L");
//      }
//      Label->Draw();      
//    }
//    if (fSetId==34)  Label = new TPaveLabel(0.2, 0.985, 0.8, 1.005, "ISET 34: incl. jets 99/00, e+","NDC");
//    if (fSetId==35)  Label = new TPaveLabel(0.2, 0.985, 0.8, 1.005, "ISET 35: H1 normalized incl. high Q2 jets 99/00, e+","NDC");
//    if (fSetId==36)  Label = new TPaveLabel(0.2, 0.985, 0.8, 1.005, "ISET 36: H1 normalized incl. high Q2 jets 99-2007","NDC");
//    if (fSetId==38)  Label = new TPaveLabel(0.2, 0.985, 0.8, 1.005, "ISET 38: H1 normalized high Q2 di-jet 99-2007","NDC");
//    if (fSetId==41)  Label = new TPaveLabel(0.2, 0.985, 0.8, 1.005, "ISET 41: ZEUS inclusive high Q2 jets 96-97","NDC");
//    if (fSetId==42)  Label = new TPaveLabel(0.2, 0.985, 0.8, 1.005, "ISET 42: ZEUS inclusive high Q2 jets 90-00","NDC");
//    Label->SetFillColor(kWhite); Label->SetBorderSize(0);
//    can->cd(0); Label->Draw();
//
//    Legend = new TLegend(0.7, 0.9, 1.0, 1.0, "                  #chi^{2}/npts","NDC"); fTrash->AddLast(Legend);
//    Legend->SetFillColor(kWhite);    Legend->SetBorderSize(1);
//    if(SetRef) {
//      temp->Form("%.01f / %d (%s)", Chi2Ref, NPoints, nameRef);
//      Legend->AddEntry(SetRef->GetTheory(0), temp->Data());
//    }
//    temp->Form("%.01f / %d (%s)", Chi2, NPoints, name);
//    Legend->AddEntry(fTheory[0], temp->Data());
//    //Legend->Draw();
//    break;
//  case 37:
//    can->Divide(3, 3, 0.0, 0.);
//    for(Int_t iq2=0; iq2<7; iq2++) {
//      for(Int_t i=0; i<fData[iq2]->GetN(); i++) {
//	Chi2 += TMath::Power((fData[iq2]->GetY()[i] - fTheory[iq2]->GetY()[i]),2.) / TMath::Power(fData[iq2]->GetEYhigh()[i],2.) ;
//	if(SetRef)
//	  Chi2Ref += TMath::Power((fData[iq2]->GetY()[i] - SetRef->GetTheory(iq2)->GetY()[i]),2.) / TMath::Power(fData[iq2]->GetEYhigh()[i],2.) ;
//	NPoints++;
//      }
//      if(fQ2Bins[iq2] == -1) {cout << "Wrong q2 bin, empty bin " << iq2<<endl; continue;}
//      if(iq2==0)       temp->Form("5 < Q^{2} < 7 GeV^{2}");
//      else if(iq2==1)  temp->Form("7 < Q^{2} < 10 GeV^{2}");
//      else if(iq2==2)  temp->Form("10 < Q^{2} < 15 GeV^{2}");
//      else if(iq2==3)  temp->Form("15 < Q^{2} < 20 GeV^{2}");
//      else if(iq2==4)  temp->Form("20 < Q^{2} < 30 GeV^{2}");
//      else if(iq2==5)  temp->Form("30 < Q^{2} < 40 GeV^{2}");
//      else if(iq2==6)  temp->Form("40 < Q^{2} < 100 GeV^{2}");
//      
//
//      can->cd(iq2+1); 
//      if(!fNormalizeToTheory) {
//	gPad->SetLogy();       gPad->SetLogx();  // Logarithmic x-y scales
//	if(iq2<3)      {fData[iq2]->SetMinimum(0.02);      fData[iq2]->SetMaximum(20.);}
//	else if(iq2<6) {fData[iq2]->SetMinimum(0.007);      fData[iq2]->SetMaximum(60.);}
//	else           {fData[iq2]->SetMinimum(0.002);      fData[iq2]->SetMaximum(9.);}
//
//	//gPad->SetLogx();  // Logarithmic x scale
//	//if(iq2<3)      {fData[iq2]->SetMinimum(0.0);      fData[iq2]->SetMaximum(6.);}
//	//else if(iq2<6) {fData[iq2]->SetMinimum(0.0);      fData[iq2]->SetMaximum(20.);}
//	//else           {fData[iq2]->SetMinimum(0.0);      fData[iq2]->SetMaximum(5.);}
//      }
//      else {
//	fData[iq2]->SetMinimum(0.5);      
//	fData[iq2]->SetMaximum(1.5);
//      }
//	
//      if(iq2==6) Label = new TPaveLabel(0.25, 0.2, 0.65, 0.29, temp->Data(), "NDC");
//      else       Label = new TPaveLabel(0.2, 0.1, 0.6, 0.2, temp->Data(), "NDC");
//      Label->SetTextSize(1.1);      fTrash->AddLast(Label);
//      Label->SetFillColor(kWhite); Label->SetBorderSize(0);
//
//      fData[iq2]->GetXaxis()->SetLimits(6., 60.);
//      fData[iq2]->SetMarkerStyle(20);
//      fData[iq2]->SetMarkerSize(0.5);
//      fData[iq2]->SetTitle("");
//      fData[iq2]->GetYaxis()->SetLabelSize(0.08); 
//      fData[iq2]->GetXaxis()->SetLabelSize(0.08); 
//      fData[iq2]->GetXaxis()->SetTitle("E_{T} [GeV]");
//      fData[iq2]->GetXaxis()->SetTitleSize(0.06);
//      fData[iq2]->GetXaxis()->SetTitleOffset(0.8);
//      fData[iq2]->Draw("AP");
//      fDataUncorrErr[iq2]->Draw("same P");
//
//      fTheoryPar[iq2]->SetLineColor(kGreen); fTheoryPar[iq2]->SetFillColor(kGreen);
//      fTheoryModel[iq2]->SetLineColor(kBlue); fTheoryModel[iq2]->SetFillColor(kBlue);
//      fTheoryMur[iq2]->SetLineColor(kRed); fTheoryMur[iq2]->SetFillColor(kRed);
//
//
//      fTheoryMur[iq2]->Draw("3");
//      fTheoryModel[iq2]->Draw("3");
//    //fTheoryPar[iq2]->Draw("3");
//  fData[iq2]->SetLineWidth(2);
//      fData[iq2]->SetMarkerSize(1.);
//    
//      fData[iq2]->Draw("same P");
//      fDataUncorrErr[iq2]->Draw("same P");
//
//      fTheory[iq2]->SetLineColor(kRed); fTheory[iq2]->SetFillColor(kWhite);
//      //fTheory[iq2]->Draw("L");
//      if(SetRef) {
//	SetRef->GetTheory(iq2)->SetLineColor(kBlue); SetRef->GetTheory(iq2)->SetFillColor(kWhite);
//	SetRef->GetTheory(iq2)->Draw("L");
//      }
//      Label->Draw();      
//    }
//    Label = new TPaveLabel(0.2, 0.985, 0.8, 1.005, "ISET 37: H1 jets at low Q2 DIS, 99-00","NDC");
//    Label->SetFillColor(kWhite); Label->SetBorderSize(0);
//    can->cd(0); Label->Draw();
//
//    Legend = new TLegend(0.5, 0.1, 0.8, 0.25, "                  #chi^{2}/npts","NDC"); fTrash->AddLast(Legend);
//    Legend->SetFillColor(kWhite);    Legend->SetBorderSize(1);
//    if(SetRef) {
//      temp->Form("%.01f / %d (%s)", Chi2Ref, NPoints, nameRef);
//      Legend->AddEntry(SetRef->GetTheory(0), temp->Data());
//    }
//    temp->Form("%.01f / %d (%s)", Chi2, NPoints, name);
//    Legend->AddEntry(fTheory[0], temp->Data());
//    //Legend->Draw();
//    break;
//
//  case 39:
//
//    n = fTheory[0]->GetN();
//    //fData[0]->SetPoint(n-1, fData[0]->GetX()[n-1], fData[0]->GetY()[n-1]);
//    //fData[0]->SetPoint(n-1, (15000.-5000.)/2., fData[0]->GetY()[n-1]);
//    fData[0]->SetPoint(0, TMath::Exp((TMath::Log(150.)+TMath::Log(200.))/2.), fData[0]->GetY()[0]);
//    fData[0]->SetPoint(1, TMath::Exp((TMath::Log(200.)+TMath::Log(270.))/2.), fData[0]->GetY()[1]);
//    fData[0]->SetPoint(2, TMath::Exp((TMath::Log(270.)+TMath::Log(400.))/2.), fData[0]->GetY()[2]);
//    fData[0]->SetPoint(3, TMath::Exp((TMath::Log(400.)+TMath::Log(700.))/2.), fData[0]->GetY()[3]);
//    fData[0]->SetPoint(4, TMath::Exp((TMath::Log(700.)+TMath::Log(5000.))/2.), fData[0]->GetY()[4]);
//    fData[0]->SetPoint(5, TMath::Exp((TMath::Log(5000.)+TMath::Log(15000.))/2.), fData[0]->GetY()[5]);
//    for(Int_t i=0; i<6; i++) {
//      fDataUncorrErr[0]->SetPoint(i, fData[0]->GetX()[i], fDataUncorrErr[0]->GetY()[i]);
//      fTheory[0]->SetPoint(i, fData[0]->GetX()[i], fTheory[0]->GetY()[i]);
//      if(SetRef) SetRef->GetTheory(0)->SetPoint(i, fData[0]->GetX()[i], SetRef->GetTheory(0)->GetY()[i]);
//    }
//
//    for(Int_t iq2=0; iq2<1; iq2++) {
//      for(Int_t i=0; i<fData[iq2]->GetN(); i++) {
//	Chi2 += TMath::Power((fData[iq2]->GetY()[i] - fTheory[iq2]->GetY()[i]),2.) / TMath::Power(fData[iq2]->GetEYhigh()[i],2.) ;
//	if(SetRef)
//	  Chi2Ref += TMath::Power((fData[iq2]->GetY()[i] - SetRef->GetTheory(iq2)->GetY()[i]),2.) / TMath::Power(fData[iq2]->GetEYhigh()[i],2.) ;
//	NPoints++;
//      }
//      if(fQ2Bins[iq2] == -1) {cout << "Wrong q2 bin, empty bin " << iq2<<endl; continue;}
//      can->cd(1); 
//
//      gPad->SetLogy();      gPad->SetLogx();
//      fData[iq2]->SetMinimum(0.0001); 
//      fData[iq2]->SetMaximum(0.001); 
//
//      //fData[iq2]->GetXaxis()->SetLimits(6., 60.);
//      fData[iq2]->SetMarkerStyle(20);
//      fData[iq2]->SetMarkerSize(1.0);
//      fData[iq2]->SetTitle("");
//      fData[iq2]->GetYaxis()->SetLabelSize(0.05); 
//      fData[iq2]->GetXaxis()->SetLabelSize(0.05); 
//      fData[iq2]->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
//      fData[iq2]->GetXaxis()->SetTitleSize(0.05);
//      fData[iq2]->GetXaxis()->SetTitleOffset(0.8);
//      fData[iq2]->Draw("AP");
//      fDataUncorrErr[iq2]->Draw("same P");
//      fTheory[iq2]->SetLineColor(kRed); fTheory[iq2]->SetFillColor(kWhite);
//      fTheory[iq2]->Draw("L");
//      if(SetRef) {
//	SetRef->GetTheory(iq2)->SetLineColor(kBlue); SetRef->GetTheory(iq2)->SetFillColor(kWhite);
//	SetRef->GetTheory(iq2)->Draw("L");
//      }
//    }
//    Label = new TPaveLabel(0.2, 0.985, 0.8, 1.005, "ISET 39: H1 normalized high Q2 3-jets, 99-2007","NDC");
//    Label->SetFillColor(kWhite); Label->SetBorderSize(0);
//    can->cd(0); Label->Draw();
//    
//    Legend = new TLegend(0.5, 0.2, 0.8, 0.35, "                  #chi^{2}/npts","NDC"); fTrash->AddLast(Legend);
//    Legend->SetFillColor(kWhite);    Legend->SetBorderSize(1);
//    if(SetRef) {
//      temp->Form("%.01f / %d (%s)", Chi2Ref, NPoints, nameRef);
//      Legend->AddEntry(SetRef->GetTheory(0), temp->Data());
//    }
//    temp->Form("%.01f / %d (%s)", Chi2, NPoints, name);
//    Legend->AddEntry(fTheory[0], temp->Data());
//    //Legend->Draw();
//    break;
//
//  case 40:
//    can->Divide(3, 2, 0.0, 0.0);
//    h = new TH1F(" ","", 1, 0., 75.);
//    h->SetTitle("");
//    h->SetStats(kFALSE);
//    h->GetYaxis()->SetLabelSize(0.06);
//    h->GetXaxis()->SetLabelSize(0.06);
//    fTrash->AddLast(h);
//    Chi2 = 0.;
//    NPoints = 0;
//    for(Int_t iq2=0; iq2<6; iq2++) {
//      for(Int_t i=0; i<fData[iq2]->GetN(); i++) {
//	Chi2 += TMath::Power((fData[iq2]->GetY()[i] - fTheory[iq2]->GetY()[i]),2.) / TMath::Power(fData[iq2]->GetEYhigh()[i],2.) ;
//	if(SetRef)
//	  Chi2Ref += TMath::Power((fData[iq2]->GetY()[i] - SetRef->GetTheory(iq2)->GetY()[i]),2.) / TMath::Power(fData[iq2]->GetEYhigh()[i],2.) ;
//	NPoints++;
//      }
//      Float_t LabelShiftX = 0.;
//      if(iq2==0)                     temp->Form("-1 < #eta_{jet1} < 0  ");
//      if(iq2==1 || iq2==2)           temp->Form(" 0 < #eta_{jet1} < 1  ");
//      if(iq2==3 || iq2==4 || iq2==5) temp->Form(" 1 < #eta_{jet1} < 2.4");
//      if(iq2==0 || iq2==3) LabelShiftX = 0.1;
//      if(iq2==0 || iq2==1 || iq2==2) Label = new TPaveLabel(0.1+LabelShiftX, 0.2, 0.4+LabelShiftX, 0.25, temp->Data(), "NDC");
//      else                           Label = new TPaveLabel(0.1+LabelShiftX, 0.25, 0.4+LabelShiftX, 0.3, temp->Data(), "NDC");
//      Label->SetTextSize(1.6);
//      if(iq2>2)      Label->SetTextSize(1.4);
//      fTrash->AddLast(Label);
//      Label->SetFillColor(kWhite); Label->SetBorderSize(0);
//       
//      if(iq2==0 || iq2==1 || iq2==3) temp->Form("-1 < #eta_{jet2} < 0  ");
//      if(iq2==2 || iq2==4)           temp->Form(" 0 < #eta_{jet2} < 1  ");
//      if(iq2==5)                     temp->Form(" 1 < #eta_{jet2} < 2.4");
//
//      if(iq2==0 || iq2==1 || iq2==2) Label2 = new TPaveLabel(0.1+LabelShiftX, 0.1, 0.4+LabelShiftX,0.15, temp->Data(), "NDC");
//      else                           Label2 = new TPaveLabel(0.1+LabelShiftX, 0.15, 0.4+LabelShiftX, 0.2, temp->Data(), "NDC");
//      Label2->SetTextSize(1.6);
//      if(iq2>2)      Label2->SetTextSize(1.4);
//      fTrash->AddLast(Label2);
//      Label2->SetFillColor(kWhite); Label2->SetBorderSize(0);
//       
//      can->cd(iq2+1); 
//
//      if(!fNormalizeToTheory) {
//	gPad->SetLogy();
//	h->SetMinimum(0.01); h->SetMaximum(50.);
//      }
//      else {
//	h->SetMinimum(0.5); h->SetMaximum(1.5);
//      }
//
//      h->Draw("");
//      fData[iq2]->SetMarkerStyle(20);
//      fData[iq2]->SetMarkerSize(0.5);
//      fData[iq2]->SetTitle("");
//      fData[iq2]->GetYaxis()->SetLabelSize(0.06);
//      fData[iq2]->Draw("P");
//      fTheory[iq2]->SetLineColor(kRed); fTheory[iq2]->SetFillColor(kWhite);
//      fTheory[iq2]->SetMarkerColor(kRed);
//      fTheory[iq2]->Draw("L");
//      if(SetRef) {
//	SetRef->GetTheory(iq2)->SetLineColor(kBlue); SetRef->GetTheory(iq2)->SetFillColor(kWhite);
//	SetRef->GetTheory(iq2)->Draw("L");
//      }
//      Label->Draw();      
//      Label2->Draw();      
//    }
//    Label = new TPaveLabel(0.2, 0.985, 0.8, 1.005, "ISET 40: ZEUS photoproduction dijet 96-97","NDC");
//    Label->SetFillColor(kWhite); Label->SetBorderSize(0);
//    can->cd(0); Label->Draw();
//
//    Legend = new TLegend(0.5, 0.55, 0.8, 0.7, "                  #chi^{2}/npts","NDC"); fTrash->AddLast(Legend);
//    Legend->SetFillColor(kWhite);    Legend->SetBorderSize(1);
//    if(SetRef) {
//      temp->Form("%.01f / %d (%s)", Chi2Ref, NPoints, nameRef);
//      Legend->AddEntry(SetRef->GetTheory(0), temp->Data());
//    }
//    temp->Form("%.01f / %d (%s)", Chi2, NPoints, name);
//    Legend->AddEntry(fTheory[0], temp->Data());
//    //Legend->Draw();
//
//    break;
//
//
//  default: Error("PlotDataSet","Set %d not implemented yet!", fSetId);
//  }
//  TPaveLabel* Name = new TPaveLabel(0.0, 0.985, 0.15, 1.005, name); Name->SetFillColor(kWhite); Name->SetBorderSize(0);
//  can->cd(0); Name->Draw();
//  delete temp;
//}
//void DataSet::SortTheory() {
//  for(Int_t iq2=0; iq2<fNMaxQ2Bins; iq2++) {
//    if(fTheory[iq2]) {
//      if(fSetId != 34 && fSetId != 35 && fSetId!=36 && fSetId != 37 && fSetId != 38 && fSetId != 39 && fSetId!=40 && fSetId!=41 && fSetId!=42) { // jets are dependent on Et, not x
//	Int_t n = fTheory[iq2]->GetN();
//	fTheory[iq2]->SetPoint(n-1, 100., 0.);
//	fTheory[iq2]->SetPointError(n, 0., 0., 0., 0.);
//      }
//      else if (fSetId == 39) { 
//	//	Int_t n = fTheory[iq2]->GetN();
//      }
//      else {  // add the 0 at Et==100
//	Int_t n = fTheory[iq2]->GetN();
//	fTheory[iq2]->SetPoint(n, 100., 0.);
//	fTheory[iq2]->SetPointError(n, 0., 0., 0., 0.);
//	if(fNormalizeToTheory) {
//	  fTheory[iq2]->SetPoint(n, 100., 1.);
//	}
//      }
//      fTheory[iq2]->Sort();
//    }
//  }
//}

