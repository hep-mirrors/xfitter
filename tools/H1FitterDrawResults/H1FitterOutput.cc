#include "H1FitterOutput.h"
#include <fstream>
#include <TObjString.h>

H1FitterOutput::H1FitterOutput(const Char_t* directory) {
  fDirectory    = new TString(directory);
  for(Int_t ipdf = 0; ipdf < fNpdfs; ipdf++) {
    fPdfs[ipdf] = new TObjArray;
    fPdfs[ipdf]->SetOwner();
  }
  fNDataSets = 0;
  fPull = new TH1F("","",20, -3., 3.);
}

H1FitterOutput::~H1FitterOutput(){
  delete fDirectory;
  for(Int_t ipdf = 0; ipdf < fNpdfs; ipdf++) {
    delete fPdfs[ipdf];
  }
  for (std::vector<DataSet*>::iterator i = fDataSets.begin(); i != fDataSets.end(); ++i)  delete *i;
  delete fPull;
}

Int_t H1FitterOutput::Prepare() {
  this->PreparePdf();
  this->PrepareDataSets();
}


Int_t H1FitterOutput::PreparePdf() {
  Double_t x, gluon, U, D, d_Ubar, d_Dbar, umin, dmin, sea, u_sea, d_sea, str, chm, bot;
  TString* filename = new TString;
  for(Int_t iq2=0; iq2<100; iq2++) {
    filename->Form("%s/pdfs_q2val_%02d.txt",fDirectory->Data(), iq2+1);
    ifstream infile(filename->Data());
    if(!infile.is_open()) break;

    for(Int_t ipdf = 0; ipdf < fNpdfs; ipdf++)
      ((TObjArray*) fPdfs[ipdf])->AddLast(new TGraph(fNpoints));
    
    for (Int_t i = 0; i < fNpoints; i++){
      infile >> x >> gluon >> U >> D >> d_Ubar >> d_Dbar >> umin >> dmin >> sea >> u_sea >> d_sea >> str >> chm >> bot;
      SetPdfPoint((Int_t)kGluon, iq2, i, x, gluon);
      SetPdfPoint((Int_t)kU    , iq2, i, x, U);
      SetPdfPoint((Int_t)kD    , iq2, i, x, D);
      SetPdfPoint((Int_t)kUv   , iq2, i, x, U - u_sea - chm);
      SetPdfPoint((Int_t)kDv   , iq2, i, x, D - d_sea - str - bot);
      SetPdfPoint((Int_t)kUb   , iq2, i, x, d_Ubar);
      SetPdfPoint((Int_t)kDb   , iq2, i, x, d_Dbar);
      SetPdfPoint((Int_t)kSea  , iq2, i, x, sea);
      SetPdfPoint((Int_t)kS    , iq2, i, x, str);
      SetPdfPoint((Int_t)kC    , iq2, i, x, chm);
      SetPdfPoint((Int_t)kB    , iq2, i, x, bot);
    }
  }  delete filename;
}

Int_t H1FitterOutput::GetNsets() {
  return fDataSets.size();
}

TGraph* H1FitterOutput::GetPdf(H1FitterOutput::pdf ipdf, Int_t Q2bin) {
  if(ipdf >= fNpdfs) {cout << "GetPdf, wrong ipdf: "<< ipdf << endl; return NULL;}
  if(Q2bin >= fPdfs[ipdf]->GetEntries()) {cout << "GetPdf, wrong iq2: "<< Q2bin << endl; return NULL;}
  return ((TGraph*)fPdfs[ipdf]->At(Q2bin));
}


void H1FitterOutput::SetPdfPoint(Int_t ipdf, Int_t iq2, Int_t ipoint, Double_t x, Double_t y) {
  if(ipdf >= fNpdfs) {cout << "SetPdfPoint, wrong ipdf: "<< ipdf << endl; return;}
  if(iq2 >= fPdfs[ipdf]->GetEntries()) {cout << "SetPdfPoint, wrong iq2: "<< iq2 << endl; return;}
  ((TGraph*)fPdfs[ipdf]->At(iq2))->SetPoint(ipoint, x, y);
}


Int_t H1FitterOutput::PrepareDataSets() {

  TString* filename = new TString;
  filename->Form("%s/fittedresults.txt",fDirectory->Data());

  Double_t v1, v2, v3, data, uncorrerr, toterr, theory, pull;
  Int_t dataset;
  Char_t buffer[120];
  TString V1, V2, V3, str, Name;

  ifstream infile(filename->Data());
  if(!infile.is_open()) { cout << "H1FitterOutput::PrepareDataSets: can not open file %s" << filename->Data()<<endl; return 1;}
  infile.getline(buffer, 120);  // number of data sets
  str.Form(buffer); 
  fNDataSets = str.Atoi();
  //fDataSets
  //  fDataSets = new DataSet[fNDataSets];

  infile.getline(buffer, 120);
  str.Form(buffer); 
  dataset = str.Atoi();
  // v1     v2    v3    data     +- uncorr.err   +-toterr      theory      pull     dataset
  while(!infile.eof()) {
    infile.getline(buffer, 120);
    Name.Form(buffer);
    infile.getline(buffer, 120);
    infile.getline(buffer, 120); 
    str.Form(buffer);
    TObjArray* array = str.Tokenize(" ");
    if(array->GetEntries() < 3) {cout << "something is wrong in fittedresults.txt" << endl; delete array; continue;}
    
    V1.Form(((TObjString*)array->At(0))->GetString().Data());
    V2.Form(((TObjString*)array->At(1))->GetString().Data());
    V3.Form(((TObjString*)array->At(2))->GetString().Data());
    delete array;

    DataSet* NewDataSet = new DataSet(dataset, Name.Data(), V1.Data(), V2.Data(), V3.Data());
    fDataSets.push_back(NewDataSet);
    while(1) {
      infile.getline(buffer, 120);
      str.Form(buffer);
      TObjArray* array = str.Tokenize(" ");
      if(array->GetEntries() == 1) {dataset = ((TObjString*)array->At(0))->GetString().Atoi(); delete array; break;}
      if(array->GetEntries() != 9) {delete array; break;}
      
      v1        = ((TObjString*)array->At(0))->GetString().Atof();
      v2        = ((TObjString*)array->At(1))->GetString().Atof();
      v3        = ((TObjString*)array->At(2))->GetString().Atof();
      data      = ((TObjString*)array->At(3))->GetString().Atof();
      uncorrerr = ((TObjString*)array->At(4))->GetString().Atof();
      toterr    = ((TObjString*)array->At(5))->GetString().Atof();
      theory    = ((TObjString*)array->At(6))->GetString().Atof();
      pull      = ((TObjString*)array->At(7))->GetString().Atof();
      dataset   = ((TObjString*)array->At(8))->GetString().Atoi();	
      delete array;
      NewDataSet->AddPoint(v1, v2, v3, data, uncorrerr, toterr, theory, pull);

      fPull->Fill(pull);
    }
  }    
  
  delete filename;
  return 0;
}






//
//Int_t H1FitterOutput::PrepareAlphas() {
//  TString* filename = new TString;
//  filename->Form("%s/minuit.out.txt",fDirectory->Data());
//  fAlphas = GetAlphas(filename->Data());
//  delete filename;
//}
//
//Int_t H1FitterOutput::PrepareErrors() {
//  AddError("FStd","FStu", kFSError, kQuadrature, kModelError);
//  AddError("MBd","MBu", kMBError, kQuadrature, kModelError);
//  AddError("MCd","MCu", kMCError, kQuadrature, kModelError);
//  AddError("Q2md","Q2mu", kQ2minError, kQuadrature, kModelError);
//  AddError("HadCd","HadCu", kHadronizationError, kQuadrature, kModelError);
//  //AddError("Pdf","Pdf", kPDF, kQuadrature, kModelError);
//  
//  //AddError("FStd","FStd", kFSErrorD);
//  //AddError("MBd","MBd", kMBErrorD);
//  //AddError("MCd","MCd", kMCErrorD);
//  //AddError("Q20d","Q20d", kQ2_0ErrorD);
//  //AddError("Q2md","Q2md", kQ2minErrorD);
//  //AddError("HadCd","HadCd", kHadronizationErrorD);
//  //
//  //AddError("FStu","FStu", kFSErrorU);
//  //AddError("MBu","MBu", kMBErrorU);
//  //AddError("MCu","MCu", kMCErrorU);
//  //AddError("Q20u","Q20u", kQ2_0ErrorU);
//  //AddError("Q2mu","Q2mu", kQ2minErrorU);
//  //AddError("HadCu","HadCu", kHadronizationErrorU);
//  //
//  //AddError("Mur1u", "Mur1u",  kMur1U);
//  //AddError("Mur2u", "Mur2u",  kMur2U);
//  //AddError("Mur3u", "Mur3u",  kMur3U);
//  //AddError("Mur23u","Mur23u", kMur23U);
//  //AddError("Mur123u","Mur123u", kMur123U);
//  //
//  //AddError("Mur1d", "Mur1d",  kMur1D);
//  //AddError("Mur2d", "Mur2d",  kMur2D);
//  //AddError("Muf3d", "Muf3d",  kMur3D);
//  //AddError("Mur23d","Mur23d", kMur23D);
//  //AddError("Mur123d","Mur123d", kMur123D);
//  //
//  //AddError("Muf1u", "Muf1u",  kMuf1U);
//  //AddError("Muf2u", "Muf2u",  kMuf2U);
//  //AddError("Muf3u", "Muf3u",  kMuf3U);
//  //AddError("Muf23u","Muf23u", kMuf23U);
//  //AddError("Muf123u","Muf123u", kMuf123U);
//  //
//  //AddError("Muf1d", "Muf1d",  kMuf1D);
//  //AddError("Muf2d", "Muf2d",  kMuf2D);
//  //AddError("Muf3d", "Muf3d",  kMuf3D);
//  //AddError("Muf23d","Muf23d", kMuf23D);
//  //AddError("Muf123d","Muf123d", kMuf123D);
//  
//  //AddError("Mur1u", "Mur1d",  kMur1);
//  //AddError("Mur1u", "Mur1d",  kMur1, kMaximal, kMurError);
//  AddError("Mur2u", "Mur2d",  kMur2, kMaximal, kMurError);
//  //AddError("Mur3u", "Mur3d",  kMur3, kMaximal, kMurError);
//  //AddError("Mur23u","Mur23d", kMur23, kMaximal, kMurError);
//  //AddError("Mur123u","Mur123d", kMur123);
//  //AddError("Mur123u","Mur123d", kMur123, kMaximal, kMurError);
//  
//  //AddError("Muf1u", "Muf1d",  kMuf1);
//  //AddError("Muf1u", "Muf1d",  kMuf1, kMaximal, kMufError);
//  AddError("Muf2u", "Muf2d",  kMuf2, kMaximal, kMufError);
//  // AddError("Muf3u", "Muf3d",  kMuf3, kMaximal, kMufError);
//  //AddError("Muf23u","Muf23d", kMuf23, kMaximal, kMufError);
//  //AddError("Muf123u","Muf123d", kMuf123);
//  //AddError("Muf123u","Muf123d", kMuf123, kMaximal, kMufError);
//  
////  AddError("Ddbar","Ddbar",k9PEuv, kMaximal, kParamError);
////  AddError("Dubar","Dubar",k9PDuv, kMaximal, kParamError);
////  AddError("Ddv","Ddv",k9PDdv, kMaximal, kParamError);
////  AddError("Dg","Dg",k9PDg, kMaximal, kParamError);
//
//
////  AddError("9P","9P",k9P);
////  AddError("9PEuv","9PEuv",k9PEuv);
////  AddError("9PDuv","9PDuv",k9PDuv);
////  AddError("9PDdv","9PDdv",k9PDdv);
////  AddError("9PDg","9PDg",k9PDg);
////  AddError("9PEuvDuv","9PEuvDuv",k9PEuvDuv, kMaximal, kParamError);
////  AddError("9PEuvDdv","9PEuvDdv",k9PEuvDdv, kMaximal, kParamError);
////  AddError("9PEuvDg","9PEuvDg",k9PEuvDg, kMaximal, kParamError);
////  //AddError("9PEuv221","9PEuv221",k9PEuv221, kMaximal, kParamError);
////  //AddError("9PEuv225","9PEuv225",k9PEuv225, kMaximal, kParamError);
////  AddError("9PEuvSqrt","9PEuvSqrt",k9PEuvSqrt);
////  AddError("9PEuvDuvSqrt","9PEuvDuvSqrt",k9PEuvDuvSqrt);
////  AddError("9PEuvDdvSqrt","9PEuvDdvSqrt",k9PEuvDdvSqrt);
////  AddError("9PEuvDgSqrt","9PEuvDgSqrt",k9PEuvDgSqrt);
////  AddError("9PEuvNg","9PEuvNg",k9PEuvNg, kMaximal, kParamError);
//  AddError("Q20d","Q20u", kQ2_0Error, kMaximal, kParamError);
//
//  AddExpError();
//
//  AddError(kFullError, kExpError);
//  AddError(kFullError, kModelError);
//  AddError(kFullError, kParamError);
//  AddError(kFullError, kMurError);
//  AddError(kFullError, kMufError);
//  //AddBandsExpError();
//}
//
//Int_t H1FitterOutput::AddMessage(const Char_t* text) {
//  for(Int_t i=0; i<fMessages->GetEntries(); i++) {
//    if(!((TObjString*) fMessages->At(i))->GetString().CompareTo(text)) return 1;
//  }
//  fMessages->AddLast(new TObjString(text));
//  return 0;
//}
//
//
//Int_t H1FitterOutput::CheckJob(const Char_t* MinuitOut) {
//
//  TString JobName(MinuitOut);
//  if(!JobName.CompareTo((TString(fDirectory->Data())+TString("/minuit.out.txt")).Data())) 
//    JobName.Form("Central Job");
//  else {
//    JobName.ReplaceAll(fDirectory->Data(), "");
//    JobName.ReplaceAll("/minuit.out.txt","");
//    JobName.ReplaceAll("_","");
//  }
//  TString temp;
//  ifstream infile(MinuitOut);
//  if(!infile.is_open()){ 
//    temp.Form("%s - can not find minuit.out.txt file", JobName.Data());
//    this->AddMessage(temp.Data());
//    return 1;
//  }
//  temp.Form("%s",MinuitOut);
//  temp.ReplaceAll("minuit.out.txt","FitPDF.out");
//  ifstream infile1(temp.Data());
//  if(!infile.is_open()){ 
//    temp.Form("%s - can not find FitPDF.out file", JobName.Data());
//    this->AddMessage(temp.Data());
//    return 1;
//  }
//  Char_t line[256];
//  Bool_t MigradConverged = kFALSE;
//  Bool_t MigradCallLimit = kFALSE;
//  Bool_t MigradNoConvergence = kFALSE;
//  Bool_t SetAlphaError = kFALSE;
//  Bool_t MigradFailsImprovement = kFALSE;
//  Bool_t MigradMachineLimit = kFALSE;
//
//  Bool_t HESSE = kFALSE;
//
//  Bool_t HesseMatrixNegative = kFALSE;
//  Bool_t HesseFails = kFALSE;
//
//  //while ( infile1.getline(line, 256)) {
//  //  temp.Form("%s",line);
//  //  if(temp.Contains("Error in SETALF")) {SetAlphaError = kTRUE;break;}
//  //}
//
//  while ( infile.getline(line, 256)) {
//    temp.Form("%s",line);
//    if(temp.Contains("HESSE")) HESSE = kTRUE;
//    if(!HESSE) {
//      if(temp.Contains("MIGRAD MINIMIZATION HAS CONVERGED")) MigradConverged = kTRUE;
//      if(temp.Contains("CALL LIMIT EXCEEDED IN MIGRAD")) MigradCallLimit = kTRUE;
//      if(temp.Contains("MIGRAD TERMINATED WITHOUT CONVERGENCE")) MigradNoConvergence = kTRUE; 
//      if(!MigradConverged && temp.Contains("MIGRAD FAILS TO FIND IMPROVEMENT")) MigradFailsImprovement = kTRUE;
//      if(temp.Contains("MACHINE ACCURACY LIMITS FURTHER IMPROVEMENT")) MigradMachineLimit = kTRUE;
//    }
//    else {
//      if(temp.Contains("MATRIX FORCED POS-DEF BY ADDING")) HesseMatrixNegative = kTRUE;
//      if(temp.Contains("MNHESS FAILS AND WILL RETURN DIAGONAL MATRIX")) HesseFails = kTRUE;
//    }
//  }
//  if(!MigradConverged) {
//    temp.Form("%s - MIGRAD CONVERGED NOT FOUND", JobName.Data());
//    this->AddMessage(temp.Data());
//  }
//  if(MigradCallLimit) {
//    temp.Form("%s - CALL LIMIT EXCEEDED IN MIGRAD", JobName.Data());
//    this->AddMessage(temp.Data());
//  }
//  if(MigradMachineLimit) {
//    temp.Form("%s - MACHINE ACCURACY LIMITS FURTHER IMPROVEMENT", JobName.Data());
//    this->AddMessage(temp.Data());
//  }
//  if(SetAlphaError) {
//    temp.Form("%s - Error in SETALF - use limits!", JobName.Data());
//    this->AddMessage(temp.Data());
//  }
//  if(MigradNoConvergence) {
//    temp.Form("%s - MIGRAD TERMINATED WITHOUT CONVERGENCE", JobName.Data());
//    this->AddMessage(temp.Data());
//  }
//  if(MigradFailsImprovement) {
//    temp.Form("%s - MIGRAD FAILS TO FIND IMPROVEMENT", JobName.Data());
//    this->AddMessage(temp.Data());
//  }
//  if(HesseMatrixNegative) {
//    temp.Form("%s - HESSE: MATRIX NEGATIVE", JobName.Data());
//    this->AddMessage(temp.Data());
//  }
//  if(HesseFails) {
//    temp.Form("%s - HESSE: MNHESSE FAILED", JobName.Data());
//    this->AddMessage(temp.Data());
//  }
//}
//
//void H1FitterOutput::AddError(H1FitterOutput::error sum, H1FitterOutput::error add) {
//  if(fAlphasUError[(Int_t)sum]<-995.) fAlphasUError[(Int_t)sum] = 0.;
//  if(fAlphasDError[(Int_t)sum]<-995.) fAlphasDError[(Int_t)sum] = 0.;
//
//  fAlphasUError[(Int_t)sum] = TMath::Sqrt(TMath::Power(fAlphasUError[(Int_t)sum], 2.0) + TMath::Power(fAlphasUError[(Int_t)add], 2.0));
//  fAlphasDError[(Int_t)sum] = TMath::Sqrt(TMath::Power(fAlphasDError[(Int_t)sum], 2.0) + TMath::Power(fAlphasDError[(Int_t)add], 2.0));
//
//  for(Int_t iq2=0; iq2<fNQ2; iq2++) {
//    AddGraphError(kGluon, iq2, sum, add);
//    AddGraphError(kU, iq2, sum, add);
//    AddGraphError(kD, iq2, sum, add);
//    AddGraphError(kUv, iq2, sum, add);
//    AddGraphError(kDv, iq2, sum, add);
//    AddGraphError(kUb, iq2, sum, add);
//    AddGraphError(kDb, iq2, sum, add);
//    AddGraphError(kSea, iq2, sum, add);
//    AddGraphError(kS, iq2, sum, add);
//    AddGraphError(kC, iq2, sum, add);
//    AddGraphError(kB, iq2, sum, add);
//  }
//}
//
//void H1FitterOutput::AddGraphError(H1FitterOutput::pdf Pdf, Int_t iq2, H1FitterOutput::error sum, H1FitterOutput::error add) {
//
//  TGraphAsymmErrors* SumGraph = GetPdf(Pdf, iq2, sum);
//  TGraphAsymmErrors* AddGraph = GetPdf(Pdf, iq2, add);
//
//  for(Int_t ipoint = 0; ipoint<fNpoints; ipoint++) {
//    Double_t x,y;
//    AddGraph->GetPoint(ipoint, x, y);
//    Double_t errh = AddGraph->GetErrorYhigh(ipoint);
//    Double_t errl = AddGraph->GetErrorYlow(ipoint);
//    AddGraphErrorPoint(SumGraph, ipoint, y+errh, y-errl, kQuadrature);
//  }
//}
//
//
//
//Double_t H1FitterOutput::GetAlphas(const Char_t* minuitOut) {
//  ifstream infile(minuitOut);
//  if(!infile.is_open()){ Warning("PrepareAlphas","can not open file %s", minuitOut); return 1;}
//  
//  Double_t ToReturn = -998.;
//  Char_t line[256];
//  Bool_t Readout = kFALSE;
//  while ( infile.getline(line, 256)) {
//    TString Line(line);
//    if(Line.Contains("MIGRAD MINIMIZATION HAS CONVERGED")) {Readout = kTRUE; continue;}
//    if(Readout && Line.Contains("alphas")) {
//      TObjArray* array = Line.Tokenize(" ");
//      ToReturn = ((TObjString*)array->At(2))->GetString().Atof();
//      delete array;
//      break;
//    }
//  }
//  return ToReturn;
//}





//
//void H1FitterOutput::AddArray(H1ArrayF* ArrayU, H1ArrayF* ArrayD, H1ArrayF* Array2, errcombscheme Scheme, H1ArrayF* ArrayRef) {
//  if(Array2->GetEntries() != ArrayU->GetEntries()) {
//    cout << "H1FitterOutput::AddArray1, arays incompatible!" << endl;
//    //exit(1);
//  }
//  if(ArrayD->GetEntries() != ArrayU->GetEntries()) {
//    cout << "H1FitterOutput::AddArray2, arays incompatible!" << endl;
//    //exit(1);
//  }
//  if(ArrayRef->GetEntries() != ArrayU->GetEntries()) {
//    cout << "H1FitterOutput::AddArray3, arays incompatible!" << endl;
//    //exit(1);
//  }
//  Int_t N = TMath::Min(ArrayU->GetEntries(), ArrayD->GetEntries());
//  N = TMath::Min(N, Array2->GetEntries());
//  N = TMath::Min(N, ArrayRef->GetEntries());
//  for(Int_t i=0; i<N; i++) {
//    if( (*ArrayD)[i] > (*ArrayRef)[i] ) {
//      cout << "ArrayD too high"<< endl; exit(1);
//    }
//    if( (*ArrayU)[i] < (*ArrayRef)[i] ) {
//      cout << "ArrayU too low"<< endl; exit(1);
//    }
//
//    (*ArrayD)[i] = (*ArrayRef)[i] - (*ArrayD)[i];
//    (*ArrayU)[i] = (*ArrayU)[i] - (*ArrayRef)[i];
//    (*Array2)[i] = (*Array2)[i] - (*ArrayRef)[i];
//
//    switch(Scheme) {
//    case kQuadrature:
//      if((*Array2)[i] > 0.)
//	(*ArrayU)[i] = TMath::Sqrt((*ArrayU)[i]*(*ArrayU)[i] + (*Array2)[i]*(*Array2)[i]);
//      else 
//	(*ArrayD)[i] = TMath::Sqrt((*ArrayD)[i]*(*ArrayD)[i] + (*Array2)[i]*(*Array2)[i]);
//      break;
//    case kMaximal:
//      if((*Array2)[i] > 0.)
//	(*ArrayU)[i] = TMath::Max((*ArrayU)[i], (*Array2)[i]);
//      else 
//	(*ArrayD)[i] = TMath::Max((*ArrayD)[i], (*Array2)[i]);
//      break;
//    }
//    (*ArrayD)[i] = (*ArrayRef)[i] - (*ArrayD)[i];
//    (*ArrayU)[i] = (*ArrayU)[i] + (*ArrayRef)[i];
//  }
//}
//
//
//void H1FitterOutput::SortArray(H1ArrayF* ArrayU, H1ArrayF* ArrayD, H1ArrayF* ArrayRef) {
//  if(ArrayU->GetEntries() != ArrayD->GetEntries()) {
//    cout << "H1FitterOutput::SortArray1, arays incompatible!" << endl;
//    return;
//  }
//  if(ArrayU->GetEntries() != ArrayRef->GetEntries()) {
//    cout << "H1FitterOutput::SortArray2, arays incompatible!" << endl;
//    return;
//  }
//  Float_t tempU;
//  Float_t tempD;
//  for(Int_t i=0; i<ArrayU->GetEntries(); i++) {
//    tempU = ArrayU->At(i);
//    tempD = ArrayD->At(i);
//    (*ArrayU)[i] = TMath::Max(tempU, tempD);
//    (*ArrayD)[i] = TMath::Min(tempU, tempD);
//    (*ArrayU)[i] = TMath::Max((*ArrayU)[i] , (*ArrayRef)[i]);
//    (*ArrayD)[i] = TMath::Min((*ArrayD)[i] , (*ArrayRef)[i]);
//  }
//}

/*
void H1FitterOutput::FillFittedResults(const Char_t* filename, H1ArrayI* SetArray, H1ArrayF* Q2Array,H1ArrayF* XArray,H1ArrayF* YArray, H1ArrayF* DataArray,H1ArrayF* UncEArray,H1ArrayF* TotEArray,H1ArrayF* TheoArray,H1ArrayF* PullArray) {

  //cout << "FillFittedResults "<< filename<<endl;

  Double_t q2, x, y, data, uncorrerr, toterr, theory, pull;
  Int_t dataset;
  ifstream infile(filename);
  if(!infile.is_open()){ Warning("FillFittedResults","can not open file %s", filename); return;}

  Char_t buffer[120];    
  // q2          x        y    data     +- uncorr.err   +-toterr      theory      pull     dataset

  Int_t idx=0;
  while(!infile.eof()) {
    infile.getline(buffer, 120);
    if(TString(buffer).Contains("uncorr.err")) { // enter the readout process
      Int_t idxDataSet = 0;
      while(1) {
	idx++;
	idxDataSet++;
	infile.getline(buffer, 120);
	TString str(buffer);
	TObjArray* array = str.Tokenize(" ");
	if(array->GetEntries() != 9) {delete array; break;}
	q2        = ((TObjString*)array->At(0))->GetString().Atof();
	x         = ((TObjString*)array->At(1))->GetString().Atof();
	y         = ((TObjString*)array->At(2))->GetString().Atof();
	data      = ((TObjString*)array->At(3))->GetString().Atof();
	uncorrerr = ((TObjString*)array->At(4))->GetString().Atof();
	toterr    = ((TObjString*)array->At(5))->GetString().Atof();
	theory    = ((TObjString*)array->At(6))->GetString().Atof();
	pull      = ((TObjString*)array->At(7))->GetString().Atof();
	dataset   = ((TObjString*)array->At(8))->GetString().Atoi();	
	delete array;

	if(dataset == 39) {// put all the q2 values to one histogram
	  x=q2;
	  q2=0.; 
	}
	if(dataset == 40) {
	  if(idxDataSet==1) {q2=1.; x=19.0;}
	  if(idxDataSet==2) {q2=2.; x=19.0;}
	  if(idxDataSet==3) {q2=2.; x=23.0;}
	  if(idxDataSet==4) {q2=2.; x=27.0;}
	  if(idxDataSet==5) {q2=3.; x=19.0;}
	  if(idxDataSet==6) {q2=3.; x=23.0;}
	  if(idxDataSet==7) {q2=3.; x=27.0;}
	  if(idxDataSet==8) {q2=3.; x=32.0;}
	  if(idxDataSet==9) {q2=3.; x=38.0;}	 
	  if(idxDataSet==10) {q2=3.; x=44.5;}
	  if(idxDataSet==11) {q2=4.; x=19.0;}
	  if(idxDataSet==12) {q2=4.; x=23.0;}
	  if(idxDataSet==13) {q2=4.; x=27.0;}
	  if(idxDataSet==14) {q2=4.; x=32.0;}
	  if(idxDataSet==15) {q2=4.; x=38.0;}
	  if(idxDataSet==16) {q2=5.; x=19.0;}
	  if(idxDataSet==17) {q2=5.; x=23.0;}
	  if(idxDataSet==18) {q2=5.; x=27.0;}
	  if(idxDataSet==19) {q2=5.; x=32.0;}
	  if(idxDataSet==20) {q2=5.; x=38.0;}
	  if(idxDataSet==21) {q2=5.; x=44.5;}
	  if(idxDataSet==22) {q2=5.; x=51.5;}
	  if(idxDataSet==23) {q2=5.; x=60.0;}
	  if(idxDataSet==24) {q2=6.; x=19.0;}
	  if(idxDataSet==25) {q2=6.; x=23.0;}
	  if(idxDataSet==26) {q2=6.; x=27.0;}
	  if(idxDataSet==27) {q2=6.; x=32.0;}
	  if(idxDataSet==28) {q2=6.; x=38.0;}
	  if(idxDataSet==29) {q2=6.; x=44.5;}
	  if(idxDataSet==30) {q2=6.; x=51.5;}
	  if(idxDataSet==31) {q2=6.; x=60.0;}
	  if(idxDataSet==32) {q2=6.; x=70.0;}
	}
	if(dataset == 41 || dataset == 42) {
	  q2 = Float_t(TMath::Floor(Float_t(idxDataSet-1)/5.))+1.;
	  if((idxDataSet-1)%5==0) x = TMath::Exp((TMath::Log(8.) +TMath::Log(10.))/2.);
	  if((idxDataSet-1)%5==1) x = TMath::Exp((TMath::Log(10.)+TMath::Log(14.))/2.);
	  if((idxDataSet-1)%5==2) x = TMath::Exp((TMath::Log(14.)+TMath::Log(18.))/2.);
	  if((idxDataSet-1)%5==3) x = TMath::Exp((TMath::Log(18.)+TMath::Log(25.))/2.);
	  if((idxDataSet-1)%5==4) x = TMath::Exp((TMath::Log(25.)+TMath::Log(100.))/2.);
	}

	if(SetArray)  SetArray->AddLast(dataset);
	if(Q2Array)   Q2Array->AddLast(q2);
	if(XArray)    XArray->AddLast(x);
	if(YArray)    YArray->AddLast(y);
	if(DataArray) DataArray->AddLast(data);
	if(UncEArray) UncEArray->AddLast(uncorrerr);
	if(TotEArray) TotEArray->AddLast(toterr);
	if(TheoArray) TheoArray->AddLast(theory);
	if(PullArray) PullArray->AddLast(pull);
      }
    }
  }
}
*/
//
//
//Int_t H1FitterOutput::AddError(const Char_t* app1, const Char_t* app2, H1FitterOutput::error err, 
//			       H1FitterOutput::errcombscheme scheme1, H1FitterOutput::error comberr1,
//			       H1FitterOutput::errcombscheme scheme2, H1FitterOutput::error comberr2) {
//  Double_t x1, gluon1, U1, D1, d_Ubar1, d_Dbar1, umin1, dmin1, sea1, u_sea1, d_sea1, str1, chm1, bot1;
//  Double_t x2, gluon2, U2, D2, d_Ubar2, d_Dbar2, umin2, dmin2, sea2, u_sea2, d_sea2, str2, chm2, bot2;
//  x1=gluon1=U1=D1=d_Ubar1=d_Dbar1=umin1=dmin1=sea1=u_sea1=d_sea1=str1=chm1=bot1=-9999.;
//  x2=gluon2=U2=D2=d_Ubar2=d_Dbar2=umin2=dmin2=sea2=u_sea2=d_sea2=str2=chm2=bot2=-9999.;
//  
//  TString* filename1 = new TString;
//  TString* filename2 = new TString;
//
//  filename1->Form("%s_%s/minuit.out.txt",fDirectory->Data(), app1);
//  filename2->Form("%s_%s/minuit.out.txt",fDirectory->Data(), app2);
//
//  this->CheckJob(filename1->Data());
//  this->CheckJob(filename2->Data());
//
//  Double_t Alphas1 = GetAlphas(filename1->Data());
//  Double_t Alphas2 = GetAlphas(filename2->Data());
//
//  if(Alphas1 >= fAlphas && Alphas2 >= fAlphas) { 
//    fAlphasUError[(Int_t)err] = TMath::Max(Alphas1, Alphas2) - fAlphas;
//    fAlphasDError[(Int_t)err] = 0.;
//  }
//  else if (Alphas1 >= fAlphas && Alphas2 <= fAlphas) {
//    fAlphasUError[(Int_t)err] = Alphas1 - fAlphas;
//    fAlphasDError[(Int_t)err] = fAlphas - Alphas2;
//  }  
//  else if (Alphas1 <= fAlphas && Alphas2 >= fAlphas) {
//    fAlphasUError[(Int_t)err] = Alphas2 - fAlphas;
//    fAlphasDError[(Int_t)err] = fAlphas - Alphas1;
//  }  
//  else if(Alphas1 <= fAlphas && Alphas2 <= fAlphas) {
//    fAlphasDError[(Int_t)err] = fAlphas - TMath::Min(Alphas1, Alphas2);
//    fAlphasUError[(Int_t)err] = 0.;
//  }
//  
//  if(comberr1 != kNoError) { // combine
//    if(fAlphasDError[(Int_t)comberr1] < 0.) fAlphasDError[(Int_t)comberr1] = 0.; 
//    if(fAlphasUError[(Int_t)comberr1] < 0.) fAlphasUError[(Int_t)comberr1] = 0.; 
//    if(scheme1==kQuadrature) {
//      fAlphasDError[(Int_t)comberr1] = TMath::Sqrt(TMath::Power(fAlphasDError[(Int_t)comberr1], 2.) + TMath::Power(fAlphasDError[(Int_t)err], 2.));
//      fAlphasUError[(Int_t)comberr1] = TMath::Sqrt(TMath::Power(fAlphasUError[(Int_t)comberr1], 2.) + TMath::Power(fAlphasUError[(Int_t)err], 2.));
//    }
//    else if(scheme1==kMaximal) {
//      fAlphasDError[(Int_t)comberr1] = TMath::Max(fAlphasDError[(Int_t)comberr1], fAlphasDError[(Int_t)err]);
//      fAlphasUError[(Int_t)comberr1] = TMath::Max(fAlphasUError[(Int_t)comberr1], fAlphasUError[(Int_t)err]);
//    }
//    else {Warning("AddError","combination scheme not supported");}
//  }
//  if(comberr2 != kNoError) {
//    if(fAlphasDError[(Int_t)comberr2] < 0.) fAlphasDError[(Int_t)comberr2] = 0.; 
//    if(fAlphasUError[(Int_t)comberr2] < 0.) fAlphasUError[(Int_t)comberr2] = 0.; 
//    if(scheme2==kQuadrature) {
//      fAlphasDError[(Int_t)comberr2] = TMath::Sqrt(TMath::Power(fAlphasDError[(Int_t)comberr2], 2.) + TMath::Power(fAlphasDError[(Int_t)err], 2.));
//      fAlphasUError[(Int_t)comberr2] = TMath::Sqrt(TMath::Power(fAlphasUError[(Int_t)comberr2], 2.) + TMath::Power(fAlphasUError[(Int_t)err], 2.));
//    }
//    else if(scheme2==kMaximal) {
//      fAlphasDError[(Int_t)comberr2] = TMath::Max(fAlphasDError[(Int_t)comberr2], fAlphasDError[(Int_t)err]);
//      fAlphasUError[(Int_t)comberr2] = TMath::Max(fAlphasUError[(Int_t)comberr2], fAlphasUError[(Int_t)err]);
//    }
//    else {Warning("AddError","combination scheme not supported");}
//  }
//
//  for(Int_t iq2=0; iq2<fNQ2; iq2++) {
//    filename1->Form("%s_%s/pdfs_q2val_%02d.txt",fDirectory->Data(), app1, iq2+1);
//    filename2->Form("%s_%s/pdfs_q2val_%02d.txt",fDirectory->Data(), app2, iq2+1);
//    ifstream infile1(filename1->Data());
//    ifstream infile2(filename2->Data());
//    if(!infile1.is_open()){ Warning("AddError","can not open file %s", filename1->Data()); continue;}
//    if(!infile2.is_open()){ Warning("AddError","can not open file %s", filename2->Data()); continue;}
//    for (Int_t i = 0; i < fNpoints; i++){
//      if(infile1.is_open())
//	infile1 >> x1 >> gluon1 >> U1 >> D1 >> d_Ubar1 >> d_Dbar1 >> umin1 >> dmin1 >> sea1 >> u_sea1 >> d_sea1 >> str1 >> chm1 >> bot1;
//      if(infile2.is_open())
//	infile2 >> x2 >> gluon2 >> U2 >> D2 >> d_Ubar2 >> d_Dbar2 >> umin2 >> dmin2 >> sea2 >> u_sea2 >> d_sea2 >> str2 >> chm2 >> bot2;
//      
//
//      this->AddGraphErrorPoint(fPdf[(Int_t)kGluon][iq2][(Int_t)err], i, gluon1, gluon2);
//      this->AddGraphErrorPoint(fPdf[(Int_t)kU][iq2][(Int_t)err],     i, U1, U2);
//      this->AddGraphErrorPoint(fPdf[(Int_t)kD][iq2][(Int_t)err],     i, D1, D2);
//      this->AddGraphErrorPoint(fPdf[(Int_t)kUv][iq2][(Int_t)err],    i, U1 - u_sea1 - chm1, U2 - u_sea2- chm2);
//      this->AddGraphErrorPoint(fPdf[(Int_t)kDv][iq2][(Int_t)err],    i, D1 - d_sea1 - str1 - bot1, D2 - d_sea2 - str2 - bot2);
//      this->AddGraphErrorPoint(fPdf[(Int_t)kUb][iq2][(Int_t)err],    i, d_Ubar1,  d_Ubar2);
//      this->AddGraphErrorPoint(fPdf[(Int_t)kDb][iq2][(Int_t)err],    i, d_Dbar1,  d_Dbar2);
//      this->AddGraphErrorPoint(fPdf[(Int_t)kSea][iq2][(Int_t)err],   i, sea1, sea2);
//      this->AddGraphErrorPoint(fPdf[(Int_t)kS][iq2][(Int_t)err],     i, str1, str2);
//      this->AddGraphErrorPoint(fPdf[(Int_t)kC][iq2][(Int_t)err],     i, chm1, chm2);
//      this->AddGraphErrorPoint(fPdf[(Int_t)kB][iq2][(Int_t)err],     i, bot1, bot2);
//      
//      if(comberr1 != kNoError) { // combine
//	this->AddGraphErrorPoint(fPdf[(Int_t)kGluon][iq2][(Int_t)comberr1], i, gluon1, gluon2, scheme1);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kU][iq2][(Int_t)comberr1],     i, U1, U2, scheme1);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kD][iq2][(Int_t)comberr1],     i, D1, D2, scheme1);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kUv][iq2][(Int_t)comberr1],    i, U1 - u_sea1 - chm1, U2 - u_sea2- chm2, scheme1);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kDv][iq2][(Int_t)comberr1],    i, D1 - d_sea1 - str1 - bot1, D2 - d_sea2 - str2 - bot2, scheme1);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kUb][iq2][(Int_t)comberr1],    i, d_Ubar1,  d_Ubar2, scheme1);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kDb][iq2][(Int_t)comberr1],    i, d_Dbar1,  d_Dbar2, scheme1);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kSea][iq2][(Int_t)comberr1],   i, sea1, sea2, scheme1);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kS][iq2][(Int_t)comberr1],     i, str1, str2, scheme1);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kC][iq2][(Int_t)comberr1],     i, chm1, chm2, scheme1);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kB][iq2][(Int_t)comberr1],     i, bot1, bot2, scheme1);
//      }
//      if(comberr2 != kNoError) { // combine
//	this->AddGraphErrorPoint(fPdf[(Int_t)kGluon][iq2][(Int_t)comberr2], i, gluon1, gluon2, scheme2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kU][iq2][(Int_t)comberr2],     i, U1, U2, scheme2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kD][iq2][(Int_t)comberr2],     i, D1, D2, scheme2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kUv][iq2][(Int_t)comberr2],    i, U1 - u_sea1 - chm1, U2 - u_sea2- chm2, scheme2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kDv][iq2][(Int_t)comberr2],    i, D1 - d_sea1 - str1 - bot1, D2 - d_sea2 - str2 - bot2, scheme2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kUb][iq2][(Int_t)comberr2],    i, d_Ubar1,  d_Ubar2, scheme2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kDb][iq2][(Int_t)comberr2],    i, d_Dbar1,  d_Dbar2, scheme2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kSea][iq2][(Int_t)comberr2],   i, sea1, sea2, scheme2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kS][iq2][(Int_t)comberr2],     i, str1, str2, scheme2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kC][iq2][(Int_t)comberr2],     i, chm1, chm2, scheme2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kB][iq2][(Int_t)comberr2],     i, bot1, bot2, scheme2);
//      }
//    }
//  }
//  delete filename1; delete filename2;
//}
//Int_t H1FitterOutput::AddBandsExpError() {
//  Double_t x1, gluon1, U1, D1, d_Ubar1, d_Dbar1, umin1, dmin1, sea1, u_sea1, d_sea1, str1, chm1, bot1;
//  Double_t x2, gluon2, U2, D2, d_Ubar2, d_Dbar2, umin2, dmin2, sea2, u_sea2, d_sea2, str2, chm2, bot2;
//  TString* filename1 = new TString;
//  TString* filename2 = new TString;
//  for(Int_t iq2=0; iq2<fNQ2; iq2++) {
//    for(Int_t ip=0; ip<11; ip++) {
//      //if(ip!=1) continue;
//      filename1->Form("%s_BANDS/pdfs_q2val_s%02dp_%02d.txt",fDirectory->Data(), ip+1, iq2+1);
//      filename2->Form("%s_BANDS/pdfs_q2val_s%02dm_%02d.txt",fDirectory->Data(), ip+1, iq2+1);
//      ifstream infile1(filename1->Data());
//      ifstream infile2(filename2->Data());
//      if(!infile1.is_open()){ Warning("AddExpError","can not open file %s", filename1->Data()); continue;}
//      if(!infile2.is_open()){ Warning("AddExpError","can not open file %s", filename2->Data()); continue;}
//      
//      for (Int_t i = 0; i < fNpoints; i++){
//	infile1 >> x1 >> gluon1 >> U1 >> D1 >> d_Ubar1 >> d_Dbar1 >> umin1 >> dmin1 >> sea1 >> u_sea1 >> d_sea1 >> str1 >> chm1 >> bot1;
//	infile2 >> x2 >> gluon2 >> U2 >> D2 >> d_Ubar2 >> d_Dbar2 >> umin2 >> dmin2 >> sea2 >> u_sea2 >> d_sea2 >> str2 >> chm2 >> bot2;
//
//	this->AddGraphErrorPoint(fPdf[(Int_t)kGluon][iq2][kExpError], i, gluon1, gluon2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kU][iq2][kExpError],     i, U1, U2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kD][iq2][kExpError],     i, D1, D2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kUv][iq2][kExpError],    i, U1 - u_sea1 - chm1, U2 - u_sea2- chm2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kDv][iq2][kExpError],    i, D1 - d_sea1 - str1 - bot1, D2 - d_sea2 - str2 - bot2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kUb][iq2][kExpError],    i, d_Ubar1,  d_Ubar2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kDb][iq2][kExpError],    i, d_Dbar1,  d_Dbar2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kSea][iq2][kExpError],   i, sea1, sea2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kS][iq2][kExpError],     i, str1, str2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kC][iq2][kExpError],     i, chm1, chm2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kB][iq2][kExpError],     i, bot1, bot2);
//
//	H1FitterOutput::error E=kExpError1;
//	switch(ip+1) {
//	case 1: E=kExpError1; break;	case 2: E=kExpError2; break;	case 3: E=kExpError3; break;
//	case 4: E=kExpError4; break;	case 5: E=kExpError5; break;	case 6: E=kExpError6; break;
//	case 7: E=kExpError7; break;	case 8: E=kExpError8; break;	case 9: E=kExpError9; break;
//	case 10: E=kExpError10; break;	case 11: E=kExpError11; break;
//	} 
//	this->AddGraphErrorPoint(fPdf[(Int_t)kGluon][iq2][E], i, gluon1, gluon2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kU][iq2][E],     i, U1, U2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kD][iq2][E],     i, D1, D2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kUv][iq2][E],    i, U1 - u_sea1 - chm1, U2 - u_sea2- chm2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kDv][iq2][E],    i, D1 - d_sea1 - str1 - bot1, D2 - d_sea2 - str2 - bot2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kUb][iq2][E],    i, d_Ubar1,  d_Ubar2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kDb][iq2][E],    i, d_Dbar1,  d_Dbar2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kSea][iq2][E],   i, sea1, sea2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kS][iq2][E],     i, str1, str2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kC][iq2][E],     i, chm1, chm2);
//	this->AddGraphErrorPoint(fPdf[(Int_t)kB][iq2][E],     i, bot1, bot2);
//
//      }
//      infile1.close();
//      infile2.close();
//
//    }
//  }
//  delete filename1; delete filename2;
//}
//
//void H1FitterOutput::AddGraphErrorPoint(TGraphAsymmErrors* graph, Int_t i, Double_t y1, Double_t y2,
//					H1FitterOutput::errcombscheme scheme) {
//  //Bool_t PRINT = kFALSE;
//  //if(scheme == kMaximal) PRINT=kTRUE;
//
//  Double_t x,y;
//  graph->GetPoint(i, x, y);
//  Double_t err1 = y-y1;
//  Double_t err2 = y-y2;	
//  if(y1 < -9000) err1 = 0.;
//  if(y2 < -9000) err2 = 0.;
//  //if(err1 == 0. && err2==0.) PRINT = kFALSE;
//  Double_t errh = graph->GetErrorYhigh(i);
//  Double_t errl = graph->GetErrorYlow(i);
//  //if(PRINT) cout <<" y: " << y << " "<< err1 <<" "<<err2 << " Old: " << errh << " " << errl;
//  if(err1 > 0.) {
//    if(err2 > 0.) {
//      if(scheme==kQuadrature)
//	errl = TMath::Sqrt(errl*errl + TMath::Max(err1, err2)*TMath::Max(err1, err2));
//      else if(scheme==kMaximal)
//	errl = TMath::Max(errl, TMath::Max(err1, err2));
//    }
//    else {
//      //if(PRINT) cout << y << " "<<y1 << " "<<y2 << endl;
//      err2 *= -1.;
//      if(scheme==kQuadrature) {
//      	errh = TMath::Sqrt(errh*errh + err2*err2);
//      	errl = TMath::Sqrt(errl*errl + err1*err1);
//      }
//      else if(scheme==kMaximal) {
//      	errh = TMath::Max(errh, err2);
//      	errl = TMath::Max(errl, err1);
//      }
//    }
//  }
//  else {
//    err1 *= -1.;
//    if(err2 > 0.) {
//      //if (PRINT) cout << "HERE" << endl;
//      if(scheme==kQuadrature) {
//	errh = TMath::Sqrt(errh*errh + err1*err1);
//	errl = TMath::Sqrt(errl*errl + err2*err2);
//      }
//      else if(scheme==kMaximal) {
//	errh = TMath::Max(errh, err1);
//	errl = TMath::Max(errl, err2);
//      }
//    }
//    else {
//      err2 *= -1.;
//      //if(PRINT) cout << err1 << " "<< err2 << " ::: " << errh;
//      if(scheme==kQuadrature) {
//      	errh = TMath::Sqrt(errh*errh + TMath::Max(err1, err2)*TMath::Max(err1, err2));
//      }
//      else if(scheme==kMaximal) {
//      	errh = TMath::Max(errh, TMath::Max(err1, err2));
//      }
//      //if(PRINT) cout << " ---> " << errh << endl;
//    }
//  }
//  //if (PRINT) cout << " New: " << errh << " "<< errl<<endl;
//  //if(PRINT) exit(0);
//  this->SetGraphPoint(graph, i, x, y, errh, errl);
//}
//



//DataSet* H1FitterOutput::FindDataSet(Int_t setId) {
//  for(Int_t iset = 0; iset<fDataSets->GetEntries(); iset++) {
//    if(((DataSet*) fDataSets->At(iset))->GetSetId() == setId) return (DataSet*) fDataSets->At(iset);
//  }
//  Error("FindDataSet","Can not find data set %d ", setId);
//  return NULL;
//}
//
//void H1FitterOutput::PlotDataSet(TCanvas* can, Int_t setId, const Char_t* name, H1FitterOutput* FitterRef, const Char_t* NameRef) {
//  DataSet* Set = FindDataSet(setId);
//  DataSet* SetRef = NULL; 
//  if(FitterRef) FitterRef->FindDataSet(setId);
//  Set->PlotDataSet(can, name, SetRef, NameRef);
//}
//
//Bool_t H1FitterOutput::UsingSet(Int_t set) {
//  for(Int_t iset = 0; iset<fDataSets->GetEntries(); iset++) 
//    if(((DataSet*) fDataSets->At(iset))->GetSetId() == set) return kTRUE;
//  return kFALSE;
//}
//Int_t H1FitterOutput::AddExpError() {
//
//  Double_t x, gluon, U, D, d_Ubar, d_Dbar, umin, dmin, sea, u_sea, d_sea, str, chm, bot;
//  Double_t gluon_rms[fNpoints], U_rms[fNpoints], D_rms[fNpoints], Uv_rms[fNpoints],Dv_rms[fNpoints],Ub_rms[fNpoints],Db_rms[fNpoints],Sea_rms[fNpoints],S_rms[fNpoints],C_rms[fNpoints],B_rms[fNpoints];
//  //Double_t gluon_mean[fNpoints], U_mean[fNpoints], D_mean[fNpoints], Uv_mean[fNpoints],	Dv_mean[fNpoints], Ub_mean[fNpoints], Db_mean[fNpoints], Sea_mean[fNpoints], S_mean[fNpoints], C_mean[fNpoints], B_mean[fNpoints];
//
//  TString* filename  = new TString;
//
//
//  filename->Form("%s/minuit.out.txt",fDirectory->Data());
//
//  ifstream infile(filename->Data());
//  if(!infile.is_open()){ Warning("PrepareAlphas","can not open file %s", filename->Data()); return 1;}
//  Char_t line[256];
//  Bool_t Converged = kFALSE;
//  Bool_t Readout = kFALSE;
////  while ( infile.getline(line, 256)) {
////    TString Line(line);
////    if(Line.Contains("MIGRAD MINIMIZATION HAS CONVERGED")) {Converged = kTRUE; continue;}
////    if(Converged && Line.Contains("HESSE")) {Readout = kTRUE; continue;}
////    if(Readout && Line.Contains("alphas")) {
////      TObjArray* array = Line.Tokenize(" ");
////      fAlphasUError[(Int_t)kExpError] = ((TObjString*)array->At(3))->GetString().Atof();
////      fAlphasDError[(Int_t)kExpError] = ((TObjString*)array->At(3))->GetString().Atof();
////      delete array;
////      break;
////    }
////  }
//  
//  // MC method for alphas
//  Int_t NErr=0;
//  Double_t Alphas_rms = 0.;  
//  for(Int_t ierr=1; ierr<=99; ierr++) {
//    filename->Form("%s_err%d/minuit.out.txt",fDirectory->Data(), ierr);
//    ifstream infile(filename->Data());
//    if(!infile.is_open()){ Warning("AddExpError","can not open file %s", filename->Data()); continue;}
//    
//    Bool_t Converged = kFALSE;
//    Bool_t Readout = kFALSE;
//    while ( infile.getline(line, 256)) {
//      TString Line(line);
//      if(Line.Contains("MIGRAD MINIMIZATION HAS CONVERGED")) {Converged = kTRUE; continue;}
//      if(Converged && Line.Contains("HESSE")) {Readout = kTRUE; continue;}
//      if(Readout && Line.Contains("alphas")) {
//	TObjArray* array = Line.Tokenize(" ");
//	Alphas_rms += TMath::Power(((TObjString*)array->At(2))->GetString().Atof() - this->GetAlphas(), 2.);
//	//cout << ((TObjString*)array->At(2))->GetString().Atof() << " - " << this->GetAlphas() << endl;
//	delete array;
//	NErr++;
//	break;
//      }
//    }
//  }
//
//  fAlphasUError[(Int_t)kExpError] = TMath::Sqrt(Alphas_rms / (Double_t) NErr);
//  fAlphasDError[(Int_t)kExpError] = TMath::Sqrt(Alphas_rms / (Double_t) NErr);
//  //exit(0);
//
//  for(Int_t iq2=0; iq2<fNQ2; iq2++) {
//    for (Int_t i=0; i<fNpoints; i++) {
//      gluon_rms[i]=0.;      U_rms[i]=0.;      D_rms[i]=0.;      Uv_rms[i]=0.;      Dv_rms[i]=0.;      Ub_rms[i]=0.;
//      Db_rms[i]=0.;      Sea_rms[i]=0.;      S_rms[i]=0.;      C_rms[i]=0.;      B_rms[i]=0.;   
//      //gluon_mean[i]=0.;      U_mean[i]=0.;      D_mean[i]=0.;      Uv_mean[i]=0.;      Dv_mean[i]=0.;      Ub_mean[i]=0.;      Db_mean[i]=0.;
//      //Sea_mean[i]=0.;      S_mean[i]=0.;      C_mean[i]=0.;      B_mean[i]=0.;
//    }
//    
//    NErr=0;
//    for(Int_t ierr=1; ierr<=99; ierr++) {
//      filename->Form("%s_err%d/pdfs_q2val_%02d.txt",fDirectory->Data(), ierr, iq2+1);
//      ifstream infile(filename->Data());
//      if(!infile.is_open()){ Warning("AddExpError","can not open file %s", filename->Data()); continue;}
//            
//      NErr++;
//      for (Int_t i = 0; i < fNpoints; i++){
//	infile >> x >> gluon >> U >> D >> d_Ubar >> d_Dbar >> umin >> dmin >> sea >> u_sea >> d_sea >> str >> chm >> bot;
//
////	gluon_mean[i] += gluon;
////	U_mean[i] += U;
////	D_mean[i] += D;
////	Uv_mean[i] += U - u_sea - chm;
////	Dv_mean[i] += D - d_sea - str - bot;
////	Ub_mean[i] += d_Ubar;
////	Db_mean[i] += d_Dbar;
////	Sea_mean[i] += sea;
////	S_mean[i] += str;
////	C_mean[i] += chm;
////	B_mean[i] += bot;
//
//	gluon_rms[i] += ((gluon-fPdf[(Int_t)kGluon][iq2][kExpError]->GetY()[i]) * (gluon-fPdf[(Int_t)kGluon][iq2][kExpError]->GetY()[i]));
//	U_rms[i] += ((U-fPdf[(Int_t)kU][iq2][kExpError]->GetY()[i])*(U-fPdf[(Int_t)kU][iq2][kExpError]->GetY()[i]));
//	D_rms[i] += ((D-fPdf[(Int_t)kD][iq2][kExpError]->GetY()[i])*(D-fPdf[(Int_t)kD][iq2][kExpError]->GetY()[i]));
//	Uv_rms[i] += (U - u_sea - chm-fPdf[(Int_t)kUv][iq2][kExpError]->GetY()[i]) * (U - u_sea - chm-fPdf[(Int_t)kUv][iq2][kExpError]->GetY()[i]);
//	Dv_rms[i] += (D - d_sea - str - bot-fPdf[(Int_t)kDv][iq2][kExpError]->GetY()[i]) * (D - d_sea - str - bot-fPdf[(Int_t)kDv][iq2][kExpError]->GetY()[i]);
//	Ub_rms[i] += (d_Ubar-fPdf[(Int_t)kUb][iq2][kExpError]->GetY()[i]) * (d_Ubar-fPdf[(Int_t)kUb][iq2][kExpError]->GetY()[i]);
//	Db_rms[i] += (d_Dbar-fPdf[(Int_t)kDb][iq2][kExpError]->GetY()[i]) * (d_Dbar-fPdf[(Int_t)kDb][iq2][kExpError]->GetY()[i]);
//	Sea_rms[i] += (sea-fPdf[(Int_t)kSea][iq2][kExpError]->GetY()[i]) * (sea-fPdf[(Int_t)kSea][iq2][kExpError]->GetY()[i]);
//	S_rms[i] += (str-fPdf[(Int_t)kS][iq2][kExpError]->GetY()[i]) * (str-fPdf[(Int_t)kS][iq2][kExpError]->GetY()[i]);
//	C_rms[i] += (chm-fPdf[(Int_t)kC][iq2][kExpError]->GetY()[i]) * (chm-fPdf[(Int_t)kC][iq2][kExpError]->GetY()[i]);
//	B_rms[i] += (bot-fPdf[(Int_t)kB][iq2][kExpError]->GetY()[i]) * (bot-fPdf[(Int_t)kB][iq2][kExpError]->GetY()[i]);
//
////	if(i==0) {
////	  cout << gluon << " " <<TMath::Sqrt(gluon_rms[i]/(Double_t)NErr)<<endl;
////	}
//
//      }
//
//      infile.close();
//    }
//    
//    for (Int_t i = 0; i < fNpoints; i++){  
//      fPdf[(Int_t)kGluon][iq2][kExpError]->SetPointEYlow(i, TMath::Sqrt(gluon_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kGluon][iq2][kExpError]->SetPointEYhigh(i, TMath::Sqrt(gluon_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kGluon][iq2][kExpError]->SetPointEXlow(i, 1e-10);
//      fPdf[(Int_t)kGluon][iq2][kExpError]->SetPointEXhigh(i, 1e-10);
//      
//      fPdf[(Int_t)kU][iq2][kExpError]->SetPointEYlow(i, TMath::Sqrt(U_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kU][iq2][kExpError]->SetPointEYhigh(i, TMath::Sqrt(U_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kU][iq2][kExpError]->SetPointEXlow(i, 1e-10);
//      fPdf[(Int_t)kU][iq2][kExpError]->SetPointEXhigh(i, 1e-10);
//      
//      fPdf[(Int_t)kD][iq2][kExpError]->SetPointEYlow(i, TMath::Sqrt(D_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kD][iq2][kExpError]->SetPointEYhigh(i, TMath::Sqrt(D_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kD][iq2][kExpError]->SetPointEXlow(i, 1e-10);
//      fPdf[(Int_t)kD][iq2][kExpError]->SetPointEXhigh(i, 1e-10);
//      
//      fPdf[(Int_t)kUv][iq2][kExpError]->SetPointEYlow(i, TMath::Sqrt(Uv_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kUv][iq2][kExpError]->SetPointEYhigh(i, TMath::Sqrt(Uv_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kUv][iq2][kExpError]->SetPointEXlow(i, 1e-10);
//      fPdf[(Int_t)kUv][iq2][kExpError]->SetPointEXhigh(i, 1e-10);
//      
//      fPdf[(Int_t)kDv][iq2][kExpError]->SetPointEYlow(i, TMath::Sqrt(Dv_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kDv][iq2][kExpError]->SetPointEYhigh(i, TMath::Sqrt(Dv_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kDv][iq2][kExpError]->SetPointEXlow(i, 1e-10);
//      fPdf[(Int_t)kDv][iq2][kExpError]->SetPointEXhigh(i, 1e-10);
//      
//      fPdf[(Int_t)kUb][iq2][kExpError]->SetPointEYlow(i, TMath::Sqrt(Ub_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kUb][iq2][kExpError]->SetPointEYhigh(i, TMath::Sqrt(Ub_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kUb][iq2][kExpError]->SetPointEXlow(i, 1e-10);
//      fPdf[(Int_t)kUb][iq2][kExpError]->SetPointEXhigh(i, 1e-10);
//      
//      fPdf[(Int_t)kDb][iq2][kExpError]->SetPointEYlow(i, TMath::Sqrt(Db_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kDb][iq2][kExpError]->SetPointEYhigh(i, TMath::Sqrt(Db_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kDb][iq2][kExpError]->SetPointEXlow(i, 1e-10);
//      fPdf[(Int_t)kDb][iq2][kExpError]->SetPointEXhigh(i, 1e-10);
//
//      fPdf[(Int_t)kSea][iq2][kExpError]->SetPointEYlow(i, TMath::Sqrt(Sea_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kSea][iq2][kExpError]->SetPointEYhigh(i, TMath::Sqrt(Sea_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kSea][iq2][kExpError]->SetPointEXlow(i, 1e-10);
//      fPdf[(Int_t)kSea][iq2][kExpError]->SetPointEXhigh(i, 1e-10);
//
//      fPdf[(Int_t)kS][iq2][kExpError]->SetPointEYlow(i, TMath::Sqrt(S_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kS][iq2][kExpError]->SetPointEYhigh(i, TMath::Sqrt(S_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kS][iq2][kExpError]->SetPointEXlow(i, 1e-10);
//      fPdf[(Int_t)kS][iq2][kExpError]->SetPointEXhigh(i, 1e-10);
//      
//      fPdf[(Int_t)kC][iq2][kExpError]->SetPointEYlow(i, TMath::Sqrt(C_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kC][iq2][kExpError]->SetPointEYhigh(i, TMath::Sqrt(C_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kC][iq2][kExpError]->SetPointEXlow(i, 1e-10);
//      fPdf[(Int_t)kC][iq2][kExpError]->SetPointEXhigh(i, 1e-10);
//      
//      fPdf[(Int_t)kB][iq2][kExpError]->SetPointEYlow(i, TMath::Sqrt(B_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kB][iq2][kExpError]->SetPointEYhigh(i, TMath::Sqrt(B_rms[i]/(Double_t)NErr));
//      fPdf[(Int_t)kB][iq2][kExpError]->SetPointEXlow(i, 1e-10);
//      fPdf[(Int_t)kB][iq2][kExpError]->SetPointEXhigh(i, 1e-10);
//    }
//  }
//  delete filename;
//}
