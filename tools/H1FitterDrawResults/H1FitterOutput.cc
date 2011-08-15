#include "H1FitterOutput.h"
#include <fstream>
#include <dirent.h>
#include <TObjString.h>
#include <stdlib.h>

H1FitterOutput::H1FitterOutput(const Char_t* directory) {
  fDirectory    = new TString(directory);
  for(Int_t ipdf = 0; ipdf < fNpdfs; ipdf++) {
    fPdfs[ipdf] = new TObjArray;
    fPdfs[ipdf]->SetOwner();
  }
  fNDataSets = 0;
  fNpoints = 0;
  fNQ2Files = 0;
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
  if(! this->CheckDirectory() ) {
    cerr << "Can not open directory " << fDirectory->Data() << endl; 
    exit(1);
  }
  this->PreparePdf();
  this->PrepareDataSets();
}

Bool_t H1FitterOutput::CheckDirectory() {
  DIR *pDir;
  bool bExists = false;
  pDir = opendir (fDirectory->Data());
  if (pDir != NULL) {
    bExists = true;    
    (void) closedir (pDir);
  }
  return bExists;
}

Int_t H1FitterOutput::PreparePdf() {
  Double_t x, gluon, U, D, d_Ubar, d_Dbar, umin, dmin, sea, u_sea, d_sea, str, chm, bot;
  TString* filename = new TString;
  for(Int_t iq2=0; iq2<100; iq2++) {
    filename->Form("%s/pdfs_q2val_%02d.txt",fDirectory->Data(), iq2+1);
    ifstream infile(filename->Data());
    if(!infile.is_open()) break;

    
    // Read Q2 value
    Int_t nx;
    Double_t xmin,xmax;
    infile >> fQ2Value[iq2] >> nx >> xmin >> xmax;
    
    fNQ2Files++;

    fNpoints = nx;

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
  if(ipdf >= fNpdfs) {cout << "GetPdf, wrong ipdf: "<< ipdf << endl; exit(1);}
  if(Q2bin >= fPdfs[ipdf]->GetEntries()) {cout << "GetPdf, wrong iq2: "<< Q2bin << endl; exit(1);}

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

  Double_t v1, v2, v3, data, uncorrerr, toterr, theory, theory_mod,pull;
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
  // v1     v2    v3    data     +- uncorr.err   +-toterr      theory  theory_mod     pull     dataset
  while(!infile.eof()) {
    infile.getline(buffer, 120);
    Name.Form(buffer);
    infile.getline(buffer, 120);
    //    infile.getline(buffer, 120); 
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
      if(array->GetEntries() != 10) {delete array; break;}

      
      v1        = ((TObjString*)array->At(0))->GetString().Atof();
      v2        = ((TObjString*)array->At(1))->GetString().Atof();
      v3        = ((TObjString*)array->At(2))->GetString().Atof();
      data      = ((TObjString*)array->At(3))->GetString().Atof();
      uncorrerr = ((TObjString*)array->At(4))->GetString().Atof();
      toterr    = ((TObjString*)array->At(5))->GetString().Atof();
      theory    = ((TObjString*)array->At(6))->GetString().Atof();
      theory_mod  = ((TObjString*)array->At(7))->GetString().Atof();
      pull      = ((TObjString*)array->At(8))->GetString().Atof();
      dataset   = ((TObjString*)array->At(9))->GetString().Atoi();	
      delete array;

      // cout << "haha" << dataset << " "<< data << endl;

      NewDataSet->AddPoint(v1, v2, v3, data, uncorrerr, toterr, theory_mod, pull);

      fPull->Fill(pull);
    }
  }    
  
  delete filename;
  return 0;
}

 const Double_t H1FitterOutput::GetQ2Value(Int_t iQ2bin) {
  if ( (iQ2bin < 0) || (iQ2bin> GetNQ2Files())) {
    cout << "Wrong Q2 bin" << iQ2bin<< " stop "<<endl;
    exit(1);
  };
  return fQ2Value[iQ2bin];
}





