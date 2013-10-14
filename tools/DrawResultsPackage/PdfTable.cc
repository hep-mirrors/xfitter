#include "PdfTable.h"
#include <fstream>
#include <iostream>
#include <math.h>
using std::cout;
using std::endl;

PdfTable::PdfTable(TString fName){
  // Initialize vars:
  fNxValues = 0;
  fQ2value  = 0;
  fNPdfs    = 0;
  fXmin     = 0;
  fXmax     = 0;
  fTable    = NULL;

  // Read the file fName
  ifstream infile(fName.Data());
  if (!infile.is_open()) {
    //std::cout << "Can not open file "<<fName.Data()<<std::endl;
    return;
  }


  // Read the file header
  infile >> fQ2value >> fNxValues >> fNPdfs >> fXmin >> fXmax;


  // Read the column names
  for (int i=0; i<=fNPdfs; i++) {
    string var;
    infile >> var;
    fColumnNames.push_back(var);
  }
  // Read the table

  fTable = new double[(fNPdfs+1)*fNxValues];

  int icount = 0;
  for (int i=0; i<fNxValues; i++) {
    for (int j=0; j<=fNPdfs; j++) {
      double val;
      infile >> val;
      fTable[icount] = val;
      icount++;
    }
  }
}

PdfTable::~PdfTable(){
  
}

const double PdfTable::GetPDF(int iX, int iPdf){
  if ( (iPdf>fNPdfs) || (iPdf<0) ) {
    cout << "Invalid PDF ID"<<iPdf << endl;
    return 0;
  }
  if ( (iX>fNxValues) || (iX<0 ) ) {
    cout << "Invalid iX = "<<iX <<endl;
    return 0;
  }

  int id = iX * (fNPdfs+1) + iPdf;
  return fTable[id];
}

const double PdfTable::GetPDF(int iX, string name){
  int idx = GetIndex(name);
  if (idx == -1) {
    cout << "Could not find column name="<<name<<endl;
    return 0;
  }
  else {
    return GetPDF(iX,idx);
  }
}

const vector<double> PdfTable::GetPDF(string name){
  int iPdf = GetIndex(name);
  vector<double> res(NULL);
  if (iPdf == -1) {
    cout << "Could not find column name="<<name<<endl;
    return res;
  }
  for ( int i=0; i<fNxValues; i++) {
    int id = i * (fNPdfs+1) + iPdf;
    res.push_back(fTable[id]);
  }
  return res;
}

TGraphAsymmErrors* PdfTable::GetPDFGraph(string name){
  int iPdf = GetIndex(name);
  if (iPdf == -1 ){
    cout << "Could not find column name="<<name<<endl;
    return NULL;
  }
  TGraphAsymmErrors* res = new TGraphAsymmErrors( GetNx() );
  for ( int i = 0; i<GetNx(); i++ ) {
    int id = i * (fNPdfs+1) + iPdf;
    int ix = i * (fNPdfs+1);
    res->SetPoint(i,fTable[ix],fTable[id]);
  }
  res->SetTitle(name.c_str());
  return res;
}

const int PdfTable::GetIndex(string name){
  for (int i=0; i<fColumnNames.size(); i++) {
    if (name == fColumnNames[i]) {
      return i;
    }
  }  
  return -1;
}

PdfTable* PdfTable::CreatePdfTable(const Char_t* filename) {
  PdfTable* tab = NULL;
  ifstream ifile(filename);
  if(ifile)  tab = new PdfTable(filename);
  return tab;
}
//-----------------------------------------------------

PdfErrorTables::PdfErrorTables(string base, int iQ2, Bool_t SymErrors, TString option) {
  
  // Read the central table:
  TString filename("");
  filename.Form("%s/pdfs_q2val_%02d.txt",base.c_str(), iQ2);


  PdfTable *t = new PdfTable(filename.Data()); 
  
  fSymmetric = SymErrors;
  
  fQ2value  = t->GetQ2();
  fXmin     = t->GetXmin();
  fXmax     = t->GetXmax();
  fNxValues = t->GetNx();
  fNPdfs    = t->GetNPdfs();
  // Copy table:

  if (fNxValues == 0) {
    return;
  }

  fTable    = new double[(fNPdfs+1)*fNxValues]; 
  for ( int i=0; i<(fNPdfs+1)*fNxValues; i++) {
    fTable[i] = t->GetTable()[i];
  }

  // Copy vars:
  for (int i=0; i<=fNPdfs; i++) {
    fColumnNames.push_back(t->GetColumnName(i));
  }
  delete t;

  if (this->GetNx() == 0) {
    cout << "Can not open file "<< filename.Data() << endl;
    return;
  }

  sErrOpt = option;

  // Read the error sets:
  if (!sErrOpt.CompareTo("mc"))
    for ( int iband = 1; iband<=1000; iband++) 
      {
	filename.Form("%s/pdfs_q2val_mc%03ds_%02d.txt",base.c_str(), iband, iQ2);
	PdfTable *eSet = CreatePdfTable(filename.Data());
	if (!eSet) break;

	if (eSet->GetNx()>0)
	  fErrorTables.push_back(eSet);
	else {delete eSet; continue;}
      }  
  else if (!sErrOpt.CompareTo("s")) 
    for (int iband = 1; iband <= 1000; iband++) 
      {
	filename.Form("%s/pdfs_q2val_s%02ds_%02d.txt",base.c_str(), iband, iQ2);
	PdfTable *eSet = CreatePdfTable(filename.Data());
	if (!eSet) break;
	
	if (eSet->GetNx()>0)
	  fErrorTables.push_back(eSet);
	else {delete eSet; continue;}
      }
  else if (!sErrOpt.CompareTo("a") || !sErrOpt.CompareTo("b"))
    for (int iband = 1; iband <= 1000; iband++) 
      {
	filename.Form("%s/pdfs_q2val_s%02dm_%02d.txt",base.c_str(), iband, iQ2);
	PdfTable *eSet = CreatePdfTable(filename.Data());
	if (!eSet) break;

	if (eSet->GetNx()>0)
	  fErrorTables.push_back(eSet);
	else {delete eSet; continue;}
    
	filename.Form("%s/pdfs_q2val_s%02dp_%02d.txt",base.c_str(), iband, iQ2);
	eSet = new PdfTable(filename.Data());
	if (!eSet) break;

	if (eSet->GetNx()>0)
	  fErrorTables.push_back(eSet);
	else {delete eSet; continue;}
      }
  else
  if ( sErrOpt.CompareTo("m") && sErrOpt.CompareTo("p") && sErrOpt.CompareTo("mc")  ) {
    for ( int iband = 1; iband<=1000; iband++) {

      filename.Form("%s/pdfs_q2val_s%02dm_%02d.txt",base.c_str(), iband, iQ2);
      PdfTable *eSet = CreatePdfTable(filename.Data());
      if (!eSet) break;

      if (eSet->GetNx()>0) {
        fErrorTables.push_back(eSet);
      }
      else {delete eSet; continue;}
    
      filename.Form("%s/pdfs_q2val_s%02dp_%02d.txt",base.c_str(), iband, iQ2);
      eSet = new PdfTable(filename.Data());
      if (!eSet) break;

      if (eSet->GetNx()>0) {
        fErrorTables.push_back(eSet);
      }
      else {delete eSet; continue;}
    }
  } else if ( !sErrOpt.CompareTo("m") ) {
    for ( int iband = 1; iband<=1000; iband++) {

      filename.Form("%s/pdfs_q2val_m%02dm_%02d.txt",base.c_str(), iband, iQ2);
      PdfTable *eSet = CreatePdfTable(filename.Data());
      if (!eSet) break;

      if (eSet->GetNx()>0) {
        fErrorTables.push_back(eSet);
      }
      else {delete eSet; continue;}
    
      filename.Form("%s/pdfs_q2val_m%02dp_%02d.txt",base.c_str(), iband, iQ2);
      eSet = new PdfTable(filename.Data());
      if (!eSet) break;

      if (eSet->GetNx()>0) {
        fErrorTables.push_back(eSet);
      }
      else {delete eSet; continue;}
    }
  } else if ( !sErrOpt.CompareTo("p") ) {
    for ( int iband = 1; iband<=1000; iband++) {

      filename.Form("%s/pdfs_q2val_p%02d_%02d.txt",base.c_str(), iband, iQ2);
      PdfTable *eSet = CreatePdfTable(filename.Data());
      if (!eSet) break;

      if (eSet->GetNx()>0) {
        fErrorTables.push_back(eSet);
      }
      else {delete eSet; continue;}
    }
  } else {
    cout << "Unknown PdfErrorTables error option = " << sErrOpt <<endl;
    return;
  }
  //cout << " Read " << fErrorTables.size() << " error sets"<<endl;
}

TGraphAsymmErrors* PdfErrorTables::GetPDFGraph(string name) {
  int nx = this->GetNx();
  if (nx == 0) {
    return NULL;
  }

  // GetPDF index
  int iPDF = this->GetIndex(name);
  if (iPDF == -1) {
    cout << "Unknown PDF = "<<name << endl;
    return NULL;
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(nx);
  for ( int ix = 0; ix<nx; ix++) {
    // x and central value
    double x = this->GetPDF(ix,0);
    double cent = this->GetPDF(ix,iPDF);
    // errors
    double elow(0),ehigh(0);
    GetPDFError(ix,iPDF,&elow,&ehigh);
    // update the graph
    graph->SetPoint(ix,x,cent);
    graph->SetPointEYlow(ix,elow);
    graph->SetPointEYhigh(ix,ehigh);
  }
  return graph;
}

void PdfErrorTables::GetPDFError(int ix, int iPDF, double* eminus, double* eplus){

  double cent = this->GetPDF(ix,iPDF);

  *eminus = 0;
  *eplus  = 0;

  //MC replica errors
  if (!sErrOpt.CompareTo("mc")) 
    {
      double sum = 0;
      double sum2 = 0;
      for (int i=0; i<fErrorTables.size(); i++) 
	{
	  sum += fErrorTables[i]->GetPDF(ix,iPDF);
	  sum2 += pow(fErrorTables[i]->GetPDF(ix,iPDF), 2);
	}
      double maxp, maxm;
      maxm = maxp = sqrt(sum2/(double)fErrorTables.size() - pow(sum/(double)fErrorTables.size(), 2));
      (*eminus) = maxm;
      (*eplus) = maxp;
    }
  else if (!sErrOpt.CompareTo("a") || !sErrOpt.CompareTo("b"))
    {  
      for (int i=0; i<fErrorTables.size(); i+=2) 
	{
	  double vm = fErrorTables[i]->GetPDF(ix,iPDF);
	  double vp = fErrorTables[i+1]->GetPDF(ix,iPDF);

	  if (fSymmetric) //option "b", symmetrise errors
	    {
	      double err = 0.5*(vp-vm);
	      (*eminus) += err*err;
	      (*eplus)  += err*err;
	    }
	  else //option "a", asymmetric errors
	    {
	      // down variation:
	      double d1 = cent - vm;
	      double d2 = cent - vp;
	      double ed = ( d2>d1) ? d2 : d1;
      
	      if (ed<0) { ed = 0;}

	      // up variation
	      d1 = -d1;
	      d2 = -d2;
	      double ep = (d2>d1) ? d2 : d1;    
	      if (ep<0) {ep = 0;}
	      (*eminus) += ed*ed;
	      (*eplus)  += ep*ep;
	    }
	}
      (*eminus) = sqrt(*eminus);
      (*eplus) = sqrt(*eplus);
    }
  else if (!sErrOpt.CompareTo("s"))
    {
      double maxp = 0;
      double maxm = 0;
    
      for (int i=0; i<fErrorTables.size(); i++) 
	{
	  double v = fErrorTables[i]->GetPDF(ix,iPDF);
	  double err = v - cent;
	  (*eminus) += err*err;
	  (*eplus)  += err*err;
	}
      (*eminus) = sqrt(*eminus);
      (*eplus) = sqrt(*eplus);
    }
  else if ( sErrOpt.CompareTo("p") ) {
  
    for (int i=0; i<fErrorTables.size(); i+=2) {
      double vm = fErrorTables[i]->GetPDF(ix,iPDF);
      double vp = fErrorTables[i+1]->GetPDF(ix,iPDF);

      if (fSymmetric) {
        double err = 0.5*(vp-vm);
        (*eminus) += err*err;
        (*eplus)  += err*err;
      }
      else {
      // down variation:
        double d1 = cent - vm;
        double d2 = cent - vp;    
        double ed = ( d2>d1) ? d2 : d1;
      
        if (ed<0) { ed = 0;}

        // up variation
        d1 = -d1;
        d2 = -d2;
        double ep = (d2>d1) ? d2 : d1;    
        if (ep<0) {ep = 0;}
   
    
        (*eminus) += ed*ed;
        (*eplus)  += ep*ep;
      }
    }
    (*eminus) = sqrt(*eminus);
    (*eplus) = sqrt(*eplus);
  
  } else {//option p, envelope of variations
    
    double maxp = 0;
    double maxm = 0;
    
    for (int i=0; i<fErrorTables.size(); i++) {
      double v = fErrorTables[i]->GetPDF(ix,iPDF);
      double d = v - cent;
      if (fSymmetric) {
	if (d<0) d = -d;
	if (d>maxp) { maxp = d; maxm = maxp; }
      } else {
	if (d>0) {
	  if (d>maxp) maxp = d;
	} else {
	  if ((-d)>maxm) maxm = -d;
	}
      }
    }
    (*eminus) = maxm;
    (*eplus) = maxp;
    
  }

  //  cout << "error: " << ix << " "<< *eminus << " "<<*eplus <<endl;
}
