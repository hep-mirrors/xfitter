/*
  @file ReactionAPFELgrid.cc
  @date 2021-11-16
*/

#include"ReactionAPFELgrid.h"
#include"TFile.h"
#include"xfitter_pars.h"
#include"xfitter_steer.h"
#include"xfitter_cpp_base.h"
#include"appl_grid/appl_grid.h"

#include<memory>
#include"BaseEvolution.h"


#include "APFELgrid/fastkernel.h"                                                                                                                                                                            
#include "APFELgrid/transform.h"                                                                                                                                                                             
// PDF in the evolution basis needed by APFELgrid
extern "C" void apfel_fnpdf_(const double& x, const double& Q, double* f);
void fkpdf (const double& x, const double& Q, const size_t& n, double* pdf)
{
  static double* lha_pdf = new double[13];
  pdf_xfxq_wrapper_(x,Q,lha_pdf);
  NNPDF::LHA2EVLN<double, double>(lha_pdf, pdf);
}
#include "APFELgridGeneration.h"

using namespace std;
using namespace xfitter;

struct GridData{
  unique_ptr<NNPDF::FKTable<double>>grid=nullptr;
  int Ndummysize=0;
  //Each grid is either a real APFELgrid or a dummy
  //For real  grids, grid is a valid pointer, and Ndummysize==0
  //For dummy grids, grid==nullptr          , and Ndummysize>0 size of grid
  //Dummy grids return a vector of zeroes of size Ndummysize
};

struct DatasetData{
  vector<GridData>grids;
  int order;
  const double*muR,*muF; // !> renormalisation and factorisation scales
  bool flagNorm;  // !> if true, multiply by bin width
};

const double ONE=1;

extern "C" ReactionAPFELgrid* create() {
  return new ReactionAPFELgrid();
}

void ReactionAPFELgrid::initTerm(TermData*td){
  DatasetData*data=new DatasetData;
  td->reactionData=(void*)data;

  NNPDF::FKTable<double> *FK;
  //Split entries in GridName and load grids
  string GridName=td->getParamS("GridName");
  try{
    std::istringstream ss(GridName);
    std::string token;
    while(std::getline(ss, token, ',')){
      data->grids.push_back(GridData());
      GridData&gd=data->grids.back();
      // Read FK table
      ifstream infile;
      infile.open(token.c_str());
      NNPDF::FKTable<double> *FK;
      if (infile) {
	FK = new NNPDF::FKTable<double>(infile);

	// Check that the relevant parameters agree with those given in the steering card
	double AsRef   = atof((FK->NNPDF::FKHeader::GetTag(NNPDF::FKHeader::THEORYINFO, "AlphasRef")).c_str());
	double QRef    = atof((FK->NNPDF::FKHeader::GetTag(NNPDF::FKHeader::THEORYINFO, "QRef")).c_str());
	double Q0      = atof((FK->NNPDF::FKHeader::GetTag(NNPDF::FKHeader::THEORYINFO, "Q0")).c_str());
	double MCharm  = atof((FK->NNPDF::FKHeader::GetTag(NNPDF::FKHeader::THEORYINFO, "MCharm")).c_str());
	double MBottom = atof((FK->NNPDF::FKHeader::GetTag(NNPDF::FKHeader::THEORYINFO, "MBottom")).c_str());
	double MTop    = atof((FK->NNPDF::FKHeader::GetTag(NNPDF::FKHeader::THEORYINFO, "MTop")).c_str());
	int    PtOrd   = atoi((FK->NNPDF::FKHeader::GetTag(NNPDF::FKHeader::THEORYINFO, "PerturbativeOrder")).c_str());
	cout << "AsRef:" << AsRef << endl;
      }
      else {
	cout << "The file '" << token << "' does not exist. Generating FK table ..." << endl;
	cout << endl;
	APFELgridGen::generateFK(token, *td->getParamD("Q0"),
				 *td->getParamD("mch"), *td->getParamD("mbt"), *td->getParamD("mtp"),
				 *td->getParamD("alphas"), *td->getParamD("Mz"),
				 OrderMap(td->getParamS("Order"))-1);
	ifstream infile1;
	infile1.open(token.c_str());
	FK = new NNPDF::FKTable<double>(infile1);
      }
      gd.grid.reset(FK);
      //      appl::grid*g=new appl::grid(token);
    }
  }
  catch ( const std::exception& e ) {
    cerr<<"[FATAL] Unhandled exception while trying to read APFELgrid file(s) \""<<GridName<<"\"; rethrowing exception"<<endl;
    throw e;
  }
  // Get Order
  data->order=OrderMap(td->getParamS("Order"));

  // Get MuR and MuF, or use 1.0 as default
  if(td->hasParam("muR"))
    data->muR=td->getParamD("muR");
  else
    data->muR=&ONE;
  
  if(td->hasParam("muF"))
    data->muF=td->getParamD("muF");
  else
    data->muF=&ONE;

  // Get if should normalize by dividing by bin width (no by default)
  data->flagNorm=false;

  if(td->hasParam("norm")){
    hf_errlog(21111601, "F: APFELgrid does not support norm (yet) " );
  }
  size_t Ngrids=data->grids.size();
}



void ReactionAPFELgrid::freeTerm(TermData*td){
  DatasetData*data=(DatasetData*)td->reactionData;
  size_t Ngrids=data->grids.size();
  delete data;
}

void ReactionAPFELgrid::compute(TermData*td,valarray<double>&val,map<string,valarray<double> >&err){

  const DatasetData&data=*(DatasetData*)td->reactionData;
  const int order=data.order;
  const double muR=*data.muR;
  const double muF=*data.muF;
  unsigned int pos = 0;

  size_t np=0;
  for(const GridData&gd:data.grids){
    np+=gd.grid->GetNData();
  }

  cout << np << endl;
  
  val.resize(np);
  
  for(const GridData&gd:data.grids){
    auto *grid=gd.grid.get();
    vector<double> gridVals;
    gridVals.resize(grid->GetNData());
    td->actualizeWrappers();
    grid->Convolute(fkpdf,1,gridVals.data());
    // insert values from this grid into output array
    copy_n(gridVals.begin(), gridVals.size(), &val[pos]);
    pos += grid->GetNData();
  }
}
