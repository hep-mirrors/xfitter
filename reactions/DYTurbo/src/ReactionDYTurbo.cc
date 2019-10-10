
/*
   @file ReactionDYTurbo.cc
   @date 2019-10-07
   @author  AddReaction.py
   Created by  AddReaction.py on 2019-10-07
*/

#include "ReactionDYTurbo.h"
#include "dyturbo/dyturbo.h"
#include "dyturbo/settings.h"
#include "dyturbo/pdf.h"

// the class factories
extern "C" ReactionDYTurbo* create() {
  return new ReactionDYTurbo();
}

// Initialize at the start of the computation
void ReactionDYTurbo::atStart()
{
  DYTurbo::init_const();
  opts.externalpdf = true;
  pdf::xfxq = pdf_xfxq_wrapper_;
  pdf::alphas = alphas_wrapper_;
}

// Main function to compute results at an iteration
void ReactionDYTurbo::compute(TermData*td,valarray<double>&val,map<string,valarray<double> >&err)
{
  td->actualizeWrappers();
  
  //read settings from input file
  string filename = "";
  if(td->hasParam("FileName"))
    filename = td->getParamS("FileName");
  opts.readfromfile(filename);
  bins.readfromfile(filename);
  //cout << "binning read from file : qt " << bins.qtbins.size() << "  m " << bins.mbins.size() << " y " << bins.ybins.size() << endl;
  
  //check settings
  opts.check_consitency();

  //Init physics parameters
  DYTurbo::init_params();

  vector <double> vals;
  vector <double> errs;
  //cout << "Calling DYTurbo::compute, vals size is " << vals.size() << endl;
  
  DYTurbo::compute(vals, errs);

  //  for (uint i = 0; i < vals.size(); i++)
  //     cout << i << "  " << vals.size() << "  " << vals[i] << "  " << errs[i] << endl;

  //check bins size
  
  //insert values into output array
  copy_n(vals.begin(), val.size(), &val[0]);
}

void ReactionDYTurbo::initTerm(TermData* td)
{
}

void ReactionDYTurbo::freeTerm(TermData* td)
{
}
