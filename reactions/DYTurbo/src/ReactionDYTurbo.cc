
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

  //check settings
  opts.check_consitency();

  //Init physics parameters
  DYTurbo::init_params();

  //Set up integration terms and bins boundaries
  DYTurbo::WarmUp();
  
  vector <double> vals;
  vector <double> errs;

  DYTurbo::compute(vals, errs);

  //insert values into output array
  copy_n(vals.begin(), vals.size(), &val[0]);
}

void ReactionDYTurbo::initTerm(TermData* td)
{
}

void ReactionDYTurbo::freeTerm(TermData* td)
{
}
