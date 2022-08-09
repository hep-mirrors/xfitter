
/*
   @file ReactionDYTurbo.cc
   @date 2019-10-07
   @author  AddReaction.py
   Created by  AddReaction.py on 2019-10-07
*/

#include "ReactionDYTurbo.h"
#include "BaseEvolution.h"
#include "xfitter_steer.h"
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
  //Init constants
  DYTurbo::init_const();

  //Attach PDFs and alphas
  opts.externalpdf = true;
  pdf::xfxq = pdf_xfxq_wrapper_;
  pdf::extalphas = alphas_wrapper_;
}

// Main function to compute results at an iteration
void ReactionDYTurbo::compute(TermData*td,valarray<double>&val,map<string,valarray<double> >&err)
{
  //Update pointers to PDFs
  td->actualizeWrappers();
  
  //read settings from input file
  string filename = "";
  if (td->hasParam("FileName"))
    filename = td->getParamS("FileName");
  opts.readfromfile(filename);
  bins.readfromfile(filename);
  //cout << "binning read from file : qt " << bins.qtbins.size() << "  m " << bins.mbins.size() << " y " << bins.ybins.size() << endl;

  //Fill PDF infos
  pdf::order = xfitter::get_evolution()->getPropertyI("OrderQCD");
  pdf::xmin = xfitter::get_evolution()->getPropertyD("XMin");
  pdf::qmin = xfitter::get_evolution()->getPropertyD("QMin");

  pdf::mc = xfitter::get_evolution()->getPropertyD("MCharm");
  pdf::mb = xfitter::get_evolution()->getPropertyD("MBottom");
  pdf::mt = xfitter::get_evolution()->getPropertyD("MTop");

  //cout << "order " << pdf::order << endl;
  //cout << "xmin " << pdf::xmin << endl;
  //cout << "qmin " << pdf::qmin << endl;
  //cout << "mc " << pdf::mc << endl;
  //cout << "mb " << pdf::mb << endl;
  //cout << "mt " << pdf::mt << endl;
  
  opts.silent      = true;
  opts.makehistos  = false;
  
  if (td->hasParam("g1"))
    opts.g1 = *(td->getParamD("g1"));

  if (td->hasParam("g2"))
    opts.g2 = *(td->getParamD("g2"));

  if (td->hasParam("g3"))
    opts.g3 = *(td->getParamD("g3"));

  //read g from LHAPDF
  pdf::g1 = xfitter::get_evolution()->getPropertyD("g", -1.);
  pdf::g1 = xfitter::get_evolution()->getPropertyD("g1", -1.);
  pdf::g2 = xfitter::get_evolution()->getPropertyD("g2", -1.);

  //read e from LHAPDF
  pdf::e = xfitter::get_evolution()->getPropertyD("e", -1.);

  //read g0 from LHAPDF
  pdf::g0 = xfitter::get_evolution()->getPropertyD("g0", -1.);

  if (td->hasParam("Q0"))
    opts.Q0 = *(td->getParamD("Q0"));

  if (td->hasParam("order"))
    opts.order = td->getParamI("order");
  
  if (td->hasParam("muR"))
    opts.kmuren = *(td->getParamD("muR"));

  if (td->hasParam("muF"))
    opts.kmufac = *(td->getParamD("muF"));

  if (td->hasParam("muRes"))
    opts.kmures = *(td->getParamD("muRes"));

  //check settings
  opts.check_consistency();

  //cout << opts.kmuren << "  " << opts.kmufac << "  " << opts.kmures << endl;
  
  //Init physics parameters
  DYTurbo::init_params();

  //opts.dumpAll();  
  //DYTurbo::PrintTable::Settings();
  
  //Compute predictions
  vector <double> vals;
  vector <double> errs;
  DYTurbo::compute(vals, errs);

  //for (uint i = 0; i < vals.size(); i++)
  //cout << i << "  " << vals.size() << "  " << vals[i] << "  " << errs[i] << endl;

  //TODO: check bins size
  
  //insert values into output array
  copy_n(vals.begin(), val.size(), &val[0]);
}
