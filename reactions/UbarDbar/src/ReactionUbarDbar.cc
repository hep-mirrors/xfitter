
/*
   @file ReactionUbarDbar.cc
   @date 2019-09-05
   @author  AddReaction.py
   Created by  AddReaction.py on 2019-09-05
*/

#include "ReactionUbarDbar.h"
#include "BaseEvolution.h"

// the class factories
extern "C" ReactionUbarDbar* create() {
  return new ReactionUbarDbar();
}


// Main function to compute results at an iteration
void ReactionUbarDbar::compute(TermData*td, valarray<double> &val, map<string, valarray<double> > &err)
{
  // get data bins:
  auto &xbins = td->getBinColumn("x");
  auto &qbins = td->getBinColumn("q");
  // get evolution
  xfitter::BaseEvolution *evol = td->getPDF();
  // compute ubar, dbar for all x values, 
  for ( size_t i = 0; i<xbins.size(); i++) {
    double x = xbins[i];
    double q = qbins[i];
    double ubar = evol->xfxQmap(x,q)[-2];
    double dbar = evol->xfxQmap(x,q)[-1];
    val[i] = ubar-dbar;
  }
}

