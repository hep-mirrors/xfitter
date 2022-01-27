#ifndef xFitter_ReactionfastNLO
#define xFitter_ReactionfastNLO

#pragma once

#include "ReactionTheory.h"
#include "fastnlotk/fastNLOReader.h"
#include <map>
#include <vector>
#include "BaseEvolution.h"
#include "xfitter_steer.h"
using xfitter::BaseEvolution;

/**
  @class' ReactionfastNLO

  @brief A wrapper class for fastNLO reaction

  Based on the ReactionTheory class

  @version 0.3
  @date 2016-12-06
  */

class fastNLOReaction : public fastNLOReader {
 public:
  fastNLOReaction(string name, BaseEvolution* ev) : fastNLOReader(name) {evolution=ev;};
 protected:
  fastNLOReaction(string name) : fastNLOReader(name) {}; // not public!
  BaseEvolution* evolution=nullptr;
  double EvolveAlphas(double q ) const {
    return evolution->getAlphaS(q);
  } //!< provide alpha_s to fastNLO
  bool InitPDF() { return true;}; //!< required by fastNLO
  vector<double> GetXFX(double xp, double muf) const {
    vector<double> pdfV(14);
    evolution->xfxQarray(xp,muf,pdfV.data());
    return pdfV;
  }//!< provide PDFs to fastNLO
};

class ReactionfastNLO : public ReactionTheory {
  public:
    ReactionfastNLO(){};
  public:
    virtual string getReactionName() const override { return  "fastNLO" ;};
    virtual void initTerm(TermData* td) override final;
    virtual void freeTerm(TermData* td) override final;
    virtual void compute(TermData* td, valarray<double> &val, map<string, valarray<double> > &err) override;
};

#endif
