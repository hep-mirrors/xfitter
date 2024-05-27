#pragma once
#include"ReactionTheory.h"
class pineappl_grid;
/**
  @class' ReactionPineAPPL

  @brief A wrapper class for PineAPPL reaction
  */
class ReactionPineAPPL:public ReactionTheory{
  std::map<std::string, pineappl_grid*> _initialized;
  std::map<std::string, std::pair<vector<double>, TermData*> > _convolved;
  std::vector<std::string> _convolved_vector_of_keys; // for parallel execuion
public:
  virtual string getReactionName() const override {return"PineAPPL";};
  virtual void initTerm(TermData*) override final;
  virtual void freeTerm(TermData*) override final;
  virtual void compute(TermData*,valarray<double>&,map<string,valarray<double> >&) override final;
  virtual void atIteration();
  int _pacount;
};
//NOTE: parameter "energy" is only read at start, once, and therefore cannot be fitted
