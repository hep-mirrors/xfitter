#pragma once

#include "ReactionBaseDISNC.h"
#include <apfel/apfelxx.h>

/*
 *   @class' ReactionN3LO_DISNC
 *
 *  @brief A wrapper class for N3LO_DISNC reaction
 *
 *  Based on the ReactionTheory class. Reads options produces 3d cross section.
 */

class ReactionN3LO_DISNC : public ReactionBaseDISNC
{
 public:
  ReactionN3LO_DISNC() {};
  //~ReactionN3LO_DISNC() {};
  //~ReactionN3LO_DISNC(const ReactionN3LO_DISNC &) {};
  //ReactionN3LO_DISNC & operator = (const ReactionAN3LO_DISNC &r) { return *(new ReactionN3LO_DISNC(r)); };

  virtual string getReactionName() const { return "N3LO_DISNC"; };
  virtual void atStart() override final;
  virtual void atIteration() override final;
  virtual void initTerm(TermData*)override final;
 
 protected:
  virtual void F2  BASE_PARS override final;
  virtual void FL  BASE_PARS override final;
  virtual void xF3 BASE_PARS override final;

 private:
  map <unsigned,valarray<double>> _f2fonll;
  map <unsigned,valarray<double>> _flfonll;
  map <unsigned,valarray<double>> _f3fonll;

  std::unique_ptr<const apfel::Grid> Grid;
  std::vector<double>                Thresholds;
  
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> F2Obj;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> FLObj;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> F3Obj;
};

