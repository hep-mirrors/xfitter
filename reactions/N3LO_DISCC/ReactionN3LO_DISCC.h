#pragma once

#include "ReactionBaseDISCC.h"
#include <apfel/apfelxx.h>

/*
 *   @class' ReactionN3LO_DISCC
 *
 *  @brief A wrapper class for N3LO_DISCC reaction
 *
 *  Based on the ReactionTheory class. Reads options produces 3d cross section.
 */

class ReactionN3LO_DISCC : public ReactionBaseDISCC
{
 public:
  ReactionN3LO_DISCC() {};
  //~ReactionN3LO_DISCC() {};
  //~ReactionN3LO_DISCC(const ReactionN3LO_DISCC &) {};
  //ReactionN3LO_DISCC & operator = (const ReactionAN3LO_DISCC &r) { return *(new ReactionN3LO_DISCC(r)); };

  virtual string getReactionName() const { return "N3LO_DISCC"; };
  virtual void atStart() override final;
  virtual void atIteration() override final;
  virtual void initTerm(TermData*)override final;
 
 protected:
  virtual valarray<double> F2(TermData*td) override;
  virtual valarray<double> FL(TermData*td) override;
  virtual valarray<double> xF3(TermData*td) override;

 private:
  map <unsigned,valarray<double>> _f2fonll;
  map <unsigned,valarray<double>> _flfonll;
  map <unsigned,valarray<double>> _f3fonll;
  map <unsigned,TermData*> _dsIDs;
  //! Allow for non-apfelxx evolution:
  bool _non_apfel_evol{false};

  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> F2PlusCCObj  ;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> F2MinusCCObj ;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> FLPlusCCObj  ;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> FLMinusCCObj ;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> F3PlusCCObj  ;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> F3MinusCCObj ;

  std::vector<double> Thresholds;
};

