#pragma once

#include "ReactionBaseDISCC.h"
#include <apfel/apfelxx.h>
#include "FKTables.h"

namespace xfitter { class EvolutionAPFELxx; }

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

  std::unique_ptr<const apfel::Grid> Grid;
  const apfel::Grid*                 _grid = nullptr;  ///< Grid, or the APFELxx evolution grid in FK mode

  /// structure functions from the evolved PDFs (default)
  void computeFromPDFs();

  // FK mode (yaml "FK: 1"), as N3LO_DISNC
  void buildFK();
  void computeFromFK();
  bool                        _fk = false;
  int                         _fkCheck = 1;
  int                         _fkVersion = -1;
  xfitter::EvolutionAPFELxx*  _apfelxx = nullptr;
  struct FKRows { std::vector<xfitter::fk::Row> f2, fl, f3; std::vector<bool> active; };
  map <unsigned, FKRows>      _fkRows;
  std::vector<double>                Thresholds;
  
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> F2PlusCCObj  ;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> F2MinusCCObj ;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> FLPlusCCObj  ;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> FLMinusCCObj ;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> F3PlusCCObj  ;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> F3MinusCCObj ;

};

