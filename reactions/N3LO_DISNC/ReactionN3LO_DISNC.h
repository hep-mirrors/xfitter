#pragma once

#include "ReactionBaseDISNC.h"
#include <apfel/apfelxx.h>
#include "FKTables.h"

namespace xfitter { class EvolutionAPFELxx; }

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
  double F2thr,F2asy,F2dmp;
  double FLthr,FLasy,FLdmp;

  map <unsigned,valarray<double>> _f2fonll;
  map <unsigned,valarray<double>> _flfonll;
  map <unsigned,valarray<double>> _f3fonll;

  std::unique_ptr<const apfel::Grid> Grid;
  const apfel::Grid*                 _grid = nullptr;  ///< Grid, or the APFELxx evolution grid in FK mode

  /// structure functions from the evolved PDFs (default)
  void computeFromPDFs();

  // FK mode (yaml "FK: 1"): rows built from the APFELxx evolution operators (FKTables.h)
  void buildFK();
  void computeFromFK();
  bool                        _fk = false;
  int                         _fkCheck = 1;      ///< compare with computeFromPDFs after each build
  int                         _fkVersion = -1;   ///< operator version of the rows
  xfitter::EvolutionAPFELxx*  _apfelxx = nullptr;
  struct FKRows { std::vector<xfitter::fk::Row> f2, fl, f3; std::vector<bool> active; };
  map <unsigned, FKRows>      _fkRows;
  std::vector<double>                Thresholds;
  std::vector<double>                Masses;
  std::vector<double>                Thresholds3F;
  std::vector<double>                Masses3F;
  std::vector<double>                Thresholds4F;
  std::vector<double>                Masses4F;
  
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> F2Obj;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> FLObj;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> F3Obj;

  // Massive coefficient functions in the 3-flavour scheme (for charm)
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> F2ObjM3F;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> FLObjM3F;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> F2ObjM03F;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> FLObjM03F;

  // Massive coefficient functions in the 4-flavour scheme (for bottom)
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> F2ObjM4F;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> FLObjM4F;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> F2ObjM04F;
  std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> FLObjM04F;
  
};

