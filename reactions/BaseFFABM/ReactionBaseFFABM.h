
#pragma once
#include "ReactionTheory.h"
#include "ABM.h"
#include "ForkPool.h"
//#include "cuba.h"

class DIS_HT;
class DIS_TMC;
class DIS_NUKE;

class ReactionBaseFFABM : public ReactionTheory
{
private:
  typedef ReactionTheory Super;

public:
  ReactionBaseFFABM(){};

public:
  //virtual string getReactionName() const override { return "BaseFFABM"; };
  //void virtual atStart() override;
  virtual void initTerm(TermData *td) override;
  virtual void atIteration() override;

protected:
  virtual abm::SFproc GetProc() = 0;
  map<int, valarray<double>> _f2abm;
  map<int, valarray<double>> _flabm;
  map<int, valarray<double>> _f3abm;
  //virtual valarray<double> F2(TermData *td) override final;
  //virtual valarray<double> FL(TermData *td) override final;
  //virtual valarray<double> xF3(TermData *td) override final;
  //virtual valarray<double> F2(TermData *td);
  //virtual valarray<double> FL(TermData *td);
  //virtual valarray<double> xF3(TermData *td);
  enum class dataFlav
  {
    incl,
    c,
    b,
    l
  }; //!< Define final state.
  virtual const dataFlav GetDataFlav(unsigned termID) = 0;

private:
  // pointers to those parameters which can change at each iteration
  const double* _mcPtr;
  const double* _mbPtr;

  double _q2mincomp;

  virtual const valarray<double> *GetBinValues(TermData *td, const string &binName) = 0;
  virtual const double GetPolarisation(unsigned termID) = 0;
  virtual const double GetCharge(unsigned termID) = 0;
  // higher twist
  //map<int, DIS_HT*> _ht;
  // target mass corrections
  //map<int, DIS_TMC*> _tmc;
  // nuclear corrections
  //map<int, DIS_NUKE*> _nuke;
  virtual DIS_HT* getHT(unsigned termID) = 0;
  virtual DIS_TMC* getTMC(unsigned termID) = 0;
  virtual DIS_NUKE* getNUKE(unsigned termID) = 0;
  struct DataPoint {
    TermData* td;
    int datasetID;
    int i;
    abm::SFproc proc_NCCC;
    //ReactionBaseDISCC::dataFlav flav;
    dataFlav flav;
    int ord;
    int ordHQ;
    int ordFL;
    int msbarmin;
    double charge;
    double polar;
    const double* sin2thetaWPtr;
    const double* mz;
    double q2;
    double x;
    DIS_NUKE* nuke;
    DIS_TMC* tmc;
    DIS_HT* ht;
    double f2;
    double fl;
    double f3;
    void calc();
  };
  ForkPool::TaskDistribution _task_distr; // 0 is chunky, 1 is cyclic
  std::map<std::string, std::vector<DataPoint> > _grouped_data_points;
};
