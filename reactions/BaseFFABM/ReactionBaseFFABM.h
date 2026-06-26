
#pragma once
#include "ABM.h"
#include <map>
#include <valarray>
#include <vector>
#include <string>
//#include "cuba.h"

using std::map;
using std::valarray;
using std::string;
using std::vector;

class DIS_HT;
class DIS_TMC;
class DIS_NUKE;
class TermData;

extern "C" {
  double numufcalflux_(const double& e); // NOMAD E(nu) flux
  double sd2_(double* acc, double (*f)(double*), void (*r)(int, double*, double*, double*)); // integrator
}

template <class BaseDIS, abm::SFproc Proc>
class ReactionBaseFFABM : public BaseDIS {
public:
  ~ReactionBaseFFABM();
  virtual void initTerm(TermData *td) override;
  virtual void atIteration() override;
  virtual void compute(TermData *td, valarray<double> &val, map<string, valarray<double>> &errors) { return BaseDIS::compute(td, val, errors); };

protected:
  map<int, valarray<double>> _f2abm;
  map<int, valarray<double>> _flabm;
  map<int, valarray<double>> _f3abm;

private:
  // pointers to those parameters which can change at each iteration
  const double* _mcPtr;
  const double* _mbPtr;
  std::vector<const double*> _ckm; // CKM matrix

  // minimum Q^2
  double _q2mincomp;

  // higher twist
  map<int, DIS_HT*> _ht;
  // target mass corrections
  map<int, DIS_TMC*> _tmc;
  // nuclear corrections
  map<int, DIS_NUKE*> _nuke;

  enum class dataFlav
  {
    incl,
    c,
    b,
    l
  }; //!< Define final state.
  virtual const dataFlav GetDataFlav(unsigned termID) { return dataFlav(BaseDIS::GetDataFlav(termID)); };
  virtual const valarray<double> *GetBinValues(TermData *td, const string &binName) { return BaseDIS::GetBinValues(td, binName); };
  virtual const double GetPolarisation(unsigned termID) { return BaseDIS::GetPolarisation(termID); };
  virtual const double GetCharge(unsigned termID) { return BaseDIS::GetCharge(termID); };
  enum class Binning {
    point_at_q2x, // requires q2, x
    point_at_e, // requires onedimvar
    point_at_x, // requires onedimvar
    point_at_sqrtshat, // requires onedimvar
    bin_q2y, // requires q2min, q2max, ymin, ymax (ymin is optional)
  };
  enum class Integrator {
    sd2,
    //cuba,
  };
  struct DataPoint;
  struct integration_params {
    DataPoint* point;
    double (*eflux)(const double& e);
    double emin;
    double emax;
    unsigned long ncalls;
  };
  struct DataPoint {
    Binning binning;
    double q2;
    double x;
    double onedimvar;
    double q2min;
    double q2max;
    double ymin;
    double ymax;
    double energy;
    int scalesemilepbr;
    Integrator integrator;
    double integrator_epsrel;
    int integrator_verbose;
    TermData* td;
    int i;
    dataFlav flav;
    bool is_beam_nu;
    int ord;
    int ordHQ;
    int ordFL;
    int msbarmin;
    double charge;
    double polar;
    const double* sin2thetaWPtr;
    const double* mz;
    const double* br0;
    const double* br1;
    double mnucl;
    double mw;
    DIS_NUKE* nuke;
    DIS_TMC* tmc;
    DIS_HT* ht;
    double f2;
    double fl;
    double f3;
    static integration_params _integration_params_static;
    static double integrand_sd2(double x[]);
    static int integrand(const int* ndim, const double* inp, const int *ncomp, double* val, void *pars_voidptr);
    static void integrand_sd2_region(int ll, double* xx, double* aa, double* bb);
    void calc();
    void calc_2d_integral();
    void calc_at_q2x();
  };
  std::map<std::string, std::vector<DataPoint> > _grouped_data_points;
  //std::map<std::string, ForkPool::SharedMemory* > _grouped_shared_memory;
};
