
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionAFB

  @brief A wrapper class for AFB reaction 

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2018-07-16
  */


class ReactionAFB : public ReactionTheory
{
  public:
    ReactionAFB(){};

//    ~ReactionAFB(){};
//    ~ReactionAFB(const ReactionAFB &){};
//    ReactionAFB & operator =(const ReactionAFB &r){return *(new ReactionAFB(r));};

  public:
    virtual string getReactionName() const { return  "AFB" ;};
    int initAtStart(const string &); 
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
  protected:
    virtual int parseOptions(){ return 0;};

  private:
    static double *propagators (double);
    static double uubarEF_funct (double *, size_t, void *);
//     double uubarEF_funct (double *, size_t, void *);
    
    double integration_uubarEF (double, double);
    double *observables (double, double);
    
    
//   public:
//     double *propagators (double Minv);
//     double uubarEF_funct (double *entries, size_t dim, void *params);
//     double integration_uubarEF (double Minv_inf, double Minv_sup);
    
    
//     double uubarEB_funct (double *entries, size_t dim, void *params);
//     double integration_uubarEB (double Minv_inf, double Minv_sup);
//     double uubarOF_funct (double *entries, size_t dim, void *params);
//     double integration_uubarOF (double Minv_inf, double Minv_sup);
//     double uubarOB_funct (double *entries, size_t dim, void *params);
//     double integration_uubarOB (double Minv_inf, double Minv_sup);
// 
//     double ubaruEF_funct (double *entries, size_t dim, void *params);
//     double integration_ubaruEF (double Minv_inf, double Minv_sup);
//     double ubaruEB_funct (double *entries, size_t dim, void *params);
//     double integration_ubaruEB (double Minv_inf, double Minv_sup);
//     double ubaruOF_funct (double *entries, size_t dim, void *params);
//     double integration_ubaruOF (double Minv_inf, double Minv_sup);
//     double ubaruOB_funct (double *entries, size_t dim, void *params);
//     double integration_ubaruOB (double Minv_inf, double Minv_sup);
// 
//     double ddbarEF_funct (double *entries, size_t dim, void *params);
//     double integration_ddbarEF (double Minv_inf, double Minv_sup);
//     double ddbarEB_funct (double *entries, size_t dim, void *params);
//     double integration_ddbarEB (double Minv_inf, double Minv_sup);
//     double ddbarOF_funct (double *entries, size_t dim, void *params);
//     double integration_ddbarOF (double Minv_inf, double Minv_sup);
//     double ddbarOB_funct (double *entries, size_t dim, void *params);
//     double integration_ddbarOB (double Minv_inf, double Minv_sup);
// 
//     double dbardEF_funct (double *entries, size_t dim, void *params);
//     double integration_dbardEF (double Minv_inf, double Minv_sup);
//     double dbardEB_funct (double *entries, size_t dim, void *params);
//     double integration_dbardEB (double Minv_inf, double Minv_sup);
//     double dbardOF_funct (double *entries, size_t dim, void *params);
//     double integration_dbardOF (double Minv_inf, double Minv_sup);
//     double dbardOB_funct (double *entries, size_t dim, void *params);
//     double integration_dbardOB (double Minv_inf, double Minv_sup);
// 
//     double observables (double Minv_inf, double Minv_sup);

};




