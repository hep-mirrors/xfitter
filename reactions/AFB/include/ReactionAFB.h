
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
    double integration_uubarEF (double, double);
    static double uubarEB_funct (double *, size_t, void *);
    double integration_uubarEB (double, double);
    static double uubarOF_funct (double *, size_t, void *);
    double integration_uubarOF (double, double);
    static double uubarOB_funct (double *, size_t, void *);
    double integration_uubarOB (double, double);

    static double ubaruEF_funct (double *, size_t, void *);
    double integration_ubaruEF (double, double);
    static double ubaruEB_funct (double *, size_t, void *);
    double integration_ubaruEB (double, double);
    static double ubaruOF_funct (double *, size_t, void *);
    double integration_ubaruOF (double, double);
    static double ubaruOB_funct (double *, size_t, void *);
    double integration_ubaruOB (double, double);

    static double ddbarEF_funct (double *, size_t, void *);
    double integration_ddbarEF (double, double);
    static double ddbarEB_funct (double *, size_t, void *);
    double integration_ddbarEB (double, double);
    static double ddbarOF_funct (double *, size_t, void *);
    double integration_ddbarOF (double, double);
    static double ddbarOB_funct (double *, size_t, void *);
    double integration_ddbarOB (double, double);

    static double dbardEF_funct (double *, size_t, void *);
    double integration_dbardEF (double, double);
    static double dbardEB_funct (double *, size_t, void *);
    double integration_dbardEB (double, double);
    static double dbardOF_funct (double *, size_t, void *);
    double integration_dbardOF (double, double);
    static double dbardOB_funct (double *, size_t, void *);
    double integration_dbardOB (double, double);

    double AFB (double, double);
};
