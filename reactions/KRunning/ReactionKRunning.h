
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionKRunning

  @brief A wrapper class for KRunning reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2019-01-16
  */

class ReactionKRunning : public ReactionTheory
{
  public:
    ReactionKRunning(){};

//    ~ReactionKRunning(){};
//    ~ReactionKRunning(const ReactionKRunning &){};
//    ReactionKRunning & operator =(const ReactionKRunning &r){return *(new ReactionKRunning(r));};

  public:
    virtual string getReactionName() const { return  "KRunning" ;};
    int atStart(const string &);
    virtual void setDatasetParameters(int dataSetID, map<string,string> pars, map<string,double> dsPars) override;
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
  protected:
    virtual int parseOptions(){ return 0;};

    std::map<int, std::string> _type;
    std::map<int, std::string> _q;
    std::map<int, double> _qValue;
    std::map<int, std::string> _q0;
    std::map<int, int> _NPoints;

    double getAlphaS(double q) { return alphaS(q); }
    double getMassMSbar(const double m0, const double q, const double as, const double as0)
    {
      const double c0 = 4.0 / 9.0;
      // m0 is m(m)
      double mMsBar = m0 * pow(as / as0, c0);
      return mMsBar;
    }
};

