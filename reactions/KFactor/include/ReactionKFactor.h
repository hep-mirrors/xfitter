
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionKFactor

  @brief A wrapper class for KFactor reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-10-28
  */

class ReactionKFactor : public ReactionTheory
{
  public:
    ReactionKFactor(){};

    //    ~ReactionKFactor(){};
    //    ~ReactionKFactor(const ReactionKFactor &){};
    //    ReactionKFactor & operator =(const ReactionAKFactor &r){return *(new ReactionKFactor(r));};

  public:
    virtual string getReactionName() const { return  "KFactor" ;};
    int atStart(const string &);
    virtual void setDatasetParameters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) override ;
    virtual void initAtIteration() override;
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
  protected:
    virtual int parseOptions(){ return 0;};
  private:
    map<int, std::vector<double> > _values;
    map<int, std::pair<std::string, double> > _parameterNames;
};

