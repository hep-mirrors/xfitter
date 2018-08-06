
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionKMatrix

  @brief A wrapper class for KMatrix reaction 

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2018-08-03
  */

class ReactionKMatrix : public ReactionTheory
{
  public:
    ReactionKMatrix(){};

//    ~ReactionKMatrix(){};
//    ~ReactionKMatrix(const ReactionKMatrix &){};
//    ReactionKMatrix & operator =(const ReactionKMatrix &r){return *(new ReactionKMatrix(r));};

  public:
    virtual string getReactionName() const { return  "KMatrix" ;};
    int initAtStart(const string &); 
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
    virtual void setDatasetParameters(int dataSetID, map<string,string> pars, map<string, double> parsDataset) override;
  protected:
    virtual int parseOptions(){ return 0;};
  private:
    map<int, std::vector<std::vector<double>> > _values2D;
    map<int, std::vector<std::vector<double>> > _values;
    map<int, std::pair<std::string, double> > _parameterNames

	
};

