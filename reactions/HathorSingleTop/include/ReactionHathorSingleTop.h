
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionHathorSingleTop

  @brief A wrapper class for HathorSingleTop reaction 

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2018-07-25
  */

class HathorSgTopT;
class HathorPdfxFitter;

class ReactionHathorSingleTop : public ReactionTheory
{
  public:
    ReactionHathorSingleTop();

    ~ReactionHathorSingleTop();

//    ~ReactionHathorSingleTop(){};
//    ~ReactionHathorSingleTop(const ReactionHathorSingleTop &){};
//    ReactionHathorSingleTop & operator =(const ReactionHathorSingleTop &r){return *(new ReactionHathorSingleTop(r));};

  public:
    virtual string getReactionName() const { return  "HathorSingleTop" ;};
    int initAtStart(const string &); 
     virtual void setDatasetParameters(int dataSetID, map<string,string> pars, map<string,double> dsPars) override;
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
    double CalculationRatio(double _mtop, double _mr, double _mf, int _scheme, HathorSgTopT* hathor);
  protected:
    virtual int parseOptions(){ return 0;};

     // this is map of key = dataset, value = pointer to Hathor instances,
    // one instance per one dataset
    std::map<int, HathorSgTopT*> _hathorArray;

    HathorPdfxFitter* _pdf;
    int* _rndStore;
    int _scheme;
    double _mtop;
    double _mr;
    double _mf;
};

