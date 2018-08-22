
#pragma once

#include "ReactionTheory.h"
#include "appl_grid/appl_grid.h"
#include <memory>

/**
  @class' ReactionAPPLgrid

  @brief A wrapper class for APPLgrid reaction 

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-03-28
  */

class ReactionAPPLgrid : public ReactionTheory
{
  public:
    ReactionAPPLgrid(){};

//    ~ReactionAPPLgrid(){};
//    ~ReactionAPPLgrid(const ReactionAPPLgrid &){};
//    ReactionAPPLgrid & operator =(const ReactionAAPPLgrid &r){return *(new ReactionAPPLgrid(r));};

 public:
    virtual string getReactionName() const { return  "APPLgrid" ;};
    int initAtStart(const string &); 
    virtual void setDatasetParameters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) override ;
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
 protected:
    virtual int parseOptions(){ return 0;};    

 private:
    enum class collision { pp, ppbar, pn};
    map<int, collision> _collType;
    map<int, std::vector<std::shared_ptr<appl::grid> > > _grids;
    map<int, int> _order;
    map<int, double> _muR, _muF; // !> renormalisation and factorisation scales
    map<int, bool> _flagNorm; // !> if true, multiply by bin width
    map<int, bool> _flagUseReference; // !> if true, prediction will be calculated from reference histogram (for tests and grids validation)
    map<int, std::vector<TH1D*> > _references;
    map<int, std::vector<double> > _eScale; // !> CMS energy
};

