
#pragma once

#include "ReactionTheory.h"
#include "appl_grid/appl_grid.h"
#include <memory>
#include "BaseEvolution.h"

/**
  @class' ReactionAPPLgrid

  @brief A wrapper class for APPLgrid reaction 

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-03-28
  */
struct DatasetData{
  std::vector<std::shared_ptr<appl::grid> >grids;
  int order;
  double muR,muF; // !> renormalisation and factorisation scales
  bool flagNorm; // !> if true, multiply by bin width
  bool flagUseReference; // !> if true, prediction will be calculated from reference histogram (for tests and grids validation)
  std::vector<TH1D*>  references;
  std::vector<double> eScale; // !> CMS energy
  double*scaleParameter=nullptr; // !> pointer to a minimization parameter which by which the predicted cross-section will be additionally multiplied. If this pointer is nullptr, no additional scaling is used.
  xfitter::BaseEvolution*evolutions[2];
};
class ReactionAPPLgrid : public ReactionTheory
{
  public:
    ReactionAPPLgrid();
    ~ReactionAPPLgrid();
    virtual string getReactionName() const { return  "APPLgrid" ;};
    int initAtStart(const string &); 
    virtual void setDatasetParameters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) override ;
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
  protected:
    virtual int parseOptions(){ return 0;};
  private:
    map<int,DatasetData>dataset_data;
};

