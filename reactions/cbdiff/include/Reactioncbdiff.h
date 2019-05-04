
#pragma once

#include "ReactionBaseHVQMNR.h"
//#include "ReactionTheory.h"

/**
  @class' Reactioncbdiff

  @brief A wrapper class for cbdiff reaction 

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2019-02-01
  */

class Reactioncbdiff : public ReactionBaseHVQMNR
//class Reactioncbdiff : public ReactionTheory
{
  public:
    Reactioncbdiff(){};

//    ~Reactioncbdiff(){};
//    ~Reactioncbdiff(const Reactioncbdiff &){};
//    Reactioncbdiff & operator =(const Reactioncbdiff &r){return *(new Reactioncbdiff(r));};

  public:
    virtual string getReactionName() const { return  "cbdiff" ;};
    virtual int initAtStart(const string &);
    virtual void initAtIteration();
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
    virtual void setDatasetParameters(int dataSetID, map<string,string> pars, map<string,double> dsPars) override;
  protected:
    virtual int parseOptions(){ return 0;};

    std::map<int, std::shared_ptr<MNR::MNR> > _mapMNR;
    std::map<int, std::shared_ptr<MNR::Grid> > _mapGrid;
    std::map<int, int> _mapMSbarMass;
    std::map<int, double> _mapMassDiff;
    std::map<int, std::shared_ptr<MNR::Grid> > _mapGridLOMassUp;
    std::map<int, std::shared_ptr<MNR::Grid> > _mapGridLOMassDown;
    std::map<int, std::shared_ptr<MNR::Grid> > _mapGridSmLOMassUp;
    std::map<int, std::shared_ptr<MNR::Grid> > _mapGridSmLOMassDown;
    std::map<int, std::shared_ptr<MNR::Grid> > _mapGridSm;
    std::map<int, std::shared_ptr<MNR::Frag> > _mapFrag;
    std::map<int, std::shared_ptr<Parameters> > _mapPar;
    std::map<int, std::vector<TH2D*> > _mapXSec;
    std::map<int, double> _mapFF;
    std::map<int, int> _mapN;
    std::map<int, double> _mapPrecision;
    std::map<int, double> _mapHadronMass;
};

