
#pragma once

#include "ReactionBaseDISNC.h"

/**
   @class' ReactionRT_DISNC

   @brief A wrapper class for RT_DISNC reaction

   Based on the ReactionTheory class. Reads options produces 3d cross section.

   @version 0.1
   @date 2017-04-10
*/

class ReactionRT_DISNC : public ReactionBaseDISNC
{
private:
   typedef ReactionBaseDISNC Super;

public:
   ReactionRT_DISNC(){};

public:
   virtual string getReactionName() const override { return "RT_DISNC"; };
   virtual void atStart() override;
   // virtual void setDatasetParameters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) override ;
   virtual void initTerm(TermData *td) override;
   virtual void atIteration() override;
   virtual void compute(TermData *, valarray<double> &val, map<string, valarray<double>> &err) override;

protected:
   virtual void F2 BASE_PARS override;
   virtual void FL BASE_PARS override;

   virtual void F2gamma_RT BASE_PARS;
   virtual void FLgamma_RT BASE_PARS;

private:
   map<unsigned, valarray<double>> _f2rt;
   map<unsigned, valarray<double>> _flrt;

   void calcF2FL(TermData *td);
};
