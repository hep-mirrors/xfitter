
#pragma once

#include "ReactionBaseDISNC.h"

/**
   @class' ReactionACOT

   @brief A wrapper class for ACOT reaction

   Based on the ReactionTheory class. Reads options produces 3d cross section.

   @version 0.1
   @date 2017-04-10
*/

class ReactionACOT : public ReactionBaseDISNC
{
private:
   typedef ReactionBaseDISNC Super;

public:
   ReactionACOT(){};

public:
   virtual string getReactionName() const override { return "ACOT"; };
   virtual void atStart() override;
   // virtual void setDatasetParameters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) override ;
   virtual void initTerm(TermData *td) override;
   virtual void atIteration() override;
   virtual void compute(TermData *, valarray<double> &val, map<string, valarray<double>> &err) override;

protected:
   virtual void F2 BASE_PARS override;
   virtual void FL BASE_PARS override;

   virtual void F2gamma_ACOT BASE_PARS;
   virtual void FLgamma_ACOT BASE_PARS;

private:
   map<unsigned, valarray<double>> _f2acot;
   map<unsigned, valarray<double>> _flacot;

   void calcF2FL(TermData *td);
};
