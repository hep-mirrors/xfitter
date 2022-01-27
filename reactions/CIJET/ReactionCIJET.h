#ifndef xFitter_ReactionCIJET
#define xFitter_ReactionCIJET

#pragma once

#include "ReactionTheory.h"
#include "CIJETReader.h"
//#include <memory> 
#include <map>
#include <vector>

/**
  @class' ReactionCIJET

  @brief A wrapper class for CIJET reaction 

  Based on the ReactionTheory class

  @version 0.3
  @date 2016-12-06
  */

class CIJETReaction : public CIJETReader {
public:
    CIJETReaction(string name, ReactionTheory* reaction) : CIJETReader(name, reaction) {};
protected:
    CIJETReaction(string name) : CIJETReader(name) {}; // not public!
};

class ReactionCIJET : public ReactionTheory {
public:
    ReactionCIJET(){};


public:
    virtual string getReactionName() const override { return  "CIJET" ;};
    virtual void initTerm(TermData* td) override final;
    virtual void atStart() override {} //< nothing todo
    virtual void compute(TermData*, valarray<double> &val, map<string, valarray<double> > &err) override;
protected:
    virtual int parseOptions(){ return 0;};
    // CIJET inherited functions 
    std::map<int,std::vector<CIJETReaction*> > ffnlo;
};

#endif
