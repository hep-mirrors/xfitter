#ifndef xFitter_ReactionfastNLO
#define xFitter_ReactionfastNLO

#pragma once

#include "ReactionTheory.h"
//#include "fastNLOReaction.h"
//#include <memory> 
#include "fastnlotk/fastNLOReader.h" 
#include <map>
#include <vector>

/**
  @class' ReactionfastNLO

  @brief A wrapper class for fastNLO reaction 

  Based on the ReactionTheory class

  @version 0.3
  @date 2016-12-06
  */

class fastNLOReaction : public fastNLOReader {
public:
   fastNLOReaction(string name, ReactionTheory* reaction) : fastNLOReader(name) {fReaction=reaction;};
protected:
   fastNLOReaction(string name) : fastNLOReader(name) {}; // not public!
   ReactionTheory* fReaction=NULL;
   double EvolveAlphas(double Q ) const { return fReaction->alphaS(Q); } //!< provide alpha_s to fastNLO      
   bool InitPDF() { return true;}; //!< required by fastNLO
   vector<double> GetXFX(double xp, double muf) const { return fReaction->xfx(xp,muf); }//!< provide PDFs to fastNLO
};

class ReactionfastNLO : public ReactionTheory {
public:
    ReactionfastNLO(){};

//    ~ReactionfastNLO(){};
//    ~ReactionfastNLO(const ReactionfastNLO &){};
//    ReactionfastNLO & operator =(const ReactionAfastNLO &r){return *(new ReactionfastNLO(r));};

public:
   virtual string getReactionName() const { return  "fastNLO" ;};
   int  initAtStart(const string &) { return 0; } //< nothing todo
   virtual void setDatasetParameters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) override ;
   virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
protected:
   virtual int parseOptions(){ return 0;};
   //std::unique_ptr<FastNLOxFitter> fnlo;
   // fastNLO inherited functions 
   std::map<int,std::vector<fastNLOReaction*> > ffnlo;
};

#endif
