
#pragma once

#include "ReactionTheory.h"
//#include <memory> 
#include "fastnlotk/fastNLOReader.h" // this is still the old interface

/**
  @class' ReactionfastNLO

  @brief A wrapper class for fastNLO reaction 

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2016-12-06
  */

class ReactionfastNLO : public ReactionTheory, public fastNLOReader {
  public:
    ReactionfastNLO(){};

//    ~ReactionfastNLO(){};
//    ~ReactionfastNLO(const ReactionfastNLO &){};
//    ReactionfastNLO & operator =(const ReactionAfastNLO &r){return *(new ReactionfastNLO(r));};

  public:
    virtual string getReactionName() const { return  "fastNLO" ;};
    int  initAtStart(const string &); 
    virtual void setDatasetParamters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) override ;
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
  protected:
   virtual int parseOptions(){ return 0;};
   //std::unique_ptr<FastNLOxFitter> fnlo;

   // fastNLO inherited functions 
   double EvolveAlphas(double Q ) const ;
   bool InitPDF();
   vector<double> GetXFX(double xp, double muf) const ;

};

