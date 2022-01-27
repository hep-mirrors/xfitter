#pragma once

#include "ReactionTheory.h"

/**
    @class' ReactionHathorMSR
  
    @brief A wrapper class for HathorMSR reaction
  
    Based on the ReactionTheory class. Reads options produces 3d cross section.
  
    @version 0.1
    @date 2019-06-04
     Author: Toni Mäkelä <toni.maekelae@desy.de>.
     Based on ReactionHathorSingleTop by Laia Parets Peris <laia.parets.peris@desy.de>, Katerina Lipka <katerina.lipka@desy.de>
     transition from pole to MSbar scheme by S. Moch (private communication), 
     transition from pole to MSR and MSR/MSBAR mass evolution by T. Mäkelä
  */

class Hathor;
class HathorPdfxFitter;
class ReactionHathorMSR : public ReactionTheory
{
public:
  ReactionHathorMSR();
  ~ReactionHathorMSR();

  virtual string getReactionName() const { return  "HathorMSR" ;};
  virtual void initTerm(TermData *td) override final;
  virtual void atStart();
  virtual void compute(TermData *td, valarray<double> &val, map<string, 
                       valarray<double> > &err);

  //Decoupling coefficients
  double d1func(double Lmu);
  double d2func(double Lmu);

  vector<double> asFactors(Hathor *XS, double muOLD, double muNEW);
  double evoInt(Hathor* XS, double (*integrand)(Hathor*,double,int), 
                double R0, double R1, int n);
  
protected:
  virtual int parseOptions(){ return 0;};

  // this is map of key = dataset, value = pointer to Hathor instances,
  // one instance per one dataset
  std::map<int, Hathor*> _hathorArray;

  HathorPdfxFitter* _pdf;
  int* _rndStore;
  map<unsigned, int> _scheme;
  map<unsigned, double> _mtop;
  map<unsigned, double> _mr;
  map<unsigned, double> _mf;
  // store term data for later access
  map<unsigned, TermData*> _tdDS;
  
  double nfl;         //#[active light flavours]
  int    orderI;      //0 LO, 1 NLO, 2 NNLO
  int    mScheme;     //0 pole, 1 MSbar, 2 MSRN, 3 MSRP
  int    mScheme_in;  //Scheme choice in input. MSR changes to MSbar if R>mt(mt)
  bool   convertMass; //0 mass input as MSR; 1 input mass MSbar, convert to MSR
  double Rscale, Rscale_in;
  double beta0, beta1, beta2, bar0, bar1;
  
  //Constants
  double pi        = 3.141592653589793;
  double const z2  = 1.644934066848226;
  double const z3  = 1.202056903159594;
  double const ln2 = 0.693147180559945;

};

