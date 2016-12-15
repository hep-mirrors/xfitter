#pragma once

#include <map>
#include <string>
#include <vector>
#include <valarray>

using std::map;
using std::string;
using std::vector;
using std::valarray;

// typedef double (*pxFx)(double*, double*);
typedef double (*pZeroParFunc)();
typedef double (*pOneParFunc)(double*);
typedef double (*pTwoParFunc)(double*, double*);
typedef void   (*pThreeParSub)(double*, double*, double*);  

// Function to emulate LHAPDF xfx behavior:
typedef void   (*pXFXlike)(double*, double*, double*);

/**
  @class ReactionTheory

  @brief A base class manages for reaction theories

  It provides an interface wich must present in the derived classes

  @author A.Sapronov <sapronov@ifh.de>

  @version 0.1
  @date 2016/01/21
  */

//class Evolution;

class ReactionTheory 
{
 public:
  ReactionTheory() {};
  ~ReactionTheory() {};

  ReactionTheory(const ReactionTheory &);
  ReactionTheory & operator =(const ReactionTheory &);

 public:
  virtual string getReactionName() const =0;  ///< Should return expected reaction name. Normally generated automatically by AddReaction.py
  virtual int  initAtStart(const string &) =0; ///< Initialization first time ReactionTheory implementation is called
  virtual void setxFitterParameters(map<string,double*> &xfitter_pars) {_xfitter_pars = xfitter_pars; }; ///< Set environment map
  virtual void setEvolFunctions(double (*palpha_S)(double *), map<string, pTwoParFunc> *func2D  ) { alpha_S = palpha_S; PDFs = func2D; }; 
				///< Set alpha_S and PDF maps
  virtual void setExtraFunctions(map<string, pZeroParFunc>, map<string, pOneParFunc>, map<string, pTwoParFunc>) { };
  virtual void initAtIteration() {};
  virtual void setXFX(pXFXlike xfx){ _xfx = xfx; };
  virtual void setBinning(int dataSetID, map<string,valarray<double> > *dsBins){ _dsIDs.push_back(dataSetID); _dsBins[dataSetID] = dsBins; } ;
  virtual void setDatasetParamters( int dataSetID, map<string,string> pars) {} ;
//  virtual void resultAt(valarray<double> *val){ _val = val; };  
  virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err) = 0;

  virtual void printInfo(){};

  // Helper functions to emmulate LHAPDF6 calls:
  void xfx(double x, double q, double* results){double q2=q*q; (*_xfx)(&x,&q2,&results[0]); };
  double xfx(double x, double q, int iPDF){ double pdfs[13]; xfx(x,q,pdfs); return pdfs[iPDF+6];};

 protected:

  virtual int parseOptions() { return 0;};
  double (*alpha_S)(double *);
  map<string, pTwoParFunc> *PDFs;

  // Check if a parameter is present on the list:
  bool checkParam(string name) 
  {
    return (_xfitter_pars.find(name) !=  _xfitter_pars.end());
  }

  // Helper function to get a parameter
  double GetParam(string name) const 
  {
    return *_xfitter_pars.at(name);
  }

  // Helper function to get bin values for a given data set, bin name. Returns null if not found
  valarray<double> *GetBinValues(int idDS, string binName)
  { 
    map<string, valarray<double> >* mapBins =  _dsBins[idDS];
    if (mapBins == NULL ) {
      return NULL;
    }
    else { 
      map<string, valarray<double> >::iterator binPair = mapBins->find(binName);
      if ( binPair == mapBins->end() ) {
	return NULL;
      }
      else {
	return &binPair->second;
      }
    }
  };

 protected:
  string _subtype;
  valarray<double> *_val;
  string _ro;
  /// dataset bins, map based on datasetID (integer)
  /// must contain 'binFlag' key
  vector<int> _dsIDs;
  map<int, map<string, valarray<double> >* > _dsBins;
  map<string, double* > _xfitter_pars;
 private:
  pXFXlike _xfx;

};


typedef ReactionTheory * create_t();
