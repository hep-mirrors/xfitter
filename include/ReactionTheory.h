#pragma once

#include <map>
#include <string>
#include <vector>
#include <valarray>
#include <functional>

#include "xfitter_cpp_base.h"

using std::map;
using std::string;
using std::vector;
using std::valarray;

// typedef double (*pxFx)(double*, double*);


typedef double (*pZeroParFunc)();
typedef double (*pOneParFunc)(const double&);
typedef double (*pTwoParFunc)(const double&, const double& );
typedef void   (*pThreeParSub)(const double& , const double&, const double&);  

// Function to emulate LHAPDF xfx behavior:
typedef void   (*pXFXlike)(const double&, const double&, double*);

//using pZeroParFunc = std::function< double() >;
//using pOneParFunc  = std::function< double(const double&) >;
//using pTwoParFunc  = std::function< double(const double&, const double&) >;
//using pXFXlike     = std::function< void(const double& , const double&, double* ) >;

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
  virtual ~ReactionTheory() {};

  ReactionTheory(const ReactionTheory &);
  ReactionTheory & operator =(const ReactionTheory &);

 public:
  virtual string getReactionName() const =0;  ///< Should return expected reaction name. Normally generated automatically by AddReaction.py
  virtual int  initAtStart(const string &) =0; ///< Initialization first time ReactionTheory implementation is called

  virtual void setxFitterParameters(map<string,double*> &xfitter_pars) {_xfitter_pars = xfitter_pars; }; ///< Set environment map
  virtual void setxFitterParametersI(map<string,int> &xfitter_pars) {_xfitter_pars_i = xfitter_pars; }; ///< Set environment map
  virtual void setxFitterParametersS(map<string,string> &xfitter_pars) {_xfitter_pars_s = xfitter_pars; }; ///< Set environment map

  virtual void setEvolFunctions(double (*palpha_S)(const double& ), map<string, pTwoParFunc> *func2D  ) { _alpha_S = palpha_S; PDFs = func2D; }; 
				///< Set alpha_S and PDF maps
  virtual void setExtraFunctions(map<string, pZeroParFunc>, map<string, pOneParFunc>, map<string, pTwoParFunc>) { };
  virtual void initAtIteration() {};
  virtual void setXFX(pXFXlike xfx, string type="p" ){ _xfx[type] = xfx; };

  
  virtual void setBinning(int dataSetID, map<string,valarray<double> > *dsBins){ _dsIDs.push_back(dataSetID); _dsBins[dataSetID] = dsBins; } ;

  //! Set dataset @param dataSetID parameters which can be term- and dataset-specific
  virtual void setDatasetParamters( int dataSetID, map<string,string> parsReaction,  map<string,double> parsDataset) {} ;

  //! Main function to compute predictions for @param dataSetID 
  virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err) = 0;  
 
  //! Provide additional optional information about the reaction
  virtual void printInfo(){};

  //! Helper function to emmulate LHAPDF6 calls to get PDFs
  void xfx(const double& x, const double& q, double* results){ (_xfx["p"])(x,q,results); };
  
  //!  Helper function to emmulate LHAPDF6 calls to get PDFs
  double xfx(double x, double q, int iPDF){ double pdfs[13]; xfx(x,q,pdfs); return pdfs[iPDF+6];};
  
  //! strong coupling at scale q [GeV]
  double alphaS(double q) { return _alpha_S(q); }

  //! Return pointer-function to XFX for external use
  const pXFXlike getXFX(string type="p") { return _xfx[type];};

  //!  Return pointer-function to alphaS for external use
  const pOneParFunc getAlphaS() { return _alpha_S;}

  //! Default helper to determine if bin is masked or not
  virtual bool notMasked(int DSID, int Bin);
 
  //! Helper function to report not implemented functionality 
  void NOT_IMPLEMENTED(const string& functionName) {
    string message = "F: Function "+functionName+" is not implemented for module" + getReactionName();
    hf_errlog_(17040701,message.c_str(),message.size());
  } 

 protected:

  virtual int parseOptions() { return 0;};
  map<string, pTwoParFunc> *PDFs;

  bool checkParam(string name)         ///< Check if a parameter is present on the list
  {
    return (_xfitter_pars.find(name) !=  _xfitter_pars.end()) 
      || (_xfitter_pars_i.find(name) !=  _xfitter_pars_i.end()) 
      || (_xfitter_pars_s.find(name) !=  _xfitter_pars_s.end()) 
      ;
  }

  // Helper function to get a parameter (double)
  double GetParam(string name) const
  {    
    if (_xfitter_pars.find(name) != _xfitter_pars.end() ) 
      return *_xfitter_pars.at(name);
    else
      return 0;
  }

  // Helper function to get a parameter (integer)
  int GetParamI(string name)  const
  {
    return _xfitter_pars_i.at(name);
  }

  // Helper function to get a parameter (string)
  string GetParamS(string name)  const
  {
    return _xfitter_pars_s.at(name);
  }


  // Helper function to get bin values for a given data set, bin name. Returns null if not found
  valarray<double> *GetBinValues(int idDS, string binName)
  { 
    map<string, valarray<double> >* mapBins =  _dsBins[idDS];
    if (mapBins == nullptr ) {
      return nullptr;
    }
    else { 
      map<string, valarray<double> >::iterator binPair = mapBins->find(binName);
      if ( binPair == mapBins->end() ) {
	return nullptr;
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
  /// GLobal double parameters:
  map<string, double* > _xfitter_pars;
  /// GLobal integer parameters:
  map<string, int > _xfitter_pars_i;
  /// GLobal string parameters:
  map<string, string > _xfitter_pars_s;

 private:
  map<string,pXFXlike> _xfx;
  pOneParFunc _alpha_S;

};


typedef ReactionTheory * create_t();
