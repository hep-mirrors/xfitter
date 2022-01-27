
#pragma once

#include "ReactionTheory.h"
#include <MNRGrid.h>
#include <MNR.h>
#include <MNRFrag.h>
#include <TH2D.h>

/**
  @class' ReactionBaseHVQMNR

  @brief A wrapper class for BaseHVQMNR reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  This is abstract class from which implementations of HVQMNR
  calculations for particular datasets should be derived.

  @version 0.1
  @date 2017-01-02
  */

class ReactionBaseHVQMNR : public ReactionTheory
{
  // allow access to alphas routine of ReactionTheory.h
  //friend class MNR::MNR;

public:
  ReactionBaseHVQMNR();
  ~ReactionBaseHVQMNR();

public:
  virtual string getReactionName() const override { return  "BaseHVQMNR" ;};
  virtual void compute(TermData *, valarray<double> &val, map<string, valarray<double>> &errors) override = 0;
  virtual void initTerm(TermData *td) override;
protected:

  // ********** common stuff for MNR calculation  **********
protected:
  // structure for particular dataset
  struct DataSet
  {
    std::string FinalState;
    double FragFraction;
    int NormY; // 0 for absolute cross sections, 1 for normalised to reference rapidity bin
    // binning
    std::valarray<double>* BinsYMin;
    std::valarray<double>* BinsYMax;
    std::valarray<double>* BinsPtMin;
    std::valarray<double>* BinsPtMax;
    std::valarray<double>* BinsYMinRef;
    std::valarray<double>* BinsYMaxRef;
  };

  // structure to store theory parameters
  struct Parameters
  {
    // heavy-quark masses
    double mc = 0.0;
    double mb = 0.0;
    // flavour (used by cbdiff reaction)
    char flav;
    // flag that mass is taken from global parameters and should be updated at each iteration
    bool flagMassIsGlobal = false;
    // scale parameters
    double mf_A_c = 0.0;
    double mf_B_c = 0.0;
    double mf_C_c = 0.0;
    double mr_A_c = 0.0;
    double mr_B_c = 0.0;
    double mr_C_c = 0.0;
    double mf_A_b = 0.0;
    double mf_B_b = 0.0;
    double mf_C_b = 0.0;
    double mr_A_b = 0.0;
    double mr_B_b = 0.0;
    double mr_C_b = 0.0;
    // fragmentation parameters
    double fragpar_c = 0.0;
    double fragpar_b = 0.0;
    // divide by bin width
    bool flagDivideBinWidth = false;
    // verbose output for debugging
    bool debug = false;
  };

  // structure to store steering parameters
  struct Steering
  {
    int    nf;
    double ptmin;
    double ptmax;
    int    npt;
    int    nptsm;
    double ymin;
    double ymax;
    int    ny;
    int    nsfnb;
    int    nx3;
    int    nx4;
    int    nbz;
    double xmin;
    double xmax;
    double mf2min;
    double mf2max;
    bool   q;
    bool   a;
    char flav; // c, b or t
  };

  // all datasets
  std::map<int, DataSet> _dataSets;
  // theory parameters
  Parameters _pars;
  // debug lebel
  int _debug;
  // store term data for later access.
  map<unsigned, TermData*> _tdDS;

  // data members for calculation
  MNR::MNR _mnr;
  MNR::Frag _frag;
  MNR::Grid _grid, _gridSmoothed;
  // histograms with calculated cross sections
  std::vector<TH2D*> _hCalculatedXSec;
  // status flags
  bool _isInitAtStart;
  //int _ifcncount_last;
  // heavy-quark mass
  //std::map

  // check if appropriate heavy-flavour scheme is used
  void CheckHFScheme();

  // read and update theory parameters
  void UpdateParameters(TermData* td = NULL);

  // print theory parameters
  void PrintParameters(Parameters const* pars = NULL) const;

  // initialise calculation with default parameters
  void DefaultInit(const Steering& steer, const double mq, MNR::MNR& mnr, MNR::Frag& frag, MNR::Grid& grid, MNR::Grid& grid_smoothed);
  void DefaultInitMNR(const Steering& steer, const double mq, MNR::MNR& mnr);
  void DefaultInitGrid(const Steering& steer, const double mq, const int npt, MNR::Grid& grid);
  void DefaultInitFrag(const Steering& steer, MNR::Frag& frag);

  // return cross section in provided pT-y bin
  double FindXSecPtYBin(const TH2* histXSec, const double ymin, const double ymax, const double ptmin, const double ptmax, const bool diff_pt, const bool diff_y);

  //private:
  // check equality of float numbers with tolerance
  bool IsEqual(const double val1, const double val2, const double eps = 1e-6);

  // TODO this old commented out code to be removed one day
  /*// read values from terminfo in format key1=value1:key2=value2:...
    int readFromTermInfo(const std::string& str, const std::string& key, int& value);
    int readFromTermInfo(const std::string& str, const std::string& key, float& value);
    int readFromTermInfo(const std::string& str, const std::string& key, double& value);
    int readFromTermInfo(const std::string& str, const std::string& key, std::string& value);*/

  // read parameters for perturbative scales from MINUIT extra parameters
  void GetMuPar(TermData* td, const char mu, const char q, double& A, double& B, double& C);

  // read fragmentation parameter from MINUIT extra parameters
  double GetFragPar(TermData* td, const char q, const map<string,string> pars = map<string,string>());

  /*// check parameter respecting priority: (1) supplied map (if supplied), (2) global
    bool checkParamInPriority(const string& name, const std::map<string,string> pars = std::map<string,string>()) const
    {
      if(pars.size() != 0)
        return (pars.find(name) != pars.end());
      else
        return checkParam(name);
    }

    // get parameter respecting priority: (1) supplied map (if supplied), (2) global
    double GetParamInPriority(const string& name, const std::map<string,string> pars = std::map<string,string>()) const
    {
      if(pars.find(name) != pars.end())
        return std::stod(pars.at(name));
      else
        return GetParam(name);
    }

    // get parameter respecting priority: (1) supplied map (if supplied), (2) global
    int GetParamIInPriority(const string& name, const std::map<string,string> pars = std::map<string,string>()) const
    {
      if(pars.find(name) != pars.end())
        return std::stod(pars.at(name));
      else
        return GetParamI(name);
    }

    // get parameter respecting priority: (1) supplied map (if supplied), (2) global
    std::string GetParamSInPriority(const string& name, const std::map<string,string> pars = std::map<string,string>()) const
    {
      if(pars.find(name) != pars.end())
        return pars.at(name);
      else
        return GetParamS(name);
    }*/
};

