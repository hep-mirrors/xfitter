
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
    virtual string getReactionName() const { return  "BaseHVQMNR" ;};
    virtual int initAtStart(const string &) = 0; 
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err) = 0;
    virtual void initAtIteration() = 0;
    virtual void setDatasetParamters(int dataSetID, map<string,string> pars, map<string,double> dsPars) override;
  protected:
    virtual int parseOptions(){ return 0;};
    
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
      double mc, mb;
      // scale parameters
      double mf_A_c, mf_B_c, mf_C_c;
      double mr_A_c, mr_B_c, mr_C_c;
      double mf_A_b, mf_B_b, mf_C_b;
      double mr_A_b, mr_B_b, mr_C_b;
      // fragmentation parameters
      double fragpar_c, fragpar_b;
    };

    // structure to store steering parameters
    struct Steering
    {
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
    };
    
    // all datasets
    std::map<int, DataSet> _dataSets;
    // theory parameters
    Parameters _pars;
    // debug lebel
    int _debug;

    // data members for calculation
    MNR::MNR _mnr;
    MNR::Frag _frag;
    MNR::Grid _grid, _gridSmoothed;
    // histograms with calculated cross sections
    std::vector<TH2D*> _hCalculatedXSec;
    // status flags
    bool _isInitAtStart;
    int _ifcncount_last;

    // check if appropriate heavy-flavour scheme is used
    void CheckHFScheme();

    // read and update theory parameters
    void UpdateParameters();
    
    // print theory parameters
    void PrintParameters() const;

    // initialise calculation with default parameters
    void DefaultInit(const Steering& steer, const double mq, MNR::MNR& mnr, MNR::Frag& frag, MNR::Grid& grid, MNR::Grid& grid_smoothed);
    
    // return cross section in provided pT-y bin
    double FindXSecPtYBin(const TH2* histXSec, const double ymin, const double ymax, const double ptmin, const double ptmax, const bool diff_pt, const bool diff_y);
    
    // error treatment
    void error(const int id, const std::string& str);

  private:    
    // check equality of float numbers with tolerance
    bool IsEqual(const double val1, const double val2, const double eps = 1e-6);
    
    // TODO this old commented out code to be removed one day
    /*// read values from terminfo in format key1=value1:key2=value2:...
    int readFromTermInfo(const std::string& str, const std::string& key, int& value);
    int readFromTermInfo(const std::string& str, const std::string& key, float& value);
    int readFromTermInfo(const std::string& str, const std::string& key, double& value);
    int readFromTermInfo(const std::string& str, const std::string& key, std::string& value);*/

    // read parameters for perturbative scales from MINUIT extra parameters
    void GetMuPar(const char mu, const char q, double& A, double& B, double& C);

    // read fragmentation parameter from MINUIT extra parameters
    double GetFragPar(const char q);    
};

