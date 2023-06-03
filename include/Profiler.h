#pragma once

#include <valarray>
#include <string>
#include "xfitter_pars.h"

namespace xfitter
{
  /**
     @class Profiler

     @brief Class to transform variables to profiled systematic sources

     Contains methods to parse profiling request, to perform evaluation of 
     the varied predictions and to add the resulting source as extra systematic
     uncertainty

     @version 0.1
     @date 2018-08-23
  */

  class Profiler
  {
  public:

    /// default constructor.
    Profiler(){};

    /// parse yaml, perform profiling
    void doProfiling();

  private:
    /// calculate predictions
    std::pair <std::valarray<double>,double> evaluatePredictions();

    /// add symmetric systematic uncertainty with the name and corresponding variations
    void addSystematics( std::string const& name, std::valarray<double> uncertainties );

    /// add asymmetric systematic uncertainty with the name and corresponding variations
    void addSystematics( std::string const& name, std::valarray<double> uncertaintiesP,  std::valarray<double> uncertaintieM );

    /// profile single parameter which is already identified to be present on gParameters list.
    void profileParameter( std::string const& name, YAML::Node const& node) ;

    /// profile PDFerrors
    void profilePDF( std::string const& evolName, YAML::Node const& node) ;


    /// Dump aux. files
    void storePdfFiles(int imember, int iPDF=0, std::string const& type="central");

    /// convert MC replicas  to eigenvectors
    void addReplicas(std::string const& pdfName,  std::vector< std::valarray<double> > const& uncertainties);

    /// write out extra files to do MC reweighting
    void addMCweightsFiles(std::string const& pdfName, std::vector<double>& chi2vals, int ndata, int nrep);

    /// Compute using fork()
    void compute_parallel(int NALL, int NPRED, int first, int iPdfSet,
			std::vector< std::valarray<double> >& preds,
			std::vector< double >& chi2vals,
			  YAML::Node gNode, BaseEvolution* evol, const std::string& errorType);
    
    /// continuous nuisance parameter number for PDFs (if several are used)
    int _ipdf{0};

    /// compute chi2 for each variation
    bool _getChi2{false};

    /// Scale factor for PDF eigenvectors
    float _scalePdfs{1.0};

    /// keep info of the number of input syst. sources
    int _nSourcesOrig{0};

    /// store pdf files or not (for plotting)
    bool _storePdfs{false};

    /// output directory name
    string _outputDir{"output"};

    /// number of processes to use by "fork"
    int _ncpu{0};

  };
  
} //namespace xfitter

