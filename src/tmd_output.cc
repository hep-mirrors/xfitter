#include<iostream>
#include<vector>
#include<algorithm>
#include<numeric>
#include<cmath>
#include<numeric>
#include<sys/stat.h>
#include<sys/types.h>
#include"xfitter_pars.h"
#include"xfitter_steer.h"
#include"xfitter_cpp_base.h"
#include"BaseEvolution.h"
#include"BaseMinimizer.h"
using namespace std;
using xfitter::BaseEvolution;



void WriteTMD(const char* f, BaseEvolution* ev, const char* const PdfType = "central")
{
  //fprintf(f, "PdfType: %s\n", PdfType);
  //fprintf(f, "Format: lhagrid1\n---\n");
  //std::cout << " in WriteTMD " << f << std::endl;
  // char* name = f ;
  ev -> Write_TMD(f) ;
  //std::cout << " after Write_TMD " << std::endl;
}

struct TMD_Options {
  BaseEvolution* pdf;
  string
  name = "TMD",
  description = "Generated using xFitter",
  setindex = "....",
  authors = "",
  reference = "",
  evolution = "",
  flavor_scheme = "",
  tmd_scheme = "" ,
  format = "",
  extrapolation_x = "",
  extrapolation_q2 = "",
  extrapolation_kt = "",
  error_type = "";
  double qmin, qmax, xmin, xmax, ktmin, ktmax;
  size_t Nx = 0, Nq = 0;
  int nmembers = 1, numflavors=5,particle=2212;

  
  bool prefer_internal_grid; //if true, try to use internal grid of evolution, if it can provide one
  static TMD_Options fromYAML(YAML::Node);
};
//Parse control block in YAML steering and fill TMD_Options
TMD_Options TMD_Options::fromYAML(YAML::Node node)
{
  TMD_Options info;

  //ranges are uninitialized before reading options
  info.xmin = info.xmax = info.qmin = info.qmax = nan("");
  info.prefer_internal_grid = false;

  //Process option "evolution" first to be able to use its name in errors
  BaseEvolution* pdf = xfitter::get_evolution(node["evolution"].as<string>(""));
  info.pdf = pdf;

  auto const& pdfMap = pdf->xfxQmap(0.1, 10.);
  cout << " tmdoutput: Options: size xfxQmap "<< pdfMap.size() << endl;

  cout << "[INFO] WriteTMD: using " << pdf->_name << endl;


  bool grid_provided = false;//true iff any of Xrange, Qrange, Xnpoints, Qnpoints is provided
  //If grid_provided==true and option "preferInternalGrid" is not given, internal grid or evolution will not be used.
  //This is useful, for example, when one wants to reduce TMD grid density, or use a smaller Q range


  auto const& pdfMap1 = pdf->xfxQmap(0.1, 10.);
  cout << " size xfxQmap "<< pdfMap1.size() << endl;


  for (const auto it : node) { //Iterate over options as key:value pairs
    string key;
    try {
      key = it.first.as<string>();
    } catch (YAML::TypedBadConversion<string>) {
      cerr << "[ERROR] WriteTMD failed to convert key to string when trying to output PDF \"" << pdf->_name <<
           "\"; trying to print key:" << endl;
      cerr << it.first << endl;
      hf_errlog(19060700, "F: Error when parsing WriteTMD, see stderr");
      abort();
    }

    try { //catch YAML::TypedBadConversion<string> exceptions
      //From now on we assume that it.second can be converted to string
      if (key == "name") {
        info.name = it.second.as<string>();
      } else if (key == "SetIndex") {
        info.setindex = it.second.as<string>();
      } else if (key == "evolution") {
        info.evolution = it.second.as<string>();
      } else if (key == "description") {
        info.description = it.second.as<string>();
      } else if (key == "authors") {
        info.authors = it.second.as<string>();
      } else if (key == "reference") {
        info.reference = it.second.as<string>();
      } else if (key == "Format") {
        info.format = it.second.as<string>();
      } else if (key == "TMDScheme") {
        info.tmd_scheme = it.second.as<string>();
      } else if (key == "FlavorScheme") {
        info.flavor_scheme = it.second.as<string>();
      } else if (key == "Extrapolation_x") {
        info.extrapolation_x = it.second.as<string>();
      } else if (key == "Extrapolation_Q2") {
        info.extrapolation_q2 = it.second.as<string>();
      } else if (key == "Extrapolation_kt") {
        info.extrapolation_kt = it.second.as<string>();        
      } else if (key == "NumFlavors") {
        info.numflavors= it.second.as<int>();        
      } else if (key == "Particle") {
        info.particle = it.second.as<int>();        
      } else if (key == "XMin") {
        info.xmin = it.second.as<double>();
      } else if (key == "XMax") {
        info.xmax = it.second.as<double>();
      } else if (key == "QMin") {
        info.qmin = it.second.as<double>();
      } else if (key == "QMax") {
        info.qmax = it.second.as<double>();
      } else if (key == "KtMin") {
        info.ktmin = it.second.as<double>();
      } else if (key == "KtMax") {
        info.ktmax = it.second.as<double>();
      } else if (key == "Xrange" or key == "Ktrange" or key == "Qrange") {
        YAML::Node n = it.second;
        if (!n.IsSequence() or n.size() != 2) {
          cerr << "[ERROR] WriteTMD: " << key << " must be given as [min, max]; error when trying to output PDF \"" << pdf->_name << "\""
               << endl;
          hf_errlog(19060700, "F: Error when parsing WriteTMD, see stderr");
          abort();
        }
        double min, max;
        try {
          min = n[0].as<double>();
          max = n[1].as<double>();
        } catch (YAML::Exception) {
          cerr << "[ERROR] WriteTMD: Failed to interpret " << key << " for PDF \"" << pdf->_name << "\"" << endl;
          hf_errlog(19060700, "F: Error when parsing WriteTMD, see stderr");
          abort();
        }
        if (not(min < max)) {
          cerr << "[ERROR] WriteTMD: In " << key << "=[min, max] min must be smaller than max; got " << key << "=[" << min << ", " << max
               << "] ; for PDF \"" << pdf->_name << "\"" << endl;
          hf_errlog(19060700, "F: Error when parsing WriteTMD, see stderr");
          abort();
        }
        if (key[0] == 'X') { //key=="Xrange"
          info.xmin = min;
          info.xmax = max; 
        } else if (key[0] == 'K') { //key=="Ktrange"
          info.ktmin = min;
          info.ktmax = max;
        } else { //key=="Qrange"
          info.qmin = min;
          info.qmax = max;
        }
        grid_provided = true;
      } else {
        cerr << "[WARN] WriteTMD: Ignoring unknown option \"" << key << ": " << it.second << "\"" << endl;
        hf_errlog(19063000, "W: WriteTMD: ignoring unknown option, see stderr");
      }
    } catch (YAML::TypedBadConversion<string>()) {
      cerr << "[ERROR] WriteTMD failed to convert value of key \"" << key << "\" to string when trying to output PDF \"" << pdf->_name
           << "\"; trying to print value" << endl;
      cerr << it.second << endl;
      hf_errlog(19060700, "F: Error when parsing WriteTMD, see stderr");
      abort();
    }
  }

  //If some options were omitted, use defaults
  if (std::isnan(info.xmin) or std::isnan(info.xmax)) {
    info.xmin = 1e-6;
    info.xmax = 1;
    if (grid_provided) {
      cerr << "[INFO] WriteTMD: using default Xrange=[" << info.xmin << ", " << info.xmax << "] for PDF \"" << pdf->_name << "\"" <<
           endl;
    }
  }
  if (std::isnan(info.qmin) or std::isnan(info.qmax)) {
    info.qmin = 1;
    info.qmax = 1e4;
    if (grid_provided) {
      cerr << "[INFO] WriteTMD: using default Qrange=[" << info.qmin << ", " << info.qmax << "] for PDF \"" << pdf->_name << "\"" <<
           endl;
    }
  }
  if (std::isnan(info.ktmin) or std::isnan(info.ktmax)) {
    info.ktmin = 0.01;
    info.ktmax = 1e4;
    if (grid_provided) {
      cerr << "[INFO] WriteTMD: using default Ktrange=[" << info.ktmin << ", " << info.ktmax << "] for PDF \"" << pdf->_name << "\"" <<
           endl;
    }
  }
  if (info.Nx == 0) {
    info.Nx = 200;
    if (grid_provided) {
      cerr << "[INFO] WriteTMD: using default Xnpoints=" << info.Nx << " for PDF \"" << pdf->_name << "\"" << endl;
    }
  }
  if (info.Nq == 0) {
    info.Nq = 120;
    if (grid_provided) {
      cerr << "[INFO] WriteTMD: using default Xnpoints=" << info.Nx << " for PDF \"" << pdf->_name << "\"" << endl;
    }
  }

  if (not grid_provided) info.prefer_internal_grid = true;

  return info;
}

void WriteTMDinfo(FILE* f, const TMD_Options& info)
{
  if (!info.description.empty()) fprintf(f, "SetDesc: \"%s\"\n"  , info.description.c_str());
  if (!info.setindex.empty()) fprintf(f, "SetIndex: \"%s\"\n"  , info.setindex.c_str());
  if (!info.authors.empty())     fprintf(f, "Authors: \"%s\"\n"  , info.authors.c_str());
  if (!info.reference.empty())   fprintf(f, "Reference: \"%s\"\n", info.reference.c_str());
  fprintf(f, "Particle: %i\n", info.particle);
  if (!info.format.empty())   fprintf(f, "Format: \"%s\"\n", info.format.c_str());
  fprintf(f, "NumMembers: %i\n", info.nmembers);
  fprintf(f, "NumFlavors: %i\n", info.numflavors);
  if (!info.flavor_scheme.empty())   fprintf(f, "FlavorScheme: \"%s\"\n", info.flavor_scheme.c_str());
  if (!info.tmd_scheme.empty())   fprintf(f, "TMDScheme: \"%s\"\n", info.tmd_scheme.c_str());
  fprintf(f, "Flavors: [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 22, 23, 24, -24, 25]\n");  

  //if (!info.evolution.empty())   fprintf(f, "evolution: \"%s\"\n", info.evolution.c_str());
  if (!info.extrapolation_x.empty())   fprintf(f, "Extrapolation_x: \"%s\"\n", info.extrapolation_x.c_str());
  if (!info.extrapolation_q2.empty())   fprintf(f, "Extrapolation_Q2: \"%s\"\n", info.extrapolation_q2.c_str());
  if (!info.extrapolation_kt.empty())   fprintf(f, "Extrapolation_kt: \"%s\"\n", info.extrapolation_kt.c_str());
  
  //fprintf(f, "DataVersion: 1\n");

  int order = OrderMap(XFITTER_PARS::getParamS("Order"));
  order -= 1; // qcdnum convention LO=1, NLO=2, ...; TMD convention LO=0, NLO=1, ...
  fprintf(f, "OrderQCD: %i\n", order);
  if (!info.error_type.empty())    fprintf(f, "ErrorType: %s\n"   , info.error_type.c_str());


  fprintf(f, "XMin: %g\n", info.xmin);
  fprintf(f, "XMax: %g\n", info.xmax);
  fprintf(f, "KtMin: %g\n", info.ktmin);
  fprintf(f, "KtMax: %g\n", info.ktmax);
  fprintf(f, "QMin: %g\n", info.qmin);
  fprintf(f, "QMax: %g\n", info.qmax);

  //TODO: different quark and MZ masses for different evolutions
  //TODO: get masses from evolutions
  using XFITTER_PARS::getParamD;
  double Mz = *getParamD("Mz");
  fprintf(f, "MZ: %g\n", Mz);
  fprintf(f, "MUp: 0\n");
  fprintf(f, "MDown: 0\n");
  fprintf(f, "MStrange: 0\n");
  fprintf(f, "MCharm: %g\n", *getParamD("mch")); 
  fprintf(f, "MBottom: %g\n", *getParamD("mbt"));
  fprintf(f, "MTop: %g\n", *getParamD("mtp"));

  BaseEvolution* pdf = info.pdf;

  auto const& pdfMap = pdf->xfxQmap(0.1, 10.);
  cout << " tmdoiutput: WriteTMDinf: size xfxQmap "<< pdfMap.size() << endl;

  //Write alphaS
  fprintf(f, "AlphaS_MZ: %g\n", pdf->getAlphaS(Mz));
  fprintf(f, "AlphaS_OrderQCD: %i\n", order);
  //fprintf(f, "AlphaS_Type: ipol\n");//I have no idea what that is --Ivan

  cout << " before calling pdfMap " << endl;
  auto const& pdfMap1 = info.pdf->xfxQmap(0.5,10.);

  cout << " tmd_output: size xfxQmap "<< pdfMap1.size() << endl;

}

bool TMDdirectoryExists(const string& path)
{
  struct stat info;
  if (stat(path.c_str(), &info) != 0) return false;
  return bool(info.st_mode & S_IFDIR);
}
size_t getTMDNmembers()
{
  //TODO this needs to be more general
  YAML::Node minuitNode = XFITTER_PARS::rootNode["MINUIT"];
  if (!minuitNode.IsMap()) return 1;
  YAML::Node doErrorsNode = minuitNode["doErrors"];
  string doErrors;
  try {
    doErrors = doErrorsNode.as<string>();
  } catch (YAML::TypedBadConversion<string>) {
    return 1;
  }
  size_t Npars = xfitter::get_minimizer()->getNpars();
  if (doErrors == "Hesse")return Npars + 1;
  if (doErrors == "Pumplin")return 2 * Npars + 1;
  return 1;
}

//FORTRAN INTERFACE
extern "C" {

  void save_data_tmd_(const int& memberID) //This is called when building bands
  {
    YAML::Node TMDnode = XFITTER_PARS::rootNode["WriteTMD"];
    if (not TMDnode) return;
    //TODO: output multiple evolutions
    //TODO: handle errors here
    TMD_Options options = TMD_Options::fromYAML(TMDnode);
   // WriteTMD(f, options.pdf, PdfType);


    options.nmembers = getTMDNmembers() - 1;
    const string& name = options.name;

    string outdir = xfitter::getOutDirName() + '/' + name;
    if (not TMDdirectoryExists(outdir)) {
      if (mkdir(outdir.c_str(), 0755) != 0) {
        cerr << "[ERROR] Failed to create directory \"" << outdir << "\" for TMD output" << endl;
        perror("mkdir error");
        hf_errlog(19060702, "F: Failed to create directory for TMD output, see stderr");
        abort();
      }
    }
    //Make file name
    string filename;
    filename.reserve(outdir.size() + name.size() + 10);
    sprintf(const_cast<char*>(filename.c_str()), "%s/%s_%04i.dat", outdir.c_str(), name.c_str(), memberID);

    FILE* f = fopen(filename.c_str(), "w");

    const char* PdfType = "central";
    const char* file = filename.c_str();
    if (memberID != 0) PdfType = "error";
    std::cout << " in save_data_TMD " << filename.c_str() << std::endl;

    if (f == nullptr) {
      cerr << "[ERROR] Failed to open file \"" << filename << "\" for TMD output" << endl;
      perror("fopen error");
      hf_errlog(19060800, "F: Failed to open file for TMD output, see stderr");
      abort();
    }
       
    WriteTMD(file, options.pdf, PdfType);
    
    if (memberID == 0) { //then write info
      filename = outdir + '/' + name + ".info";
      f = fopen(filename.c_str(), "w");
      if (f == nullptr) {
        cerr << "[ERROR] Failed to open file \"" << filename << "\" for TMD output" << endl;
        perror("fopen error");
        hf_errlog(19060800, "F: Failed to open file for TMD output, see stderr");
        abort();
      }
      WriteTMDinfo(f, options);
      fclose(f);
    }

    
    //WriteTMD(f, options.pdf, PdfType);
    //fclose(f);

   
  }

  void print_tmd_()
  {
    save_data_tmd_(0);
  }

//TODO: replace all calls to this function with calls to print_tmd and delete this one
  void print_tmd_opt_()
  {
    save_data_tmd_(0);
  }

}
