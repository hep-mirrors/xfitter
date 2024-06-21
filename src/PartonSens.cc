#include "PartonSens.h"
#include "xfitter_cpp.h"
#include "xfitter_pars.h"
#include <valarray>

extern "C" {
  void update_theory_iteration_();
}

namespace xfitter
{
  std::vector<int> PartonSens::_global_partons = std::vector<int>();
  std::vector<int> PartonSens::_global_partons1 = std::vector<int>();
  double PartonSens::_global_x1 = -1.0;
  double PartonSens::_global_x2 = -1.0;

  void PartonSens::doPartonSens() {
    using namespace std;
    const YAML::Node node=XFITTER_PARS::rootNode["PartonSens"];
    if(!node)return;
    if(!node.IsMap()) hf_errlog(2024062001,"F: Cannot do parton sensitivity: PartonSens node is not a YAML map");
    if (node["Status"]) {
      if (node["Status"].as<string>() == "Off") return;
      if (node["Status"].as<string>() != "On") {
        const auto msg = "F: Unknown PartonSens Status = " + node["Status"].as<string>();
        hf_errlog(20240620, msg);
      }
    }
    YAML::Node const xrange = node["xrange"];
    if(!xrange.IsSequence() || xrange.size() != 3){
      hf_errlog(2024062002,"F: PartonSens: xrange must be sequence of xmin, xmax, nx");
    }
    const double xmin = xrange[0].as<double>();
    const double xmax = xrange[1].as<double>();
    const int nx = xrange[2].as<int>();
    bool xlog = true;
    YAML::Node const node_xlog = node["xlog"];
    if(node_xlog) {
      xlog = node_xlog.as<bool>();
    }
    std::valarray<double> xbins(nx);
    if(xlog) {
      double step = (log(xmax) - log(xmin)) / (nx - 1);
      for(int ix = 0; ix < nx; ix++) 
        xbins[ix] = exp(log(xmin) + step * ix);
    }
    else {
      double step = (xmax - xmin) / (nx - 1);
      for(int ix = 0; ix < nx; ix++) 
        xbins[ix] = xmin + step * ix;
    }
    int input_parton = node["input_parton"] ? node["input_parton"].as<int>() : 1;
    if(input_parton != 1 && input_parton != 2) {
      hf_errlog(2024062003,"F: PartonSens: input_parton must be 1 or 2");
    }
    std::vector<int>& _global_partons_ref = (input_parton==1) ? _global_partons : _global_partons1;
    YAML::Node const partons = node["partons"];
    if(!partons.IsSequence()){
      hf_errlog(2024062004,"F: PartonSens: partons must be sequence");
    }
    std::valarray<std::valarray<std::valarray<double> > > pred(partons.size());
    for(size_t i = 0; i < partons.size(); i++) {
      pred[i].resize(xbins.size());
      _global_partons_ref.resize(partons[i].size());
      for(size_t ip = 0; ip < partons[i].size(); ip++) {
        _global_partons_ref[ip] = partons[i][ip].as<int>();
      }
      for(size_t ix = 0; ix < xbins.size() - 1; ix++) {
        printf("xbins[%ld] = %f, %f\n", ix, xbins[ix], xbins[ix+1]);
        _global_x1 = xbins[ix];
        _global_x2 = xbins[ix+1];
        update_theory_iteration_();
        //Return theory predictions
        pred[i][ix] = valarray<double>(c_theo_.theo, cndatapoints_.npoints);
        //printf("pred.size() = %ld\n", pred.size());
        //for(size_t ip = 0; ip < pred.size(); ip++) {
        //  printf("pred[%ld] = %f\n", ip, pred[ip]);
        //}
      }
    }
    _global_partons_ref.clear();
    // cross check
    if (1 == 1) {
      update_theory_iteration_();
      //Return theory predictions
      const auto pred_sum = valarray<double>(c_theo_.theo, cndatapoints_.npoints);
      valarray<double> pred_sum_check(pred_sum.size());
      for(size_t i = 0; i < partons.size(); i++) {
        for(size_t ix = 0; ix < xbins.size() - 1; ix++) {
          for(size_t id = 0; id < pred[i][ix].size(); id++) {
            pred_sum_check[id] += pred[i][ix][id];
          }
        }
      }
      double maxdiff = 0.;
      for(size_t id = 0; id < pred_sum.size(); id++) {
        double diff = std::abs(pred_sum_check[id] - pred_sum[id])/pred_sum[id];
        printf("diff[%ld] = %f\n", id, diff);
        maxdiff = std::max(maxdiff, diff);
      }
      printf("maxdiff = %f\n", maxdiff);
    }
    string outputdir = node["OutputDirectory"] ? node["OutputDirectory"].as<string>() : "output";
    for(size_t i = 0; i < partons.size(); i++) {
      FILE* fout = fopen((outputdir+"/partonsens_"+to_string(i)+".dat").c_str(), "w");
      //fprintf(fout, "%13s%13s%14s%14s\n", "xmin", "xmax", "data1", "...");
      fprintf(fout, "%13s%13s", "xmin", "xmax");
      for(size_t id = 0; id < pred[i][0].size(); id++) {
        fprintf(fout, "%14s", ("data"+to_string(id+1)).c_str());
      }
      fprintf(fout, "\n");
      for(size_t ix = 0; ix < xbins.size() - 1; ix++) {
        fprintf(fout, "%13.5e%13.5e", xbins[ix], xbins[ix+1]);
        for(size_t id = 0; id < pred[i][ix].size(); id++) {
          fprintf(fout, "%14.5e", pred[i][ix][id]);
        }
        fprintf(fout, "\n");
      }
      fclose(fout);
    }
    //chi2data_theory_(1);
    //writefittedpoints_();
    //chi2data_theory_(2);
    //bool cp = system(((string)"cp " + _outputDir + "/fittedresults.txt "
		// + _outputDir + "/fittedresults.txt_partonsens" + tag).c_str());
  }
} //namespace xfitter
