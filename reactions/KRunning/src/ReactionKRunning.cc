 
/*
   @file ReactionKRunning.cc
   @date 2019-01-16
   @author  AddReaction.py
   Created by  AddReaction.py on 2019-01-16
*/

#include "ReactionKRunning.h"

// the class factories
extern "C" ReactionKRunning* create() {
  return new ReactionKRunning();
}


// Initialize at the start of the computation
int ReactionKRunning::atStart(const string &s)
{
  return 0;
}

void ReactionKRunning::setDatasetParameters(int dataSetID, map<std::string, std::string> pars, map<std::string, double> dsPars)
{
  // check if dataset with provided ID already exists
  if(_type.find(dataSetID) != _type.end())
    hf_errlog(19011501, "F: dataset with id " + std::to_string(dataSetID) + " already exists");

  // read type of running
  // presently only alphaS (accesses whatever running implemented in xFitter) and mass MSbar running (at NLO) is implemented
  auto it = pars.find("type");
  if(it == pars.end())
    hf_errlog(19011502, "F: no type of running for dataset with id " + std::to_string(dataSetID));
  if(it->second != "as" && it->second != "massMSbarNLO")
    hf_errlog(19011503, "F: unsupported running type = " + it->second);
  _type[dataSetID] = it->second;
  
  // read scale
  it = pars.find("q");
  if(!checkParam(it->second)) // value provided
    _qValue[dataSetID] = stod(it->second);
  else // parameter name is provided
    _q[dataSetID] = it->second;

  // for type=massMSbarNLO read q0: scale at which m(m) is quoted
  if(_type[dataSetID] == "massMSbarNLO")
    _q0[dataSetID] = pars.find("q0")->second;

  // read optional number of points (if not provided, use the number of data points)
  it = pars.find("N");
  if(it == pars.end())
    _NPoints[dataSetID] = _dsBins[dataSetID]->begin()->second.size();
  else
    _NPoints[dataSetID] = atoi(pars["N"].c_str());
  //printf("npoints: %d\n", _NPoints[dataSetID]);
}

// Main function to compute results at an iteration
int ReactionKRunning::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
  double q = (_qValue.find(dataSetID) != _qValue.end()) ? _qValue[dataSetID] : GetParam(_q[dataSetID]);
  if(_type[dataSetID] == "as")
    val = valarray<double>(getAlphaS(q), _NPoints[dataSetID]);
  else if(_type[dataSetID] == "massMSbarNLO")
  {
    double q0 = GetParam(_q0[dataSetID]);
    val = valarray<double>(getMassMSbar(q0, q, getAlphaS(q0), getAlphaS(q)), _NPoints[dataSetID]);
  }
  //for(int i = 0; i < val.size(); i++)
  //  printf("val[%d] = %f\n", i, val[i]);
  return 0;
}
