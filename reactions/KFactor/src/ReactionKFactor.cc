
/*
   @file ReactionKFactor.cc
   @date 2017-10-28
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-10-28
*/

#include "ReactionKFactor.h"
#include <iostream>
#include <fstream>
#include <sstream>

// the class factories
extern "C" ReactionKFactor* create()
{
  return new ReactionKFactor();
}


// Initialize at the start of the computation
int ReactionKFactor::initAtStart(const string &s)
{
  return 0;
}

// Initialisze for a given dataset:
void ReactionKFactor::setDatasetParamters(int dataSetID, map<string,string> pars, map<string, double> parsDataset)
{
  // check if kfactors should be read from separate file
  if (pars.find("FileName") != pars.end())
  {
    // file name to read kfactors
    std::string fileName = pars["FileName"];

    // requested column from file (by default 1st)
    int column = 1;
    if(pars.find("FileColumn") != pars.end())
      column = atoi(pars["FileColumn"].c_str());

    // check that the column is reasonable
    if(column < 1)
      hf_errlog(17102800, "F: wrong column = " + std::to_string(column));

    // open file
    std::ifstream file(fileName.c_str());
    string line;
    if (!file.is_open())
      hf_errlog(17102802, "F: error opening kfactor file = " + fileName);

    while (1)
    {
      // read new line
      getline(file, line);
      if (true == file.eof()) break;
      if (line.at(0) == '#' ) continue; //ignore comments
      line.erase(line.find_last_not_of(" \n\r\t")+1); // trim trailing whitespaces
      std::stringstream sline(line);

      // count columns in line
      int nColumns = 0;
      while (sline.good())
      {
        string dummy;
        sline >> dummy;
        nColumns++;
      }

      // check that the number of columns is not smaller than the requested column
      if(column > nColumns)
        hf_errlog(17102801, "F: no column = " + std::to_string(column) + " in file " + fileName);

      // read kfactor value
      sline.clear();
      sline.seekg(0);
      sline.str(line);
      double val = 0.0;
      for (int col = 0; col < column; col++)
        sline >> val;

      // store kfactor value
      _values[dataSetID].push_back(val);
    }
    file.close();
  }
  // check if kfactors should be read from data file (soecifying column is mandatory)
  else if (pars.find("DataColumn") != pars.end())
  {
    std::string columnName = pars["DataColumn"];
    valarray<double>* vals = this->GetBinValues(dataSetID, columnName);
    if(!vals)
      hf_errlog(17102803, "F: no column = " + columnName + " in data file");
    // store kfactor values
    _values[dataSetID].resize(vals->size());
    _values[dataSetID].assign(std::begin(*vals), std::end(*vals));
  }
  else
    hf_errlog(17102804, "F: either FileName or DataColumn must be provided for KFactor");
}
// Main function to compute results at an iteration
int ReactionKFactor::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
  val = std::valarray<double>(_values[dataSetID].data(), _values[dataSetID].size());
  return 0;
}

