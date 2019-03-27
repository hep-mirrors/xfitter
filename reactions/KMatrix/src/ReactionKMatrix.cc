
/*
   @file ReactionKMatrix.cc
   @date 2018-08-03
   @author  AddReaction.py
   Created by  AddReaction.py on 2018-08-03
*/

#include "ReactionKMatrix.h"
#include <iostream>
#include <fstream>
#include <sstream>

// the class factories
extern "C" ReactionKMatrix* create()
{
  return new ReactionKMatrix();
}


// Initialize at the start of the computation
int ReactionKMatrix::initAtStart(const string &s)
{
  return 0;
}

// Initialisze for a given dataset:
void ReactionKMatrix::setDatasetParameters(int dataSetID, map<string,string> pars, map<string, double> parsDataset){
  std::vector<double> temp;
  double var;
  // check if KMatrixs should be read from separate file
  if (pars.find("FileName") != pars.end())
  {
    // file name to read KMatrixs
    std::string fileName = pars["FileName"];

    // TODO: how to get number of data points?
    // not very elegant way below
    int np = _dsBins[dataSetID]->begin()->second.size();

    // requested starting column from file (by default 1st)
    int column_start = 1;
    if(pars.find("FileColumnStart") != pars.end()){
      column_start = atoi(pars["FileColumnStart"].c_str());
    }
    // requested finishing column from file (by default 1st)
    int column_finish = np;
    if(pars.find("FileColumnFinish") != pars.end()){
      column_finish = atoi(pars["FileColumnFinish"].c_str());
    }

    // check that the column is reasonable
    if(column_start < 1){
      hf_errlog(18080700, "F: wrong starting column = " + std::to_string(column_start));
    }
    if(column_start > column_finish){
      hf_errlog(18080701, "F:  starting column greater than finishing column ");
    }

    // requested starting line from file (by default 1st) and last line (by default last line)
    int lineStart = 1;
    if(pars.find("FileLineStart") != pars.end())
      lineStart = atoi(pars["FileLineStart"].c_str());
    int lineFinish = -1;
    if(pars.find("FileLineFinish") != pars.end())
      lineFinish = atoi(pars["FileLineFinish"].c_str());

    // open file
    std::ifstream file(fileName.c_str());
    string line;
    if (!file.is_open()){
      hf_errlog(18080702, "F: error opening KMatrix file = " + fileName);
    }

    // skip lineStart lines
    int readLines = 0;
    for(int l = 1; l < lineStart; l++){
      readLines++;
      getline(file, line);
    }

    while(getline(file,line))
    {
      readLines++;
      if(lineFinish != -1 && readLines > lineFinish)
        break;
      if(line.at(0) == '#') continue; //ignore comments
      line.erase(line.find_last_not_of(" \n\r\t")+1); // trim trailing whitespaces

      std::stringstream sline(line);

      int current_col = 1;
      while(sline.good()){
        sline >> var;
        if(current_col >= column_start && current_col <= column_finish){ //Only using range between spezified columns
          temp.push_back(var);
        }
        current_col++;
      }
      _values2D[dataSetID].push_back(temp);
      temp.clear();
      sline.clear();

    }
    file.close();

    //mapping 2d matrix (m x n) to 1d vector (m*n): list of column vectors, mapping with vec(i*n + j) = mat(j,i)
    int m = _values2D[dataSetID].size();
    int n = _values2D[dataSetID].at(0).size();
    _values[dataSetID].resize(n*m);
    for(int i = 0; i < m; i++){
      for(int j = 0; j < n; j++){
        _values[dataSetID].at(i*n+j)=_values2D[dataSetID].at(i).at(j);
      }
    }
  }else{
    hf_errlog(18080703, "F: FileName must be provided for KMatrix");
  }
}

// Main function to compute results at an iteration
int ReactionKMatrix::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
  // kmatrix is constant value read in setDatasetParameters()
  val = std::valarray<double>(_values[dataSetID].data(), _values[dataSetID].size());

  return 0;
}

