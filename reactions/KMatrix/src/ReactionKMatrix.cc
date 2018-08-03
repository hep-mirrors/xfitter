 
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

      // check that the column is reasonable//TODO ERROR-log numbering
      if(column_start < 1){
        hf_errlog(17102800, "F: wrong starting column = " + std::to_string(column));
      }
      if(column_start > column_finish){
        hf_errlog(17102800, "F:  starting column greater than finishing column ");
      }

      // requested starting line from file (by default 1st)
      int lineStart = 1;
      if(pars.find("FileLine") != pars.end()){
        lineStart = atoi(pars["FileLine"].c_str());
      }

      // open file
      std::ifstream file(fileName.c_str());
      string line;
      if (!file.is_open()){
        hf_errlog(17102802, "F: error opening KMatrix file = " + fileName);
      }

      // skip lineStart lines
      for(int l = 1; l < lineStart; l++){
        getline(file, line);
      }

      while(getline(file,line))
      {
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
          _values[dataSetID].push_back(temp);
          temp.clear();
          sline.clear();

          }
      file.close();
    }


    /*// check if KMatrixs should be read from data file (soecifying column is mandatory)
    else if (pars.find("DataColumn") != pars.end())
    {
      std::string columnName = pars["DataColumn"];
      valarray<double>* vals = this->GetBinValues(dataSetID, columnName);
      if(!vals)
        hf_errlog(17102803, "F: no column = " + columnName + " in data file");
      // store KMatrix values
      _values[dataSetID].resize(vals->size());
      _values[dataSetID].assign(std::begin(*vals), std::end(*vals));
    }
    else
      hf_errlog(17102804, "F: FileName or DataColumn or Parameter must be provided for KMatrix");
  }*/
}

// Main function to compute results at an iteration
int ReactionKMatrix::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
  return 0;
}

