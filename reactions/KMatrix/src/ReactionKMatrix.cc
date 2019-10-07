
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
#include"hf_errlog.h"

// the class factories
extern "C" ReactionKMatrix* create()
{
  return new ReactionKMatrix();
}
struct ReactionData{
  valarray<double>_values;
};
// Initialisze for a given dataset:
void ReactionKMatrix::initTerm(TermData*td){
  ReactionData*rd=new ReactionData;
  td->reactionData=(void*)rd;
  int np=td->getNbins();
  // requested starting column from file (by default 1st)
  int column_start = 1;
  if(td->hasParam("FileColumnStart")){
    column_start=td->getParamI("FileColumnStart");
    if(column_start<1){
      hf_errlog(18080700, "F: wrong starting column = "+std::to_string(column_start));
    }
  }
  // requested finishing column from file (by default 1st)
  int column_finish = np;
  if(td->hasParam("FileColumnFinish")){
    column_finish=td->getParamI("FileColumnFinish");
  }
  // check that the column is reasonable
  if(column_start > column_finish){
    hf_errlog(18080701, "F:  starting column greater than finishing column ");
  }
  // requested starting line from file (by default 1st)
  int lineStart = 1;
  if(td->hasParam("FileLineStart")){
    lineStart=td->getParamI("FileLineStart");
    if(lineStart<0)hf_errlog(19050720,"F: Negative FileLineStart");
  }
  int lineFinish=-1;
  if(td->hasParam("FileLineFinish")){
    lineFinish=td->getParamI("FileLineFinish");
    if(lineFinish<lineStart)hf_errlog(19050721,"F: FileLineFinish<FileLineStart");
  }

  // file name to read KMatrix
  string fileName=td->getParamS("FileName");
  // open file
  std::ifstream file(fileName.c_str());
  if (!file.is_open()){
    hf_errlog(18080702, "F: error opening KMatrix file = " + fileName);
  }

  string line;
  // skip lineStart lines
  int readLines=0;//number of next line
  for(int l = 1; l < lineStart; l++){
    readLines++;
    getline(file, line);
  }

  double var;
  vector<vector<double>>_values2D;
  while(getline(file,line)){
    readLines++;
    if(lineFinish!=-1&&readLines>lineFinish)break;
    if(line.at(0) == '#') continue; //ignore comments
    line.erase(line.find_last_not_of(" \n\r\t")+1); // trim trailing whitespaces
    _values2D.push_back(vector<double>());
    vector<double>&col=_values2D.back();

    std::stringstream sline(line);

    int current_col = 1;
    while(sline.good()){
      sline >> var;
      if(current_col >= column_start && current_col <= column_finish){ //Only using range between spezified columns
        col.push_back(var);
      }
      current_col++;
    }
    sline.clear();
  }
  file.close();

  //mapping 2d matrix (m x n) to 1d vector (m*n): list of column vectors, mapping with vec(i*n + j) = mat(j,i)
  int m = _values2D.size();
  int n = _values2D.at(0).size();
  valarray<double>&_values=rd->_values;
  _values.resize(n*m);
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      _values[i*n+j]=_values2D[i][j];
    }
  }
}
void ReactionKMatrix::freeTerm(TermData*td){
  delete(ReactionData*)td->reactionData;
}
// Main function to compute results at an iteration
void ReactionKMatrix::compute(TermData*td,valarray<double>&val,map<string,valarray<double> >&err){
  val=((ReactionData*)td->reactionData)->_values;
}

