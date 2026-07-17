
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
#include"hf_errlog.h"
#include "xfitter_cpp_base.h"
#ifdef CRUNDEC_FOUND
#include "CRunDec.h"
#else
struct CRunDec {
  ~CRunDec(){};
};
#endif
using namespace std;

// the class factories
extern "C" ReactionKFactor* create()
{
  return new ReactionKFactor();
}
struct ReactionData{
  ~ReactionData() {
    if(crd) {
      delete crd;
    }
  }
  valarray<double>_values;
  const double*parameter=nullptr;
  /* If parameter "Parameter" is given, KFactor returns column with all values equal to that parameter
   * This parameter can be fitted
   * If FileName is given, values are loaded from that file, and KFactor returns this constant column
   * If DataColumn is given, the constant column is read from the datafile
  */
  size_t Npoints=0;//if parameter!=nullptr, this is the size of returned array. Usually it is the number of datapoints in the dataset, but can be overridden using option "N"
  ReactionKFactor::Convert convert;
  CRunDec* crd;
  int nf;
  int nl;
};
// Initialisze for a given dataset:
void ReactionKFactor::initTerm(TermData*td){
  ReactionData*rd=new ReactionData;
  td->reactionData=(void*)rd;
  rd->convert=ReactionKFactor::Convert::None;
  rd->crd = nullptr;
  // check if kfactors should be read from separate file
  if(td->hasParam("FileName")){
    // file name to read kfactors
    string fileName=td->getParamS("FileName");
    // requested column from file (by default 1st)
    int column=1;
    if(td->hasParam("FileColumn")){
      column=td->getParamI("FileColumn");
      if(column<1)hf_errlog(17102800, "F: wrong column = " + std::to_string(column));
    }
    // requested starting line from file (by default 1st)
    int lineStart=1;
    if(td->hasParam("FileLine"))lineStart=td->getParamI("FileLine");
    // requested last line -> how many values to read (by default the same as number of data points)
    int np=td->getNbins();
    int lineFinish=lineStart+np-1;
    if(td->hasParam("FileLineFinish")){
      lineFinish=td->getParamI("FileLineFinish");
      np = lineFinish - lineStart + 1;
    }
    // open file
    std::ifstream file(fileName.c_str());
    string line;
    if (!file.is_open())
      hf_errlog(17102802, "F: error opening kfactor file = " + fileName);

    // skip lineStart lines
    for(int l = 1; l < lineStart; l++)
      getline(file, line);

    vector<double>_values;
    for(int p = 0; p < np; p++)
    //while (1)
    {
      // read new line
      getline(file, line);
      if(file.eof())break;
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
      _values.push_back(val);
    }
    file.close();
    rd->_values=valarray<double>(_values.data(),_values.size());
  }
  // check if kfactors should be read from data file (specifying column is mandatory)
  else if(td->hasParam("DataColumn"))
  {
    rd->_values=td->getBinColumn(td->getParamS("DataColumn"));
  }
  // check if kfactor is a parameter (possibly free); the value of this parameter will be used for all bins in data set
  else if(td->hasParam("Parameter"))
  {
    std::string parameterName = td->getParamS("Parameter");
    rd->parameter=td->getParamD(parameterName);
    if(td->hasParam("N")){
      int Npoints=td->getParamI("N");
      if(Npoints<=0){
        cerr<<"[ERROR] Requested nonpositive array size in reaction KFactor, N="<<Npoints<<"; termID="<<td->id<<endl;
        hf_errlog(19050710, "F: Requested nonpositive array size in reaction KFactor");
      }
      rd->Npoints=size_t(Npoints);
    }else{
      rd->Npoints=td->getNbins();
    }
    if(td->hasParam("MSBAR_to_pole") && td->getParamS("MSBAR_to_pole") == parameterName)
    {
      rd->convert = ReactionKFactor::Convert::MSBAR_to_Pole;
      rd->nf = td->getParamI("NFlavour");
      rd->nl = OrderMap(td->getParamS("Order"));
#ifdef CRUNDEC_FOUND
      rd->crd = new CRunDec(rd->nf);
#else
      hf_errlog(2026071701,"F:CRunDec not found, please install");
#endif
    }
  }
  else
    hf_errlog(17102804, "F: FileName or DataColumn or Parameter must be provided for KFactor");
}
void ReactionKFactor::freeTerm(TermData*td){delete(ReactionData*)td->reactionData;}
// Main function to compute results at an iteration
void ReactionKFactor::compute(TermData*td, valarray<double> &val, map<string, valarray<double> > &err)
{
  ReactionData*rd=(ReactionData*)td->reactionData;
  if(rd->parameter){
    if(rd->convert == ReactionKFactor::Convert::None) {
      double parval=*rd->parameter;
      val=valarray<double>(parval,rd->Npoints);
    }
    if(rd->convert == ReactionKFactor::Convert::MSBAR_to_Pole) {
#ifdef CRUNDEC_FOUND
      double m_msbar=*rd->parameter;
      td->actualizeWrappers();
      double as = alphas_wrapper_(m_msbar);
      double m_pole = rd->crd->mMS2mOS(m_msbar, nullptr, as, m_msbar, rd->nf, rd->nl);
      //printf("m_msbar -> m_pole: %f -> %f\n", m_msbar, m_pole);
      val=valarray<double>(m_pole,rd->Npoints);
#else
      hf_errlog(2026071701,"F:CRunDec not found, please install");
#endif
    }
  }else{ // kfactor is constant value read in setDatasetParameters()
    val=rd->_values;
  }
}

