/*!
 @file TheorEval.cc
 @date Tue Aug 12 2013
 @author Andrey Sapronov <Andrey.Sapronov@cern.ch>

 Contains TheorEval class member function implementations.
 */

#include <fstream>
#include <list>
#include <sstream>
#include <stack>
#include <float.h>
#include <valarray>
#include <dlfcn.h>

#include "TheorEval.h"
#include "ReactionTheory.h"
#include"TermData.h"
#include "xfitter_cpp.h"
#include <string.h>

#include <yaml-cpp/yaml.h>
#include "xfitter_pars.h"
#include "xfitter_steer.h"

// ROOT spline can be uncommented (here and below in the code) and used e.g. for cross checks (obviously requires ROOT)
//#include <TSpline.h>
#include <spline.h>

using namespace std;
TheorEval::TheorEval(const int dsId, const int nTerms, const std::vector<string> stn, const std::vector<string> stt,
                     const std::vector<string> sti, const std::vector<string> sts, const string& expr) : _dsId(dsId), _nTerms(nTerms)
{
  for (int it= 0 ; it<nTerms; it++ ){
    _termNames.push_back(stn[it]);
    _termTypes.push_back(stt[it]);
    _termInfos.push_back(sti[it]);
    _termSources.push_back(sts[it]);
  }
  _expr.assign(expr);
}

TheorEval::~TheorEval()
{
  for(auto it:_exprRPN){
    if(it.val and it.ownsVal){
      delete it.val;
      it.val=nullptr;
    }
  }

  // OZ delete reactions
  for (tNameReactionmap::iterator itt = gNameReaction.begin(); itt!=gNameReaction.end(); itt++)
  {
    if(itt->second)
    {
      delete itt->second;
      itt->second = NULL;
    }
  }
  //Delete all instances of TermData
  for(const auto it:term_datas)delete it;
  term_datas.clear();
}

void
TheorEval::initTheory()
{
  list<tToken> sl;
  this->assignTokens(sl);
  this->convertToRPN(sl);
}

void
TheorEval::initReactionToken(tToken& t,const string& name){
  const vector<string>::iterator found_term = find(_termNames.begin(), _termNames.end(), name);
  if ( found_term == _termNames.end() ) {
    cerr<<"[ERROR] Undeclared reaction term \""<<name<<"\" in expression \""<<_expr<<'\"'<<endl;
    hf_errlog(19051430,"F: Undeclared reaction term, see stderr");
  }
  int iterm = int(found_term-_termNames.begin());
  t.opr  = 0;
  t.name = name;
  const auto it=_mapInitdTerms.find(name);
  if(it!=_mapInitdTerms.end()){
    t.val = it->second;
    t.ownsVal = false;
  }else{
    t.val = new valarray<double>(0.,getNbins());
    t.ownsVal = true;
    initTerm(iterm,t.val);
    _mapInitdTerms[name] = t.val;
  }
}

void
TheorEval::assignTokens(list<tToken> &sl)
{
  stringstream strexpr(_expr);
  int it = 0;
  const int nb = this->getNbins();
  char c;
  string term;
  tToken t;
  while (1){
    strexpr.get(c);
    if ( strexpr.eof() ) break;
    if ( isspace(c) ) continue; // skip whitespaces.
    // Oh noes! doesn't work after fortran reading expression with spaces :(.
    if ( isdigit(c) ) {  // process numbers
      term.assign(1,c);
      do {
        strexpr.get(c);
        if ( strexpr.eof() ) break;
        if ( isdigit(c) || c=='.' )  {
          term.append(1,c);
        }  else if ( c=='E' || c=='e' ) { // read mantissa including sign in scientific notation
          term.append(1,c);
          strexpr.get(c);
          if ( strexpr.eof() ) break;
          if ( isdigit(c) || c == '-' ){
            term.append(1,c);
          } else {
            cerr << "Theory expression syntax error: " << _expr << endl;
            hf_errlog(19072200, "F: Syntax error in theory expression, see stderr");
            return;
          }
        } else {
          strexpr.putback(c);
          break;
        }
      } while (1);
      double dterm = atof(term.c_str());

      t.opr = 0;
      t.name = term;
      t.val = new valarray<double>(dterm, nb);
      sl.push_back(t);
      continue;
    } else if ( isalpha(c) ) {  // process literal terms
      term.assign(1,c);
      while (strexpr.get(c) ) {
        if ( isalnum(c) ) term.append(1,c);
        else {
          strexpr.putback(c);
          break;
        }
      }
      if ( term == string("sum") ) { // special case for sum() function
        t.opr = 4;
        t.name = "sum";
        t.val = new valarray<double>(0., nb);
        sl.push_back(t);
        continue;
      }
      if ( term == string("spline") || term == string("splinederivative") )
      {
        // special case for natural cubic spline interpolation
        if ( term == string("spline"))
        {
          t.opr = 6;
          t.name = "spline";
        }
        else if ( term == string("splinederivative"))
        {
          t.opr = 7;
          t.name = "splinederivative";
        }
        t.ownsVal = false;
        // push spline
        sl.push_back(t);
        int& narg_spline = sl.back().narg;

        // process arguments
        t.narg = 0;
        t.opr = 0;
        // format: spline[x1,y1,x2,y2,x3,y3,x4,y4,...,x]
        strexpr.get(c);
        if(c != '[')
          hf_errlog(18090900, "F: Theory expression syntax error: expected [");
        narg_spline = 0;
        bool flagDone = false;
        while(true)
        {
          strexpr.get(c);
          int nsymbols = 0;
          term.assign(1,c);
          while(strexpr.get(c))
          {
            if(c == ',' || c == ']')
            {
              if(nsymbols == 0)
                hf_errlog(18090903, "F: Theory expression syntax error: error reading arguments");
              if(c == ']')
                flagDone = true;
              break;
            }
            if (!isalnum(c))
              hf_errlog(18090904, "F: Theory expression syntax error: error reading arguments");
            term.append(1,c);
            nsymbols++;
          }

          // have read new argument: push it
          if(nsymbols > 0)
          {
            initReactionToken(t,term);
            sl.push_back(t);
            narg_spline++;
            // finish reading spline arguments
            if(flagDone)
              break;
          }
          else
          {
            if(!flagDone)
              // should not be here
              //assert(0);
              hf_errlog(18090901, "F: Theory expression syntax error reading spline arguments");
            if(narg_spline % 2 == 0)
              hf_errlog(18090901, "F: Theory expression syntax error: spline expected odd number of arguments");
            if(narg_spline < 9)
              hf_errlog(18090902, "F: Theory expression syntax error: spline expected at least 9 arguments");
            break;
          }
        }
        continue;
      }
      if ( term == string("norm") )
      {
        // special case for normalised expression: norm(A)=A/sum(A)
        // (A is coomputed once)
        t.opr = 8;
        t.name = "norm";
        t.val = new valarray<double>(0., nb);
        sl.push_back(t);
        continue;
      }
      /*
      if ( term == string("avg") ) { // special case for avg() function
        t.opr = 4;
        t.name = "avg";
        t.val = new valarray<double>(0., nb);
        sl.push_back(t);
        continue;
      }
      */

      initReactionToken(t,term);
      sl.push_back(t);
      term.clear();
      continue;
    } else {
      switch(c){
        case '(': t.opr = -1; break;
        case ')': t.opr = -2; break;
        case '+': t.opr = 1; break;
        case '-': t.opr = 1; break;
        case '*': t.opr = 3; break;
        case '/': t.opr = 3; break;

        case '.': t.opr = 5; break;

        default: cout << "Unknown operator "<< c << " in expression " << _expr << endl;
      }
      t.name.assign(1,c);
      t.val = new valarray<double>(0., nb);
      sl.push_back(t);
    }
  }
}

void
TheorEval::convertToRPN(list<tToken> &sl)
{
  stack<tToken> tknstk;
  // convert to RPN
  while ( 0!=sl.size()){
    tToken t = sl.front();
    sl.pop_front();
    if ( 0 == t.opr ) {
      _exprRPN.push_back(t);
    }
    //if ( 4 == t.opr ){ // push functions
    //  tknstk.push(t);
    //}
    if ( t.opr >0 ) {
      while ( tknstk.size() > 0 && t.opr <= tknstk.top().opr ) {
        _exprRPN.push_back(tknstk.top());
        tknstk.pop();
      }

      tknstk.push(t);
    }
    if ( t.opr == -1 ){ tknstk.push(t); delete t.val;} // left parenthesis
    if ( t.opr == -2 ){                   // right parenthesis
      while ( tknstk.top().opr != -1 ) {
        if ( tknstk.size() == 0 ) cout << "ERROR: Wrong syntax in theoretical expression: "<< _expr << endl;
        _exprRPN.push_back(tknstk.top());
        tknstk.pop();
      }
      delete t.val;
      tknstk.pop();
    }
  }
  while ( tknstk.size() != 0 ){
    if (tknstk.top().opr == -1 ) cout << "ERROR: Wrong syntax in theoretical expression: "<< _expr << endl;
    _exprRPN.push_back(tknstk.top());
    tknstk.pop();
  }


  /*
  vector<tToken>::iterator it= _exprRPN.begin();
  for (;it!=_exprRPN.end(); it++){
    cout << it->name << " " ;
  }
  cout << endl;
  */
}

void TheorEval::initTerm(int iterm, valarray<double> *val)
{
  string term_type = _termTypes.at(iterm);
  if(term_type != string("reaction")) {
    std::cerr<<"[ERROR] Unknown term_type=\""<<term_type<<"\" in expression for term \""<<_termNames[iterm]<<'\"'<<std::endl;
    hf_errlog(15102301,"S: Unknown term type, see stderr");
  }
  this->initReactionTerm(iterm, val);
}

const string GetParamDS(const string&parName,const std::string&dsName,int dsIndex){
  using XFITTER_PARS::rootNode;
  //Get option referred to by "use:" from YAML
  YAML::Node parNode=rootNode[parName];
  if(!parNode.IsDefined()){
    cerr<<"[ERROR] Undefined key \""<<parName<<"\" in \"use:\" in dataset \""<<dsName<<"\" (index="<<dsIndex<<"). This key must be defined in YAML steering."<<endl;
    hf_errlog(19042000,"F: Undefined key in \"use:\", see stderr");
    std::abort();//unreachable
  }
  if(parNode.IsScalar())return parNode.as<string>();
  {
  YAML::Node nameNode=parNode[dsName];
  YAML::Node indexNode=parNode[dsIndex];
  if(nameNode.IsDefined()){
    if(indexNode.IsDefined()){
      cerr<<"[WARN] Value for key \""<<parName<<"\" referenced in \"use:\" is given both by dataset name (\""<<dsName<<"\") and by its index ("<<dsIndex<<"). Using the value provided by name and ignoring the value provided by index."<<endl;
      hf_errlog(19042001,"W: Value for \"use:\" key given both by name and by index, see stderr");
    }
    return nameNode.as<string>();
  }else if(indexNode.IsDefined()){
    return indexNode.as<string>();
  }
  }
  YAML::Node defaultNode=parNode["defaultValue"];//I think "default" would be nicer --Ivan
  if(defaultNode.IsDefined())return defaultNode.as<string>();
  cerr<<"[ERROR] No value fiven for key \""<<parName<<"\" referenced in \"use:\" in dataset \""<<dsName<<"\" (index="<<dsIndex<<")"<<endl;
  hf_errlog(19042002,"F: Key in \"use:\" has no value, see stderr");
  std::abort();//unreachable
}
void TheorEval::initReactionTerm(int iterm, valarray<double> *val, bool change_source){
  string reactionName = _termSources.at(iterm);
  string term_info =  _termInfos.at(iterm);
  if(beginsWith(reactionName,"use:")){//then redefine term source
    //I am not sure this works correctly right now --Ivan
    //Replace dsPars
    reactionName = GetParamDS(reactionName.substr(4),_ds_name,_dsIndex);
  }
  ReactionTheory* rt = xfitter::getReaction(reactionName);
  size_t termID=_dsId*1000+iterm;
  TermData*term_data=new TermData(termID,rt,this,term_info.c_str());
  term_data->val=val;
  if (change_source)
    {
      delete term_datas[iterm];
      term_datas[iterm] = term_data;
    }
  else
    term_datas.push_back(term_data);
  rt->initTerm(term_data);
}

void TheorEval::setBins(int nBinDim, int nPoints, int *binFlags, double *allBins,map<string,size_t>&_columnNameMap){
  //Copy data to _binFlags and _dsBins
  _binFlags.resize(nPoints);
  for(int i=0;i<nPoints;i++){
    _binFlags[i]=binFlags[i];
  }

  _dsBins.resize(nBinDim);
  for(int ibd = 0; ibd < nBinDim; ibd++){
    valarray<double>&bins=_dsBins[ibd];
    bins.resize(nPoints);
    for(int i=0;i<nPoints;i++){
      bins[i]=allBins[i*10+ibd];
    }
  }
  columnNameMap=_columnNameMap;
}

void TheorEval::Evaluate(valarray<double> &vte )
{
  // get values from grids
  this->updateReactionValues();

  // resize valarray for safe arithmetics (+,-,*,/), e.g. for APPLgrid when grid is longer than data
  auto resize_to_min = [](std::valarray<double>& v1, std::valarray<double>& v2) {
    if (v1.size() > v2.size()) {
      v1 = std::valarray<double>(&v1[0], v2.size());
    }
  };

  // calculate expression result
  stack<valarray<double> > stk;
  vector<tToken>::iterator it = _exprRPN.begin();
  while(it!= _exprRPN.end()){
    if ( it->opr < 0 ){
      cout << "ERROR: Expression RPN is wrong" << endl;
      return;
    }
    if ( it->opr == 0 ){
      stk.push(*(it->val));
    } else if ( it->name == string("sum") ){
      double sum = stk.top().sum();
      stk.top() = sum;
/*    } else if ( it->name == string("avg") ){
      if (0 == stk.top().size()) {
        cout << "ERROR: avg() argument dimension is 0." << endl;
      }
      double avg = stk.top().sum()/stk.top().size();
      stk.top() = avg;*/
    } else if ( it->name == string("spline") || it->name == string("splinederivative") )
    {
      // load all arguments
      int narg = it->narg;
      std::valarray<double> x0 = stk.top();
      stk.pop();
      int nsections = (it->narg - 1) / 2;
      std::valarray<std::valarray<double> > x(nsections);
      std::valarray<std::valarray<double> > y(nsections);
      for(int sect = nsections - 1; sect >= 0; sect--)
      {
        y[sect] = stk.top();
        stk.pop();
        x[sect] = stk.top();
        stk.pop();
      }
      auto result = x0;
      for(int p = 0; p < x0.size(); p++)
      {
        std::vector<double> xSpline(nsections);
        std::vector<double> ySpline(nsections);
        for(int sect = 0; sect < nsections; sect++)
        {
          xSpline[sect] = x[sect][p];
          ySpline[sect] = y[sect][p];
        }
        //TSpline3 spline("", &xSpline[0], &ySpline[0], ySpline.size());
        tk::spline spline;
        spline.set_points(xSpline, ySpline);
        if(it->name == string("spline"))
        {
          //result[p] = spline.Eval(x0[p]);
          result[p] = spline(x0[0]);
        }
        else if(it->name == string("splinederivative"))
        {
          //result[p] = spline.Derivative(x0[p]);
          result[p] = spline(x0[0], 1);
        }
      }
      stk.push(result);
    }
    else if ( it->name == string("norm") )
    {
      double sum = stk.top().sum();
      stk.top() = stk.top() / sum;
    } else if ( it->name == string("+") ){
      valarray<double> a(stk.top());
      stk.pop();
      resize_to_min(stk.top(), a);
      stk.top() += a;
    } else if ( it->name == string("-") ){
      valarray<double> a(stk.top());
      stk.pop();
      resize_to_min(stk.top(), a);
      stk.top() -= a;
    } else if ( it->name == string("*") ){
      valarray<double> a(stk.top());
      stk.pop();
      resize_to_min(stk.top(), a);
      stk.top() *= a;
    } else if ( it->name == string("/") ){
      valarray<double> a(stk.top());
      stk.pop();
      resize_to_min(stk.top(), a);
      stk.top() /= a;
    }
    else if ( it->name == string(".") ){
      valarray<double> temp;
      valarray<double> result;

      valarray<double> a(stk.top());
      int size_a = a.size();
      stk.pop();
      valarray<double> b(stk.top());
      int size_b = b.size();

      if(size_a % size_b == 0){  // Matrix * Vector
        int size_return = size_a / size_b;
        result.resize(size_return);
        for ( int n = 0; n < size_b; n++){
          temp.resize(size_return);
          temp = a[std::slice(n*size_return, size_return, 1)]; //creating nth colum vector
          temp *= b[n];
          result += temp;
        }
        stk.top() = result;
      }else if(size_b % size_a == 0){  //  Transposed(Vector)*Matrix -> Transposed(Matrix) vector
        int size_return = size_b / size_a;
        result.resize(size_return);
        for ( int n = 0; n < size_a; n++){
          temp.resize(size_return);
          temp = b[std::slice(n, size_return, size_a)]; // creating nth row vector -> nth colum vector
          temp *= a[n];
          result += temp;
        }
        stk.top() = result;
      }else{
        cout<<"ERROR: Dimensions do not match"<<endl;
      }
      /*if(it + 1 ->name == string("kmatrix")){//possible matrix matrix multiplication
          int nb1 = ?;//TODO find dimensions of matrices for check and multiplication
          int mb1 = ?;
          int nb2 = ?;
          int mb2 = ?;
          result.resize(mb1*nb2);
          for(int m = 0; m < mb1; m++){
              for(int n = 0; n < nb2; n++){
                  temp.resize(nb1);
                  temp = M.slize(m*nb1,1, nb);
                  temp *= M2.slize(n, mb2, nb2);
                  result[m*nb1 + n] = temp.sum();
              }
          }
      }*/
    }
    it++;
  }

  if (stk.size() != 1 ) {
    cout << "ERROR: Expression RPN calculation error." << endl;
    return;
  } else {
    vte = stk.top();
    //Normalised cross section
    if (_normalised)
      {
        double integral = 0;
        for (int bin = 0; bin < _binFlags.size(); bin++)
          if (!(vte[bin] != vte[bin])) //protection against nan
            integral += (_dsBins[1][bin] - _dsBins[0][bin]) * vte[bin];
        if (integral != 0)
          for (int bin = 0; bin < _binFlags.size(); bin++)
            vte[bin] *= _normalisation/integral;
      }
  }
}

void TheorEval::updateReactionValues(){
  map<string,valarray<double> > errors;//The errors returned by reaction are ignored
  //TODO: actually use the errors reported by ReactionTheory in chi2 calculation
  for(const auto td:term_datas){
    td->reaction->compute(td,*td->val,errors);
  }
}
const valarray<double>*TheorEval::getBinColumn(const string&n)const{
  auto it=columnNameMap.find(n);
  if(it==columnNameMap.end())return nullptr;
  return &_dsBins[it->second];
}

//This method is used by chi2scan to change the theory input file
void TheorEval::ChangeTheorySource(string term, string source)
{
  vector<string>::iterator found_term = find(_termNames.begin(), _termNames.end(), term);
  if ( found_term == _termNames.end())
    {
      string msg = (string) "S: Undeclared term " + term;
      hf_errlog_(14020603, msg.c_str(), msg.size());
    }
  int iterm = int(found_term-_termNames.begin());
  //  cout << "switch " << _termInfos[iterm] << " to " << source << endl;
  _termInfos[iterm] = source;

  initReactionTerm(iterm, _mapInitdTerms[term], true);
}

//This method is used by chi2scan to get the theory input file for the central prediction
string TheorEval::GetTheorySource(string term)
{
  vector<string>::iterator found_term = find(_termNames.begin(), _termNames.end(), term);
  if ( found_term == _termNames.end())
    {
      hf_errlog(14020603,(string) "S: Undeclared term " + term);
    }
  int iterm = int(found_term-_termNames.begin());
  return _termInfos.at(iterm);
}
