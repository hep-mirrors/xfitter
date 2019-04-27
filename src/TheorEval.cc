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
    if(it.val){
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

int
TheorEval::initTheory()
{
  list<tToken> sl;
  this->assignTokens(sl);
  this->convertToRPN(sl);
}

int
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
            cout << "Theory expression syntax error: " << _expr << endl;
            return -1;
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

      /*
      if ( term == string("avg") ) { // special case for avg() function
        t.opr = 4;
        t.name = "avg";
        t.val = new valarray<double>(0., nb);
        sl.push_back(t);
        continue;
      }
      */

      vector<string>::iterator found_term = find(_termNames.begin(), _termNames.end(), term);
      if ( found_term == _termNames.end() ) {
        cout << "Undeclared term " << term << " in expression " << _expr << endl;
        return -1;
      } else {
        t.opr = 0;
        t.name = term;
        if ( _mapInitdTerms.find(term) != _mapInitdTerms.end()){
          t.val = _mapInitdTerms[term];
        } else {
          t.val = new valarray<double>(0.,nb);
          this->initTerm(int(found_term-_termNames.begin()), t.val);
          _mapInitdTerms[term] = t.val;
        }
        sl.push_back(t);
      }
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

        case '.': t.opr = 5; break; //change

        default: cout << "Unknown operator "<< c << " in expression " << _expr << endl;
      }
      t.name.assign(1,c);
      t.val = new valarray<double>(0., nb);
      sl.push_back(t);
    }
  }
}

int
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
  string term_type=_termTypes.at(iterm);
  if(term_type!=string("reaction")) {
    std::cerr<<"[ERROR] Unknown term_type=\""<<term_type<<"\" in expression for term \""<<_termNames[iterm]<<'\"'<<std::endl;
    hf_errlog(15102301,"S: Unknown term type, see stderr");
  }
  this->initReactionTerm(iterm, val);
}
ReactionTheory*getReaction(const string&name){
  string libname = gReactionLibs[name];
  if (libname == "") {
    hf_errlog(16120501,"F: Reaction " +name+ " not present in Reactions.txt file");
  }
  if ( gNameReaction.find(name) == gNameReaction.end()) {
    string path_to_lib=PREFIX+string("/lib/")+libname;
    void *theory_handler = dlopen(path_to_lib.c_str(), RTLD_NOW);
    if (theory_handler == NULL)  {
      std::cerr<<"Failed to open shared library "<<path_to_lib<<" for "<<name<<"; error:\n"
               <<dlerror()<<"\n Check that the correct library is given in Reactions.txt"<<std::endl;
      hf_errlog(16120502,"F: Failed to open reaction shared library, see stderr for details");
    }

    // reset errors
    dlerror();

    void*dispatch_theory=dlsym(theory_handler, "create");
    ReactionTheory*rt=(ReactionTheory*)((void*(*)())dispatch_theory)();//Create the ReactionTheory
    gNameReaction[name] = rt;
    // First make sure the name matches:
    if ( rt->getReactionName() == name) {
      hf_errlog(17041610,"I: Loaded reaction "+name);
    }else{
      hf_errlog(16120801,"F: Reaction "+name+" does not match library: "+rt->getReactionName());
    }
    // initialize
    rt->atStart();
    return rt;
  } else {
    return gNameReaction.at(name);
  }
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
  }else{
    return indexNode.as<string>();
  }
  }
  YAML::Node defaultNode=parNode["defaultValue"];//I think "default" would be nicer --Ivan
  if(defaultNode.IsDefined())return defaultNode.as<string>();
  cerr<<"[ERROR] No value fiven for key \""<<parName<<"\" referenced in \"use:\" in dataset \""<<dsName<<"\" (index="<<dsIndex<<")"<<endl;
  hf_errlog(19042002,"F: Key in \"use:\" has no value, see stderr");
  std::abort();//unreachable
}
void TheorEval::initReactionTerm(int iterm, valarray<double> *val){
  string term_source = _termSources.at(iterm);
  string term_info =  _termInfos.at(iterm);
  if(beginsWith(term_source,"use:")){//then redefine term source
    //I am not sure this works correctly right now --Ivan
    //Replace dsPars
    term_source=GetParamDS(term_source.substr(4),_ds_name,_dsIndex);
  }
  ReactionTheory*rt=getReaction(term_source);
  size_t termID=_dsId*1000+iterm;
  TermData*term_data=new TermData(termID,rt,this,term_info.c_str());
  term_data->val=val;
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
    } else if ( it->name == string("+") ){
      valarray<double> a(stk.top());
      stk.pop();
      stk.top() += a;
    } else if ( it->name == string("-") ){
      valarray<double> a(stk.top());
      stk.pop();
      stk.top() -= a;
    } else if ( it->name == string("*") ){
      valarray<double> a(stk.top());
      stk.pop();
      stk.top() *= a;
    } else if ( it->name == string("/") ){
      valarray<double> a(stk.top());
      stk.pop();
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
            vte[bin] /= integral;
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

/* What are those? They are currently unused, and I am not sure they work correctly now, with TermData. --Ivan
void TheorEval::ChangeTheorySource(string term, string source)
{
  vector<string>::iterator found_term = find(_termNames.begin(), _termNames.end(), term);
  if ( found_term == _termNames.end())
    {
      string msg = (string) "S: Undeclared term " + term;
      hf_errlog_(14020603, msg.c_str(), msg.size());
    }
  int iterm = int(found_term-_termNames.begin());
  //  cout << "switch " << _termSources[iterm] << " to " << source << endl;
  _termSources[iterm] = source;

  initTerm(int(found_term-_termNames.begin()), _mapInitdTerms[term]);
}

string TheorEval::GetTheorySource(string term)
{
  vector<string>::iterator found_term = find(_termNames.begin(), _termNames.end(), term);
  if ( found_term == _termNames.end())
    {
      hf_errlog(14020603,(string) "S: Undeclared term " + term);
    }
  int iterm = int(found_term-_termNames.begin());
  return _termSources[iterm];
}
*/
