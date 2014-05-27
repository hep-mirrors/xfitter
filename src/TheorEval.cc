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

#include "TheorEval.h"

#include "appl_grid/appl_grid.h"

using namespace std;
using namespace appl;

#include "herafitter_cpp.h"

void appl_fnpdf_bar(const double& x, const double& Q, double* f)
{
  appl_fnpdf_( x, Q, f);    
  double fbar[13];
  for (int fi = 0; fi < 13; fi++)
    fbar[fi] = f[12-fi];
  for (int fi = 0; fi < 13; fi++)
    f[fi] = fbar[fi];
  return; 
}


TheorEval::TheorEval(const int dsId, const int nTerms, const string* stn, const string* stt, 
                     const string* sts, const string& expr)
{
  _dsId = dsId;
  _nTerms = nTerms;
  for (int it= 0 ; it<nTerms; it++ ){
    _termNames.push_back(stn[it]);
    _termTypes.push_back(stt[it]);
    _termSources.push_back(sts[it]);
  }
  _expr.assign(expr);
}

TheorEval::~TheorEval()
{
  map<appl::grid*, valarray<double>* >::iterator itm = _mapGridToken.begin();
  for (; itm!= _mapGridToken.end(); itm++){
    delete itm->first;
  }

  vector<tToken>::iterator it = _exprRPN.begin();
  for (; it!=_exprRPN.end(); it++){
    if ( ! it->val ) { delete it->val; it->val = NULL; }
  }
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
        
      if ( term == string("avg") ) { // special case for avg() function
        t.opr = 4;
        t.name = "avg";
	t.val = new valarray<double>(0., nb);
	sl.push_back(t);
	continue;
      }
        
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
      if ( tknstk.size() > 0 && t.opr <= tknstk.top().opr ) {
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
    if (it->opr == 0 ){
      cout << it->name << "\t" << (*it->val)[0] << endl;
    }
  }
  */
}

int
TheorEval::initTerm(int iterm, valarray<double> *val)
{
  string term_type =  _termTypes.at(iterm);
  if ( term_type == string("applgrid") ){
    this->initGridTerm(iterm, val);
  } else if (term_type == string("kfactor")) {
    this->initKfTerm(iterm, val);
  } else {
    cout << "ERROR: Unknown term type \"" << term_type << "\""<< endl;
    return -1;
  }
}

int
TheorEval::initGridTerm(int iterm, valarray<double> *val)
{
  string term_source = _termSources.at(iterm);
  appl::grid *g = new appl::grid(term_source);
  if (_dynamicscale != 0)
    {
#ifdef APPLGRID_DYNSCALE
      g->setDynamicScale( _dynamicscale );
#else
      int id = 2204201401;
      char text[] = "S: Cannot use dynamic scale emulation in Applgrid, use v1.4.43 or higher";
      int textlen = strlen(text);
      hf_errlog_(id, text, textlen);
#endif
    }

  g->trim();

  // read binning information from grid and compare it to that of data
  int n_agbins = g->Nobs();
  /*
  if (n_agbins != _dsBins.at(0).size() ) {
    cout << "ERROR: number of bins doesn't match for " << term_source << " in dataset " << _dsId << endl;
    return -1;
  }
  */
  
  // I assume that we have only 1d binning in the dataset.
  // Didn't find a good way to deal with multidimentional binning,
  // since applgrids are always 1d.  -- AS
  for (int igb = 0; igb <n_agbins-1; igb++){
    if ( igb >= _binFlags.size() ) break;
    if ( _binFlags.at(igb) == 0  ) continue;
    if ( 0 == (_binFlags.at(igb) & 2) ) {
      if (fabs(g->obslow(igb) - _dsBins.at(0).at(igb)) > 100*DBL_MIN ||
          fabs(g->obslow(igb+1) - _dsBins.at(1).at(igb)) > 100*DBL_MIN) { 
        cout << "ERROR: Bin boundaries don't match for bin" << igb << " in dataset " << _dsId << endl;
        return -1;
      }
    }
  }

  // associate grid and valarray pointers in token
  _mapGridToken[g] = val;
}

int
TheorEval::initKfTerm(int iterm, valarray<double> *val)
{
  string term_source(_termSources.at(iterm));
  // read k-Factor table and compare it's binning to the data
  cout << "reading k-factor table from " << term_source << endl;
  ifstream kff(term_source.c_str());
  vector<double> vbl, vbu, vkf;
  if (kff.is_open()) {
    while (1){
      double bl(0.), bu(0.), kf(0.);
      kff >> bl >> bu >> kf;
      if (true == kff.eof()) break;
      vbl.push_back(bl);
      vbu.push_back(bu);
      vkf.push_back(kf);
    }
    kff.close();
  } else {
    cout << "ERROR: problem opening k-factor file " << term_source << endl;
    return -1;
  }

  // check that k-factor file binnig is compatible with data
  if ( vkf.size() != this->getNbins()) {
    cout << "ERROR: number of bins doesn't match for " << term_source << " in dataset " << _dsId << endl;
    return -1;
  }
  vector<double>::iterator ikf = vkf.begin();
  for (; ikf < vkf.end(); ikf++){
    int ind = int(ikf-vkf.begin());
    if ( _binFlags.at(ind) == 0 ) continue;
    if ( 0 == (_binFlags.at(ind) & 2) ) {
      if (fabs(vbl.at(ind) - _dsBins.at(0).at(ind)) > 100*DBL_MIN ||
          fabs(vbu.at(ind) - _dsBins.at(1).at(ind)) > 100*DBL_MIN) { 
        cout << "ERROR: Bin boundaries don't match for bin" << ind << " in dataset " << _dsId << endl;
        return -1;
      }
    }
  }

  // write k-factor array to the token valarray
  *val = valarray<double>(vkf.data(), vkf.size());
}  

int
TheorEval::setBins(int nBinDim, int nPoints, int *binFlags, double *allBins)
{
  for(int ip = 0; ip<nPoints; ip++){
    _binFlags.push_back(binFlags[ip]);
  }

  for(int ibd = 0; ibd < nBinDim; ibd++){ 
    vector<double> bins;
    bins.clear();
    for(int ip = 0; ip<nPoints; ip++){
      bins.push_back(allBins[ip*10 + ibd]);
    }
    _dsBins.push_back(bins);
  }
}

int 
TheorEval::setCKM(const vector<double> &v_ckm)
{
#ifdef APPLGRID_CKM
  map<appl::grid*, valarray<double>* >::iterator itm = _mapGridToken.begin();
  for(; itm != _mapGridToken.end(); itm++){
    appl::grid* g = itm->first;
    g->setckm(v_ckm);
  }
#else
   int id = 611201320;
   char text[] = "S: Cannot set CKM in Applgrid, use v1.4.33 or higher";
   int textlen = strlen(text);
   hf_errlog_(id, text, textlen);
#endif
}

int
TheorEval::Evaluate(const int iorder, const double mur, const double muf, valarray<double> &vte )
{
  // get values from grids
  this->getGridValues(iorder, mur, muf);

  // calculate expression result
  stack<valarray<double> > stk;
  vector<tToken>::iterator it = _exprRPN.begin();
  while(it!= _exprRPN.end()){
    if ( it->opr < 0 ){
      cout << "ERROR: Expression RPN is wrong" << endl;
      return -1;
    }
    if ( it->opr == 0 ){
      stk.push(*(it->val));
    } else if ( it->name == string("sum") ){
      double sum = stk.top().sum();
      stk.top() = sum;
    } else if ( it->name == string("avg") ){
      if (0 == stk.top().size()) {
        cout << "ERROR: avg() argument dimension is 0." << endl;
      }
      double avg = stk.top().sum()/stk.top().size();
      stk.top() = avg;
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

    it++;
  }

  if (stk.size() != 1 ) {
    cout << "ERROR: Expression RPN calculation error." << endl;
    return -1;
  } else {
    vte = stk.top();
    //Normalised cross section
    if (_normalised)
      {
	double integral = 0;
	for (int bin = 0; bin < _binFlags.size(); bin++)
	  if (!(vte[bin] != vte[bin])) //protection against nan
	    integral += (_dsBins.at(1).at(bin) - _dsBins.at(0).at(bin)) * vte[bin];
	for (int bin = 0; bin < _binFlags.size(); bin++)
	  vte[bin] /= integral;
      }
    vte /= _units;
  }
}

int
TheorEval::getGridValues(const int iorder, const double mur, const double muf)
{
  map<appl::grid*, valarray<double>* >::iterator itm = _mapGridToken.begin();
  for(; itm != _mapGridToken.end(); itm++){
    appl::grid* g = itm->first;
    vector<double> xs;
    if (_ppbar)
      xs = g->vconvolute(appl_fnpdf_, appl_fnpdf_bar, appl_fnalphas_, iorder, mur, muf);
    else
      xs = g->vconvolute(appl_fnpdf_, appl_fnalphas_, iorder, mur, muf);

    *(itm->second) = valarray<double>(xs.data(), xs.size());
    /*
    for (int i = 0; i<xs.size(); i++){
      cout << xs[i] << endl;
    }
    */
  }
}

int
TheorEval::getNbins()
{
  return _dsBins[0].size();
}
