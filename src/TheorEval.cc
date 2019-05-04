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
#include "CommonGrid.h"
#include "ReactionTheory.h"
#include "xfitter_cpp.h"
#include "get_pdfs.h"
#include <string.h> 

#include <yaml-cpp/yaml.h>
#include "xfitter_pars.h"

// ROOT spline can be uncommented (here and below in the code) and used e.g. for cross checks (obviously requires ROOT)
//#include <TSpline.h>
#include <spline.h>

using namespace std;

// extern struct ord_scales {
//    double datasetmur[150];
//    double datasetmuf[150];
//    int datasetiorder[150];
// } cscales_;



TheorEval::TheorEval(const int dsId, const int nTerms, const std::vector<string> stn, const std::vector<string> stt, 
                     const std::vector<string> sti, const std::vector<string> sts, const string& expr) : _dsId(dsId), _nTerms(nTerms)
{
  // _iOrd = cscales_.datasetiorder[_dsId-1];
  // _xmur = cscales_.datasetmur[_dsId-1];
  // _xmuf = cscales_.datasetmuf[_dsId-1];
  for (int it= 0 ; it<nTerms; it++ ){
    _termNames.push_back(stn[it]);
    _termTypes.push_back(stt[it]);
    _termInfos.push_back(sti[it]);
    _termSources.push_back(sts[it]);
  }
  _expr.assign(expr);

  _ppbar = false;
}

TheorEval::~TheorEval()
{
  map<CommonGrid*, valarray<double>* >::iterator itm = _mapGridToken.begin();
  for (; itm!= _mapGridToken.end(); itm++){
    delete itm->first;
  }

  vector<tToken>::iterator it = _exprRPN.begin();
  for (; it!=_exprRPN.end(); it++){
    if ( ! it->val ) { delete it->val; it->val = NULL; }
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
        // push spline
        sl.push_back(t);
        int& narg_spline = sl.back().narg;

        // process arguments
        t.val = new valarray<double>(0., nb);
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
            }
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

        case '.': t.opr = 5; break;

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

int
TheorEval::initTerm(int iterm, valarray<double> *val)
{
   
  string term_type =  _termTypes.at(iterm);
  if ( term_type.find("grid") != string::npos || term_type.find("ast") != string::npos ){ //appl'grid' or f'ast'NLO
    this->initGridTerm(iterm, val);
  } else if ( term_type == string("reaction")) {
    this->initReactionTerm(iterm, val);
  } else if ( term_type == string("kfactor")) {
    this->initKfTerm(iterm, val);
  } else {
    int id = 15102301;
    char text[] = "S: Unknown term type in expression for term";
    std::cout << "Unknown term type in expression for term " << _termNames[iterm] << std::endl;
    int textlen = strlen(text);
    hf_errlog_(id, text, textlen);
    return -1;
  }
}

int
TheorEval::initGridTerm(int iterm, valarray<double> *val)
{
  string term_source = _termSources.at(iterm);
  string term_type =  _termTypes.at(iterm);
  string term_info =  _termInfos.at(iterm);
  CommonGrid *g = new CommonGrid(term_type, term_source); 
  if (  term_type.find("grid") != string::npos && !(term_type.find("apfel") != string::npos) ) {
     // set the collision for the grid term
     string collision ("pp"); // default is pp

     // this is to have backward-compatibility with Tevatron datasets
     if ( _ppbar ) collision.assign(string("ppbar"));
     // otherwise we check beams in the TermInfo lines
     else {
       size_t beams_pos = term_info.find(string("beams"));
       if ( beams_pos != string::npos ){
         size_t semicol_pos = term_info.find(';', beams_pos);
         size_t eq_pos = term_info.find('=', beams_pos);
	 collision.assign(term_info.substr(eq_pos+1, semicol_pos - eq_pos-1));
       }
     }

     // strip blanks
     collision.erase(std::remove(collision.begin(), collision.end(), ' '), collision.end());

     // and set the collision
     g->SetCollisions(collision);
     g->SetDynamicScale( _dynamicscale );
     
     // check the binning with the grids, will be ignored for normalisation grids
     g->checkBins(_binFlags, _dsBins);
  }
  else if ( term_type.find("ast") != string::npos ){
     bool PublicationUnits = true; // todo: take from new steering flag 'TermNorm'
     //FastNLOReader* fnlo = g->getHBins().back().f;
     FastNLOxFitter* fnlo = g->getHBins().back().f; 
     if(PublicationUnits)
	fnlo->SetUnits(fastNLO::kPublicationUnits);
     else 
	fnlo->SetUnits(fastNLO::kAbsoluteUnits);
     
  // OZ 9.05.2017
  if(term_type.find("norm") != string::npos)
    fnlo->SetUnits(fastNLO::kAbsoluteUnits);

     // --- set scales
     if(_MurDef>=0)
	fnlo->SetMuRFunctionalForm((fastNLO::EScaleFunctionalForm) ((int) (_MurDef)));
     if(_MufDef>=0)
	fnlo->SetMuFFunctionalForm((fastNLO::EScaleFunctionalForm) ((int) (_MufDef)));
     if ( _xmur!=1 || _xmuf!=1 )
	fnlo->SetScaleFactorsMuRMuF(_xmur, _xmuf);

     // --- set order  
     if ( _iOrd == 1 ) {
	fnlo->SetContributionON(fastNLO::kFixedOrder,1,false); // switch 'off' NLO
     }
     else if (_iOrd==2) {
	// that's fastNLO default
     }
     else if (_iOrd==3) {
	fnlo->SetContributionON(fastNLO::kFixedOrder,2,true); // switch 'on' NNLO
     }
     else {
	printf("fastNLO pert. order is not defined, ordercalc = %d:\n",_iOrd);
	exit(1);
     }
  }
  else if ( term_type.find("apfel") != string::npos ){
    // set the collision for the grid term                                                                                                                                                  
    string collision ("pp"); // default is pp                                                                                                                                               

    // this is to have backward-compatibility with Tevatron datasets                                                                                                                        
    if ( _ppbar ) collision.assign(string("ppbar"));
    // otherwise we check beams in the TermInfo lines                                                                                                                                       
    else {
      size_t beams_pos = term_info.find(string("beams"));
      if ( beams_pos != string::npos ){
	size_t semicol_pos = term_info.find(';', beams_pos);
	size_t eq_pos = term_info.find('=', beams_pos);
	collision.assign(term_info.substr(eq_pos+1, semicol_pos - eq_pos-1));
      }
    }

    // strip blanks                                                                                                                                                                         
    collision.erase(std::remove(collision.begin(), collision.end(), ' '), collision.end());

    // and set the collision                                                                                                                                                                
    g->SetCollisions(collision);

    //g->SetCollisions("pp");
  }
  /*
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

  */

  // associate grid and valarray pointers in token
  _mapGridToken[g] = val;
}

int
TheorEval::initReactionTerm(int iterm, valarray<double> *val)
{
  string term_source = _termSources.at(iterm);
  string term_type =  _termTypes.at(iterm);
  string term_info =  _termInfos.at(iterm);
//  ReactionTheory *rt = ReactionTheoryDispatcher::getInstance().getReactionTheory(_termSources.at(iterm)); 
  
  // Re-define term-source if "use:" string is found:
  if ( term_source.find("use:") != std::string::npos ) {
    term_source =  GetParamDS(term_source.substr(4),GetDSname(),_dsPars["FileIndex"]);
  }

  string libname = gReactionLibs[term_source];
  if (libname == "") {
    string text = "F: Reaction " +term_source + " not present in Reactions.txt file";
    hf_errlog_(16120501,text.c_str(),text.size());
  }

  ReactionTheory * rt;
  if ( gNameReaction.find(term_source) == gNameReaction.end()) {
    void *theory_handler = dlopen((PREFIX+string("/lib/")+libname).c_str(), RTLD_NOW);
    if (theory_handler == NULL)  { 
      std::cout  << dlerror() << std::endl;
      string text = "F: Reaction shared library ./lib/"  + libname  +  " not present for " +term_source + ". Check Reactions.txt file" ;
      hf_errlog_(16120502,text.c_str(),text.size());
    }
    
    // reset errors
    dlerror();
 
    create_t *dispatch_theory = (create_t*) dlsym(theory_handler, "create");
    rt = dispatch_theory();
    gNameReaction[term_source] = rt;


  // First make sure the name matches:
    if ( rt->getReactionName() == term_source) {
      string msg =  "I: Use reaction "+ rt->getReactionName();
      hf_errlog_(17041610+_dsId,msg.c_str(),msg.size());
    }
    else {
      string text = "F: Reaction " +term_source + " does not match with library: "+rt->getReactionName();
      hf_errlog_(16120801,text.c_str(),text.size());
    }


    // Some initial stuff:

    // transfer the parameters:
    rt->setxFitterParameters(XFITTER_PARS::gParameters);
    rt->setxFitterParametersI(XFITTER_PARS::gParametersI);
    rt->setxFitterParametersS(XFITTER_PARS::gParametersS);
    rt->setxFitterparametersVec(XFITTER_PARS::gParametersV);
    rt->setxFitterparametersYaml(XFITTER_PARS::gParametersY);
  
    // Override some global pars for reaction specific:
    if ( XFITTER_PARS::gParametersY[term_source] ) {
      rt->resetParameters(XFITTER_PARS::gParametersY[term_source]);
    }

    // set alpha_S, pdfs:
    rt->setEvolFunctions( &HF_GET_ALPHASQ_WRAP, &g2Dfunctions);

    // simplify interfaces to LHAPDF:
    rt->setXFX(&HF_GET_PDFSQ_WRAP);           // proton
    rt->setXFX(&HF_GET_PDFSQ_BAR_WRAP,"pbar"); // anti-proton
    rt->setXFX(&HF_GET_PDFSQ_N_WRAP,"n");   // neutron

    // initialize
    if (rt->initAtStart("") != 0) {
      // failed to init, somehow ...
      string text = "F:Failed to init reaction " +term_source  ;
      hf_errlog_(16120803,text.c_str(),text.size());
    };
 
  } else {
    rt = gNameReaction[term_source];
  }


  /// Reaction-term / dataset specific:

  // Set bins
  rt->setBinning(_dsId*1000+iterm, &gDataBins[_dsId]);
  
  // split term_info into map<string, string> according to key1=value1:key2=value2:key=value3...
  map<string, string> pars = SplitTermInfo(term_info);

  // and transfer to the module
    //printf("pars\n");
    //for(map<string,string>::iterator it = pars.begin(); it != pars.end(); it++)
    //  printf("%s = %s\n", it->first.c_str(), it->second.c_str());
  rt->setDatasetParameters(_dsId*1000+iterm, pars, _dsPars);

  _mapReactionToken[ std::pair<ReactionTheory*,int>(rt,iterm) ] = val;
}

int
TheorEval::initKfTerm(int iterm, valarray<double> *val)
{
  string term_source(_termSources.at(iterm));
  // read k-Factor table and compare it's binning to the data
  cout << "reading k-factor table from " << term_source << endl;
  vector<double> tv;
  vector<vector<double> > bkf(_dsBins.size(),tv);
  vector<double> vkf;
  ifstream kff(term_source.c_str());
  string line;
  if (kff.is_open()){
    while (1) {
      getline(kff,line);
      if (true == kff.eof()) break;
      if (line.at(0) == '#' ) continue; //ignore comments
      line.erase(line.find_last_not_of(" \n\r\t")+1); // trim trailing whitespaces
      stringstream sl(line);
      // first count words
      int nw(0);
      while (sl.good()) {
        string ts;
	sl >> ts;
	nw++;
      }
      // check that we have even number of bins (low and high columns)
      if (0!=(nw-1)%2) {
        int id = 14040340;
        char text[] = "S: Bad number of bins in k-factor file. Each bin must have low and high value.";
        int textlen = strlen(text);
        hf_errlog_(id, text, textlen);
      }
      // check that the number of bins is equal to data binning dimension
      if ((nw-1) != _dsBins.size()) {
        int id = 14040341;
        char text[] = "S: Bad number of bins in k-factor file. Must be equal to data binning dimension.";
        int textlen = strlen(text);
        hf_errlog_(id, text, textlen);
      }

      // now read bins
      sl.clear();
      sl.seekg(0);
      sl.str(line);
      double tb(0);
      for (int iw=0; iw<nw-1; iw++) {
        sl >> tb;
	bkf.at(iw).push_back(tb);
      }
      
      // and k-factor
      sl>>tb;
      vkf.push_back(tb);
    }
    kff.close();
  } else {
    int id = 14040339;
    char text[] = "S: Error reading k-factor file.";
    int textlen = strlen(text);
    hf_errlog_(id, text, textlen);
  }

  // check that k-factor file binning is compatible with data
  for (int iv = 0; iv<_dsBins.size(); iv++){
    for (int ib = 0; ib<_dsBins.at(iv).size(); ib++){
      if ( _binFlags.at(ib) == 0 ) continue;
      if ( 0 == (_binFlags.at(ib) & 2) ) {
        if (fabs(bkf[iv][ib] - _dsBins[iv][ib]) > 100*DBL_MIN) { 
          int id = 14040338;
          char text[] = "S: Data and grid bins don't match.";
          int textlen = strlen(text);
          hf_errlog_(id, text, textlen);
          return -1;
        }
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

  return _dsBins.size();
}

int 
TheorEval::setCKM(const vector<double> &v_ckm)
{
#ifdef APPLGRID_CKM
  map<CommonGrid*, valarray<double>* >::iterator itm = _mapGridToken.begin();
  for(; itm != _mapGridToken.end(); itm++){
    itm->first->setCKM(v_ckm);
  }
#else
   int id = 611201320;
   char text[] = "S: Cannot set CKM in Applgrid, use v1.4.33 or higher";
   int textlen = strlen(text);
   hf_errlog_(id, text, textlen);
#endif
}

int
TheorEval::Evaluate(valarray<double> &vte )
{
  // get values from grids
  this->getGridValues();
  this->getReactionValues();

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
      for(int arg = 0; arg < narg; arg++)
      {
        //it++;
        //stk.push(*(it->val));
      }
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
    }
    else if ( it->name == string("+") ){
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
		char error[] = "ERROR: Dimensions do not match ";
		cout<<error<<endl;}


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
	if (integral != 0)
	  for (int bin = 0; bin < _binFlags.size(); bin++)
	    vte[bin] /= integral;
      }
    //vte /= _units;
  }
}

int
TheorEval::getGridValues()
{
  map<CommonGrid*, valarray<double>*>::iterator itm;
  for(itm = _mapGridToken.begin(); itm != _mapGridToken.end(); itm++){
    CommonGrid* g = itm->first;
    vector<double> xs;
    std::vector< std::vector<double> > result = g->vconvolute(_iOrd, _xmur, _xmuf);        
    for(int i = 0; i < result.size(); i++)
      for(int j = 0; j < result[i].size(); j++)
        xs.push_back(result[i][j]);

    (itm->second)->resize(xs.size());
    *(itm->second) = valarray<double>(xs.data(), xs.size());
    /*
    for (int i = 0; i<xs.size(); i++){
      cout << xs[i] << endl;
    }
    */
    
    
  }
}

int
TheorEval::getReactionValues()
{
  //  map<ReactionTheory*, valarray<double>*>::iterator itm;
  for(auto itm = _mapReactionToken.begin(); itm != _mapReactionToken.end(); itm++){
    ReactionTheory* rt = (itm->first).first;
    int idTerm =  (itm->first).second;
    map<string, valarray<double> > errors;
     
    int result =  rt->compute(_dsId*1000+idTerm, *(itm->second), errors);
     
    if (result != 0) {
      string text = "F:(from TheorEval::getReactionValues)  Failed to compute theory";
      hf_errlog_(16081202,text.c_str(),text.size());
    }
  }
  
  return 1;
}


int
TheorEval::getNbins()
{
  return _dsBins[0].size();
}

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

  //delete old applgrid
  map<CommonGrid*, valarray<double>* >::iterator itm = _mapGridToken.begin();
  for (; itm!= _mapGridToken.end(); itm++)
    {
      if (itm->second == _mapInitdTerms[term])
	{
	  delete itm->first;
	  _mapGridToken.erase(itm);
	  break;
	}
    }

  initTerm(int(found_term-_termNames.begin()), _mapInitdTerms[term]);
}

string TheorEval::GetTheorySource(string term)
{
  vector<string>::iterator found_term = find(_termNames.begin(), _termNames.end(), term);
  if ( found_term == _termNames.end())
    {
      string msg = (string) "S: Undeclared term " + term;
      hf_errlog_(14020603, msg.c_str(), msg.size());
    }
  int iterm = int(found_term-_termNames.begin());
  return _termSources[iterm];
}

map<string, string> TheorEval::SplitTermInfo(const string& term_info)
{
  // split term_info into map<string, string> according to key1=value1:key2=value2:key=value3...
  map<string, string> pars;
  std::size_t pos0 = 0;
  while(pos0 < term_info.size())
  {
    std::size_t pos1 = term_info.find("=", pos0);
    std::size_t pos2 = term_info.find(":", pos0);
    if(pos2 == std::string::npos) // last key=value does not have trailing :
      pos2 = term_info.size();
    // check for possible wrong format
    if(pos0 == 0 && pos1 == std::string::npos)
    { // no = in non empty term_info
      string text = "W: Wrong TermInfo format. The correct format is key1=value1:key2=value2:...";
      hf_errlog_(17020101, text.c_str(), text.size());
    }
    if(pos2 < pos1)
    { // two : : without = between them
      string text = "W: Wrong TermInfo format. The correct format is key1=value1:key2=value2:...";
      hf_errlog_(17020101, text.c_str(), text.size());
    }
    std::size_t pos11 = term_info.find("=", pos1 + 1);
    if(pos11 != std::string::npos && pos11 < pos2)
    { // two = = without : between them
      string text = "W: Wrong TermInfo format. The correct format is key1=value1:key2=value2:...";
      hf_errlog_(17020101, text.c_str(), text.size());
    }
    // split into key value pair
    std::string key = std::string(term_info, pos0, pos1 - pos0);
    std::string value = std::string(term_info, pos1 + 1, pos2 - pos1 - 1);
    // check if this key already exists
    if(pars.find(key) != pars.end())
    {
      string text = "W: Replacing existing key when reading TermInfo";
      hf_errlog_(17020102, text.c_str(), text.size());
    }
    // add to map
    pars[key] = value;
    // start next search iteration after current :
    pos0 = pos2 + 1;
  }
  //printf("read term_info %s\n", term_info.c_str());
  //for(map<string, string>::iterator it = pars.begin(); it != pars.end(); it++)
  //  printf("  %s=%s\n", (it->first).c_str(), (it->second).c_str());
  return pars;
}

const std::string GetParamDS(const std::string& ParName, const std::string& DSname, int DSindex) {
  // First check the list of strings, if present there. If yes, just return
  if ( XFITTER_PARS::gParametersS.find(ParName) != XFITTER_PARS::gParametersS.end() ) {
    return XFITTER_PARS::gParametersS[ParName];
  }
  // Now check the complex list:
  if ( XFITTER_PARS::gParametersY.find(ParName) != XFITTER_PARS::gParametersY.end() ) {
    YAML::Node Node = XFITTER_PARS::gParametersY[ParName];

    // Default:
    if ( Node["defaultValue"]) {
      std::string Val = Node["defaultValue"].as<string>();

      if (Node[DSname]) {
	Val = Node[DSname].as<string>();
      }
      if (Node[DSindex]) {
	Val = Node[DSindex].as<string>();
      }

      return Val;
    }
    else {
      string text = "F: missing value field for parameter " + ParName;
      hf_errlog_(17041101,text.c_str(),text.size());            
      return "";
    }
  }
  else {
    return "";
  }
}
