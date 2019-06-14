#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>
#include<sys/stat.h>
#include<sys/types.h>
#include"xfitter_pars.h"
#include"xfitter_steer.h"
#include"xfitter_cpp_base.h"
#include"BaseEvolution.h"
#include"BaseMinimizer.h"
using namespace std;
using xfitter::BaseEvolution;

/*!
\brief Clip points to a given range
\param a      Lower edge of range
\param b      Upper edge of range
\param points Sorted array of points
\return Sorted array containing a, b and all points from input array that are >a and <b
*/
vector<double> clipPoints(double a, double b, const vector<double>& points){
  vector<double> ret;
  ret.push_back(a);
  auto it  = points.begin();
  auto end = points.end();
  for (; it<end; it++) {
    if (*it>a) break;
  }
  for( ; it<end; it++){
    if (*it>=b) break;
    ret.push_back(*it);
  }
  ret.push_back(b);
  return ret;
}

//Webster's method
/*
\brief Proportionally distribute a given number of points between bins
\details
  Divide a given number of points between some bins, proportionally to weights assigned to the bins, using Webster's method.
  Let N --- total number of points to distribute; W --- array of weights assigned to each bin; then, according to Webster's method, bin i gets
  R[i]=round(x*W[i])
  where x --- some value common for all bins, which is determined from condition that the total number of points assigned to all bins is equal to N.
  This function is used to distribute grid points between subintervals.
\param[in] W Array of weights assigned to bins, len(W) is number of bins, W[i] is weight assigned to bin i
\param[in] N Number of points to distribute
\return Integer array R, where R[i] is number of points assigned to bin i, which is proportional to its weight W[i]. It is guaranteed that sum(R)=N.
*/
vector<size_t> apportionWebsters(const vector<double>& W, const size_t N){
  size_t size = W.size();
  double Wsum = accumulate(W.begin(), W.end(), double(0));//Wsum=sum of elements of W
  if (std::isnan(Wsum)) {
    cerr<<"[FATAL ERROR] nan in apportionWebsters"<<endl;
    abort();
  }
  //start with an approximate solution
  vector<size_t> R(size);
  for (size_t i=0; i<size; ++i) R[i] = round(W[i]/Wsum*N);
  size_t Rsum = accumulate(R.begin(), R.end(), size_t(0));
  //if the rounded number of points is as requested, we are lucky: return
  if (Rsum == N) return R;
  //else add or remove some seates so that N==sum(R)
  vector<double> P(size);//P[i] is priority
  if (Rsum<N) {
    for (size_t i=0; i<size; ++i) P[i] = W[i] / (2*double(R[i])+1);
    while (true) {
      size_t i = max_element(P.begin(), P.end()) - P.begin();
      R[i]++;
      P[i] = W[i] / (2*double(R[i])+1);
      if (++Rsum == N) return R;
    }
  }else{
    for (size_t i=0; i<size; ++i) P[i] = (2*double(R[i])-1) / W[i];
    while (true) {
      size_t i = max_element(P.begin(), P.end()) - P.begin();
      R[i]--;
      P[i] = (2*double(R[i])-1) / W[i];
      if (--Rsum == N) return R;
    }
  }
}

/*!
\brief Make a grid of Q points
\details
  The points are picked according to the spacing function f(Q).
  That means that if one applies the spacing function to each point in Q grid, the resulting numbers would be (roughly) equally spaced.
  In the regular case the spacing function is f(Q)=ln(ln(Q/Lambda)).
  In case Qmin<=Lambda this spacing function cannot be used, so makeQgrid will fall back to logarithmic spacing function f(Q)=ln(Q) and issue a warning.
  The algorithm makes sure that quark mass thresholds and Z-boson mass are included as points in the returned grid, if these points are between Qmin and Qmax.
\param min Lower edge of grid, included in the grid
\param max Upper edge of grid, included in the grid
\param N total number of points in the grid
\return Array of size N, containing the Q points
*/
vector<double> makeQgrid(const double min, const double max, const size_t N){
  const static size_t MAX_POINTS = 1000000;//Protection against overflow
  if (N > MAX_POINTS) abort();//TODO make error message
  //TODO: mass thresholds are not always equal to masses: take kmuc etc into account
  //TODO: or, even better, get mass thresholds from evolution
  using XFITTER_PARS::getParamD;
  double mch = *getParamD("mch");
  double mbt = *getParamD("mbt");
  double mtp = *getParamD("mtp");
  double Mz  = *getParamD("Mz");
  vector<double> P = clipPoints(min, max, vector<double>{mch, mbt, Mz, mtp});
  //apply spacing transformation
  //Spacing function is ln(ln(Q/Lambda)
  const static double Lambda = 0.25;//GeV
  if (min > Lambda) {
    for (double&q:P) q = log( log(q) - log(Lambda) );
  }else{//use a simple logarithmic spacing
    cerr<<"[WARN] In makeQgrid: qmin<Lambda, normal spacing function cannot be used, falling back to plain log spacing; Q sampling points may be suboptimal"<<endl;
    hf_errlog(19060702, "W: makeQgrid: qmin<Lambda, falling back to log spacing");
    for (double&q:P) q = log(q);
  }
  //distribute points
  //calculate sizes of each subrange
  vector<double> W(P.size()-1);
  const size_t Wsize = W.size();
  for (size_t i=0; i<Wsize; ++i) W[i] = P[i+1] - P[i];
  if (N < P.size()) abort();//TODO: Error: not enough points
  vector<size_t> NP = apportionWebsters(W, N-P.size());
  vector<double> Q;
  Q.reserve(N);
  for (size_t i=0; i<Wsize; ++i){//for each subgrid
    //linear spread of points
    size_t endj = NP[i]+1;
    double q = P[i];
    double dq = (P[i+1]-q) / endj;
    for (size_t j=0; j<endj; ++j, q+=dq ) Q.push_back(q);
  }
  //apply reverse spacing transformation
  if (min > Lambda) {
    for (double&q:Q) q = exp( exp(q) + log(Lambda) );
  }else{//using simple logarithmic spacing
    for (double&q:Q) q = exp(q);
  }
  Q.push_back(max);
  return Q;
}

//X spacing function
const double X_SPACING_CONSTANT = 5;

double Xspacing(double x){
  return log(x) + X_SPACING_CONSTANT*x;
}

//The inverse of Xspacing
//let a=X_SPACING_CONSTANT
//Then X spacing function is y=ln(x)+ax
//Its inverse is found by solving transcendental equation exp(y)=x*exp(a*x) with respect to x
//The equation is solved numerically using Newton's method
double invXspacing(double y){
  const double EPSILON = 1e-10;
  const double a = X_SPACING_CONSTANT;
  const size_t MAX_ITERATIONS = 100;//usually 1--10 iterations is enough
  if (std::isnan(y)) return NAN;
  double x;
  if (y>0.1) x = y/a;//initial approximation
  else x = exp(y);
  y = exp(y)/a;
  for (size_t i=0; i<MAX_ITERATIONS; ++i) {
    double nx = ( x*x + y*exp(-a*x) ) / (x + 1./a);
    if (fabs(nx-x) < EPSILON) return nx;
    x = nx;
  }
  abort();//TODO proper error: too many iterations
}

/*!
\brief Make a grid of X points
\details
  The points are picked according to the spacing function f(x)=ln(x)+5*x
  That means that if one applies the spacing function to each point in X grid, the resulting numbers would be equally spaced.
\param min Lower edge of grid, included in the grid
\param max Upper edge of grid, included in the grid
\param N total number of points in the grid
\return Array of size N, containing the X points
\see Xspacing invXspacing
*/
vector<double> makeXgrid(const double xmin, const double xmax, const size_t N){
  const size_t MAX_POINTS = 1000000;//Protection against overflow
  if (N > MAX_POINTS) abort();//TODO make error message
  if (N < 2) abort();//TODO error message
  vector<double> X(N);
  double ymin = Xspacing(xmin);
  double ymax = Xspacing(xmax);
  double h = (ymax-ymin) / (N-1);
  X[0] = xmin;
  double y = ymin+h;
  for (size_t i=1; i<N-1; ++i, y+=h) X[i] = invXspacing(y);
  X[N-1] = xmax;
  return X;
}

struct Q_Subgrid{
  double* begin;//pointer to start of array of Q points
  double* end;  //pointer to just after the end of array of Q points
  int nflavors; //number of active flavors in this Q-subrange, usually between 3 (u,d,s) and 6 (u,d,s,c,b,t)
  //Offsets are applied to first and last points in the Q array
  //This is used to sample Q points just above or just below mass thresholds
  double offset_first = 0;
  double offset_last  = 0;
};

void WriteLHAPDF6subgrid(FILE* f, BaseEvolution* pdf, const vector<double>& X, const Q_Subgrid& subgrid){
  //,double*const qbegin,double*const qend,const int nflavors,double offset_first=0,double offset_last=0
  double* qbegin = subgrid.begin;
  double* qend   = subgrid.end;
  //Write x and Q points
  for (double x:X) fprintf(f, "%e ", x);
  fprintf(f, "\n");
  for (const double* q=qbegin; q!=qend; ++q) fprintf(f, "%e ", *q);
  fprintf(f, "\n");

  //        Flavor       tb bb cb sb ub db  g  d  u  s  c  b  t
  //        C++ QCDNUM    0  1  2  3  4  5  6  7  8  9 10 11 12     <-used by BaseEvolution::xfxQarray()
  const int PDG_codes[]={-6,-5,-4,-3,-2,-1,21, 1, 2, 3, 4, 5, 6};// <-used by LHAPDF
  const size_t MAX_FLAVORS = sizeof(PDG_codes) / sizeof(int);

  //TODO: handle photon

  //Determine active flavors in C++ QCDNUM convention
  const int nflavors = subgrid.nflavors;
  size_t ibegin = 6-nflavors;
  size_t iend   = 7+nflavors;

  //Write flavors
  for (size_t i=ibegin; i<iend; ++i) fprintf(f, "%i ", PDG_codes[i]);
  fprintf(f, "\n");

  //Shift first and last q-values a little bit
  //This is needed to sample PDF above or below mass threshold
  double saved_qfirst = *qbegin;
  double saved_qlast  = *(qend-1);
  *qbegin  += subgrid.offset_first;
  *(qend-1)+= subgrid.offset_last;

  //Write the PDF values
  double pdfs[MAX_FLAVORS];
  for (const double x:X) {
    for (const double* q=qbegin; q!=qend; ++q) {
      pdf->xfxQarray(x, *q, pdfs);//fills pdfs
      for (size_t i=ibegin; i<iend; ++i) fprintf(f, "%e ", pdfs[i]);
      fprintf(f, "\n");
    }
  }

  //Restore the q grid to its original state
  *qbegin   = saved_qfirst;
  *(qend-1) = saved_qlast;

  //Mark end of subgrid
  fprintf(f, "---\n");
}

//return Q_Subgrid with all points between min and max, including min and max
Q_Subgrid makeQ_Subgrid(const vector<double>& Q, double min, double max, int nflavors){
  const static double EPSILON = 1e-5;
  auto itmin = lower_bound(Q.begin(), Q.end(), min-EPSILON);
  auto itmax = lower_bound(itmin    , Q.end(), max-EPSILON);
  Q_Subgrid ret;
  ret.begin = const_cast<double*>(&*itmin);
  ret.end   = const_cast<double*>(&*(itmax+1));
  ret.nflavors = nflavors;
  return ret;
}

vector<Q_Subgrid> makeQ_Subgrids(const vector<double>& Q){
  const static double EPSILON = 1e-5;
  using XFITTER_PARS::getParamD;
  double mch = *getParamD("mch");
  double mbt = *getParamD("mbt");
  double mtp = *getParamD("mtp");
  const static size_t NTHRESHOLDS = 3;
  double thresholds[NTHRESHOLDS]={mch,mbt,mtp};
  const size_t N = Q.size();
  const double min = Q[0];
  const double max = Q[N-1];
  vector<double> subrangeEdges = {min};
  size_t i = 0;
  for (; i<NTHRESHOLDS; ++i) {
    if (thresholds[i] > min) break;
  }
  int first_nf = 3 + int(i);
  for (; i<NTHRESHOLDS; ++i){
    if (thresholds[i] >= max) break;
    subrangeEdges.push_back(thresholds[i]);
  }
  subrangeEdges.push_back(max);
  size_t endi = subrangeEdges.size()-1;
  int nflavors = first_nf + endi - 1;//Even though different subranges have different number of flavors, lhapdf requires flavors to be the same for all subgrids
  vector<Q_Subgrid> R;
  for (size_t i = 0; i<endi; ++i) {
    double min = subrangeEdges[i];
    double max = subrangeEdges[i+1];
    Q_Subgrid subgrid = makeQ_Subgrid(Q, min, max, nflavors);
    if (i != 0)      subgrid.offset_first = +EPSILON;
    if (i != endi-1) subgrid.offset_last  = -EPSILON;
    if (subgrid.begin != subgrid.end) R.push_back(subgrid);
  }
  return R;
}

struct QX_Grid{
  vector<double> X;
  vector<double> Q;
  vector<Q_Subgrid> subgrids;
};

QX_Grid makeQX_Grid(double xmin, double xmax, size_t Nx, double qmin, double qmax, size_t Nq) {
  QX_Grid ret;
  ret.X = makeXgrid(xmin, xmax, Nx);
  ret.Q = makeQgrid(qmin, qmax, Nq);
  ret.subgrids = makeQ_Subgrids(ret.Q);
  return ret;
}

void WriteLHAPDF6(FILE* f, BaseEvolution* ev, const QX_Grid& qx_grid, const char*const PdfType = "central"){
  fprintf(f, "PdfType: %s\n", PdfType);
  fprintf(f, "Format: lhagrid1\n---\n");
  for (const Q_Subgrid& subgrid:qx_grid.subgrids) {
    WriteLHAPDF6subgrid(f, ev, qx_grid.X, subgrid);
  }
}

struct LHAPDF6_Options{
  BaseEvolution* pdf;
  string
    name = "xfitter_pdf",
    description = "Generated using xFitter",
    authors = "",
    reference = "",
    flavor_scheme = "",
    error_type = "";
  double qmin, qmax, xmin, xmax;
  size_t Nx, Nq;
  int nmembers = 1;
  static LHAPDF6_Options fromYAML(YAML::Node);
};

//Returns "hessian" or "symmhessian" for PDF error type
const char* getErrorType(){
  //PLACEHOLDER: I am not sure right now how to get error type in the general case
  //for example when CERES is used instead of MINUIT
  YAML::Node minuitNode = XFITTER_PARS::rootNode["MINUIT"];
  string doErrors;
  YAML::Node doErrorsNode;
  if ( !minuitNode.IsMap() ) goto failed;
  doErrorsNode = minuitNode["doErrors"];
  doErrors = doErrorsNode.as<string>("");
  if ( !doErrorsNode.IsScalar() ) goto failed;
  if      (doErrors == "Hesse"  ) return "symmhessian";
  else if (doErrors == "Pumplin") return "hessian";
  else goto failed;
  failed:
  hf_errlog(19052700, "W: LHAPDF6 output: failed to determine error type, assuming symmhessian");
  return "symmhessian";
}

const char* getFlavorScheme(){
  int nf = XFITTER_PARS::getParamI("NFlavour");
  //PLACEHOLDER: I am not sure right now how to get flavor scheme in the general case
  //Maybe we should add BaseEvolution::getFlavorScheme ?
  if (nf<5) return "fixed";
  else      return "variable";
}

size_t getNmembers(){
  //TODO this needs to be more general
  YAML::Node minuitNode = XFITTER_PARS::rootNode["MINUIT"];
  if ( !minuitNode.IsMap() ) return 1;
  YAML::Node doErrorsNode = minuitNode["doErrors"];
  string doErrors;
  try{
    doErrors = doErrorsNode.as<string>();
  }catch(YAML::TypedBadConversion<string>){
    return 1;
  }
  size_t Npars=xfitter::get_minimizer()->getNpars();
  if (doErrors == "Hesse")return Npars+1;
  if (doErrors == "Pumplin")return 2*Npars+1;
  return 1;
}

//Parse control block in YAML steering and fill LHAPDF6_Options
LHAPDF6_Options LHAPDF6_Options::fromYAML(YAML::Node node){
  LHAPDF6_Options info;
  BaseEvolution*pdf = xfitter::get_evolution(node["evolution"].as<string>(""));
  info.pdf = pdf;
  {
  YAML::Node n = node["name"];
  if (n.IsDefined()) {
    try{
      info.name = n.as<string>();
    }catch(YAML::TypedBadConversion<string>()){
      cerr<<"[ERROR] WriteLHAPDF6 failed to convert name to string when trying to output PDF \""<<pdf->_name<<"\""<<endl;
      hf_errlog(19060700, "F: Error when parsing WriteLHAPDF6, see stderr");
      abort();
    }
  }
  }
  {
  YAML::Node n = node["description"];
  if (n.IsDefined()) {
    try{
      info.description = n.as<string>();
    }catch(YAML::TypedBadConversion<string>()){
      cerr<<"[ERROR] WriteLHAPDF6 failed to convert description to string when trying to output PDF \""<<pdf->_name<<"\""<<endl;
      hf_errlog(19060700, "F: Error when parsing WriteLHAPDF6, see stderr");
      abort();
    }
  }
  }
  {
  YAML::Node n = node["authors"];
  if (n.IsDefined()) {
    try{
      info.authors = n.as<string>();
    }catch(YAML::TypedBadConversion<string>()){
      cerr<<"[ERROR] WriteLHAPDF6 failed to convert authors to string when trying to output PDF \""<<pdf->_name<<"\""<<endl;
      hf_errlog(19060700, "F: Error when parsing WriteLHAPDF6, see stderr");
      abort();
    }
  }
  }
  {
  YAML::Node n = node["reference"];
  if(n.IsDefined()){
    try{
      info.reference = n.as<string>();
    }catch(YAML::TypedBadConversion<string>()){
      cerr<<"[ERROR] WriteLHAPDF6 failed to convert reference to string when trying to output PDF \""<<pdf->_name<<"\""<<endl;
      hf_errlog(19060700, "F: Error when parsing WriteLHAPDF6, see stderr");
      abort();
    }
  }
  }
  info.error_type = getErrorType();
  info.flavor_scheme = getFlavorScheme();
  {
  YAML::Node n = node["Xrange"];
  if (n.IsDefined()) {
    if (!n.IsSequence() or n.size()!=2) {
      cerr<<"[ERROR] WriteLHAPDF6: Xrange must be given as [min, max]; error when trying to output PDF \""<<pdf->_name<<"\""<<endl;
      hf_errlog(19060700, "F: Error when parsing WriteLHAPDF6, see stderr");
      abort();
    }
    try{
      info.xmin = n[0].as<double>();
      info.xmax = n[1].as<double>();
    }catch(YAML::Exception){
      cerr<<"[ERROR] WriteLHAPDF6: Failed to interpret Xrange for PDF \""<<pdf->_name<<"\""<<endl;
      hf_errlog(19060700, "F: Error when parsing WriteLHAPDF6, see stderr");
      abort();
    }
    if ( not ( info.xmin < info.xmax ) ) {
      cerr<<"[ERROR] WriteLHAPDF6: In Xrange=[min, max] min must be smaller than max; for PDF \""<<pdf->_name<<"\""<<endl;
      hf_errlog(19060700, "F: Error when parsing WriteLHAPDF6, see stderr");
      abort();
    }
    YAML::Node nodeNx = node["Xnpoints"];
    if (nodeNx.IsDefined()) {
      int Nx;
      try{
        Nx = nodeNx.as<int>();
      }catch(YAML::TypedBadConversion<int>){
        cerr<<"[ERROR] WriteLHAPDF6: Failed to interpret Xnpoints for PDF \""<<pdf->_name<<"\", expected number of points"<<endl;
        hf_errlog(19060700, "F: Error when parsing WriteLHAPDF6, see stderr");
        abort();
      }
      if ( Nx <= 0) {
        cerr<<"[ERROR] WriteLHAPDF6: Nonpositive number of X points for PDF \""<<pdf->_name<<"\""<<endl;
        hf_errlog(19060700, "F: Error when parsing WriteLHAPDF6, see stderr");
        abort();
      }
      info.Nx = Nx;
    }else{
      cerr<<"[WARN] WriteLHAPDF6: for PDF \""<<pdf->_name<<"\" Xrange is provided, but not Xnpoints; using Xnpoints=200"<<endl;
      hf_errlog(19060701, "W: Xnpoints not provided, using Xnpoints=200, see stderr");
      info.Nx = 200;
    }
  }else{//Use default X parameters
    info.xmin = 1e-6;
    info.xmax = 1;
    info.Nx = 200;
  }
  }
  {
  YAML::Node n = node["Qrange"];
  if(n.IsDefined()){
    if (!n.IsSequence() or n.size()!=2) {
      cerr<<"[ERROR] WriteLHAPDF6: Qrange must be given as [min, max]; error when trying to output PDF \""<<pdf->_name<<"\""<<endl;
      hf_errlog(19060700, "F: Error when parsing WriteLHAPDF6, see stderr");
      abort();
    }
    try{
      info.qmin = n[0].as<double>();
      info.qmax = n[1].as<double>();
    }catch(YAML::Exception){
      cerr<<"[ERROR] WriteLHAPDF6: Failed to interpret Qrange for PDF \""<<pdf->_name<<"\""<<endl;
      hf_errlog(19060700, "F: Error when parsing WriteLHAPDF6, see stderr");
      abort();
    }
    if ( not ( info.qmin < info.qmax ) ) {
      cerr<<"[ERROR] WriteLHAPDF6: In Qrange=[min, max] min must be smaller than max; for PDF \""<<pdf->_name<<"\""<<endl;
      hf_errlog(19060700, "F: Error when parsing WriteLHAPDF6, see stderr");
      abort();
    }
    YAML::Node nodeNq = node["Qnpoints"];
    if (n.IsDefined()) {
      int Nq;
      try{
        Nq = nodeNq.as<int>();
      }catch(YAML::TypedBadConversion<int>){
        cerr<<"[ERROR] WriteLHAPDF6: Failed to interpret Qnpoints for PDF \""<<pdf->_name<<"\", expected number of points"<<endl;
        hf_errlog(19060700, "F: Error when parsing WriteLHAPDF6, see stderr");
        abort();
      }
      if ( Nq <= 0 ) {
        cerr<<"[ERROR] WriteLHAPDF6: Nonpositive number of Q points for PDF \""<<pdf->_name<<"\""<<endl;
        hf_errlog(19060700, "F: Error when parsing WriteLHAPDF6, see stderr");
        abort();
      }
      info.Nq = Nq;
    }else{
      info.Nq = 120;
      cerr<<"[WARN] WriteLHAPDF6: for PDF \""<<pdf->_name<<"\" Qrange is provided, but not Qnpoints; using Qnpoints="<<info.Nq<<endl;
      hf_errlog(19060703, "W: Qnpoints not provided, using Qnpoints=120, see stderr");
    }
  }else{//Use default Q parameters
    info.qmin = 1;
    info.qmax = 1e4;
    info.Nq = 120;
  }
  }
  return info;
}

void WriteLHAPDF6info(FILE* f, const LHAPDF6_Options& info, const QX_Grid& qx_grid){
  if (!info.description.empty()) fprintf(f, "SetDesc: %s\n"  , info.description.c_str());
  if (!info.authors.empty())     fprintf(f, "Authors: %s\n"  , info.authors.c_str());
  if (!info.reference.empty())   fprintf(f, "Reference: %s\n", info.reference.c_str());
  fprintf(f, "Format: lhagrid1\n");
  fprintf(f, "DataVersion: 1\n");
  fprintf(f, "NumMembers: %i\n", info.nmembers);

  int nflavors = qx_grid.subgrids.back().nflavors;
  fprintf(f, "Flavors: [");
  for (int i=-nflavors; i<0; ++i) fprintf(f, "%i, ", i);
  fprintf(f, "21, ");
  for (int i=1; i<nflavors; ++i) fprintf(f, "%i, ", i);
  fprintf(f, "%i]\n", nflavors);
  //TODO: photon

  int order = OrderMap(XFITTER_PARS::getParamS("Order"));
  fprintf(f, "OrderQCD: %i\n", order); // qcdnum notation LO=1, ...; LHAPDF6 LO=0, ...
  if (!info.flavor_scheme.empty()) fprintf(f, "FlavorScheme: %s\n", info.flavor_scheme.c_str());
  if (!info.error_type.empty())    fprintf(f, "ErrorType: %s\n"   , info.error_type.c_str());

  fprintf(f, "XMin: %g\n", qx_grid.X.front());
  fprintf(f, "XMax: %g\n", qx_grid.X.back());

  fprintf(f, "QMin: %g\n", qx_grid.Q.front());
  fprintf(f, "QMax: %g\n", qx_grid.Q.back());

  //TODO: different quark and MZ masses for different evolutions
  //TODO: get masses from evolutions
  using XFITTER_PARS::getParamD;
  double Mz = *getParamD("Mz");
  fprintf(f, "MZ: %g\n", Mz);
  fprintf(f, "MUp: 0\n");
  fprintf(f, "MDown: 0\n");
  fprintf(f, "MStrange: 0\n");
  fprintf(f, "MCharm: %g\n", *getParamD("mch"));
  fprintf(f, "MBottom: %g\n", *getParamD("mbt"));
  fprintf(f, "MTop: %g\n", *getParamD("mtp"));

  BaseEvolution* pdf = info.pdf;
  //Write alphaS
  fprintf(f, "AlphaS_MZ: %g\n", pdf->getAlphaS(Mz));
  fprintf(f, "AlphaS_OrderQCD: %i\n", order);//How is this different from OrderQCD? --Ivan
  fprintf(f, "AlphaS_Type: ipol\n");//I have no idea what that is --Ivan

  //Tabulate alphaS
  //Note that mass thresholds  are printed twice: just below and just above the threshold
  fprintf(f, "AlphaS_Qs: [");
  bool first = true;
  for (const Q_Subgrid& subgrid: qx_grid.subgrids) {
    double* qend = subgrid.end;
    for (double* q=subgrid.begin; q!=qend; ++q) {
      if (first) {
        fprintf(f, "%g", *q);
        first = false;
      }else{
        fprintf(f, ", %g", *q);
      }
    }
  }
  fprintf(f, "]\nAlphaS_Vals: [");
  first = true;
  for (const Q_Subgrid& subgrid: qx_grid.subgrids) {
    double* qbegin = subgrid.begin;
    double* qend = subgrid.end;
    for(double*qp=subgrid.begin;qp!=qend;++qp){
      double q = *qp;
      if      (qp == qbegin) q += subgrid.offset_first;
      else if (qp == qend-1) q += subgrid.offset_last;
      double alpha = pdf->getAlphaS(q);
      if (first) {
        fprintf(f, "%g", alpha);
        first = false;
      } else {
        fprintf(f, ", %g", alpha);
      }
    }
  }
  fprintf(f, "]\n");
}

map<string,QX_Grid> cached_QXgrids;

bool directoryExists(const string& path){
  struct stat info;
  if (stat(path.c_str(), &info)!=0) return false;
  return bool(info.st_mode & S_IFDIR);
}

//FORTRAN INTERFACE
extern "C"{

void save_data_lhapdf6_(const int& memberID){//This is called when building bands
  YAML::Node lhapdf6node=XFITTER_PARS::rootNode["WriteLHAPDF6"];
  if (not lhapdf6node) return;
  //TODO: output multiple evolutions
  //TODO: handle errors here
  LHAPDF6_Options options = LHAPDF6_Options::fromYAML(lhapdf6node);
  options.nmembers = getNmembers();
  const string& name = options.name;
  const QX_Grid* qx_grid;
  const auto it = cached_QXgrids.find(name);
  if (it != cached_QXgrids.end()) {
    qx_grid = &(it->second);
  }else{
    //TODO: If options.pdf->getClassName()=="QCDNUM" use QCDNUM's grid to avoid double interpolation
    cached_QXgrids[name] = makeQX_Grid(options.xmin, options.xmax, options.Nx, options.qmin, options.qmax, options.Nq);
    qx_grid = &(cached_QXgrids.at(name));
  }
  string outdir = xfitter::getOutDirName() + '/' + name;
  if (not directoryExists(outdir)){
    if (mkdir(outdir.c_str(),0755) != 0) {
      cerr<<"[ERROR] Failed to create directory \""<<outdir<<"\" for lhapdf6 output"<<endl;
      perror("mkdir error");
      hf_errlog(19060702, "F: Failed to create directory for lhapdf6 output, see stderr");
      abort();
    }
  }
  //Make file name
  string filename;
  filename.reserve( outdir.size()+name.size()+10 );
  sprintf(const_cast<char*>(filename.c_str()), "%s/%s_%04i.dat", outdir.c_str(), name.c_str(), memberID);

  FILE* f = fopen(filename.c_str(),"w");
  if (f==nullptr) {
    cerr<<"[ERROR] Failed to open file \""<<filename<<"\" for lhapdf6 output"<<endl;
    perror("fopen error");
    hf_errlog(19060800, "F: Failed to open file for lhapdf6 output, see stderr");
    abort();
  }
  const char* PdfType = "central";
  if (memberID!=0) PdfType = "error";
  WriteLHAPDF6(f, options.pdf, *qx_grid, PdfType);
  fclose(f);

  if (memberID==0) {//then write info
    filename = outdir+'/'+name+".info";
    f = fopen(filename.c_str(),"w");
    if (f==nullptr) {
      cerr<<"[ERROR] Failed to open file \""<<filename<<"\" for lhapdf6 output"<<endl;
      perror("fopen error");
      hf_errlog(19060800, "F: Failed to open file for lhapdf6 output, see stderr");
      abort();
    }
    WriteLHAPDF6info(f, options, *qx_grid);
    fclose(f);
  }
}

void print_lhapdf6_(){
  save_data_lhapdf6_(0);
}

//TODO: replace all calls to this function with calls to print_lhapdf6 and delete this one
void print_lhapdf6_opt_(){
  save_data_lhapdf6_(0);
}

}
