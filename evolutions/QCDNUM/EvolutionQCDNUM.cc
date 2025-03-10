
/*
   @file EvolutionQCDNUM.cc
   @date 2018-08-14
   @author  AddEvolution.py
   Created by  AddEvolution.py on 2018-08-14
*/

#include "EvolutionQCDNUM.h"
#include "QCDNUM/QCDNUM.h"
#include "QCDNUM_Manager.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"
#include "xfitter_steer.h"
#include"BasePdfDecomposition.h"
#include <algorithm>
#include <cmath>

// Global var to hold current pdfDecomposition
xfitter::BasePdfDecomposition* gPdfDecomp;
xfitter::BaseEvolution* gExternalEvolution;

// Wrapper for QCDNUM
double funcPDF(int *ipdf, double *x) {
  if (*ipdf<0) return 0.;
  const std::map <int,int> ip = { {1,-3}, {2,-2}, {3,-1}, {4,1}, {5,2}, {6,3}, {7,-4}, {8,4}, {9,-5}, {10,5}, {11,-6}, {12,6}, {0,21} } ;
  return gPdfDecomp->xfxMap(*x)[ip.at(*ipdf)];
}

// Wrapper for QCDNUM
double funcPDFext(int* ipdf, double* x,  double* q2, bool* first){
  const std::map <int,int> ip =
  {
    {-6,-6}, {-5,-5}, {-4,-4}, {-3,-3}, {-2,-2}, {-1,-1}, {1,1}, {2,2}, {3,3}, {4,4}, {5,5}, {6,6},
    {0,21},{7,22}
  };
  double q = sqrt(*q2);
  // avoid crash if copying external PDFs which do not have photon
  try {
    return gExternalEvolution->xfxQmap(*x, q)[ip.at(*ipdf)];
  }
  catch(std::out_of_range&ex) {
    if (*ipdf == 7) {
      return 0.;
    }
    return 0.; // avoid warning
  }
}


// helper to parse yaml sequences of uniform type
template <class T>
vector<T> getSeq(const YAML::Node node) {
  if(!node.IsSequence()){
    std::cerr<<"[ERROR]getSeq: node=\n"<<node<<std::endl;
    hf_errlog(180829150,"F: In QCDNUM in function getSeq: wrong node type, expected sequence");
  }
  size_t len = node.size();
  vector<T> v(len);
  for (size_t i=0; i<len; i++) {
    v[i] = node[i].as<T>();
  }
  return v;
}

//
double static   qcdnumDef[] = {
// tb bb  cb  sb ub   db  g   d   u   s   c   b   t
  0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  // sb
  0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,  // ub
  0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,  // db
  0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,  // d
  0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,  // u
  0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,  // s
  0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  // cb
  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,  // c
  0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  // bb
  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,  // b
  1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  // tb
  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.  // t
};


namespace xfitter {

// the class factories
  extern "C" EvolutionQCDNUM* create(const char*name) {
    return new EvolutionQCDNUM(name);
  }
  const char*EvolutionQCDNUM::getClassName()const{return "QCDNUM";}

  // Initialize at the start of the computation
  void EvolutionQCDNUM::atStart()
  {
    initQCDNUM();

    // check version:
    int ver;
    QCDNUM::getint("vers",ver);
    if (ver < 170114) {
      hf_errlog(231020151, "F: Obsolete version of QCDNUM. Install version 17.01.14 or later.");
    }

    // check memory limits for different pars:
    std::string memT[] = {"mxg0","mxx0","mqq0","mst0","mce0","mbf0","mky0","nwf0"};
    int         siz[]  = { 5,     300,    150,   30,   20,     10,   50,    1200000};
    for (int i=0; i<8; i++) {
      int itest;
      QCDNUM::getint(memT[i].c_str(),itest);

      if (itest < siz[i]) {
        std::string mess = "F: QCDNUM memory allocation insufficient. Recompile QCDNUM with at least "+memT[i] + " = " + std::to_string(siz[i]) + " in qcdnum.inc.";
        hf_errlog(231020152,mess.c_str());
      }
    }

    // Evolution order:
    const int     PtOrder    = OrderMap(XFITTER_PARS::getParamS("Order")) ;

    const double* Q0         = XFITTER_PARS::getParamD("Q0");
    double q20 = (*Q0) * (*Q0);

    Mz = XFITTER_PARS::getParamD("Mz");
    alphas = XFITTER_PARS::getParamD("alphas");

    const double* mch     = XFITTER_PARS::getParamD("mch");
    const double* mbt     = XFITTER_PARS::getParamD("mbt");

    //  const double* mtp     = XFITTER_PARS::getParamD("mtp");   // no top PDF treatment yet XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    YAML::Node yQCDNUM=XFITTER_PARS::getEvolutionNode(_name);
    _icheck      = yQCDNUM["ICheck"].as<int>();
    _splineOrder = yQCDNUM["SplineOrder"].as<int>();
    _readTables  = yQCDNUM["Read_QCDNUM_Tables"].as<int>();

    // Get grids
    vector<double> xGrid   = getSeq<double>(yQCDNUM["xGrid"]);
    vector<int>    xGridW  = getSeq<int>   (yQCDNUM["xGridW"]);
    vector<double> Q2Grid  = getSeq<double>(yQCDNUM["Q2Grid"]);
    vector<double> Q2GridW = getSeq<double>(yQCDNUM["Q2GridW"]);

    int nxgrid             = yQCDNUM["NXbins"].as<int>();
    int nxout = 0;

    QCDNUM::setord(PtOrder);
    std::cout << "Set evolution order "<<PtOrder <<"\n";
    std::cout << "xGrid[0]: " << xGrid[0] << std::endl;
    QCDNUM::gxmake(xGrid.data(),xGridW.data(),xGrid.size(),nxgrid,nxout,_splineOrder); // x-grid
    std::cout << "Requested (actual) number of x  grid points: "<<nxgrid << "(" << nxout << ")\n";

    // for Q2 grid, add matching points and starting scale
    Q2Grid.push_back(q20);
    Q2GridW.push_back(4.0);

    Q2Grid.push_back( (*mch)*(*mch) );
    Q2GridW.push_back(2.0);

    Q2Grid.push_back( (*mbt)*(*mbt) );
    Q2GridW.push_back(1.3);

    Q2Grid.push_back( (*Mz) * (*Mz) );
    Q2GridW.push_back(1.1);

    //    for (size_t i = 0; i<Q2Grid.size(); i++) {
    //  std::cout << " Q2= " << Q2Grid[i] << " weight = " << Q2GridW[i] << "\n";
    //}
    // sort
    std::vector < pair<double,double> > tmp;
    for ( size_t i = 0; i<Q2Grid.size(); i++) {
      tmp.push_back({Q2Grid[i],Q2GridW[i]});
    }
    std::sort(tmp.begin(),tmp.end(), []( pair<double,double> a, pair<double,double> b) {return a.first<b.first;} );

    for ( size_t i = 0; i<Q2Grid.size(); i++) {
      Q2Grid[i]  = tmp[i].first;
      Q2GridW[i] = tmp[i].second;
    }

    //for (size_t i = 0; i<Q2Grid.size(); i++) {
    //  std::cout << " Q2= " << Q2Grid[i] << " weight = " << Q2GridW[i] << "\n";
    //}

    int nq2grid             = yQCDNUM["NQ2bins"].as<int>();
    int nq2out = 0;

    QCDNUM::gqmake(Q2Grid.data(), Q2GridW.data(), Q2Grid.size(), nq2grid, nq2out);
    std::cout << "Requested (actual) number of Q2 grid points: "<<nq2grid << "(" << nq2out << ")\n";

    // set VFNS thresholds
    int iqc = QCDNUM::iqfrmq( (*mch)*(*mch) + 1.e-6 );
    int iqb = QCDNUM::iqfrmq( (*mbt)*(*mbt) + 1.e-6 );
    int iqt = 0;  // top off for now

    // For now VFNS only and NFlavour = 3 only
    int nflavour = XFITTER_PARS::gParametersI.at("NFlavour");
    if(nflavour < 3 || nflavour > 6)
      hf_errlog(110520201, "F: Unsupported NFlavour = " + std::to_string(nflavour));
    int isFFNS = 0; // VFNS by default
    if(XFITTER_PARS::gParametersI.find("isFFNS") != XFITTER_PARS::gParametersI.end())
      isFFNS = XFITTER_PARS::gParametersI.at("isFFNS");
    if(isFFNS == 1)
    {
      std::cout << "Fixed Flavour Number Scheme set with nf=" << nflavour << std::endl;
      QCDNUM::setcbt(nflavour,0,0,0);
    }
    else if(isFFNS == 0)
    {
      std::cout << "Variable Flavour Number Scheme set with nf=" << nflavour << std::endl;
      if(nflavour == 3)
        QCDNUM::setcbt(0,0,0,0);
      else if(nflavour == 4)
        QCDNUM::setcbt(0,iqc,0,0);
      else if(nflavour == 5)
        QCDNUM::setcbt(0,iqc,iqb,0);
      else if(nflavour == 6)
        QCDNUM::setcbt(0,iqc,iqb,iqt);
    }
    else
      hf_errlog(110520202, "F: Unsupported isFFNS = " + std::to_string(isFFNS));

    // Init SF
    int id1=0;      int id2=0;      int nw=0;      int ierr=1;
    if (_readTables>0 ) {
      QCDNUM::readwt(22,"unpolarised.wgt",id1,id2,nw,ierr);
    }
    // Fill the tables if did not read correctly
    if (ierr != 0) {
      QCDNUM::fillwt(0,id1,id2,nw);
      QCDNUM::dmpwgt(1,22,"unpolarised.wgt");
    }

    //Special case for QCDNUM, allow for external PDFs
    if (yQCDNUM["EvolutionCopy"]) {
      auto evolCopyName = yQCDNUM["EvolutionCopy"].as<string>();
      gExternalEvolution = xfitter::get_evolution(evolCopyName);
      _itype = 5; // external
    }
    else {
    //Evolution gets its decomposition from YAML
      gPdfDecomp=XFITTER_PARS::getInputDecomposition(yQCDNUM);
    }
  }

  void EvolutionQCDNUM::atConfigurationChange(){}

  void EvolutionQCDNUM::atIteration(){
    if (_itype <5 ) {
      const double* q0 = XFITTER_PARS::getParamD("Q0");
      int iq0  = QCDNUM::iqfrmq( (*q0) * (*q0) );
      double epsi = 0;

      //Re-read Mz and alphas at each iteration so that they can be fitted
      double alS = *alphas;
      if ( std::isnan(alS) ) {
        alS = 0.0001; //this should make chi2 bad and force minimizer to an earlier, valid, value of alphas
        cerr<<"[WARN] QCDNUM got alphas = NaN; using alphas = "<<alS<<" instead"<<endl;
        hf_errlog(19070500, "W: alphas = NaN in QCDNUM, see stderr");
      }

      const double MZ = *Mz;
      QCDNUM::setalf( alS, MZ * MZ );

      QCDNUM::evolfg(_itype,funcPDF,qcdnumDef,iq0,epsi);
    }
  }

  void EvolutionQCDNUM::afterIteration(){
    if (_itype != 5) return;
    // External evolution
    double epsi = 0;
    QCDNUM::extpdf(funcPDFext,_itype,_nExt,0.001,epsi);
  }

  std::map<int,double>EvolutionQCDNUM::xfxQmap(double x,double Q){
    std::map<int, double> res;
    for (int ipdf =-6; ipdf<7; ipdf++) {
      int ii = ( ipdf == 0 ) ? 21 : ipdf ;
      res[ii] = QCDNUM::fvalxq(_itype,ipdf,x,Q*Q,_icheck);
    }
    return res;
  }

  void EvolutionQCDNUM::xfxQarray(double x,double Q,double*pdfs){
    QCDNUM::allfxq(_itype,x,Q*Q,pdfs,_nExt,_icheck);
  }

  double EvolutionQCDNUM::xfxQ(int i,double x,double Q){
    return QCDNUM::fvalxq(_itype,i,x,Q*Q,_icheck);
  }

  double EvolutionQCDNUM::getAlphaS(double Q){
    if (_itype == 5 ) { //if external evolution is used
      return gExternalEvolution->getAlphaS(Q);
    } else {
      int nfout,ierr;
      return QCDNUM::asfunc(Q*Q,nfout,ierr);
    }
  }

vector<double> EvolutionQCDNUM::getXgrid(){
  int Nx;
  {// get Nx
  int i_ignore;
  double d_ignore;
  QCDNUM::grpars(Nx,d_ignore,d_ignore,i_ignore,d_ignore,d_ignore,i_ignore);
  }
  vector<double> xgrid(Nx);
  for (int i=0;i<Nx;++i) xgrid[i] = QCDNUM::xfrmix(i+1);
  return xgrid;
}

vector<double> EvolutionQCDNUM::getQgrid(){
  int Nq;
  {// get Nq
  int i_ignore;
  double d_ignore;
  QCDNUM::grpars(i_ignore,d_ignore,d_ignore,Nq,d_ignore,d_ignore,i_ignore);
  }
  vector<double> qgrid(Nq);
  for (int i=0;i<Nq;++i) qgrid[i] = sqrt( QCDNUM::qfrmiq(i+1) );
  return qgrid;
}
}
