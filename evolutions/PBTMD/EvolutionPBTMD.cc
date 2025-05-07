
/*
   @file EvolutionPBTMD.cc
   @date 2023-07-22
   @author  H. Jung
*/

#include "EvolutionPBTMD.h"
#include "QCDNUM/QCDNUM.h"
#include "QCDNUM_Manager.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"
#include "xfitter_steer.h"
#include"BasePdfDecomposition.h"
#include <algorithm>
#include <cmath>
#include <stdio.h>
#include <string.h>

// Global var to hold current pdfDecomposition

xfitter::BasePdfDecomposition *gPdfDecomp = nullptr;
xfitter::BaseEvolution* gExternalEvolution;

// Wrapper for QCDNUM
double funcPDF(int *ipdf, double *x) {
  if (*ipdf<0) return 0.;
  const std::map <int,int> ip = { {1,-3}, {2,-2}, {3,-1}, {4,1}, {5,2}, {6,3}, {7,-4}, {8,4}, {9,-5}, {10,5}, {11,-6}, {12,6}, {0,21} } ;
  return gPdfDecomp->xfxMap(*x)[ip.at(*ipdf)];
}

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


extern "C" {



double funcstart_(int *ipdf, double *x) {
//double funcPDF(int *ipdf, double *x) {
  if (gPdfDecomp == nullptr) {
    std::cout<< "PBTMD PDF decomposition is not set, cannot compute PDFs" << std::endl;
  }

    /*Parton codes are LHAPDF convention:
     *      i  -6 -5 -4 -3 -2 -1 21 1  2  3  4  5  6  22 11 13 15 -11 -13 -15
     * pdfs[i] tb bb cb sb ub db gl d  u  s  c  b  t  ga e  mu tau ae amu atau
     */

    // std::map<int,double>xfxMap(double x)const=0;

  // if (*ipdf<0) return 0.;
  //const std::map <int,int> ip = { {1,-3}, {2,-2}, {3,-1}, {4,1}, {5,2}, {6,3}, {7,-4}, {8,4}, {9,-5}, {10,5}, {11,-6}, {12,6}, {0,21} } ;
  const std::map <int,int> ip =
  {
    {-6,-6}, {-5,-5}, {-4,-4}, {-3,-3}, {-2,-2}, {-1,-1}, {1,1}, {2,2}, {3,3}, {4,4}, {5,5}, {6,6},
    {0,21},{7,22}, {8,23}, {9,24}, {10,-24}, {11,25}
  };  
  // std::cout << " in PBTMD in funcSTART " << *ipdf << " " << *x << std::endl;
  //std::cout << " in PBTMD in funcSTART " << gPdfDecomp->xfxMap(*x)[ip.at(*ipdf)]<< std::endl;
  return gPdfDecomp->xfxMap(*x)[ip.at(*ipdf)];
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


namespace xfitter {

// the class factories
  extern "C" EvolutionPBTMD* create(const char*name) 
  {
    return new EvolutionPBTMD(name);
  }
  
  EvolutionPBTMD::EvolutionPBTMD(const char* name): BaseEvolution(name)
  {
   // _pdf = nullptr;
  }  
  const char*EvolutionPBTMD::getClassName()const{return "PBTMD";}


/// Global initialization
  void EvolutionPBTMD::atStart()
  {
  std::cout<< " using PBTMD evolution  " << _name << std::endl;
    YAML::Node yPBTMD=XFITTER_PARS::getEvolutionNode(_name);

    //Evolution gets its decomposition from YAML
    gPdfDecomp=XFITTER_PARS::getInputDecomposition(yPBTMD);
    Mz = XFITTER_PARS::getParamD("Mz");
    alphas = XFITTER_PARS::getParamD("alphas");
    const double* Q0         = XFITTER_PARS::getParamD("Q0");
    double q20 = (*Q0) * (*Q0);
   }

  void EvolutionPBTMD::atConfigurationChange(){
    std::cout<< " in PBTMD atConfigurationChange " << std::endl;

  }

  void EvolutionPBTMD::atIteration(){
//      std::cout<< " in PBTMD atIteration " << std::endl;
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
//        cout << " PBTMD at iteration: as = " << alS << " " << MZ << endl;
      QCDNUM::setalf( alS, MZ * MZ );

      QCDNUM::evolfg(_itype,funcPDF,qcdnumDef,iq0,epsi);

  }

  void EvolutionPBTMD::afterIteration(){
  // std::cout<< " in PBTMD afterIteration " << std::endl;
  }

  /// Return PDFs as a map <int,double> where int is PDF ID (-6, ... 6, 21)
  std::map<int,double>EvolutionPBTMD::xfxQmap(double x,double Q)
  {
     //std::map<int, double> res;
     double q2 = Q*Q ;
     //double test(-6,11);
     double* mc   = XFITTER_PARS::getParamD("mch");
     double* mb   = XFITTER_PARS::getParamD("mbt");
     double* mt   = XFITTER_PARS::getParamD("mtp");
     // cout << " in PBTMD xfxQmap: mc = " << *mc << " mb = " << *mb << endl;
     
     //std::cout<< " in PBTMD xfxQmap " << x << " " << q2 << std::endl;
     int npdfMax = 11 ;
//     int npdfMax = 9 ;
     
     string name = "TMDgrids" ;
     //char* file = name.data();
     // do this to be on safe side with linux and Mac
     char *file = new char[name.length() + 1];
     strcpy(file, name.c_str());
     //
     
     //const char* file = name.c_str();
     int length = strlen(file);
     
     //cout << " callin xfxQmap from PBTMD "<< endl;
     // convert s to dest stringToFortran(dest, maxlen,s)
     //stringToFortran(file, 132, name);
     // std::cout << " file = "<< file << " length " << length<< std::endl;
     
   // if (  x != xSav || q2 != Q2Sav) {  
        for (int ipdf = -6; ipdf <= npdfMax; ipdf++){
            int ii = (ipdf == 0) ? 21 : ipdf;
            // photon PDF:
            if (ipdf == 7) ii = 22;
            if (ipdf == 8) ii = 23;
            if (ipdf == 9) ii = 24;
            if (ipdf == 10) ii = -24;
            if (ipdf == 11) ii = 25;
            // std::cout<< " ipdf = " << ipdf << " ii = " << ii << " x = " << x << " q2 = " << q2 << std::endl;
            res[ii] = pbtmdsubr(ipdf, x, q2, mc, mb, file, length);
            
        }
        xSav = x; Q2Sav = q2;        
//        cout << " xfxQmap [22,23} = " << res[22] << " " << res[23] << endl;
    // }
    
  return res;
  }

  /// Returns PDFs as a function of i, x, Q
  double EvolutionPBTMD::xfxQ(int i, double x, double Q)
  {
    std::cout<< " in PBTMD xfxQ i=  " << i << std::endl;
    double pdf=0; 
   // return _pdf->xfxQ(i, x, Q);
   return pdf;
  }

  void EvolutionPBTMD::xfxQarray(double x,double Q,double*pdfs){
  std::cout<< " in PBTMD xfxQarray " << std::endl;
   }

  double EvolutionPBTMD::getAlphaS(double Q){
  // std::cout<< " in PBTMD getAlphaS " << std::endl;
    if (_itype == 5 ) { //if external evolution is used
      return gExternalEvolution->getAlphaS(Q);
    } else {
      int nfout,ierr;
      return QCDNUM::asfunc(Q*Q,nfout,ierr);
    }
  }

  std::vector<double> EvolutionPBTMD::getXgrid(){
  double test=0 ;
  std::cout<< " in PBTMD getXgrid - adapted from QCDnum" << std::endl;
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

  std::vector<double> EvolutionPBTMD::getQgrid(){
  double test=0;
  std::cout<< " in PBTMD getQgrid - adapted from QCDnum" << std::endl;
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

  void EvolutionPBTMD::Write_TMD(const char* name){
     std::cout<< " in evolutionPBTMD Write_PBTMD " << name << std::endl;
     int length = strlen(name);
//     std::cout << " file = "<< name << " length " << length<< std::endl;
     write_pbtmd(name, length);
  }

}
