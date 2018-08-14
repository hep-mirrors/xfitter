 
/*
   @file EvolutionQCDNUM.cc
   @date 2018-08-14
   @author  AddEvolution.py
   Created by  AddEvolution.py on 2018-08-14
*/

#include "EvolutionQCDNUM.h"
#include "QCDNUM/QCDNUM.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"

// Global var to hold current pdfDecomposition
std::function<std::map<int,double>(double const& x)> gPdfDecomp; 

// Wrapper for QCDNUM
double funcPDF(int *ipdf, double *x) {
  const std::map <int,int> ip = { {1,-3}, {2,-2}, {3,-1}, {4,1}, {5,2}, {6,3}, {0,21} } ;
  return gPdfDecomp(*x)[ip.at(*ipdf)];
}

// PDF type

int const itype = 1;

// Check or not outside boundaries
int icheck = 0;

// Number of extra pdfs ( e.g. photon, etc.
int nEXT = 0;

//
double epsi = 1e-5;


//
double static   qcdnumDef[] = {
// tb bb  cb  sb ub   db  g   d   u   s   c   b   t
  0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  // sb
  0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,  // ub
  0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,  // db
  0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,  // d
  0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,  // u
  0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,  // s
  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  // 
  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  // 
  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  // 
  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  // 
  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  // 
  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.  // 
};


namespace xfitter
{

// the class factories
  extern "C" EvolutionQCDNUM* create() {
    return new EvolutionQCDNUM();
  }



  // Initialize at the start of the computation
  void EvolutionQCDNUM::initAtStart()
  {
    QCDNUM::qcinit(6," ");

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
    const int     PtOrder    = OrderMap(XFITTER_PARS::gParametersS.at("Order")) ;
    QCDNUM::setord(PtOrder);
    
    const double* Q0         = XFITTER_PARS::gParameters.at("Q0");
    const double* Q_ref      = XFITTER_PARS::gParameters.at("Mz");
    const double* Alphas_ref = XFITTER_PARS::gParameters.at("alphas");

    return ;
  }

  // Initialize at 
  void EvolutionQCDNUM::initAtIteration()
  {
    gPdfDecomp = _inPDFs;
    // XXXXXXXXXXXXXX
    const double* q0 = XFITTER_PARS::gParameters.at("Q0");
    int iq0  = QCDNUM::iqfrmq( (*q0) * (*q0) );  

    QCDNUM::evolfg(itype,funcPDF,qcdnumDef,iq0,epsi);
    return ;
  }

  std::function<std::map<int,double>(double const& x, double const& Q)> EvolutionQCDNUM::xfxQMap() {
    
    const auto _f0 =  [=] (double const& x, double const& Q) -> std::map<int, double> {
      std::map<int, double> res;
      for (int ipdf =-6; ipdf<7; ipdf++) {
	int ii = ( ipdf == 0 ) ? 21 : ipdf ;
	res[ii] = QCDNUM::fvalxq(itype,ipdf,x,Q*Q,icheck);
      }
      return res;
    };
    return _f0;
  }

  std::function<void(double const& x, double const& Q, double* pdfs)> EvolutionQCDNUM::xfxQArray() {
    const auto _f0 = [=] (double const& x, double const& Q, double* pdfs) {
      QCDNUM::allfxq(itype,x,Q*Q,pdfs,nEXT,icheck);
    };
    return _f0;
  }

  std::function<double(int const& i, double const& x, double const& Q)> EvolutionQCDNUM::xfxQDouble() {
    const auto _f0 = [=] (int const& i, double const& x, double const& Q) -> double {
      return  QCDNUM::fvalxq(itype,i,x,Q*Q,icheck);
    };
    return _f0;
  }

  std::function<double(double const& Q)> EvolutionQCDNUM::AlphaQCD() {
  }
}
