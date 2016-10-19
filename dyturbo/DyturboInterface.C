#include "DyturboInterface.h"

#include "xfitter_cpp.h"

#include "dyturbo/dyturbo.h"
#include "dyturbo/phasespace.h"
#include "dyturbo/settings.h"
#include "dyturbo/cubacall.h"
//#include "dyturbo/integr.h"
#include "dyturbo/ctintegr.h"
#include "dyturbo/pdf.h"
#include "dyturbo/pdfevol.h"
#include "dyturbo/pegasus.h"
#include "dyturbo/resint.h"
#include "dyturbo/coupling.h"
#include "dyturbo/rapint.h"
#include "dyturbo/mcfm_interface.h"
#include "dyturbo/dyres_interface.h"

#include "dyturbo/mesq.h"
#include "dyturbo/anomalous.h"
#include "dyturbo/switch.h"
#include "dyturbo/gaussrules.h"
#include "dyturbo/resconst.h"

#ifdef LHAPDF_ENABLED
#include <LHAPDF/LHAPDF.h>
#endif

#include <vector>
#include <string.h>
#include <algorithm>
#include <iomanip>

using namespace std;

//Constructor
Dyturbo::Dyturbo(string file): infile(file)
{
  //int argc = 2;
  //char** argv = new char*[10] ;
  //std::strcpy(argv[0], file.c_str());
  //DYTurbo::Init(argc,argv);
  gaussinit_();             //initialisation of fortran gaussian quadrature nodes and weights
  coupling::SMparameters(); //initialisation of unused MCFM parameters
  // parsing options from input file
  opts.readfromfile      (file);
  DYTurbo::init_params();


  //DYTurbo::WarmUp();
  //DYTurbo::PrintTable::Header();
  //DYTurbo::PrintTable::Settings();
  
  proc = opts.nproc;
  
  //set up the right PDF in LHAPDF, to override DYTURBO pdf
  string lhapdfset = string(clhapdf_.lhapdfset_, 128);
  lhapdfset = lhapdfset.erase(lhapdfset.find_last_not_of(" ")+1, string::npos);
  int member = clhapdf_.ilhapdfset_;
#ifdef LHAPDF_ENABLED
  LHAPDF::initPDFSet(lhapdfset.c_str());
  LHAPDF::initPDF(member);
  c_alphas_.alphas_ = LHAPDF::alphasPDF(boson_masses_.mz_);
#endif

  int iord;
  double mur0;
  double muf0;
}

void Dyturbo::SetBins(vector <double> ledge, vector <double> uedge, double ylow, double yhigh, double mlow, double mhigh)
{
  lowedge = ledge;
  upedge = uedge;
  vector<double>::iterator itl = lowedge.begin();
  vector<double>::iterator itu = upedge.begin();
  values.resize(lowedge.size());
  yl = ylow;
  yh = yhigh;
  ml = mlow;
  mh = mhigh;
}

void Dyturbo::SetOrdScales(int iord, double kmuren, double kmufac, double kmures, double muC3)
{
  order = iord-1;
  
  opts.order = order;
  nnlo_.order_ = opts.order;
  opts.kmuren = kmuren;
  opts.kmufac = kmufac;
  if (kmures > 0)
    opts.kmures = kmures;
  coupling::initscales();
  opts.C3 = muC3;
}

void Dyturbo::Calculate(const double muren, const double mufac, const double mures, const double muC3)
{
  //re-read input file
  opts.readfromfile(infile.c_str());

  //restore overwritten settings
  opts.order = order;
  opts.nproc = proc;
  //opts.mlow = ml;
  //opts.mhigh = mh;
  string lhapdfset = string(clhapdf_.lhapdfset_, 128);
  lhapdfset = lhapdfset.erase(lhapdfset.find_last_not_of(" ")+1, string::npos);
  int member = clhapdf_.ilhapdfset_;
  opts.LHAPDFset = lhapdfset;
  opts.LHAPDFmember = member;

  //set scales
  opts.kmuren = muren;
  opts.kmufac = mufac;
  opts.kmures = mures;
  opts.C3 = muC3;
  
  //Minimal reinitialisation
  /*
  dyres::init();
  mcfm::init();
  iniflavreduce_();
  coupling::initscales();
  */

  //Full reinitialisation
    DYTurbo::init_params();
    /*
    dofill_.doFill_ = 0;
    dyres::init();
    mcfm::init();
    iniflavreduce_(); //need to call this after nproc_.nproc_ is set
    coupling::initscales();
    //C++ resum
    //initialise all the C modules
    gr::init(); //nodes and weights of gaussian quadrature rules
    mellinint::initgauss(); //gaussian quadrature for mellin inversion
    mesq::init(); //EW couplings for born amplitudes
    rapint::init(); //allocate memory for the rapidity quadrature
    resconst::init(); //calculate beta, A and B coefficients
    anomalous::init(); //calculate anomalous dimensions, C1, C2 and gamma coefficients
    pdfevol::init(); //transform the PDF from x- to N-space at the factorisation scale
    pegasus::init(); //initialise Pegasus QCD and transform the PDF from x- to N-space at the starting scale
    resint::init(); //initialise dequad integration for the bessel integral
    //end C++ resum
    //V+j fixed order initialisation
    vjint::init();
    vjloint::init();
    //
    loint::init(); //Born term initialisation
    switching::init(); //switching function initialisation
    rescinit_();
    */
  //End reinitialisation

  setg();
  setalphas();

  int mode = 0;
  double costh = 0.1; double m = 91; double qt = 5; double y = 0.2;
  resumm_(costh,m,qt,y,mode);
  double f[opts.totpdf];
  ctint_(costh,m,qt,y,mode,f);

  /*
#ifdef LHAPDF_ENABLED
  LHAPDF::PDFInfo info(lhapdfset, member);
  double gformfactor = info.get_entry_as<double>("g", g_param_.g_param_);
  cout << "g for this PDF set is " << gformfactor << endl;
  g_param_.g_param_ = gformfactor;
  cout << "alphas for this set is " << info.get_entry_as<double>("AlphaS_MZ") << endl;
  setalphas();
#endif
  */
  
  //recompute PDF moments
  if (!opts.fixedorder)
    if (opts.resumcpp)
      {
	pdfevol::init();
	pegasus::init();
	resint::init();
      }
    else
      initmoments_();

  vector<double>::iterator itl = lowedge.begin();
  vector<double>::iterator itu = upedge.begin();
  //cout << "xw " << opts.xw << endl;  
  //  cout << "scale_.scale_ " << scale_.scale_ << endl;  
  for (vector<double>::iterator it = values.begin(); it != values.end(); it++, itl++, itu++)
    {
      //Setbounds
      //cout << yl << "  " << yh << "  " << ml << "  " << mh << "  " << opts.nproc << endl;
      phasespace::setbounds(ml, mh, *itl, *itu, yl, yh);
      //get cross section
      double error;
      vector <double> vals;
      if (!opts.fixedorder)
	{
	  if (opts.resumcpp)
	    //C++ resum
	    rapint::cache(phasespace::ymin, phasespace::ymax);
	  //end C++ resum
	  else	    
	    cacheyrapint_(phasespace::ymin, phasespace::ymax);
	  resintegr2d(vals, error);
	  cout << "resummation result " << vals[0]/(*itu - *itl) << "  " << error/(*itu - *itl) << endl;
	  *it = vals[0];
	  ctintegr2d(vals, error);
	  cout << "counterterm result " << vals[0]/(*itu - *itl) << "  " << error/(*itu - *itl) << endl;
	  *it += vals[0];
	  *it /= (*itu - *itl);
	  cout << *itl << "  " << *itu << "  " << *it << endl;
	}
      else
	{
	  vjlointegr5d(vals, error);
	  //cout << "V+J LO result " << vals[0]/(*itu - *itl) << "  " << error/(*itu - *itl) << endl;
	  //cout << "V+J LO result " << *itl << "  " << *itu << "  " << vals[0] << "  " << error << endl;
	  *it = vals[0];
	  //*it /= (*itu - *itl);
	}
    }
  return;
}
