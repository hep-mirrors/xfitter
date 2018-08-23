#include "DyturboInterface.h"

#include "xfitter_cpp.h"

#include "dyturbo/dyturbo.h"
#include "dyturbo/vjint.h"
#include "dyturbo/vjloint.h"
#include "dyturbo/loint.h"
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
#include "dyturbo/abint.h"

#ifdef LHAPDF_ENABLED
#include <LHAPDF/LHAPDF.h>
#endif

#include <vector>
#include <string.h>
#include <algorithm>
#include <iomanip>

using namespace std;

string Dyturbo::pdfname;
int Dyturbo::pdfmember;

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
  opts.readfromfile(file);
  bins.readfromfile(file);
  opts.check_consitency();
  DYTurbo::init_params();


  //DYTurbo::WarmUp();
  //DYTurbo::PrintTable::Header();
  
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

  //DYTurbo::PrintTable::Settings();
  
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
  //order = iord-1;
  order = (iord-1)+1; //when running with applgrid V+jet, the order of applgrid is one order less
  
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

  //restore settings which are set in the data file INFO
  opts.order = order;
  opts.nproc = proc;
  //opts.mlow = ml;
  //opts.mhigh = mh;

  //restore LHAPDF settings (needed to read g from LHAPDF)
  opts.LHAPDFset = pdfname;
  opts.LHAPDFmember = pdfmember;
  
  // string lhapdfset = string(clhapdf_.lhapdfset_, 128);
  //issue with deleting last character from the PDF name
  //if (RunningMode == "LHAPDF Analysis")
  //lhapdfset = lhapdfset.erase(lhapdfset.find_last_not_of(" ")+1, string::npos); //--> use this for LHAPDF analysis
  //else if (RunningMode == "Chi2 Scan")
  //lhapdfset = lhapdfset.erase(lhapdfset.find_last_not_of(" "), string::npos); //--> use this for chi2scan analysis
  //int member = clhapdf_.ilhapdfset_;
  //opts.LHAPDFset = lhapdfset;
  //opts.LHAPDFmember = member;

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

  //recompute PDF moments
  pdfevol::init();
  pegasus::init();

  resint::init();
  */

  //Full reinitialisation
  //cout << "Start reinit" << endl;
  //DYTurbo::init_params();

    // init filling
    dofill_.doFill_ = 0;
    dyres::init();
    mcfm::init();               //This functions calls coupling::init()
    //pdf::init();                //Set up PDFs from LHAPDF, set alphas, and read g from the PDF
    iniflavreduce_();           //need to call this after nproc_.nproc_ is set
    coupling::initscales();
    //cc::init();                 //nodes and weights of Clenshaw-Curtis quadrature rules        --> skip this to speed up
    //gr::init();                 //nodes and weights of gaussian quadrature rules               --> skip this to speed up
    mellinint::initgauss();     //gaussian quadrature for mellin inversion
    mesq::init();               //EW couplings for born amplitudes
    rapint::init();             //allocate memory for the rapidity quadrature
    resconst::init();           //calculate beta, A and B coefficients
    anomalous::init();          //calculate anomalous dimensions, C1, C2 and gamma coefficients
    pdfevol::init();            //transform the PDF from x- to N-space at the factorisation scale
    if (!opts.fixedorder)
      {
	pegasus::init();        //initialise Pegasus QCD and transform the PDF from x- to N-space at the starting scale
	resint::init();         //initialise dequad integration for the bessel integral
      }
    vjint::init();
    vjloint::init();
    abint::init();              //alfa beta integration initialisation
    switching::init();          //switching function initialisation
    //itilde::init();             //itilde dequad initialisation
    rescinit_();

  //cout << "End reinit" << endl;
  //End reinitialisation

  //setup g and alphas from the current PDF set
  pdf::setg();
  pdf::setalphas();
    
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

  vector<double>::iterator itl = lowedge.begin();
  vector<double>::iterator itu = upedge.begin();
  //cout << "xw " << opts.xw << endl;  
  //cout << "scale_.scale_ " << scale_.scale_ << endl;  
  for (vector<double>::iterator it = values.begin(); it != values.end(); it++, itl++, itu++)
    {
      //Setbounds
      //cout << "bounds" << endl;
      //cout << yl << "  " << yh << "  " << ml << "  " << mh << "  " << opts.nproc << endl;
      //cout << yl << "  " << yh << "  " << *itl << "  " << *itu << "  " << opts.costhmin << "  " << opts.costhmax << "  " << opts.nproc << endl;
      //cout << "order: " << opts.order << endl;

      phasespace::setbounds(ml, mh, *itl, *itu, yl, yh); //pt bins
      //phasespace::setbounds(*itl, *itu, 0, 4000, yl, yh); //m bins
      //phasespace::setbounds(ml, mh, 0, 4000, *itl, *itu); //y bins
      phasespace::setcthbounds(opts.costhmin, opts.costhmax);

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
	  //resintegr2d(vals, error);
	  resintegr1d(vals, error);
	  //if (it - values.begin() < 10)
	  //cout << "RES result " << vals[0]/(*itu - *itl) << "  " << error << endl;
	  *it = vals[0];
	  ctintegr1d(vals, error);
	  //if (it - values.begin() < 10)
	  //cout << "CT result " << vals[0]/(*itu - *itl) << "  " << error << endl;
	  *it += vals[0];
	  //	  vjlointegr5d(vals, error);
	  //	  vjintegr3d(vals, error);
	  //	  cout << "V+J result " << vals[0]/(*itu - *itl) << "  " << error << endl;
	  //	  *it += vals[0];
	}
      else
	{
	  //vjlointegr5d(vals, error);
	  bornintegr2d(vals, error);
	  *it = vals[0];
	  cout << "BORN result " << vals[0] << "  " << error << endl;
	  if (opts.order > 0)
	    {
	      ctintegr2d(vals, error);
	      cout << "CT result " << vals[0] << "  " << error << endl;
	      *it += vals[0];
	      vjlointegr5d(vals, error);
	      cout << "V+J result " << vals[0] << "  " << error << endl;
	      *it += vals[0];
	    }
	}
      *it /= (*itu - *itl);
      //if ((it-values.begin()) < 10)
      //cout << "TOT result " << *itl << "  " << *itu << "  " << *it *(*itu - *itl) << endl;
    }
  return;
}
