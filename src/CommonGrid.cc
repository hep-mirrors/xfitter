/*!
 @file CommonGrid.cc
 @date Tue March 25 2014
 @author Andrey Sapronov <Andrey.Sapronov@cern.ch>

 Contains CommonGrid class member function implementations.
 */

#include <algorithm>
#include <fstream>
#include <list>
#include <sstream>
#include <float.h>

#include "xfitter_cpp.h"
#include "CommonGrid.h"

#ifdef APPLGRID_ENABLED
#include "appl_grid/appl_grid.h"
#endif

#include <FastNLOxFitter.h>
#include <string.h>

using namespace std;
// using namespace appl;

// external pdf functions
extern "C" void appl_fnpdf_(const double& x, const double& Q, double* f);
extern "C" void appl_fnpdf_bar_(const double& x, const double& Q, double* f);
extern "C" void appl_fnpdf_neut_(const double& x, const double& Q, double* f);
extern "C" double appl_fnalphas_(const double& Q);
extern "C" void hf_errlog_(const int &id, const char *text, int); 

#ifdef APFELGRID_ENABLED
#include "APFELgrid/fastkernel.h"
#include "APFELgrid/transform.h"
// PDF in the evolution basis needed by APFELgrid
extern "C" void apfel_fnpdf_(const double& x, const double& Q, double* f);
void fkpdf (const double& x, const double& Q, const size_t& n, double* pdf)
{
  static double* lha_pdf = new double[13];
  apfel_fnpdf_(x,Q,lha_pdf);
  NNPDF::LHA2EVLN<double, double>(lha_pdf, pdf);
}
#include "APFELgridGeneration.h"
#endif

/*
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
*/

CommonGrid::CommonGrid(const string & grid_type, const string &grid_source): _dynamicscale(0)
{
  if ( grid_type == "applgrid" ) { 
    _flag = 0;
    this->readAPPLgrid(grid_source);
  } else if ( grid_type == "applgrid_norm" ) { 
    _flag = 1;
    this->readAPPLgrid(grid_source);
  } else if ( grid_type == "virtgrid" ) { 
    _flag = 2; 
    this->readVirtGrid(grid_source);
  } else if ( grid_type == "virtgrid_norm" ) { 
    _flag = 3; 
    this->readVirtGrid(grid_source);
  } else if ( grid_type == "applgrid_normhyperbin" ) { 
    _flag = 4;
    this->readAPPLgrid(grid_source);
  } else if ( grid_type == "virtgrid_normhyperbin" ) { 
    _flag = 5;
    this->readVirtGrid(grid_source);
  } else if ( grid_type.find("ast") != string::npos ) { // fastNLO
     this->initfastNLO(grid_source);
  } else if ( grid_type == "apfelgrid" ) { // APFELgrid
    this->readAPFELgrid(grid_source);
  } else {
    int id = 14032542;
    char text[] = "S: Unknown grid type in theory expression.";
    int textlen = strlen(text);
    hf_errlog_(id, text, textlen);
  }
}

CommonGrid::~CommonGrid(){
    vector<tHyperBin>::iterator ihb;
    for (ihb = _hbins.begin(); ihb != _hbins.end(); ihb++){
        delete[] ihb->b;
        if ( ihb->f  )  delete ihb->f;
        
#ifdef APPLGRID_ENABLED
        if ( ihb->g  )  delete ihb->g;
#endif
    
#ifdef APFELGRID_ENABLED
        if ( ihb->fk )  delete ihb->fk;
#endif
    }
}

int
CommonGrid::readAPPLgrid(const string &grid_source)
{

#ifdef APPLGRID_ENABLED
  double *b = new double;
  appl::grid *g = new appl::grid(grid_source);
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

  tHyperBin hb;
  hb.b   = b;
  hb.g   = g;
  hb.f   = NULL;
  hb.fk  = NULL;
  hb.ngb = g->Nobs();
  _hbins.push_back(hb);
  _ndim = 2;
  return _hbins.size();
#else
  int id = 16102401;
  char text[] = "S: APPLgrid must be present. Recompile with --enable-applgrid to use this option.";
  int textlen = strlen(text);
  hf_errlog_(id, text, textlen);
#endif
}


int
CommonGrid::readAPFELgrid(const string &grid_source)
{
#ifdef APFELGRID_ENABLED
  // Read FK table
  ifstream infile;
  infile.open(grid_source.c_str());
  NNPDF::FKTable<double> *FK;

  // If an FK table exists, read it
  if(infile) {
    FK = new NNPDF::FKTable<double>(infile);

    // Check that the relevant parameters agree with those given in the steering card
    double AsRef   = atof((FK->NNPDF::FKHeader::GetTag(NNPDF::FKHeader::THEORYINFO, "AlphasRef")).c_str());
    double QRef    = atof((FK->NNPDF::FKHeader::GetTag(NNPDF::FKHeader::THEORYINFO, "QRef")).c_str());
    double Q0      = atof((FK->NNPDF::FKHeader::GetTag(NNPDF::FKHeader::THEORYINFO, "Q0")).c_str());
    double MCharm  = atof((FK->NNPDF::FKHeader::GetTag(NNPDF::FKHeader::THEORYINFO, "MCharm")).c_str());
    double MBottom = atof((FK->NNPDF::FKHeader::GetTag(NNPDF::FKHeader::THEORYINFO, "MBottom")).c_str());
    double MTop    = atof((FK->NNPDF::FKHeader::GetTag(NNPDF::FKHeader::THEORYINFO, "MTop")).c_str());
    int    PtOrd   = atoi((FK->NNPDF::FKHeader::GetTag(NNPDF::FKHeader::THEORYINFO, "PerturbativeOrder")).c_str());

    double dAsRef   = abs( AsRef / c_alphas_.alphas_ - 1 );
    double dQRef    = abs( QRef / boson_masses_.mz_ - 1 );
    double dQ0      = abs( Q0 / sqrt(steering_.starting_scale_) - 1 );
    double dMCharm  = abs( MCharm / fermion_masses_.mch_ - 1);
    double dMBottom = abs( MBottom / fermion_masses_.mbt_ - 1);
    double dMTop    = abs( MTop / fermion_masses_.mtp_ - 1);
    int    dPtOrd   = abs( PtOrd + 1 - steering_.i_fit_order_);

    double toll = 1e-5;
    if( dAsRef > toll || dQRef > toll || dQ0 > toll || dMCharm > toll || dMBottom > toll || dMTop > toll || dPtOrd > 0 ) {
      cout << "The evolution parameters in '" << grid_source << "' do not correspond to those in the steering card. Regenerating FK table ..." << endl;
      cout << endl;

      // Move original file
      cout << "Moving '" << grid_source << "' into '" << grid_source << "-old'" << endl;
      cout << endl;
      rename(grid_source.c_str(), (grid_source + "-old").c_str());

      // Generate FK table is absent
      APFELgridGen::generateFK(grid_source, sqrt(steering_.starting_scale_),
      			       fermion_masses_.mch_, fermion_masses_.mbt_, fermion_masses_.mtp_,
			       c_alphas_.alphas_, boson_masses_.mz_,
			       steering_.i_fit_order_ - 1);

      // Reread FK table
      ifstream infile1;
      infile1.open(grid_source.c_str());
      FK = new NNPDF::FKTable<double>(infile1);
    }
  }
  // If no FK table exists, generate it
  else {
    cout << "The file '" << grid_source << "' does not exist. Generating FK table ..." << endl;
    cout << endl;

    // Generate FK table is absent
    APFELgridGen::generateFK(grid_source, sqrt(steering_.starting_scale_),
			     fermion_masses_.mch_, fermion_masses_.mbt_, fermion_masses_.mtp_,
			     c_alphas_.alphas_, boson_masses_.mz_,
			     steering_.i_fit_order_ - 1);

    // Read FK table
    ifstream infile1;
    infile1.open(grid_source.c_str());
    FK = new NNPDF::FKTable<double>(infile1);
  }

  tHyperBin hb;
  hb.b   = NULL;
  hb.g   = NULL;
  hb.f   = NULL;
  hb.fk  = FK;
  hb.ngb = FK->GetNData();
  _hbins.push_back(hb);
  _ndim = 2;
  return _hbins.size();
#else
  int id = 16051601;
  char text[] = "S: APFELgrid must be present. Recompile with --enable-apfelgrid to use this option.";
  int textlen = strlen(text);
  hf_errlog_(id, text, textlen);
#endif
}

int
CommonGrid::initfastNLO(const string &grid_source)
{
  FastNLOxFitter* fnlo = new FastNLOxFitter(grid_source);
   
  tHyperBin hb;
  hb.b   = NULL;
  hb.g   = NULL;
  hb.f   = fnlo;
  hb.fk  = NULL;
  hb.ngb = fnlo->GetNObsBin();
  _hbins.push_back(hb);
  _ndim = fnlo->GetNumDiffBin(); // _ndim not needed for fastNLO
  for ( int i = 0 ; i<_ndim ; i++ ) 
    if ( fnlo->GetIDiffBin(i) != 1 ) _ndim++;
  return _hbins.size();
}

int
CommonGrid::readVirtGrid(const string &grid_source)
{

#ifdef APPLGRID_ENABLED

  cout << "reading virtual grid from " << grid_source << endl;
  ifstream vg(grid_source.c_str());
  string line;
  if (vg.is_open()){
    while (1) {
      getline(vg,line);
      if (true == vg.eof()) break;
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
      if (0!=(nw-2)%2) {
        int id = 14040142;
        char text[] = "S: Bad number of bins in virtual grid. Each bin must have low and high value.";
        int textlen = strlen(text);
        hf_errlog_(id, text, textlen);
      }
      // set the grid dimension
      _ndim = nw;

      // now read bins
      sl.clear();
      sl.seekg(0);
      sl.str(line);
      double *b = new double[nw-2];
      for (int iw=0; iw<nw-2; iw++) sl >> b[iw];
      
      // and grid source
      string ag_source;
      sl>>ag_source;
      int n_grid_bins(0);
      sl>>n_grid_bins;
      appl::grid *g = new  appl::grid(ag_source);
      g->trim();

      // add the hyperbin
      tHyperBin hb;
      hb.b = b;
      hb.g = g;
      hb.ngb = n_grid_bins;
      _hbins.push_back(hb);
    }
    vg.close();
  }
  return _hbins.size();
#else
  int id = 16102401;
  char text[] = "S: APPLgrid must be present. Recompile with --enable-applgrid to use this option.";
  int textlen = strlen(text);
  return 0;
#endif

}

std::vector< std::vector<double> >
CommonGrid::vconvolute(const int iorder, const double mur, const double muf)
{
   // calculate cross sections with fastNLO or applgrid
   vector<tHyperBin>::iterator ihb;
   std::vector< std::vector<double> > result;
   for (ihb = _hbins.begin(); ihb != _hbins.end(); ihb++){
      if ( ihb->g ) 
	 // have to decrement iorder, since for appl_grid 0 -> LO, 1 -> NLO, etc.
	result.push_back(vconvolute_appl(iorder-1,mur,muf,&(*ihb)));
      else if ( ihb->f )
	result.push_back(vconvolute_fastnlo(iorder,mur,muf,ihb->f));
      else if ( ihb->fk )
	result.push_back(vconvolute_apfelg(&(*ihb)));
      else {
	 int id = 15010262;
	 char text[] = "S: Either applgrid or fastNLO or APFELgrid must be present.";
	 int textlen = strlen(text);
	 hf_errlog_(id, text, textlen);	 
      }
   }
   return result;
}

std::vector<double> 
CommonGrid::vconvolute_fastnlo(const int iorder, const double mur, const double muf, FastNLOxFitter* fnlo)
{
   //! calculate fastNLO cross sections
   if ( mur != fnlo->GetScaleFactorMuR() || mur != fnlo->GetScaleFactorMuF()) 
      fnlo->SetScaleFactorsMuRMuF(mur, muf);
   // iorder is already set earlier
   fnlo->FillAlphasCache();
   fnlo->FillPDFCache();
   fnlo->CalcCrossSection();
   return fnlo->GetCrossSection();
}

std::vector<double> 
CommonGrid::vconvolute_appl(const int iorder, const double mur, const double muf, tHyperBin* ihb)
{
#ifdef APPLGRID_ENABLED
  vector<double> xs;
  vector<double> gxs;
  appl::grid *g = ihb->g;
  // extract convoluted cross sections
  switch (_collision) {
    case PP : gxs = g->vconvolute(appl_fnpdf_, appl_fnalphas_, iorder, mur, muf); break;
    case PPBAR : gxs = g->vconvolute(appl_fnpdf_, appl_fnpdf_bar_, appl_fnalphas_, iorder, mur, muf); break;
    case PN : gxs = g->vconvolute(appl_fnpdf_, appl_fnpdf_neut_, appl_fnalphas_, iorder, mur, muf); break;
    default: {
      int id = 16020542;
      char text[] = "S: Unknown collision type for applgrid selected.";
      int textlen = strlen(text);
      hf_errlog_(id, text, textlen);
      break;
    }
  }

  // compute virtual bin width if available,
  double bw(1.);
  for (int ibb = 0; ibb<(_ndim-2)/2; ibb++){
     bw *= ihb->b[2*ibb+1] - ihb->b[2*ibb];
  }
  if ( 0!= bw){
     for(vector<double>::iterator ixs = gxs.begin(); ixs!=gxs.end(); ixs++) (*ixs)/=bw;
  } else {
     int id = 14040842;
     char text[] = "S: Bin width for virtual grid is 0, cannot scale.";
     int textlen = strlen(text);
     hf_errlog_(id, text, textlen);
  }

  // scale by hyperbin width if requested
  if(0 != (_flag & 1) && _flag <= 3) 
    for(vector<double>::iterator ixs = gxs.begin(); ixs!=gxs.end(); ixs++)
      (*ixs) *= bw * g->deltaobs(int(ixs - gxs.begin()));
  // scale by bin width if normalization is requested
  else if(_flag > 3) 
    for(vector<double>::iterator ixs = gxs.begin(); ixs!=gxs.end(); ixs++)
      (*ixs) *= bw;

  xs.insert(xs.end(), gxs.begin(), gxs.begin()+ihb->ngb);
  return xs;
#else
  int id = 16102401;
  char text[] = "S: APPLgrid must be present. Recompile with --enable-applgrid to use this option.";
  int textlen = strlen(text);
#endif
}

std::vector<double> 
CommonGrid::vconvolute_apfelg(tHyperBin* ihb)
{
#ifdef APFELGRID_ENABLED
  NNPDF::FKTable<double> *FK = ihb->fk;
  double* gxs = new double[FK->GetNData()];

  // extract convoluted cross sections
  switch (_collision) {
    case PP : FK->Convolute(fkpdf, 1, gxs); break;
    default: {
      int id = 1605121;
      char text[] = "S: Only PP collisions currently available in APFELgrid.";
      int textlen = strlen(text);
      hf_errlog_(id, text, textlen);
      break;
    }
  }
  vector<double> xs(gxs, gxs + FK->GetNData());
  return xs;
#else
  cout << "ERROR: calling dummy procedure for theory evaluation. Recompile with --enable-apfelgrid to use this option." << endl;
#endif
}

int
CommonGrid::checkBins(vector<int> &bin_flags, vector<vector<double> > &data_bins)
{
#ifdef APPLGRID_ENABLED
  // do not check bins for normalization grids
  if ( 0 != (_flag & 1) || _flag > 3 ) return 0;

  // compute number of bins in the virtual grid. Compare to the number of data bins
  int n_gbins(0);
  vector<tHyperBin>::iterator ihb;
  for (ihb = _hbins.begin(); ihb!=_hbins.end(); ihb++){
    n_gbins += ihb->ngb;
  }

  // Not checking number of bins, because grids can have wider range
  /*
  if ( n_gbins != data_bins[0].size() ) {
    int id = 14040242;
    char text[] = "S: Number of bins in grid doesn't match for data one.";
    int textlen = strlen(text);
    hf_errlog_(id, text, textlen);
  }
  */
    
  // compare the bins if flag is set
  vector<double> tv;
  vector<vector<double> > gb(data_bins.size(),tv);
  for (ihb = _hbins.begin(); ihb!=_hbins.end(); ihb++){
    appl::grid *g = ihb->g;
    
    for (int igb = 0; igb <ihb->ngb; igb++){  // cycle over applgrid bins
      for (int id = 0; id <_ndim-2; id++ ){ // cycle over hyperbins
        gb.at(id).push_back(ihb->b[id]);
      }
      gb.at(_ndim-2).push_back(g->obslow(igb));  // fill last dimension bins
      gb.at(_ndim-1).push_back(g->obslow(igb+1));  // with applgrid values
    }
  }

  for (int iv = 0; iv<_ndim; iv++){
    for (int ib = 0; ib<data_bins.at(iv).size(); ib++){
      if ( bin_flags.at(ib) == 0 ) continue;
      if ( 0 == (bin_flags.at(ib) & 2) ) {
        bool found_bin = false;
        for (int ibg = 0; ibg < gb[iv].size(); ibg++){
	  //          if (fabs(gb[iv][ibg] - data_bins[data_bins.size()-_ndim+iv][ib]) < 100*DBL_MIN) 
          if (fabs(*(short*)& gb[iv][ibg] - *(short*)& data_bins[data_bins.size()-_ndim+iv][ib]) < 2 ) 
	    found_bin = true;
	}
	if ( !found_bin ){ 
            int id = 14040342;
            char text[] = "S: Data and grid bins don't match.";
            int textlen = strlen(text);
            hf_errlog_(id, text, textlen);
            return -1;
	}
      }
    }
  }
  
  return 0;
#else
  int id = 16102401;
  char text[] = "S: APPLgrid must be present. Recompile with --enable-applgrid to use this option.";
  int textlen = strlen(text);
  hf_errlog_(id, text, textlen);
  return 0;
#endif

}

void
CommonGrid::SetCollisions(const string &collision)
{
  if (collision == string("pp")) _collision = PP;
  else if (collision == string("ppbar")) _collision = PPBAR;
  else if (collision == string("pn")) _collision = PN;
  else if (collision == string("pd")) _collision = PD;
  else if (collision == string("pnuc")) _collision = PNUC;
  else if (collision == string("nucnuc")) _collision = NUCNUC;
  else {
    int id = 16020401;
    char text[] = "S: Unknown collision type.";
    std::cerr << "Unknown collision type: " << collision << endl;
    int textlen = strlen(text);
    hf_errlog_(id, text, textlen);
  }
}

int 
CommonGrid::setCKM(const vector<double> &v_ckm)
{
#ifdef APPLGRID_ENABLED
  _hbins.at(0).g->setckm(v_ckm);
  return 0;
#else
  int id = 16102401;
  char text[] = "S: APPLgrid must be present. Recompile with --enable-applgrid to use this option.";
  int textlen = strlen(text);
  hf_errlog_(id, text, textlen);
#endif

}
