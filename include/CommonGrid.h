/**
  @class CommonGrid

  @brief A class for general grid representation

  Eh, what does it do?

  @author A.Sapronov <sapronov@ifh.de>

  @version 1.0

  @date 2014/03/25
 */
#ifndef CommonGrid_h
#define CommonGrid_h 1

#include <list>
#include <map>
#include <string>
#include <valarray>
#include <vector>

#ifdef APPLGRID_ENABLED
#include "appl_grid/appl_grid.h"
#endif
#include <FastNLOxFitter.h>
#ifdef APFELGRID_ENABLED
#include "APFELgrid/fastkernel.h"
#endif

using namespace std;
#ifdef APPLGRID_ENABLED
using namespace appl;
#endif

struct tHyperBin {
  double *b;
  int ngb;
#ifdef APPLGRID_ENABLED
  appl::grid *g;
#else
  void* g;
#endif
  FastNLOxFitter* f;
#ifdef APFELGRID_ENABLED
  NNPDF::FKTable<double> *fk;
#else
  void* fk;
#endif
  // default constructor
tHyperBin(): b(NULL), ngb(0), g(NULL), f(NULL), fk(NULL) {}
};

enum tCollisions { PP, PPBAR, PN, PD, PNUC, NUCNUC };

class CommonGrid {
 private:
  //! Hidden default constructor
  CommonGrid(){};
 public:
  ~CommonGrid();
  //! Proper constructor
  /*!
   \param grid_type the type of grid
   \param grid_source shows from where the grid should be read
  */
  CommonGrid(const string & grid_type, const string &grid_source);

  //! Selects if we have a proton-antiproton collision
  void SetCollisions(const string &collisions);// {_ppbar = (ppbar == 1);};
  void SetDynamicScale(double dynscale) {_dynamicscale = dynscale;};
  //! Check that the data and grid bins are consistent
  int checkBins(vector<int> &bin_flags, vector<vector<double> > &data_bins);

  //! Convolute the PDFs and return the cross section
  std::vector< std::vector<double> > vconvolute(const int iorder, const double mur, const double muf);

  //! Set custom CKM matrix for APPLgrid
  /*!
    \param v_ckm the 3x3 CMK matrix to be set in the APPLgrid

    The CKM matrix values used in APPLgrids calculations can be updated
    here.
   */
  int setCKM(const vector<double> &v_ckm);
   vector<tHyperBin>& getHBins() { return _hbins;} //!< get _hbins
   
 private:
  //! Read the applgrid
  int readAPPLgrid(const string &grid_source);
  //! Read the virtual grid
  int readVirtGrid(const string &grid_source);
  //! Read the virtual grid
  int initfastNLO(const string &grid_source);
  //! Read the APFELgrid
  int readAPFELgrid(const string &grid_source);

   //! convolute applgrid
   std::vector<double> vconvolute_appl(const int iorder, const double mur, const double muf, tHyperBin* ihb);
   //! calculate fastNLO
   std::vector<double> vconvolute_fastnlo(const int iorder, const double mur, const double muf, FastNLOxFitter* fnlo);
   //! convolute APFELgrid
   std::vector<double> vconvolute_apfelg(tHyperBin* ihb);

 private:
  int _ndim;
  /*!
   The _flag is a binary value:
    - first bin: 0/1 - standard/normalization grid
    - second bin: 0/1 - applgrid/virtgrid selection.
   For example, if datafile requests 'virtgrid_norm' option, the flag
   will be '3'
  */
  unsigned short _flag;
  vector<tHyperBin> _hbins;

  /// ppbar PDF
  tCollisions _collision;

  /// bin-by-bin dynamic scale
  double _dynamicscale;
};

#endif
