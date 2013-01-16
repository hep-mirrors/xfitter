#ifndef IntSteps_h
#define IntSteps_h 1

#include <string>

//! Base class for all others in DY. 
/*! 
  Maintains integration grid for 
  invariant mass and boson rapidity. Stores binning information. 
 */
class IntSteps
{
 public:
  IntSteps(){};
  ~IntSteps();

  IntSteps(const IntSteps&);
  //! Proper constructor
  /*! 
    Contains tests against boson and variables names and 
    creates corresponding grids.
    \param boz       Boson name, "W" or "Z"
    \param ranges    7 element array with integration ranges 
    { m_low, m_up, y_low, y_up, eta_low, eta_up, pt_el} 
    \param var_name  Binning variable name, "y" or "eta"
    \param n_bins    Number of bins
    \param var_bins  Bin edges array for the variable 
   */
  IntSteps(const std::string boz, const double *ranges, 
    const std::string var_name, const int n_bins, const double*var_bins);

 protected:
  //! Boson name, "W" or "Z"
  std::string _boz;
  //! Binning variable name, "y" or "eta"
  std::string _var_name;
  //! Number of bins
  int _nbins;
  //! Array of bins
  double *_bins;
  //! Array of bin starting steps
  int *_ssb; 

  //! 2-element Arrays of ranges for mass, rapidity and eta
  double *_mr;
  double *_yr;
  double *_etar;

  //! Number of Simpson integration steps for mass and rapidity
  int _nms;
  int _nys;
  //! Arrays with Simpson integration steps for mass and rapidity
  double *_msteps;
  double *_ysteps;

  //! Dirty lepton pt cut storage
  double _leptPtCut;

 private:
  //! Create integration steps for masses and rapidities
  /*!
    Create inv mass integration steps for Z
   */
  void _makeMstepsZ();
  /*!
    Create inv mass integration steps for W
   */
  void _makeMstepsW();
  /*!
   * Create rapidity integration steps.
   */
  void _makeYsteps();
  /*!
   * Create rapidity integration steps if we have y binning.
   */
  void _makeYstepsBinned();

 public:
  //! Returns boson name
  std::string getBozName(){ return _boz;}

  //! Returns grid dimension info
  void getDims(int &nms, int &nys){
    nms = _nms;
    nys = _nys;
  }

  //! Returns m steps array
  void getMsteps(int &nms, double *m_steps){
    nms=_nms;
    m_steps=_msteps;
  }

  //! Returns y steps array
  void getYsteps(int &nys, double *y_steps){
    nys=_nys;
    y_steps=_ysteps;
  }

  //! dirty lepton pt cut storage
  void setLeptPtCut(const double leptPtCut){
    _leptPtCut = leptPtCut;
  }

  double getLeptPtCut() { return _leptPtCut; }
};


#endif
