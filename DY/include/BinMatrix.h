#ifndef BinMatrix_h
#define BinMatrix_h 1

#include <string>

#include "IntSteps.h"

typedef double**** double4d;

//!Class maintaining so called BinMatrix(BM).  
/*!
  BM is a PDF independent
  component of Drell-Yan cross section organized in a 4-dim matrix
  for each integration points of mass, rapidity, direction/flavor,
  eta bin. The arrangement is done to speed optimize PDF fitting 
  - the component matrix is calculated once for all iterations.
  Inherits from IntSteps.
 */

class BinMatrix : public IntSteps
{
 public:
  BinMatrix(){};
  ~BinMatrix();

  //! Proper constructor
  /*!
    Contains test on boson and variable names and
    calls corresponding bin matrix calculation routines.
    \param be   Beam energy
    \param ist  Integration steps to initialize the bin matrix
   */
  BinMatrix(const double *be, const IntSteps* ist);

 public:
  //! The bin matrix is accessible directly.
  double4d BM;

 protected:
  //! Methods for bin matrix construction for various cases.
  /*!
    Builds bin matrix for W process with eta binning.
    Contains nested loops in mass, rapidity, eta bin. 
    Direction index is (0 - fw, 1 - bw).
   */
  void BuildBM_W_eta();
  /*!
    Builds bin matrix for Z process with eta binning.
    Contains nested loops in mass, rapidity, eta bin, and
    component/flavor/direction index. The last one is
    a number from 0 to 11, with values corresponding to the following:
    - first 4 are for photon component
    - second 4 are for interference component
    - third 4 are for Z component
    Flavors and direction are accessed like this:
    - int comp = cdf/4; // component: 0 - gamma, 1 - gZ, 2 - Z
    - double dir = pow(-1.,cdf%4/2); // 0, 1 -> 1 (forw) ; 2,3 -> -1 (bw);
    - int flav = cdf%4%2; // 0 -> 0 (d), 1 -> 1(u), 2 -> 0(d), 3 -> 1(u)
   */
  void BuildBM_Z_eta();
  /*!
    Builds bin matrix for Z process with y binning.
    \sa BuildBM_Z_eta()
   */
  void BuildBM_Z_y();

  //! Analytical integration in eta
  /*! W cross section 
    \param m   Invariant mass
    \param y   Rapidity
    \param c_low_eta_lab Lower cosine(Theta) boundary in LAB frame
           corresponding to eta ranges.
    \param c_up_eta_lab Upper cosine(Theta) boundary in LAB frame
           corresponding to eta ranges.
   */
  double CosthAnIntW(const double &m,const double &y,
    const double &c_low_eta_lab, const double &c_up_eta_lab);
  /*! Z cross section 
    \param m   Invariant mass
    \param ry   Rapidity
    \param cdf Component/direction/flavor index
     - int comp = cdf/4; // component: 0 - gamma, 1 - gZ, 2 - Z
     - double dir = pow(-1.,cdf%4/2); // 0, 1 -> 1 (forw) ; 2,3 -> -1 (bw);
     - int flav = cdf%4%2; // 0 -> 0 (d), 1 -> 1(u), 2 -> 0(d), 3 -> 1(u)
    \param rc_low_eta_lab Lower cosine(Theta) boundary in LAB frame
           corresponding to eta ranges.
    \param rc_up_eta_lab Upper cosine(Theta) boundary in LAB frame
           corresponding to eta ranges.
   */
  double CosthAnIntZ(const double &m,const double &ry,
    const int cdf, 
    const double &rc_low_eta_lab, const double &rc_up_eta_lab);

  double _beam_en;

 public:
  //! Get beam energy
  double getBeamEn() { return _beam_en;}
  //! Get number of bins
  int getNbins() { return _nbins; }

 private:
  //! Lorentz transformation for cosine(Theta)
  /*!
   * \param b  Boost velocity
   * \param c  Cosine(Theta) to boost
   * \return   Boosted cosine
   */
  double costhLT(const double &b, const double &c);
};



#endif
