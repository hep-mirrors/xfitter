#ifndef PDFconv_h
#define PDFconv_h 1

#include <string>

#include "IntSteps.h"

//! Manges PDF convolutions
/*!
  PDFconv inherits from IntSteps. The class maintains 
  parton function grids for optimized interpolation in QCDNUM.
  Provides convolutions of PDFS.
 */
class PDFconv:public IntSteps
{
 public:
  PDFconv(){};
  ~PDFconv();

  //! Proper constructor
  /*!
    Sets data. Initialiezes fdef matrix of interpolation
    linear coefficients for QCDNUM. Assigns function pointers
    to convolution routines based on the boson name. 
    Calls init() for further initialisation.
    \sa init()
    \param chg_prod  Beam charge product to specify LHC or Tevatron
    \param beam_en   Beam CM energy
    \param ist       IntSteps object to initialize with
   */
  PDFconv(const int chg_prod, const double * beam_en,
         const IntSteps *ist);
  PDFconv(const PDFconv&);

  int interpPDF();

  //! Function pointer to convolution routine.
  int (PDFconv::*getPDFconv)(const int, const int, double, 
      double &, double &);
  //! Convolution routine for Z process
  /*!
   \param imp   Mass integration point index
   \param iyp   Rapidity integration point index
   \param dir   Direction, -1 or 1
   \param xfxcD Convolution for down flavors
   \param xfxcU Convolution for up flavors
   */
  int getPDFconvZ(const int imp, const int iyp, double dir, 
    double &xfxcD, double &xfxcU);
  //! Convolution routine for W process
  /*!
   \param imp    Mass integration point index
   \param iyp    Rapidity integration point index
   \param dir    Direction, -1 or 1
   \param xfxcWm Convolution for Wminus
   \param xfxcWp Convolution for Wplus
   */
  int getPDFconvW(const int imp, const int iyp, double dir, 
    double &xfxcWm, double &xfxcWp);

 private:
  //! Number of interpolation points for QCDNUM grid
  int _NINT;
  //! Index array
  /*! 
   An array of indeces of integration points to interpolation
   elements association. The dimensions are number of mass 
   points by number of y points. The values are indeces of 
   corresponding elements in QCDNUM grid.
   */
  int **_IDX;

  //! Charge product -1 or 1
  int _chg_prod;
  //! Beam energy
  double _beam_en;

  //! Number of flavors, 13
  /*!
   Flavor numbering:
   \bar{t   b   c   s   u   d}  g   d   u   s   c   b   t
        0   1   2   3   4   5   6   7   8   9  10  11  12
   */
  static const int _nfl;
  //! fdef matrix of QCDNUM linear coefficients of size _nfl*_nfl
  double *_fdef;

 protected:
  //! Array of interpolation points in x
  double *_XINT;
  //! Array of interpolation points in q2
  double *_Q2INT;
  //! Matrix of fast interpolated PDF values
  /*!
   Dimensions are _nfl*_NINT
   */
  double **_XFINT;

 private:
  //! Initialiser
  /*!
   Initialises interpalotion arrays _XINT, _Q2INT. Sets correspondence
   between Simpson integration points and QCDNUM interpolation elements.
   */
  void init();

 public:
  //! Checks that we have same beam in this class
  bool isSameBeam(const int chg_prod, const double *beam_en){
    if ( chg_prod == _chg_prod && *beam_en == _beam_en ) return true;
    else return false;
  }
};

#endif
