/**
  @file ext_pdfs.h

  @brief Generic external PDFs and alphaS calls for ACOT codes

  Provides an interface to call PDFs.

  @version 0.1
  @date 2017/04/16
 */

// Function to emulate LHAPDF xfx behavior:
typedef void   (*pXFXlike)(const double&, const double&, double*);
// and also PDF-like
typedef double (*pOneParFunc)(const double&);

extern "C" {
  void acot_get_pdfs_(const double&, const double&, double*);  //!< Return PDFs
  double acot_get_alphas_(const double&);                      //!< Return alphaS
  void acot_set_pdfs_alphaS( pXFXlike xfx, pOneParFunc aS);    //!< Set PDFs and alphaS
  
}
