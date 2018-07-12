/**
  @file ext_pdfs.h

  @brief Generic external PDFs and alphaS calls for fortran APFEL codes

  Provides an interface to call PDFs.

  @version 0.1
  @date 2018/07/10
 */

// Function to emulate LHAPDF xfx behavior:
typedef void (*pXFXlike)(const double&, const double&, double*);

extern "C" {
  void externalsetapfel1_(const double&, const double&, double*);  //!< Return PDFs
  void APFEL_set_pdfs(pXFXlike xfx);  //! Set PDFs
}

