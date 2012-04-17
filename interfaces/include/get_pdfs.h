#ifndef get_pdfs_h
#define get_pdfs_h

// Define interfaces to PDF calls

//
// Function to obtain PDFs
// 
#define HF_GET_PDFS_WRAP F77_FUNC_ (hf_get_pdfs, HF_GET_PDFS)
#ifdef __cplusplus
extern "C"  /* prevent C++ name mangling */
#endif
void HF_GET_PDFS_WRAP(const double *x, const double *q2, double *pdfs);


//
// Function to get alphaS
//
#define HF_GET_ALPHAS_WRAP F77_FUNC_ (hf_get_alphas, HF_GET_ALPHAS)
#ifdef __cplusplus
extern "C"  /* prevent C++ name mangling */
#endif
double HF_GET_ALPHAS_WRAP(const double *q2);


//
// Function to get muR scale
// 
#define HF_GET_MUR_WRAP F77_FUNC_ (hf_get_mur, HF_GET_MUR)
#ifdef __cplusplus
extern "C"  /* prevent C++ name mangling */
#endif
double HF_GET_MUR_WRAP(const int *idataset);


//
// Function to get muF scale
// 
#define HF_GET_MUF_WRAP F77_FUNC_ (hf_get_muf, HF_GET_MUF)
#ifdef __cplusplus
extern "C"  /* prevent C++ name mangling */
#endif
double HF_GET_MUF_WRAP(const int *idataset);


//
// Function to get PDFs fast 
//
#define HF_PDFFAST_WRAP F77_FUNC_ (hf_pdffast, HF_PDFFAST)
#ifdef __cplusplus
extern "C"  /* prevent C++ name mangling */
#endif
void HF_PDFFAST_WRAP(const double *pdffdef, 
		     const double *x, 
		     const double *q2,
		     double *PDFout,
		     const int   *Npoints);

//
//--------------- Internal -------------------------
//

#define HF_GET_PDFS_UNCACHED_WRAP F77_FUNC_ (hf_get_pdfs_uncached, HF_GET_PDFS_UNCACHED)
#ifdef __cplusplus
extern "C"  /* prevent C++ name mangling */
#endif
void HF_GET_PDFS_UNCACHED_WRAP(const double *x, const double *q2, double *pdfs);


#endif
