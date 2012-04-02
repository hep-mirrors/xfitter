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
//--------------- Internal -------------------------
//

#define HF_GET_PDFS_UNCACHED_WRAP F77_FUNC_ (hf_get_pdfs_uncached, HF_GET_PDFS_UNCACHED)
#ifdef __cplusplus
extern "C"  /* prevent C++ name mangling */
#endif
void HF_GET_PDFS_UNCACHED_WRAP(const double *x, const double *q2, double *pdfs);


#endif
