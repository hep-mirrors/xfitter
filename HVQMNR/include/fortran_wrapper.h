#ifndef HVQMNR_FORTRAN_WRAPPER_H
#define HVQMNR_FORTRAN_WRAPPER_H


extern "C"
{
  // Interface to xFitter FORTRAN routines
	void hf_get_pdfs_(double *x, double *q2, double* pdf);
	double hf_get_alphas_(double* q2);
  void hf_stop_();
}


#endif // HVQMNR_FORTRAN_WRAPPER_H
