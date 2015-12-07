#ifndef rotate_h
#define rotate_h

#include "pdf2yaml.h"

int symmhessian_rotate(PdfSet *pdf_set, double **rot_matrix, PdfSet *rotated);
int hessian_quadapprox_rotate(PdfSet *pdf_set_in, double **rot_matrix, PdfSet *rotated);
int hessian_by_M_direction_rotate(PdfSet *pdf_set_in, double **rot_matrix, PdfSet *rotated);

#endif
