#ifndef utils_h
#define utils_h
#include"pdf2yaml.h"
void pdf_set_error_sets(PdfSet *dest, PdfSet *src);
void pdf_set_join_central_errors(Pdf *central, PdfSet *errors, Info *set_info, PdfSet *dest);
int pdf_set_split_hessian_errors(const PdfSet *pdf_set, PdfSet *pdf_set_m, PdfSet *pdf_set_p);
int pdf_set_join_hessian_errors(const PdfSet *original, const PdfSet *pdf_set_m, const PdfSet *pdf_set_p, PdfSet *dest);

typedef struct rot_matrix_s {
        double **val;
        int n;
        char *pdf_name;
} RotMatrix;

//int load_rotation_matrix(double ***rot_matrix, int *n_matrix, char *path, char *pdf_name);
int load_rotation_matrix(RotMatrix *rot_matrix, char *file_path);

typedef struct shifts_s {
        double *val;
        double *err;
        int n;
        char *pdf_name;
} Shifts;

int load_shifts(Shifts *shifts, char *file_path);
#endif
