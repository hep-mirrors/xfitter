#include <stdio.h>
#include "pdf2yaml.h"
#include "utils.h"

void pdf_set_error_sets(PdfSet *dest, PdfSet *src) {
        int im;
        pdf_set_initialize(dest, src->n_members-1, NULL);

        for(im=1; im<src->n_members; im++)
                pdf_cpy(&dest->members[im-1], &src->members[im]);
}

void pdf_set_join_central_errors(Pdf *central, PdfSet *errors, Info *set_info, PdfSet *dest) {
        int im;

        pdf_set_initialize(dest, errors->n_members+1, set_info);
        for(im=0; im<errors->n_members; im++)
                pdf_cpy(&dest->members[im+1], &errors->members[im]);

                pdf_cpy(&dest->members[0], central);
}

int pdf_set_split_hessian_errors(const PdfSet *pdf_set, PdfSet *pdf_set_m, PdfSet *pdf_set_p) {
        int im, n_split_members;
        Info_Node *node=info_node_where(pdf_set->info, "ErrorType");
        if(strcmp(node->value.string, "hessian")) {
                fprintf(stderr, "Wrong ErrorType, need hessian. At %s, line %d.\n", __FILE__, __LINE__);
                return EXIT_FAILURE;
        }

        n_split_members=(pdf_set->n_members-1)/2;

        pdf_set_initialize(pdf_set_m, n_split_members, NULL);
        pdf_set_initialize(pdf_set_p, n_split_members, NULL);
        for(im=0; im<n_split_members; im++) {
                pdf_cpy(&pdf_set_p->members[im], &pdf_set->members[im*2+1]); 
                pdf_cpy(&pdf_set_m->members[im], &pdf_set->members[im*2+2]); 
        }

        return EXIT_SUCCESS;
}

int pdf_set_join_hessian_errors(const PdfSet *original, const PdfSet *pdf_set_m, const PdfSet *pdf_set_p, PdfSet *dest) {

        int im;

        Info_Node *node=info_node_where(original->info, "ErrorType");
        if(strcmp(node->value.string, "hessian")) {
                fprintf(stderr, "Wrong ErrorType, need hessian. At %s, line %d.\n", __FILE__, __LINE__);
                return EXIT_FAILURE;
        }

        if(pdf_set_p->n_members!=pdf_set_m->n_members) {
                fprintf(stderr, "pdf_set_m and pdf_set_p are inconsistent. At %s, line %d.\n", __FILE__, __LINE__);
                return EXIT_FAILURE;
        }

        if(pdf_set_p->n_members!=(original->n_members-1)/2) {
                fprintf(stderr, "pdf_set_m, pdf_set_p and original are inconsistent. At %s, line %d.\n", __FILE__, __LINE__);
                return EXIT_FAILURE;
        }

        pdf_set_initialize(dest, original->n_members, original->info);
        for(im=0; im<pdf_set_m->n_members; im++) {
                pdf_cpy(&dest->members[im*2+1], &pdf_set_p->members[im]); 
                pdf_cpy(&dest->members[im*2+2], &pdf_set_m->members[im]); 
        }
        pdf_cpy(&dest->members[0], &original->members[0]);


        return EXIT_SUCCESS;
}

int load_rotation_matrix(RotMatrix *rot_matrix, char *file_path) {
        int i_tmp, i, j;
        char *line=NULL;
        size_t len=0;

        FILE *fp=fopen(file_path, "r");
        if(!fp) {
                fputs("cant open rotation matrix file!", stderr);
                return(EXIT_FAILURE);
        }

        if(getline(&line, &len, fp)==-1) {
                fputs("cant open rotation matrix file!", stderr);
                return(EXIT_FAILURE);
        }
        sscanf(line, "LHAPDF set=%ms", &rot_matrix->pdf_name);

//        int n_matrix;
        i_tmp=fscanf(fp, "%d", &rot_matrix->n );
        if (!i_tmp || i_tmp==EOF){ 
                fputs("Bad rotation matrix", stderr);
                return(EXIT_FAILURE);
        }

        rot_matrix->val=malloc(sizeof(double*)*rot_matrix->n); //FREE
        for(i=0; i< rot_matrix->n; i++) {
                i_tmp=fscanf(fp, "%*d");
                rot_matrix->val[i]= malloc(sizeof(double)*rot_matrix->n); //FREE
                for(j=0; j<rot_matrix->n; j++) { 
                        i_tmp=fscanf(fp, "%lg", &rot_matrix->val[i][j]); 
                        if(!i_tmp || i_tmp==EOF) {
                                fprintf(stderr, "bad rotation matrix format. At %s, line %d.\n", __FILE__, __LINE__);
                                return(EXIT_FAILURE);
                        }
                }
        }
        return EXIT_SUCCESS;
}


static void shifts_file_error(void) {
                fputs("wrong format of shifts file!", stderr);
                exit(EXIT_FAILURE);
}

int load_shifts(Shifts *shifts, char *file_path) {
        int i, i_err;
        char *line=NULL;
        size_t len;
        FILE *fp=fopen(file_path, "r");
        if(!fp) {
                fputs("cant open shifts file!", stderr);
                exit(1);
        }


        if(getline(&line, &len, fp)==-1) shifts_file_error();
        sscanf(line, "LHAPDF set=%ms", &shifts->pdf_name);

        i_err=fscanf(fp, "%d", &shifts->n);
        if(!i_err || i_err==EOF) shifts_file_error();

        shifts->val=malloc(sizeof(double)*shifts->n);
        shifts->err=malloc(sizeof(double)*shifts->n);
        for(i=0; i< shifts->n; i++) {
                i_err=fscanf(fp, "%*d%lg%lg", &shifts->val[i], &shifts->err[i]);
                if(!i_err || i_err==EOF) shifts_file_error();
        }
        fclose(fp);
        return EXIT_SUCCESS;
}
