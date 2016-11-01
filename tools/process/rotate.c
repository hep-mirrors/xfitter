#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <sys/stat.h>
#include "pdf2yaml.h"
#include "rotate.h"
#include "utils.h"





static void help() {
        puts("usage: xfitter-process rotate [--piecewise-linear] matrix_file pdf_dir_in pdf_dir_out\n");
	puts("Options are");
        puts("  --piecewise-linear: use piecewise linear approximation (default: quadratic approximation )");
        exit(0);
}




int rotate(int argc,char* argv[]) {
        int i;
        int flagc=0;
        int quad_approx=1;
        for(i=0; i<argc; i++) {
                if(!strcmp(argv[i], "--piecewise-linear")) { 
                        quad_approx=0;
                        flagc++;
                }
        }

        if((argc-flagc)!=3 || !strcmp(argv[0],"--help")) { help(); exit(0);}

        argv+=flagc;

        char *mfile=argv[0];
        char *in_path=argv[1];
        char *out_path=argv[2];

        char *in_path_tmp=strdup(in_path); //FREE
        char *pdf_in_name=basename(in_path_tmp); 

        PdfSet pdf_set_in, rotated;
        load_lhapdf6_set(&pdf_set_in, in_path);

        Info_Node *node=info_node_where(pdf_set_in.info, "ErrorType");
        char *error_type=strdup(node->value.string);
        printf("ErrorType: %s\n", error_type);
        if(!strcmp(error_type,"replicas")) { 
                fprintf(stderr, "Wrong ErrorType, need symmhessian or hessian\n");
                exit(EXIT_FAILURE);
        }

        RotMatrix rot_matrix;
        load_rotation_matrix(&rot_matrix, mfile);

        if(strcmp(rot_matrix.pdf_name,pdf_in_name)) {
                fprintf(stderr, "input PDF set and rotation matrix file are inconsistent:\n"
                        "\t\"%s\" in rotation matrix file, \"%s\" in input PDF set\n", rot_matrix.pdf_name, pdf_in_name);
                return(EXIT_FAILURE);
        }

        if(!strcmp(error_type,"hessian")) { 

                if(rot_matrix.n!=(pdf_set_in.n_members-1)/2.0) {
                        fprintf(stderr, "Number of pdf members and matrix dimension are inconsistent\n");
                        exit(EXIT_FAILURE);
                }

                if(quad_approx) hessian_quadapprox_rotate(&pdf_set_in, rot_matrix.val, &rotated);
                else hessian_by_M_direction_rotate(&pdf_set_in, rot_matrix.val, &rotated);

        } else if(!strcmp(error_type,"symmhessian")) {

                if(rot_matrix.n!=pdf_set_in.n_members-1) {
                        fprintf(stderr, "Number of pdf members and matrix dimension are inconsistent\n");
                        exit(EXIT_FAILURE);
                }

                symmhessian_rotate(&pdf_set_in, rot_matrix.val, &rotated);
        } else {
                fprintf(stderr, "Wrong ErrorType, need symmhessian or hessian\n");
                exit(1);
        }

        Info_Node* FP = info_node_where(rotated.info, "ForcePositive");
        if(FP == NULL) rotated.info = info_add_node_str(rotated.info, "ForcePositive", "0");
        else info_node_update_str(FP, "0");

        save_lhapdf6_set(&rotated, out_path);
        
        return EXIT_SUCCESS;
}




int hessian_quadapprox_rotate(PdfSet *pdf_set_in, double **rot_matrix, PdfSet *rotated){

        int j, im, ig, ix, iq, ifl; // counters and dummy vars
        double ome, gamma;
        PdfSet pdf_set_m, pdf_set_p, *rot_m, *rot_p;

        printf("Use hessian quadratic approximation\n");

        pdf_set_split_hessian_errors(pdf_set_in, &pdf_set_m, &pdf_set_p);

        rot_m=pdf_set_dup(&pdf_set_m);
        rot_p=pdf_set_dup(&pdf_set_p);

        EACH_IN_SET(&pdf_set_m, im, ig, ix, iq, ifl) 
                pdf_set_m.members[im].val[ig][ix][iq][ifl]-=pdf_set_in->members[0].val[ig][ix][iq][ifl];
        EACH_IN_SET(&pdf_set_p, im, ig, ix, iq, ifl) 
                pdf_set_p.members[im].val[ig][ix][iq][ifl]-=pdf_set_in->members[0].val[ig][ix][iq][ifl];

        EACH_IN_SET(&pdf_set_m, im, ig, ix, iq, ifl) {
                rot_m->members[im].val[ig][ix][iq][ifl]=0;
                rot_p->members[im].val[ig][ix][iq][ifl]=0;
                for(j=0; j<pdf_set_m.n_members;j++) {
                        gamma=0.5*(pdf_set_p.members[j].val[ig][ix][iq][ifl]
                                        -pdf_set_m.members[j].val[ig][ix][iq][ifl]);
                        ome=0.5*(pdf_set_m.members[j].val[ig][ix][iq][ifl]
                                        +pdf_set_p.members[j].val[ig][ix][iq][ifl]);
                        rot_p->members[im].val[ig][ix][iq][ifl]+=
                                rot_matrix[im][j]*(gamma+ome*rot_matrix[im][j]);

                        rot_m->members[im].val[ig][ix][iq][ifl]-=
                                rot_matrix[im][j]*(gamma-ome*rot_matrix[im][j]);
                }

        }

        EACH_IN_SET(&pdf_set_m, im, ig, ix, iq, ifl) {
                        rot_m->members[im].val[ig][ix][iq][ifl]+=
                                pdf_set_in->members[0].val[ig][ix][iq][ifl];
                        rot_p->members[im].val[ig][ix][iq][ifl]+=
                                pdf_set_in->members[0].val[ig][ix][iq][ifl];

        }
        pdf_set_join_hessian_errors(pdf_set_in, rot_m, rot_p, rotated);
        return EXIT_SUCCESS;
}




int hessian_by_M_direction_rotate(PdfSet *pdf_set_in, double **rot_matrix, PdfSet *rotated) {
        int j, im, ig, ix, iq, ifl; // counters and dummy vars
        PdfSet pdf_set_m, pdf_set_p, *rot_m, *rot_p;
        pdf_set_split_hessian_errors(pdf_set_in, &pdf_set_m, &pdf_set_p);

        printf("Use piecewise-linear approximation\n");
        rot_m=pdf_set_dup(&pdf_set_m);
        rot_p=pdf_set_dup(&pdf_set_p);

        EACH_IN_SET(&pdf_set_m, im, ig, ix, iq, ifl) 
                pdf_set_m.members[im].val[ig][ix][iq][ifl]-=pdf_set_in->members[0].val[ig][ix][iq][ifl];
        EACH_IN_SET(&pdf_set_p, im, ig, ix, iq, ifl) 
                pdf_set_p.members[im].val[ig][ix][iq][ifl]-=pdf_set_in->members[0].val[ig][ix][iq][ifl];

        EACH_IN_SET(&pdf_set_m, im, ig, ix, iq, ifl) {
                rot_m->members[im].val[ig][ix][iq][ifl]=0;
                for(j=0; j<pdf_set_m.n_members;j++) {
                        if(rot_matrix[im][j]>0) {
                                rot_m->members[im].val[ig][ix][iq][ifl]+=
                                        rot_matrix[im][j]*pdf_set_m.members[j].val[ig][ix][iq][ifl]; 
                        } else {
                                rot_m->members[im].val[ig][ix][iq][ifl]-=
                                        rot_matrix[im][j]*pdf_set_p.members[j].val[ig][ix][iq][ifl]; 
                        }

                }
        }

        EACH_IN_SET(&pdf_set_p, im, ig, ix, iq, ifl) {
                rot_p->members[im].val[ig][ix][iq][ifl]=0;
                for(j=0; j<pdf_set_p.n_members;j++) {
                        if(rot_matrix[im][j]>0) {
                                rot_p->members[im].val[ig][ix][iq][ifl]+=
                                        rot_matrix[im][j]*pdf_set_p.members[j].val[ig][ix][iq][ifl]; 
                        } else {
                                rot_p->members[im].val[ig][ix][iq][ifl]-=
                                        rot_matrix[im][j]*pdf_set_m.members[j].val[ig][ix][iq][ifl]; 
                        }

                }
        }

        EACH_IN_SET(&pdf_set_m, im, ig, ix, iq, ifl) {
                        rot_m->members[im].val[ig][ix][iq][ifl]+=
                                pdf_set_in->members[0].val[ig][ix][iq][ifl];
                        rot_p->members[im].val[ig][ix][iq][ifl]+=
                                pdf_set_in->members[0].val[ig][ix][iq][ifl];

        }
        pdf_set_join_hessian_errors(pdf_set_in, rot_m, rot_p, rotated);
        return EXIT_SUCCESS;
}




int symmhessian_rotate(PdfSet *pdf_set, double **rot_matrix, PdfSet *rotated) {
        int j, im, ig, ix, iq, ifl; // counters and dummy vars
        PdfSet errs, *rot_errs;
        pdf_set_error_sets(&errs, pdf_set);

        rot_errs=pdf_set_dup(&errs);


        EACH_IN_SET(&errs, im, ig, ix, iq, ifl) 
                                errs.members[im].val[ig][ix][iq][ifl]-=
                                pdf_set->members[0].val[ig][ix][iq][ifl];


        EACH_IN_SET(&errs, im, ig, ix, iq, ifl) {
                        rot_errs->members[im].val[ig][ix][iq][ifl]=0;
                        for(j=0; j<errs.n_members; j++) {
                        rot_errs->members[im].val[ig][ix][iq][ifl]+=
                                errs.members[j].val[ig][ix][iq][ifl]*rot_matrix[im][j];
                }
        }

        EACH_IN_SET(&errs, im, ig, ix, iq, ifl) 
                                rot_errs->members[im].val[ig][ix][iq][ifl]+=
                                pdf_set->members[0].val[ig][ix][iq][ifl];

        pdf_set_join_central_errors(&pdf_set->members[0], rot_errs, pdf_set->info, rotated); 
        
        return EXIT_SUCCESS;
}
