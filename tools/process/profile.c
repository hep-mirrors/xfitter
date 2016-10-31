#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include "pdf2yaml.h"
#include "utils.h"
#include "rotate.h"
#include <math.h>


static void help(){
        puts("xfitter-process profile [--piecewise-linear] pdf_shifts pdf_rotation pdf_dir_in pdf_dir_out");
	puts("Options are");
        puts("  --piecewise-linear: use piecewise linear approximation (default: quadratic approximation )");
        exit(0);
}

int profile(int argc, char* argv[]) {

        int i, ig, ix, iq, ifl;

        int flagc=0;
        int quad_approx=1;
        for(i=0; i<argc; i++) {
                if(!strcmp(argv[i], "--piecewise-linear")) { 
                        quad_approx=0;
                        flagc++;
                }
        }
       
        if((argc-flagc)!=4 || !strcmp(argv[0],"--help")) { help(); exit(0);}

        argv+=flagc;
        char *shifts_path=argv[0];
        char *rot_path=argv[1];
        char *in_path=argv[2];
        char *out_path=argv[3];
        char *line=NULL;
        Shifts shifts;

        

        char *in_path_tmp=strdup(in_path); //FREE
        char *pdf_in_name=basename(in_path_tmp); 

        //--------------- load shifts

        if( load_shifts(&shifts, shifts_path)) return EXIT_FAILURE;

        if(strcmp(shifts.pdf_name,pdf_in_name)) {
                fprintf(stderr, "input PDF set and shifts file are inconsistent:\n"
                                "\t\"%s\" in shifts, \"%s\" in input PDF set\n", shifts.pdf_name, pdf_in_name);
                exit(1);
        }


        //--------------- load rotation matrix

        RotMatrix rot_matrix;
        if(load_rotation_matrix(&rot_matrix, rot_path)) return EXIT_FAILURE;
        if(strcmp(rot_matrix.pdf_name,pdf_in_name)) {
                fprintf(stderr, "input PDF set and rotation matrix file are inconsistent:\n"
                        "\t\"%s\" in rotation matrix file, \"%s\" in input PDF set\n", rot_matrix.pdf_name, pdf_in_name);
                return(EXIT_FAILURE);
        }



        PdfSet pdf_set, shifted;
        if(load_lhapdf6_set(&pdf_set, in_path)) return 1;

        Info_Node *node=info_node_where(pdf_set.info, "ErrorType");
        char *error_type=strdup(node->value.string);
        printf("ErrorType: %s\n", error_type);

        Pdf  *Center;
        Center=pdf_dup(&pdf_set.members[0]);



        // prepare pdf delta, {M+, M-} - M0
        for(i=1; i<pdf_set.n_members; i++ ) 
                EACH_IN_PDF(Center, ig, ix, iq, ifl)
                pdf_set.members[i].val[ig][ix][iq][ifl]-=Center->val[ig][ix][iq][ifl]; 
        



        if(!strcmp(error_type,"hessian")) { 
	  if( 1  ) {

	    EACH_IN_PDF(&pdf_set.members[0], ig, ix, iq, ifl)
	      for(i=0; i<shifts.n; i++ ) {
		pdf_set.members[0].val[ig][ix][iq][ifl]+= 
		  shifts.val[i]*
		  (pdf_set.members[i*2+1].val[ig][ix][iq][ifl]
		   -pdf_set.members[i*2+2].val[ig][ix][iq][ifl])/2.0;
		
		pdf_set.members[0].val[ig][ix][iq][ifl] += 
		  shifts.val[i]*shifts.val[i]*
		  (pdf_set.members[i*2+1].val[ig][ix][iq][ifl]
		   +pdf_set.members[i*2+2].val[ig][ix][iq][ifl])/2.0;
	      }
	  } 
	  else {
	    EACH_IN_PDF(&pdf_set.members[0], ig, ix, iq, ifl)
	      for(i=0; i<shifts.n; i++ ) 
		pdf_set.members[0].val[ig][ix][iq][ifl]-= 
		  fabs(shifts.val[i])*( shifts.val[i] >0 ?
					pdf_set.members[i*2+2].val[ig][ix][iq][ifl]:
					pdf_set.members[i*2+1].val[ig][ix][iq][ifl]);
	  }
        } else if(!strcmp(error_type,"symmhessian")) {
	  
	  EACH_IN_PDF(&pdf_set.members[0], ig, ix, iq, ifl)
	    for(i=0; i<shifts.n; i++ ) 
	      pdf_set.members[0].val[ig][ix][iq][ifl] += 
		shifts.val[i]*pdf_set.members[i+1].val[ig][ix][iq][ifl];
        }
	



        for(i=1; i<pdf_set.n_members; i++ ) 
                EACH_IN_PDF(Center, ig, ix, iq, ifl)
                pdf_set.members[i].val[ig][ix][iq][ifl]+=pdf_set.members[0].val[ig][ix][iq][ifl]; 


        //rotation

        if(!strcmp(error_type,"hessian")) { 

                if(rot_matrix.n!=(pdf_set.n_members-1)/2.0) {
                        fprintf(stderr, "Number of pdf members and matrix dimension are inconsistent\n");
                        exit(EXIT_FAILURE);
                }

                if(quad_approx) hessian_quadapprox_rotate(&pdf_set, rot_matrix.val, &shifted);
                else hessian_by_M_direction_rotate(&pdf_set, rot_matrix.val, &shifted);

        } else if(!strcmp(error_type,"symmhessian")) {

                if(rot_matrix.n!=pdf_set.n_members-1) {
                        fprintf(stderr, "Number of pdf members and matrix dimension are inconsistent\n");
                        exit(EXIT_FAILURE);
                }

                symmhessian_rotate(&pdf_set, rot_matrix.val, &shifted);
        } else {
                fprintf(stderr, "Wrong ErrorType, need symmhessian or hessian\n");
                exit(1);
        }

	Info_Node* FP = info_node_where(shifted.info, "ForcePositive");
	if(FP == NULL) shifted.info = info_add_node_str(shifted.info, "ForcePositive", "0");
	else info_node_update_str(FP, "0");

        save_lhapdf6_set(&shifted, out_path);       
        puts("profiled\n");
        free(in_path_tmp);
        free(shifts.val);
        free(shifts.err);
        free(line);
        
        pdf_set_free(&pdf_set);
        pdf_free(Center);
        return 0;
}

