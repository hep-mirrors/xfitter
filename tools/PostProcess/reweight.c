#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <sys/stat.h>
#include "utils.h"
#include "pdf2yaml.h"

typedef struct weights_s {
        double *val;
        int n;
} Weights;

static void help(void) {
        puts("usage: postproc reweight number_output_replicas pdf_weights pdf_dir_in pdf_dir_out");
        exit(0);
}

static void weights_file_error(void) {
                fputs("wrong format of weights file!", stderr);
                exit(1);
}

int reweight(int argc,char* argv[]) {

        int flagc=0;
        if((argc-flagc)!=4 || !strcmp(argv[0],"--help")) { help(); exit(0);}

        int i_err, i, j, im, ig, ix, iq, ifl;

        argv+=flagc;
        int n_outreps=atoi(argv[0]);
        char *weights_path=argv[1];
        char *in_path=argv[2];
        char *out_path=argv[3];
        char *pdf_name;
        char *line=NULL;
        size_t len;
        FILE *fp;
        Weights weights;       

        char *in_path_tmp=strdup(in_path); //FREE
        char *pdf_in_name=basename(in_path_tmp); 

        //--------------- parse weights

        fp=fopen(weights_path, "r");
        if(!fp) {
                fputs("cant open weights file!", stderr);
                exit(1);
        }

        if(getline(&line, &len, fp)==-1) weights_file_error();
        sscanf(line, "LHAPDF set=%ms", &pdf_name);
        if(strcmp(pdf_name,pdf_in_name)) {
                fprintf(stderr, "input PDF set and weights file are inconsistent:\n"
                                "\t\"%s\" in weights, \"%s\" in input PDF set\n", pdf_name, pdf_in_name);
                exit(1);
        }
        free(pdf_name);

        i_err=fscanf(fp, "%d", &weights.n);
        if(!i_err || i_err==EOF) weights_file_error();

        weights.val=malloc(sizeof(double)*weights.n);
        for(i=0; i< weights.n; i++) {
                i_err=fscanf(fp, "%*d%*g%lg", &weights.val[i]);
                if(!i_err || i_err==EOF) weights_file_error();
        }
        fclose(fp);

        PdfSet pdf_set_in;
        if(load_lhapdf6_set(&pdf_set_in, in_path)) return 1;

        double step=(1.0/(double)(n_outreps+1.0));
        double nrep_old = (double)pdf_set_in.n_members;

        int *repnums=malloc(sizeof(int)*n_outreps);
        double *wC=malloc(sizeof(double)*nrep_old);//weight cumulants

        wC[0]=weights.val[0]/nrep_old;

        for(i=1; i< nrep_old; i++) {
                wC[i]=(wC[i-1]+weights.val[i]/(double)nrep_old);
        }

        // if ((wC.back()-1)>1e-12)?????
        /* if ((wC[i]-1)>1e-12) */
        /*   {   */
        /*     fputs("Error with weight normalisation", stderr); */
        /*     exit(1); */
        /*   } */


        j=0;
        for(i=1; i< n_outreps; i++) {

                double xpos = (double)i*step;

                while (xpos>wC[j])
                        j++;

                repnums[i]=j;//(repnums_old[j]); // vector containing new replicas in terms of member numbers of input set 
        }		
        /// 

        // printout kept replicas?

        //create new set (DIFFERENT number of members!)
        PdfSet pdf_set;
        Info *info=info_dup(pdf_set_in.info);
        Info_Node *desc=info_node_where(info, "SetDesc");
        info_node_update_str(desc, "created by the HeraFitter package");
        Info_Node *nmem=info_node_where(info, "NumMembers");
        char buf[20];
        info_node_update_str(nmem, n2str(buf, n_outreps));

        pdf_set_initialize(&pdf_set, n_outreps, info);
        for(im=0; im<pdf_set.n_members; im++) {
                pdf_initialize_as(&pdf_set.members[im], &pdf_set_in.members[im]);
                pdf_set.members[im].info=
                        info_add_node_str(pdf_set.members[im].info, "SetDesc", "reweighted");
        }

        EACH_IN_SET_ERRORS(&pdf_set, im, ig, ix, iq, ifl) {
                // here modify set (unweighting!!!!!!!!)
                pdf_set.members[im].val[ig][ix][iq][ifl]=pdf_set_in.members[repnums[im]].val[ig][ix][iq][ifl];
                pdf_set.members[0].val[ig][ix][iq][ifl]+=pdf_set_in.members[repnums[im]].val[ig][ix][iq][ifl]/(double)pdf_set.n_members;
        }


        save_lhapdf6_set(&pdf_set, out_path);
        pdf_set_free(&pdf_set);
        pdf_set_free(&pdf_set_in);
        //        free(desc);

        return 0;
}
