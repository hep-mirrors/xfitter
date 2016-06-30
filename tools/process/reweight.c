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
        puts("usage: xfitter-process reweight number_output_replicas pdf_weights pdf_dir_in pdf_dir_out");
        exit(0);
}

static void weights_file_error(void) {
                fputs("wrong format of weights file!", stderr);
                exit(1);
}

int reweight(int argc,char* argv[]) {

        int flagc=0;
        if((argc-flagc)!=4 || !strcmp(argv[0],"--help")) { help(); exit(0);}

        int i_err, i, im, ig, ix, iq, ifl;
	size_t j;
        argv+=flagc;
        int n_outreps=atoi(argv[0]);
        char *weights_path=argv[1];
        char *in_path=argv[2];
        char *out_path=argv[3];
        char *pdf_name;
        char *method;
	char* methodname="GIELE-KELLER";
        char *line=NULL;
        size_t len;
        FILE *fp;
        Weights weights;       
        Weights chi2;       

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
	int GK_method=0;

        if(getline(&line, &len, fp)==-1) weights_file_error();
	sscanf(line, "Reweight method=%ms", &method);
	if(strcmp(method,methodname)) GK_method=1;
	free(method);

	int ndata;
        if(getline(&line, &len, fp)==-1) weights_file_error();
	sscanf(line, "ndata=  %u", &ndata);

        i_err=fscanf(fp, "%d", &weights.n);
	i_err=i_err+1;

        if(!i_err || i_err==EOF) weights_file_error();

        weights.val=malloc(sizeof(double)*weights.n);
        chi2.val=malloc(sizeof(double)*weights.n);
        for(i=0; i< weights.n; i++) {
	  i_err=fscanf(fp, "%*d%lg%lg", &chi2.val[i], &weights.val[i]);  
             if(!i_err || i_err==EOF) weights_file_error();
        }

      fclose(fp);

        PdfSet pdf_set_in;
        if(load_lhapdf6_set(&pdf_set_in, in_path)) return 1;

        double step=(1.0/(double)(n_outreps+1.0));
        double nrep_old = (double)pdf_set_in.n_members-1;

        int *repnums=malloc(sizeof(int)*(n_outreps+1));
        double *wC=malloc(sizeof(double)*nrep_old);//weight cumulants
	double *vec_weight=malloc(sizeof(double)*nrep_old);//
	double *vec_chi2=malloc(sizeof(double)*nrep_old);//

        wC[0]=weights.val[0]/nrep_old;

        for(i=0; i< nrep_old; i++) {
	  wC[i+1]=(wC[i]+weights.val[i]/(double)nrep_old);
	  vec_weight[i]=((weights.val[i]));
	  vec_chi2[i]=((chi2.val[i]));
        }

#ifdef ROOT_ENABLED
	FillwHist(vec_weight, vec_chi2, ndata, nrep_old, GK_method);
#endif

        j=0;
        repnums[0]=0; // to copy info from original central set
	fprintf(stdout, "These are the PDF sets that are written out from the original set to the new reweighted PDF set. Please note, that some sets can occur twice or even more often. If there is only one or very few sets to be written out, then the reweighted set will not be so reliable (or the reweighting can be deemed a failure). Please make sure you were using consistently NNLO predictions with NNLO sets etc.\n");

	int ndiffsets = 0;
	int lastset = 0;

        for(i=1; i< n_outreps+1; i++) {
	 
	  double xpos = (double)i*step;

	  while (xpos>wC[j])
	    j++;

	  if (lastset!=j) {lastset=j;ndiffsets++;};

	  repnums[i]=j;//(repnums_old[j]); // vector containing new replicas in terms of member numbers of input set 

	  fprintf(stdout, "%i / %d\n", (int) j , n_outreps);

        }		
        /// 
	fprintf(stdout, "kept %i replicas out of %i input replicas \n", ndiffsets , (int)nrep_old);

        //create new set (DIFFERENT number of members!)
        PdfSet pdf_set;
        Info *info=info_dup(pdf_set_in.info);
        Info_Node *desc=info_node_where(info, "SetDesc");
        info_node_update_str(desc, "created by the xFitter package");
        Info_Node *nmem=info_node_where(info, "NumMembers");
        char buf[20];
        info_node_update_str(nmem, n2str(buf, n_outreps+1));

        pdf_set_initialize(&pdf_set, n_outreps+1, info);
        for(im=0; im<pdf_set.n_members; im++) {
                pdf_initialize_as(&pdf_set.members[im], &pdf_set_in.members[repnums[im]]);
                pdf_set.members[im].info=
                        info_add_node_str(pdf_set.members[im].info, "SetDesc", "reweighted");
        }

	EACH_IN_SET(&pdf_set, im, ig, ix, iq, ifl) {
	  pdf_set.members[0].val[ig][ix][iq][ifl]=0;
	}

       EACH_IN_SET_ERRORS(&pdf_set, im, ig, ix, iq, ifl) {
                pdf_set.members[im].val[ig][ix][iq][ifl]=pdf_set_in.members[repnums[im]].val[ig][ix][iq][ifl];
                pdf_set.members[0].val[ig][ix][iq][ifl]+=pdf_set_in.members[repnums[im]].val[ig][ix][iq][ifl]/(double)pdf_set.n_members;
        }


        save_lhapdf6_set(&pdf_set, out_path);
        pdf_set_free(&pdf_set);
        pdf_set_free(&pdf_set_in);
        //        free(desc);

        return 0;
}

