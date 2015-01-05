#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <sys/stat.h>
#include "pdf2yaml.h"

static void help() {
        puts("usage: postproc rotate matrix_file pdf_dir_in pdf_dir_out\n");
        exit(0);
}
int rotate(int argc,char* argv[]) {
        if(argc!=3) help();
        int i, j, ix, iq, ifl, i_tmp; // counters and dummy vars
        char *mfile=argv[0];
        char *in_path=argv[1];
        char *out_path=argv[2];

        char *in_path_tmp=strdup(in_path); //FREE
        char *pdf_in_name=basename(in_path_tmp); 

        char *out_path_tmp=strdup(out_path); //FREE
        char *pdf_out_name=basename(out_path_tmp);
        
        size_t len=0;
        char *info_path=NULL; 
        FILE *ss=open_memstream(&info_path, &len);
        fprintf(ss,"%s/%s.info",in_path,pdf_in_name);
        fflush(ss);
        FILE *info_fp=fopen(info_path, "r");
        fclose(ss);
        free(info_path);
        
        Info *info=info_load(info_fp);
        fclose(info_fp);


        Info_Node *node=info_node_where(info, "NumMembers");
        int n_members=atoi(node->value.string);
        printf("NumMembers: %d\n", n_members);
        node=info_node_where(info, "ErrorType");
        char *error_type=strdup(node->value.string);

        printf("ErrorType: %s\n", error_type);
        if(!strcmp(error_type,"replicas")) { 
                fprintf(stderr, "Wrong ErrorType, need symmhessian or hessian\n");
                exit(1);
        }

        Pdf *members=malloc(sizeof(Pdf)*(n_members)); 
        char *member_path;          //FREE
        for(i=0; i<n_members; i++) { 
                member_path=NULL; 
                ss=open_memstream(&member_path, &len);
                fprintf(ss,"%s/%s_%04d.dat",in_path,pdf_in_name,i);
                fflush(ss);
                printf("load: %s\n", member_path);
                load_lhapdf6_member(&members[i], member_path); //FREE
                fclose(ss);
                free(member_path);
        }
        
        FILE *matrix_fp=fopen(mfile, "r");
        if(!matrix_fp) { fprintf(stderr, "Error loading rotation matrix: %s\n", mfile); exit(1);}
        int n_matrix;
        if( fscanf(matrix_fp, "%d", &n_matrix ));

        printf("n_matrix: %d\n", n_matrix);

        double **rot_matrix=malloc(sizeof(double*)*n_matrix); //FREE
        for(i=0; i< n_matrix; i++) {
                if(fscanf(matrix_fp, "%d", &i_tmp));
                rot_matrix[i]= malloc(sizeof(double)*n_matrix); //FREE
                for(j=0; j<n_matrix; j++) if(!fscanf(matrix_fp, "%lg", &rot_matrix[i][j])) 
                        fprintf(stderr, "bad rotation matrix format\n");
        }

        
        Pdf *rot_members=malloc(sizeof(Pdf)*(n_members-1)); //FREE
        rot_members--;
        for(i=1; i<n_members;i++) {
                pdf_initialize(&rot_members[i], 
                                members[i].nx, members[i].nq, members[i].n_pdf_flavours, members[i].info);
                memcpy(rot_members[i].x, members[i].x, sizeof(double)*members[i].nx);
                memcpy(rot_members[i].q, members[i].q, sizeof(double)*members[i].nq);
                memcpy(rot_members[i].pdf_flavours, members[i].pdf_flavours, sizeof(int)*members[i].n_pdf_flavours);
        }

        // prepare pdf delta, {M+, M-} - M0
        for(i=1; i<n_members; i++ ) 
        for(ix=0; ix<members[i].nx; ix++ ) 
        for(iq=0; iq<members[i].nq; iq++ ) 
        for(ifl=0; ifl< members[i].n_pdf_flavours; ifl++ ) 
                members[i].val[ix][iq][ifl]-=members[0].val[ix][iq][ifl]; 
                
        // matrix product ( M- and M+ in parallel )
        if(!strcmp(error_type,"hessian")) { 
        if(n_matrix!=(n_members-1)/2) {
                fprintf(stderr, "Number of pdf members and matrix dimension are inconsistent\n");
        }
        for(i=1; i<n_members; i++ ) 
        for(ix=0; ix<rot_members[i].nx; ix++ ) 
        for(iq=0; iq<rot_members[i].nq; iq++ ) 
        for(ifl=0; ifl< rot_members[i].n_pdf_flavours; ifl++ ) {
                        rot_members[i].val[ix][iq][ifl]=0;
                        for(j=(i+1)%2+1; j<n_members; j+=2) {
                        rot_members[i].val[ix][iq][ifl]+=members[j].val[ix][iq][ifl]*rot_matrix[(i-1)/2][(j-1)/2];
                }
        }
        
        } else if(!strcmp(error_type,"symmhessian")) {
        if(n_matrix!=n_members-1) {
                fprintf(stderr, "Number of pdf members and matrix dimension are inconsistent\n");
        }
        for(i=1; i<n_members; i++ ) 
        for(ix=0; ix<rot_members[i].nx; ix++ ) 
        for(iq=0; iq<rot_members[i].nq; iq++ ) 
        for(ifl=0; ifl< rot_members[i].n_pdf_flavours; ifl++ ) {
                        rot_members[i].val[ix][iq][ifl]=0;
                        for(j=1; j<n_members; j++) {
                        rot_members[i].val[ix][iq][ifl]+=members[j].val[ix][iq][ifl]*rot_matrix[i-1][j-1];
                }
        }
        } else {
                fprintf(stderr, "Wrong ErrorType, need symmhessian or hessian\n");
        }
                
        // {M'-, M'+} + M0 
        for(i=1; i<n_members; i++ ) 
        for(ix=0; ix<rot_members[i].nx; ix++ ) 
        for(iq=0; iq<rot_members[i].nq; iq++ ) 
        for(ifl=0; ifl< rot_members[i].n_pdf_flavours; ifl++ ) 
                rot_members[i].val[ix][iq][ifl]+=members[0].val[ix][iq][ifl]; 

        // save result
        mkdir(out_path, 0755);

        info_path=NULL; 
        ss=open_memstream(&info_path, &len);
        fprintf(ss,"%s/%s.info", out_path, pdf_out_name);
        fflush(ss);
        info_fp=fopen(info_path,"w");
        if(!info_fp) fprintf(stderr, "Error: can't open file %s\n", info_path);
        info_save(info, info_fp);
        fclose(ss);
        free(info_path);

        for(i=0; i<n_members; i++) { 
                member_path=NULL; 
                ss=open_memstream(&member_path, &len);
                fprintf(ss, "%s/%s_%04d.dat", out_path, pdf_out_name, i);
                fflush(ss);
                printf("save: %s\n", member_path);
                if(i==0) save_lhapdf6_member(&members[i], member_path); 
                else save_lhapdf6_member(&rot_members[i], member_path); 
                fclose(ss);
                free(member_path);
        }

        //clean up
        pdf_free(&members[0]);
        for(i=1; i<n_members; i++) {
                pdf_free(&members[i]);
//                rot_members[i].info=NULL;
                pdf_free(&rot_members[i]);
        }
        free(members);
        free(++rot_members);

        info_free(info);
        for(i=0; i<n_matrix; i++) free(rot_matrix[i]);
        free(rot_matrix);
        free(in_path_tmp);
        free(out_path_tmp);
        free(error_type);


        return 0;
}

