#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include "pdf2yaml.h"

typedef struct shifts_s {
        double *val;
        double *err;
        int n;
} Shifts;

static void help(){
        puts("postproc profile [--no-omega] [--quad-appr] ");
        exit(0);
}

static void shifts_file_error(void) {
                fputs("wrong format of shifts file!", stderr);
                exit(1);
}

//void shifts_free(Shifts *shifts);

int profile(int argc, char* argv[]) {

        int i_err, i, j, ix, iq, ifl;

        int flagc=0;
        int omega=1;
        int quad_approx=0;
        for(i=0; i<argc; i++) {
                if(!strcmp(argv[i], "--no-omega")) { 
                        omega=0;
                        puts("no omega");
                        flagc++;
                }
                if(!strcmp(argv[i], "--quad-appr")) { 
                        quad_approx=1;
                        puts("quadratic approximation");
                        flagc++;
                }
        }
        if((argc-flagc)!=4 || !strcmp(argv[0],"--help")) { help(); exit(0);}
        argv+=flagc;
        char *shifts_path=argv[0];
        char *rot_path=argv[1];
        char *in_path=argv[2];
        char *out_path=argv[3];
        char *pdf_name;
        char *line=NULL;
        size_t len;
        FILE *fp;
        Shifts shifts;

        

        char *in_path_tmp=strdup(in_path); //FREE
        char *pdf_in_name=basename(in_path_tmp); 

        //--------------- parse shifts

        fp=fopen(shifts_path, "r");
        if(!fp) {
                fputs("cant open shifts file!", stderr);
                exit(1);
        }


        if(getline(&line, &len, fp)==-1) shifts_file_error();
        sscanf(line, "LHAPDF set=%ms", &pdf_name);
        if(strcmp(pdf_name,pdf_in_name)) {
                fprintf(stderr, "input PDF set and shifts file are inconsistent:\n"
                                "\t\"%s\" in shifts, \"%s\" in input PDF set\n", pdf_name, pdf_in_name);
                exit(1);
        }
        free(pdf_name);

        i_err=fscanf(fp, "%d", &shifts.n);
        if(!i_err || i_err==EOF) shifts_file_error();

        shifts.val=malloc(sizeof(double)*shifts.n);
        shifts.err=malloc(sizeof(double)*shifts.n);
        for(i=0; i< shifts.n; i++) {
                i_err=fscanf(fp, "%*d%lg%lg", &shifts.val[i], &shifts.err[i]);
                if(!i_err || i_err==EOF) shifts_file_error();
        }
        fclose(fp);


        //--------------- parse correlation matrix

        fp=fopen(rot_path, "r");
        if(!fp) {
                fputs("cant open rotation matrix file!", stderr);
                exit(1);
        }

        if(getline(&line, &len, fp)==-1) shifts_file_error();
        sscanf(line, "LHAPDF set=%ms", &pdf_name);
        if(strcmp(pdf_name,pdf_in_name)) {
                fprintf(stderr, "input PDF set and rotation matrix file are inconsistent:\n"
                        "\t\"%s\" in rotation matrix file, \"%s\" in input PDF set\n", pdf_name, pdf_in_name);
                exit(1);
        }
        free(pdf_name);

        int n_matrix;
        i_err=fscanf(fp, "%d", &n_matrix );
        if (!i_err || i_err==EOF){ 
                fputs("Bad rotation matrix", stderr);
                exit(1);
        }

//        printf("n_matrix: %d\n", n_matrix);

        double **rot_matrix=malloc(sizeof(double*)*n_matrix); //FREE
        for(i=0; i< n_matrix; i++) {
                i_err=fscanf(fp, "%*d");
                rot_matrix[i]= malloc(sizeof(double)*n_matrix); //FREE
                for(j=0; j<n_matrix; j++) { 
                        i_err=fscanf(fp, "%lg", &rot_matrix[i][j]); 
                        if(!i_err || i_err==EOF) {
                                fprintf(stderr, "bad rotation matrix format\n");
                                exit(1);
                        }
                }
        }

        PdfSet pdf_set;
        if(load_lhapdf6_set(&pdf_set, in_path)) return 1;

        Pdf  *Center;
        Center=pdf_dup(&pdf_set.members[0]);

        // prepare pdf delta, {M+, M-} - M0
        for(i=1; i<pdf_set.n_members; i++ ) 
        for(ix=0; ix<pdf_set.members[i].nx; ix++ ) 
        for(iq=0; iq<pdf_set.members[i].nq; iq++ ) 
        for(ifl=0; ifl< pdf_set.members[i].n_pdf_flavours; ifl++ ) 
                pdf_set.members[i].val[ix][iq][ifl]-=Center->val[ix][iq][ifl]; 
        
        for(ix=0; ix<pdf_set.members[0].nx; ix++ ) 
        for(iq=0; iq<pdf_set.members[0].nq; iq++ ) 
        for(ifl=0; ifl< pdf_set.members[0].n_pdf_flavours; ifl++ ) 
        for(i=0; i<shifts.n; i++ ) {
                pdf_set.members[0].val[ix][iq][ifl]+= 
                        shifts.val[i]*
                        (pdf_set.members[i*2+1].val[ix][iq][ifl]
                         -pdf_set.members[i*2+2].val[ix][iq][ifl])/2.0;

                         if(omega) pdf_set.members[0].val[ix][iq][ifl]-= 
                         shifts.val[i]*shifts.val[i]*
                         (pdf_set.members[i*2+1].val[ix][iq][ifl]
                          +pdf_set.members[i*2+2].val[ix][iq][ifl])/2.0;
                         }

        PdfSet *rot_set=pdf_set_dup(&pdf_set);

        
        if(quad_approx) {

        double ome, gamma;

        for(i=1; i<rot_set->n_members; i+=2 ) 
        for(ix=0; ix<rot_set->members[i].nx; ix++ ) 
        for(iq=0; iq<rot_set->members[i].nq; iq++ ) 
        for(ifl=0; ifl< rot_set->members[i].n_pdf_flavours; ifl++ ) {
                        rot_set->members[i].val[ix][iq][ifl]=0;
                        for(j=(i+1)%2+1; j<rot_set->n_members; j+=2) {
                        gamma=0.5*(pdf_set.members[j].val[ix][iq][ifl]-pdf_set.members[j+1].val[ix][iq][ifl]);
                        ome=0.5*(pdf_set.members[j].val[ix][iq][ifl]-pdf_set.members[j+1].val[ix][iq][ifl]);
                        rot_set->members[i].val[ix][iq][ifl]+=
                                rot_matrix[(i-1)/2][(j-1)/2]*(gamma+ome*rot_matrix[(i-1)/2][(j-1)/2]);
                        rot_set->members[i+1].val[ix][iq][ifl]-=
                                rot_matrix[(i-1)/2][(j-1)/2]*(gamma-ome*rot_matrix[(i-1)/2][(j-1)/2]);
                }
        }

        } else {
        for(i=1; i<rot_set->n_members; i+=2 ) 
        for(ix=0; ix<rot_set->members[i].nx; ix++ ) 
        for(iq=0; iq<rot_set->members[i].nq; iq++ ) 
        for(ifl=0; ifl< rot_set->members[i].n_pdf_flavours; ifl++ ) {
                        rot_set->members[i].val[ix][iq][ifl]=0;
                        for(j=(i+1)%2+1; j<rot_set->n_members; j+=2) {
                        if(rot_matrix[(i-1)/2][(j-1)/2]>0)
                        rot_set->members[i].val[ix][iq][ifl]+=
                                pdf_set.members[j].val[ix][iq][ifl]*rot_matrix[(i-1)/2][(j-1)/2];
                        else
                        rot_set->members[i].val[ix][iq][ifl]-=
                                pdf_set.members[j+1].val[ix][iq][ifl]*rot_matrix[(i-1)/2][(j-1)/2];
                }
        }

        for(i=2; i<rot_set->n_members; i+=2 ) 
        for(ix=0; ix<rot_set->members[i].nx; ix++ ) 
        for(iq=0; iq<rot_set->members[i].nq; iq++ ) 
        for(ifl=0; ifl< rot_set->members[i].n_pdf_flavours; ifl++ ) {
                        rot_set->members[i].val[ix][iq][ifl]=0;
                        for(j=(i+1)%2+1; j<rot_set->n_members; j+=2) {
                        if(rot_matrix[(i-1)/2][(j-1)/2]>0)
                        rot_set->members[i].val[ix][iq][ifl]+=
                                pdf_set.members[j].val[ix][iq][ifl]*rot_matrix[(i-1)/2][(j-1)/2];
                        else
                        rot_set->members[i].val[ix][iq][ifl]-=
                                pdf_set.members[j-1].val[ix][iq][ifl]*rot_matrix[(i-1)/2][(j-1)/2];
                }
        }
        }
        
        
        // {M'-, M'+} + M0 
        for(i=1; i<rot_set->n_members; i++ ) 
        for(ix=0; ix<rot_set->members[i].nx; ix++ ) 
        for(iq=0; iq<rot_set->members[i].nq; iq++ ) 
        for(ifl=0; ifl< rot_set->members[i].n_pdf_flavours; ifl++ ) 
                rot_set->members[i].val[ix][iq][ifl]=rot_set->members[0].val[ix][iq][ifl]+rot_set->members[i].val[ix][iq][ifl];
        
        save_lhapdf6_set(rot_set, out_path);       
        puts("profiled\n");
        free(in_path_tmp);
        free(shifts.val);
        free(shifts.err);
        free(line);
        
        for(i=0; i<n_matrix; i++) { 
                free(rot_matrix[i]);
        }

        pdf_set_free(&pdf_set);
        pdf_free(Center);
        pdf_set_free(rot_set);
        return 0;
}

