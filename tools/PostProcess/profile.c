#include <stdio.h>
#include <stdlib.h>
#include "pdf2yaml.h"

typedef struct shifts_s {
        double *val;
        double *err;
        int n;
} Shifts;

static void help(){
        puts("Not implemented!!!");
        exit(0);
}

static Shifts parse_results(char *file) {

        FILE *fp=fopen(file, "r");
        char *line=NULL, *token=NULL, *(tokens[4]);
        size_t len;
        int i=0,j=0;
        Shifts shifts={NULL, NULL, 0};
        while(!feof(fp)) { 
                if(getline(&line, &len, fp));
                if(strstr(line, "PDF_nuisance_param_")) shifts.n++;
                free(line); line=NULL; len=0;
        }
        shifts.val=malloc(sizeof(double)*shifts.n);
        shifts.err=malloc(sizeof(double)*shifts.n);
        rewind(fp);
        while(!feof(fp)) {
                if(getline(&line, &len, fp));
                if(strstr(line, "PDF_nuisance_param_")) {
                        puts(line);
                        for(i=0, token=strtok(line, " \n\t"); token=strtok(NULL, " \n\t"); i++) tokens[i]=token;
                        shifts.val[j]=atof(tokens[1]);
                        shifts.err[j]=atof(tokens[3]);
                        j++;
                }
                free(line); line=NULL; len=0;
        }
        close(fp);
        return shifts;
}

//void shifts_free(Shifts *shifts);

int profile(int argc, char* argv[]) {
        if(!strcmp(argv[0],"--help")) help();

        Pdf pdf;
        load_lhapdf6_member(&pdf , argv[0]);
        printf(">> %d   %d\n", pdf.nx, pdf.nq);
//        Shifts shifts=parse_results("../../output/Results.txt");
//        int i;
//        for(i=0; i< shifts.n; i++) { 
//                printf("> %g ", shifts.val[i]);
//                printf("%g\n", shifts.err[i]);
//        }
        save_lhapdf6_member(&pdf , "/tmp/member");
//        pdf_free(&pdf);
        puts("profiled\n");
        return 0;
}

