#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <sys/stat.h>
#include "pdf2yaml.h"

static void help() {
        puts("template for user custom function");
        exit(0);
}

int custom(int argc,char* argv[]) {
        if(!strcmp(argv[0],"--help")) help();

        int i, j, ix, iq, ifl, i_tmp; // counters and dummy vars
        char *in_path1=argv[0];
        char *in_path2=argv[1];

        /***                      CUSTOM CODE                      ***/
        PdfSet pdf_set1;
        PdfSet pdf_set2;
        if(load_lhapdf6_set(&pdf_set1, in_path1)) return 1;
        if(load_lhapdf6_set(&pdf_set2, in_path2)) return 1;

        if(pdf_set_cmp(&pdf_set1, &pdf_set2)) puts("ok");

        /***                    END CUSTOM CODE                    ***/
        pdf_set_free(&pdf_set2);
        pdf_set_free(&pdf_set1);

        return 0;
}
