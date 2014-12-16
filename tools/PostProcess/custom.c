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
        char *in_path=argv[0];

        /***                      CUSTOM CODE                      ***/
        PdfSet pdf_set;
        if(load_lhapdf6_set(&pdf_set, in_path)) return 1;

        Pdf *copy=pdf_dup(&pdf_set.members[0]);

        info_save(copy->info, stdout);
        
        /***                    END CUSTOM CODE                    ***/
        pdf_set_free(&pdf_set);
        pdf_free(copy);

        return 0;
}
