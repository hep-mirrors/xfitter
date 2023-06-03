#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <sys/stat.h>
#include "pdf2yaml.h"
#include "utils.h"

static void help() {
        puts("xfitter-process scale90to68 pdf_dir_in pdf_dir_out");
        exit(0);
}

int scale90to68(int argc,char* argv[]) {
        if((argc) != 2 || !strcmp(argv[0],"--help")) help();

        int im, ig, ix, iq, ifl; // counters and dummy vars
        char *in_path=argv[0];
        char *out_path=argv[1];

        double scale_factor = 1./1.645;

        //sscanf(argv[2], "%lg", &scale_factor);

        PdfSet in_set;
        PdfSet err_set;
        PdfSet out_set;
        if(load_lhapdf6_set(&in_set, in_path)) return 1;

        pdf_set_error_sets(&err_set, &in_set);
        Pdf *central=in_set.members;

        EACH_IN_SET(&err_set, im, ig, ix, iq, ifl) {
                err_set.members[im].val[ig][ix][iq][ifl]-=central->val[ig][ix][iq][ifl];
                err_set.members[im].val[ig][ix][iq][ifl]*=scale_factor;
                err_set.members[im].val[ig][ix][iq][ifl]+=central->val[ig][ix][iq][ifl];
        }

        pdf_set_join_central_errors(central, &err_set, in_set.info, &out_set);

        save_lhapdf6_set(&out_set, out_path);



        pdf_set_free(&in_set);
        pdf_set_free(&err_set);
        pdf_set_free(&out_set);

        return 0;
}
