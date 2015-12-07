#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <sys/stat.h>
#include "pdf2yaml.h"
#include "utils.h"

static void help() {
        puts("usage: xfitter-process symmetrize pdf_dir_in pdf_dir_out");
        exit(EXIT_SUCCESS);
}

int symmetrize(int argc,char* argv[]) {
        if(argc!=2) help();

        int im, ig, ix, iq, ifl; // counters and dummy vars
        char *in_path=argv[0];

        PdfSet pdf_set, pdf_set_m, pdf_set_p, *errs, sym_set;
        if(load_lhapdf6_set(&pdf_set, in_path)) return EXIT_FAILURE;

        pdf_set_split_hessian_errors(&pdf_set, &pdf_set_m, &pdf_set_p);

        errs=pdf_set_dup(&pdf_set_m);
        EACH_IN_SET(errs, im, ig, ix, iq, ifl)
                errs->members[im].val[ig][ix][iq][ifl]=
                (pdf_set_p.members[im].val[ig][ix][iq][ifl]-pdf_set_m.members[im].val[ig][ix][iq][ifl])/2.0
                        +pdf_set.members[0].val[ig][ix][iq][ifl];

        pdf_set_join_central_errors(&pdf_set.members[0], errs, NULL, &sym_set);

        Info *sym_info=info_dup(pdf_set.info);

        char buf[30];
        Info_Node *node=info_node_where(sym_info, "NumMembers");
        n2str(buf, sym_set.n_members);
        info_node_update_str(node, buf);

        node=info_node_where(sym_info, "ErrorType");
        info_node_update_str(node, "symmhessian");

        sym_set.info=sym_info;
        save_lhapdf6_set(&sym_set, argv[1]);
        
        return 0;
}
