#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <sys/stat.h>
#include "pdf2yaml.h"

static void help() {
        puts("usage: postproc symmetrize pdf_dir_in pdf_dir_out");
        exit(0);
}

int symmetrize(int argc,char* argv[]) {
        if(argc!=2) help();

        int i, j, ix, iq, ifl, i_tmp; // counters and dummy vars
        char *in_path=argv[0];

        PdfSet pdf_set;
        if(load_lhapdf6_set(&pdf_set, in_path)) return 1;

        for(i=1; i<pdf_set.n_members; i+=2) 
        for(ix=0; ix<pdf_set.members[i].nx; ix++) 
        for(iq=0; iq<pdf_set.members[i].nq; iq++) 
        for(ifl=0; ifl<pdf_set.members[i].n_pdf_flavours; ifl++) {
                pdf_set.members[i].val[ix][iq][ifl]=(pdf_set.members[i].val[ix][iq][ifl] //TODO
                                -pdf_set.members[i+1].val[ix][iq][ifl])/2.0;
                        //+pdf_set.members[0].val[ix][iq][ifl];
        }
        for(i=1; i<pdf_set.n_members; i+=2) pdf_free(&pdf_set.members[i+1]);

        Pdf *asym_members=pdf_set.members;
        Pdf *sym_members=malloc(sizeof(Pdf)*((pdf_set.n_members-1)/2+1));
        sym_members[0]=asym_members[0];
        for(i=1; i< (pdf_set.n_members-1)/2+1; i++)
                sym_members[i]=asym_members[i*2-1];
        pdf_set.n_members=(pdf_set.n_members-1)/2+1;

        char buf[30];
        Info_Node *node_n=info_node_where(pdf_set.info, "NumMembers");
        free(node_n->value.string);
        node_n->value.string=n2str(buf, pdf_set.n_members);
        pdf_set.members=sym_members;

        Info_Node *node_err=info_node_where(pdf_set.info, "ErrorType");
        free(node_err->value.string);
        node_err->value.string="symmhessian";

        save_lhapdf6_set(&pdf_set, argv[1]);
        node_n->value.string=NULL;
        node_err->value.string=NULL;
        
        pdf_set_free(&pdf_set);
        free(asym_members);

        return 0;
}
