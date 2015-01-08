#include "pdf2yaml.h"
#include <yaml.h>
#include <string.h>
#include "c2yaml.h"
#include "list.h"
#define CHECK(FUNC) \
        do { printf("checking %s: ",#FUNC); \
        if(!FUNC()) printf("failed\n"); else printf("ok\n"); } while(0);

int TEST_info_add_node_str(){
        int res;
        Info *info=NULL;
        info=info_add_node_str(info, "KEY", "VALUE");
        Info_Node *node=info->data;
        char *str=node->value.string;
        res = (strcmp("VALUE", str)==0 && strcmp("KEY", node->key)==0 );
        info_free(info);
        return res;
}

int TEST_info_add_node_darray(){
        int res;
        Info *info=NULL;
        double darray[]={0.0001,3.14,30};
        info=info_add_node_darray(info, "KEY", darray, 3);
        Info_Node *node=info->data;
        res = (node->value.darray.size==3 && 
                        node->value.darray.vals[1]==3.14 && 
                        strcmp("KEY", node->key)==0 );
        info_free(info);
        return res;
}

int TEST_info_save() {
        int res;
        Info *info=NULL;
        info=info_add_node_str(info, "KEY", "VALUE");
        FILE *fp=tmpfile();
        info_save(info, fp);
        fp=freopen(NULL, "r", fp);
        size_t len;
        char *str=NULL;
        if(!getline(&str, &len, fp)) return 0;
        res = (strcmp("KEY: VALUE\n", str)==0);
        close(fp);
        info_free(info);
        free(str);
        return res;
}

int TEST_info_dup() {
        int res;
        Info *info=NULL;
        info=info_add_node_str(info, "KEY1", "VALUE1");
        double darray[]={10,20,30};
        info=info_add_node_darray(info, "KEY2", darray, 3);

        Info *copy=info_dup(info);
        res = info_cmp(copy, info);
        info_free(info);
        info_free(copy);
        return res;
}

int TEST_info_node_where() {
        int res;
        double darray[]={10,20,30};
        Info *info=NULL;
        info=info_add_node_str(info, "KEY1", "VALUE1");
        info=info_add_node_darray(info, "KEY2", darray, 3);
        info=info_add_node_str(info, "KEY3", "VALUE3");
        info=info_add_node_str(info, "KEY4", "VALUE4");
        Info_Node *node=info_node_where(info, "KEY3");
        res =(strcmp("VALUE3", node->value.string)==0);
        return res;
}

int TEST_info_node_update_str() {
        int res;
        double darray[]={10,20,30};
        Info *info=NULL;
        info=info_add_node_str(info, "KEY1", "VALUE1");
        info=info_add_node_darray(info, "KEY2", darray, 3);
        info=info_add_node_str(info, "KEY3", "VALUE3");
        info=info_add_node_str(info, "KEY4", "VALUE4");
        Info_Node *node=info_node_where(info, "KEY3");
        info_node_update_str(node, "UPDATED");
        res =(strcmp("UPDATED", node->value.string)==0);
        return res;
}

int TEST_info_node_dup() {
        int res;
        Info *info;
        info=info_add_node_str(info, "KEY1", "VALUE1");
        Info_Node *node1=info_node_where(info, "KEY1");
        Info_Node *node2=info_node_dup(node1);

        res = info_node_cmp(node1, node2 );
        return res;
}

int TEST_info_load() {
        int res;
        Info *info=NULL;
        info=info_add_node_str(info, "KEY1", "VALUE1");
        info=info_add_node_str(info, "KEY2", "VALUE2");
        FILE *fp=tmpfile();
        info_save(info, fp);
        fp=freopen(NULL, "r", fp);
        Info *loaded=info_load(fp);
        res = info_cmp(info, loaded );
        return res;
}

int TEST_pdf_initialize() {
        int res;
        Pdf pdf;
        pdf_initialize(&pdf, 100, 100, 6, NULL);
        pdf.val[10][10][1]=1.5;
        res = (
                        pdf.nx==100 && 
                        pdf.nq==100 && 
                        pdf.n_pdf_flavours==6 &&
                        pdf.val[10][10][1]==1.5
                        );
        
        return res;
}

int TEST_save_load_lhapdf6_member() {
        int res;
        int ix, iq, ifl;
        Pdf pdf;
        Info *info=info_add_node_str(NULL, "KEY", "VALUE");
        pdf_initialize(&pdf, 2, 2, 2, info);
        for(ix=0; ix< pdf.nx ; ix++)
        for(iq=0; iq< pdf.nq ; iq++)
        for(ifl=0; ifl< pdf.n_pdf_flavours ; ifl++) pdf.val[ix][iq][ifl]=0.0;
        pdf.val[1][1][1]=1.5;
        pdf.x[0]=0;
        pdf.x[1]=1;
        pdf.q[0]=1.9;
        pdf.q[1]=100.0;
        pdf.pdf_flavours[0]=-1;
        pdf.pdf_flavours[1]=1;

        char *path=tmpnam(NULL);
        save_lhapdf6_member(&pdf, path);
        Pdf loaded;
        load_lhapdf6_member(&loaded, path);
        res =( pdf_cmp(&loaded, &pdf));
        return res;
}

int TEST_pdf_dup() {
        int res;
        int ix, iq, ifl;
        Pdf pdf;
        Info *info=info_add_node_str(NULL, "KEY", "VALUE");
        pdf_initialize(&pdf, 2, 2, 2, info);
        for(ix=0; ix< pdf.nx ; ix++)
        for(iq=0; iq< pdf.nq ; iq++)
        for(ifl=0; ifl< pdf.n_pdf_flavours ; ifl++) pdf.val[ix][iq][ifl]=0.0;
        pdf.val[1][1][1]=1.5;
        pdf.x[0]=0;
        pdf.x[1]=1;
        pdf.q[0]=1.9;
        pdf.q[1]=100.0;
        pdf.pdf_flavours[0]=-1;
        pdf.pdf_flavours[1]=1;

        Pdf *copy=pdf_dup(&pdf);

        res = (pdf_cmp(&pdf, copy));

        return res;
}
int TEST_pdf_cpy() {
        int res;
        int ix, iq, ifl;
        Pdf pdf;
        Info *info=info_add_node_str(NULL, "KEY", "VALUE");
        pdf_initialize(&pdf, 2, 2, 2, info);
        for(ix=0; ix< pdf.nx ; ix++)
        for(iq=0; iq< pdf.nq ; iq++)
        for(ifl=0; ifl< pdf.n_pdf_flavours ; ifl++) pdf.val[ix][iq][ifl]=0.0;
        pdf.val[1][1][1]=1.5;
        pdf.x[0]=0;
        pdf.x[1]=1;
        pdf.q[0]=1.9;
        pdf.q[1]=100.0;
        pdf.pdf_flavours[0]=-1;
        pdf.pdf_flavours[1]=1;

        Pdf copy;
        pdf_initialize(&copy, pdf.nx, pdf.nq, pdf.n_pdf_flavours, pdf.info);
        pdf_cpy(&copy, &pdf);

        res = pdf_cmp(&copy, &pdf);
        return res;
}

int TEST_pdf_set_initialize() {
        int res;
        int ix, iq, ifl;
        PdfSet pdf_set;
        Pdf pdf, pdf2;
        Info *info=NULL, *minfo=NULL;
        info=info_add_node_str(NULL, "KEY", "VALUE");
        minfo=info_add_node_str(NULL, "MKEY", "MVALUE");
        pdf_initialize(&pdf, 2, 2, 2, minfo);
        for(ix=0; ix< pdf.nx ; ix++)
        for(iq=0; iq< pdf.nq ; iq++)
        for(ifl=0; ifl< pdf.n_pdf_flavours ; ifl++) pdf.val[ix][iq][ifl]=0.0;
        pdf.val[1][1][1]=1.5;
        pdf.x[0]=0;
        pdf.x[1]=1;
        pdf.q[0]=1.9;
        pdf.q[1]=100.0;
        pdf.pdf_flavours[0]=-1;
        pdf.pdf_flavours[1]=1;

        pdf_initialize(&pdf2, 2, 2, 2, minfo);
        pdf_cpy(&pdf2, &pdf);
        pdf2.q[0]=4.9;

        pdf_set_initialize(&pdf_set, 2, 2, 2, 2, info, minfo);
        pdf_set.members[0]=pdf;
        pdf_set.members[1]=pdf2;
        res = (pdf_set.members[0].q[0]==1.9 && pdf_set.members[1].q[0]==4.9);
        return res;
}

int TEST_save_load_lhapdf6_set() {
        int res;
        puts("");
        int ix, iq, ifl;
        PdfSet pdf_set;
        Pdf pdf, pdf2;
        Info *info=NULL, *minfo=NULL;
        info=info_add_node_str(NULL, "NumMembers", "2");
        minfo=info_add_node_str(NULL, "MKEY", "MVALUE");
        pdf_initialize(&pdf, 2, 2, 2, minfo);
        for(ix=0; ix< pdf.nx ; ix++)
        for(iq=0; iq< pdf.nq ; iq++)
        for(ifl=0; ifl< pdf.n_pdf_flavours ; ifl++) pdf.val[ix][iq][ifl]=0.0;
        pdf.val[1][1][1]=1.5;
        pdf.x[0]=0;
        pdf.x[1]=1;
        pdf.q[0]=1.9;
        pdf.q[1]=100.0;
        pdf.pdf_flavours[0]=-1;
        pdf.pdf_flavours[1]=1;

        pdf_initialize(&pdf2, pdf.nx, pdf.nq, pdf.n_pdf_flavours, info);
        pdf_cpy(&pdf2, &pdf);
        pdf2.q[0]=4.9;

        pdf_set_initialize(&pdf_set, 2, 2, 2, 2, info, minfo);
        pdf_set.members[0]=pdf;
        pdf_set.members[1]=pdf2;
        char *path=tmpnam(NULL);
        save_lhapdf6_set(&pdf_set, path);
        PdfSet loaded;
        load_lhapdf6_set(&loaded, path);
        
        res =(pdf_set_cmp(&loaded, &pdf_set));
        return res;
}

int TEST_pdf_set_dup() {
        int res;
        int ix, iq, ifl;
        PdfSet pdf_set;
        Pdf pdf, pdf2;
        Info *info=NULL, *minfo=NULL;
        info=info_add_node_str(NULL, "NumMembers", "2");
        minfo=info_add_node_str(NULL, "MKEY", "MVALUE");
        pdf_initialize(&pdf, 2, 2, 2, minfo);
        for(ix=0; ix< pdf.nx ; ix++)
        for(iq=0; iq< pdf.nq ; iq++)
        for(ifl=0; ifl< pdf.n_pdf_flavours ; ifl++) pdf.val[ix][iq][ifl]=0.0;
        pdf.val[1][1][1]=1.5;
        pdf.x[0]=0;
        pdf.x[1]=1;
        pdf.q[0]=1.9;
        pdf.q[1]=100.0;
        pdf.pdf_flavours[0]=-1;
        pdf.pdf_flavours[1]=1;

        pdf_initialize(&pdf2, pdf.nx, pdf.nq, pdf.n_pdf_flavours, info);
        pdf_cpy(&pdf2, &pdf);
        pdf2.q[0]=4.9;

        pdf_set_initialize(&pdf_set, 2, 2, 2, 2, info, minfo);
        pdf_set.members[0]=pdf;
        pdf_set.members[1]=pdf2;

        PdfSet *copy=pdf_set_dup(&pdf_set);


        res = (pdf_set_cmp(copy, &pdf_set));
        return res;
}

//utility function for int/double convertion to char*
int TEST_n2str();

main() {

CHECK(TEST_info_add_node_str);
CHECK(TEST_info_add_node_darray);
CHECK(TEST_info_save); 
CHECK(TEST_info_dup);

CHECK(TEST_info_node_where);
CHECK(TEST_info_node_update_str);
CHECK(TEST_info_node_dup);
CHECK(TEST_info_load);

CHECK(TEST_pdf_initialize);
CHECK(TEST_pdf_dup);
CHECK(TEST_pdf_cpy);
CHECK(TEST_save_load_lhapdf6_member);

CHECK(TEST_pdf_set_initialize);
CHECK(TEST_pdf_set_dup);
CHECK(TEST_save_load_lhapdf6_set);
//utility function for int/double convertion to char*
/* TEST_n2str(); */
return 0;
}
