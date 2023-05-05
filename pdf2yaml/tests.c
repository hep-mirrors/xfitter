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
        fclose(fp);
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
        pdf_initialize(&pdf, NULL);
        pdf_add_grid(&pdf, 200, 100, 6);
        pdf_add_grid(&pdf, 200, 100, 6);
        pdf_add_grid(&pdf, 100, 200, 10);
        pdf.val[2][10][10][1]=1.5;
        res = (
                        pdf.nx[2]==100 && 
                        pdf.nq[2]==200 && 
                        pdf.n_pdf_flavours[2]==10 &&
                        pdf.val[2][10][10][1]==1.5
                        );
        
        return res;
}

int TEST_pdf_initialize_as() {
        int res;
        Pdf pdf;
        Pdf pdf_test;
        pdf_initialize(&pdf, NULL);
        pdf_add_grid(&pdf, 200, 100, 6);
        pdf_add_grid(&pdf, 200, 100, 6);
        pdf_add_grid(&pdf, 100, 200, 10);
        pdf.val[2][10][10][1]=1.5;
        
        pdf_initialize_as(&pdf_test, &pdf);
        pdf.val[2][10][10][1]=1.5;

        res = (
                        pdf.nx[2]==100 && 
                        pdf.nq[2]==200 && 
                        pdf.n_pdf_flavours[2]==10 &&
                        pdf.val[2][10][10][1]==1.5
                        );
        return res;
}

int TEST_pdf_dup() {
        int res;
        int ix, iq, ifl;
        Pdf pdf;
        Info *info=info_add_node_str(NULL, "KEY", "VALUE");
        pdf_initialize(&pdf, info);
        pdf_add_grid(&pdf, 2, 2, 2);
        for(ix=0; ix< pdf.nx[0] ; ix++)
        for(iq=0; iq< pdf.nq[0] ; iq++)
        for(ifl=0; ifl< pdf.n_pdf_flavours[0] ; ifl++) pdf.val[0][ix][iq][ifl]=0.0;
        pdf.val[0][1][1][1]=1.5;
        pdf.x[0][0]=0;
        pdf.x[0][1]=1;
        pdf.q[0][0]=1.9;
        pdf.q[0][1]=100.0;
        pdf.pdf_flavours[0][0]=-1;
        pdf.pdf_flavours[0][1]=1;

        pdf_add_grid(&pdf, 2, 2, 2);
        for(ix=0; ix< pdf.nx[1] ; ix++)
        for(iq=0; iq< pdf.nq[1] ; iq++)
        for(ifl=0; ifl< pdf.n_pdf_flavours[1] ; ifl++) pdf.val[1][ix][iq][ifl]=0.0;
        pdf.val[1][1][1][1]=1.8;
        pdf.x[1][0]=9;
        pdf.x[1][1]=1;
        pdf.q[1][0]=1.3;
        pdf.q[1][1]=10.0;
        pdf.pdf_flavours[1][0]=-2;
        pdf.pdf_flavours[1][1]=2;

        Pdf *copy=pdf_dup(&pdf);

        res = (pdf_cmp(&pdf, copy));

        return res;
}
int TEST_save_load_lhapdf6_member() {
        int res;
        int ix, iq, ifl;
        Pdf pdf;
        Info *info=info_add_node_str(NULL, "KEY", "VALUE");
        pdf_initialize(&pdf, info);
        pdf_add_grid(&pdf,2,2,2);
        for(ix=0; ix< pdf.nx[0] ; ix++)
        for(iq=0; iq< pdf.nq[0] ; iq++)
        for(ifl=0; ifl< pdf.n_pdf_flavours[0] ; ifl++) pdf.val[0][ix][iq][ifl]=0.0;
        pdf.val[0][1][1][1]=1.5;
        pdf.x[0][0]=0;
        pdf.x[0][1]=1;
        pdf.q[0][0]=1.9;
        pdf.q[0][1]=100.0;
        pdf.pdf_flavours[0][0]=-1;
        pdf.pdf_flavours[0][1]=1;

        pdf_add_grid(&pdf,2,2,2);
        for(ix=0; ix< pdf.nx[1] ; ix++)
        for(iq=0; iq< pdf.nq[1] ; iq++)
        for(ifl=0; ifl< pdf.n_pdf_flavours[1] ; ifl++) pdf.val[1][ix][iq][ifl]=1.0;
        pdf.val[1][1][1][1]=1.5;
        pdf.x[1][0]=0;
        pdf.x[1][1]=1;
        pdf.q[1][0]=1.9;
        pdf.q[1][1]=110.0;
        pdf.pdf_flavours[1][0]=-2;
        pdf.pdf_flavours[1][1]=2;

	//        char *path=tmpnam(NULL);
	char path[] = "/tmp/tempfileXXXXXX";
	int fd = mkstemp(path);
	if (fd == -1) {
	  perror("mkstemp");
	  return 1;
	}
        save_lhapdf6_member(&pdf, path);
        Pdf loaded;
        load_lhapdf6_member(&loaded, path);
        res =( pdf_cmp(&loaded, &pdf));

        pdf_add_grid(&loaded,2,3,3);
        for(ix=0; ix< loaded.nx[2] ; ix++)
        for(iq=0; iq< loaded.nq[2] ; iq++)
        for(ifl=0; ifl< loaded.n_pdf_flavours[2] ; ifl++) loaded.val[2][ix][iq][ifl]=1.0;
        loaded.val[2][1][1][1]=1.5;
        loaded.x[2][0]=0;
        loaded.x[2][1]=1;
        loaded.q[2][0]=1.9;
        loaded.q[2][1]=110.0;
        loaded.q[2][2]=150.0;
        loaded.pdf_flavours[2][0]=-2;
        loaded.pdf_flavours[2][1]=2;
        loaded.pdf_flavours[2][2]=4;

        save_lhapdf6_member(&loaded, "/tmp/PDF");
        return res;
}

int TEST_pdf_cpy() {
        int res;
        int ix, iq, ifl;
        Pdf pdf;
        Pdf copy;

        Info *info=info_add_node_str(NULL, "KEY", "VALUE");
        pdf_initialize(&pdf, info);
        pdf_add_grid(&pdf,2,2,2);
        for(ix=0; ix< pdf.nx[0] ; ix++)
        for(iq=0; iq< pdf.nq[0] ; iq++)
        for(ifl=0; ifl< pdf.n_pdf_flavours[0] ; ifl++) pdf.val[0][ix][iq][ifl]=0.0;
        pdf.val[0][1][1][1]=1.5;
        pdf.x[0][0]=0;
        pdf.x[0][1]=1;
        pdf.q[0][0]=1.9;
        pdf.q[0][1]=100.0;
        pdf.pdf_flavours[0][0]=-1;
        pdf.pdf_flavours[0][1]=1;

        pdf_add_grid(&pdf,2,2,2);
        for(ix=0; ix< pdf.nx[1] ; ix++)
        for(iq=0; iq< pdf.nq[1] ; iq++)
        for(ifl=0; ifl< pdf.n_pdf_flavours[1] ; ifl++) pdf.val[1][ix][iq][ifl]=1.0;
        pdf.val[1][1][1][1]=1.5;
        pdf.x[1][0]=0;
        pdf.x[1][1]=1;
        pdf.q[1][0]=1.9;
        pdf.q[1][1]=110.0;
        pdf.pdf_flavours[1][0]=-2;
        pdf.pdf_flavours[1][1]=2;
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

        pdf_set_initialize(&pdf_set, 2, info);

        pdf_initialize(&pdf_set.members[0], minfo);
        pdf_add_grid(&pdf_set.members[0], 2, 3, 3);
        for(ix=0; ix< pdf_set.members[0].nx[0] ; ix++)
        for(iq=0; iq< pdf_set.members[0].nq[0] ; iq++)
        for(ifl=0; ifl< pdf_set.members[0].n_pdf_flavours[0] ; ifl++) pdf_set.members[0].val[0][ix][iq][ifl]=1.0;
        pdf_set.members[0].val[0][1][1][1]=1.5;
        pdf_set.members[0].x[0][0]=0;
        pdf_set.members[0].x[0][1]=1;
        pdf_set.members[0].q[0][0]=1.9;
        pdf_set.members[0].q[0][1]=110.0;
        pdf_set.members[0].q[0][2]=150.0;
        pdf_set.members[0].pdf_flavours[0][0]=-2;
        pdf_set.members[0].pdf_flavours[0][1]=2;
        pdf_set.members[0].pdf_flavours[0][2]=4;

        pdf_cpy(&pdf_set.members[1],&pdf_set.members[0]);

        res=(pdf_set.members[0].val[0][1][1][1]==1.5);
        return res;
}

int TEST_pdf_set_initialize_as() {
        int res;
        int ix, iq, ifl;
        PdfSet pdf_set;
        Pdf pdf, pdf2;
        Info *info=NULL, *minfo=NULL;
        info=info_add_node_str(NULL, "KEY", "VALUE");
        minfo=info_add_node_str(NULL, "MKEY", "MVALUE");

        pdf_set_initialize(&pdf_set, 2, info);

        pdf_initialize(&pdf_set.members[0], minfo);
        pdf_add_grid(&pdf_set.members[0], 2, 3, 3);
        pdf_cpy(&pdf_set.members[1],&pdf_set.members[0]);
        PdfSet new_set;
        pdf_set_initialize_as(&new_set, &pdf_set);
        for(ix=0; ix< new_set.members[0].nx[0] ; ix++)
        for(iq=0; iq< new_set.members[0].nq[0] ; iq++)
        for(ifl=0; ifl< new_set.members[0].n_pdf_flavours[0] ; ifl++) new_set.members[0].val[0][ix][iq][ifl]=1.0;
        new_set.members[0].val[0][1][1][1]=1.5;
        new_set.members[0].x[0][0]=0;
        new_set.members[0].x[0][1]=1;
        new_set.members[0].q[0][0]=1.9;
        new_set.members[0].q[0][1]=110.0;
        new_set.members[0].q[0][2]=150.0;
        new_set.members[0].pdf_flavours[0][0]=-2;
        new_set.members[0].pdf_flavours[0][1]=2;
        new_set.members[0].pdf_flavours[0][2]=4;


        res=(new_set.members[0].val[0][1][1][1]==1.5);
        return res;
}

int TEST_save_load_lhapdf6_set() {
        int res;
        puts("");
        int ix, iq, ifl;
        PdfSet pdf_set;
        Info *info=NULL, *minfo=NULL;
        info=info_add_node_str(NULL, "NumMembers", "2");
        minfo=info_add_node_str(NULL, "MKEY", "MVALUE");

        pdf_set_initialize(&pdf_set, 2, info);

        pdf_initialize(&pdf_set.members[0], minfo);
        pdf_add_grid(&pdf_set.members[0], 2, 3, 3);
        for(ix=0; ix< pdf_set.members[0].nx[0] ; ix++)
        for(iq=0; iq< pdf_set.members[0].nq[0] ; iq++)
        for(ifl=0; ifl< pdf_set.members[0].n_pdf_flavours[0] ; ifl++) pdf_set.members[0].val[0][ix][iq][ifl]=1.0;
        pdf_set.members[0].val[0][1][1][1]=1.5;
        pdf_set.members[0].x[0][0]=0;
        pdf_set.members[0].x[0][1]=1;
        pdf_set.members[0].q[0][0]=1.9;
        pdf_set.members[0].q[0][1]=110.0;
        pdf_set.members[0].q[0][2]=150.0;
        pdf_set.members[0].pdf_flavours[0][0]=-2;
        pdf_set.members[0].pdf_flavours[0][1]=2;
        pdf_set.members[0].pdf_flavours[0][2]=4;

        pdf_add_grid(&pdf_set.members[0], 2, 2, 2);
        for(ix=0; ix< pdf_set.members[0].nx[1] ; ix++)
        for(iq=0; iq< pdf_set.members[0].nq[1] ; iq++)
        for(ifl=0; ifl< pdf_set.members[0].n_pdf_flavours[1] ; ifl++) pdf_set.members[0].val[1][ix][iq][ifl]=1.0;
        pdf_set.members[0].val[1][1][1][1]=5.5;
        pdf_set.members[0].x[1][0]=0;
        pdf_set.members[0].x[1][1]=1;
        pdf_set.members[0].q[1][0]=2.9;
        pdf_set.members[0].q[1][1]=150.0;
        pdf_set.members[0].pdf_flavours[1][0]=-1;
        pdf_set.members[0].pdf_flavours[1][1]=1;

        pdf_cpy(&pdf_set.members[1],&pdf_set.members[0]);
        for(ix=0; ix< pdf_set.members[0].nx[0] ; ix++)
        for(iq=0; iq< pdf_set.members[0].nq[0] ; iq++)
        for(ifl=0; ifl< pdf_set.members[0].n_pdf_flavours[0] ; ifl++) pdf_set.members[1].val[0][ix][iq][ifl]=3.3;

        //char *path=tmpnam(NULL);
	char path[] = "/tmp/tempfileXXXXXX";
        int fd = mkstemp(path);
        if (fd == -1) {
          perror("mkstemp");
          return 1;
        }
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

        Info *info=NULL, *minfo=NULL;
        info=info_add_node_str(NULL, "NumMembers", "2");
        minfo=info_add_node_str(NULL, "MKEY", "MVALUE");

        pdf_set_initialize(&pdf_set, 2, info);

        pdf_initialize(&pdf_set.members[0], minfo);
        pdf_add_grid(&pdf_set.members[0], 2, 3, 3);
        for(ix=0; ix< pdf_set.members[0].nx[0] ; ix++)
        for(iq=0; iq< pdf_set.members[0].nq[0] ; iq++)
        for(ifl=0; ifl< pdf_set.members[0].n_pdf_flavours[0] ; ifl++) pdf_set.members[0].val[0][ix][iq][ifl]=1.0;
        pdf_set.members[0].val[0][1][1][1]=1.5;
        pdf_set.members[0].x[0][0]=0;
        pdf_set.members[0].x[0][1]=1;
        pdf_set.members[0].q[0][0]=1.9;
        pdf_set.members[0].q[0][1]=110.0;
        pdf_set.members[0].q[0][2]=150.0;
        pdf_set.members[0].pdf_flavours[0][0]=-2;
        pdf_set.members[0].pdf_flavours[0][1]=2;
        pdf_set.members[0].pdf_flavours[0][2]=4;

        pdf_add_grid(&pdf_set.members[0], 2, 2, 2);
        for(ix=0; ix< pdf_set.members[0].nx[1] ; ix++)
        for(iq=0; iq< pdf_set.members[0].nq[1] ; iq++)
        for(ifl=0; ifl< pdf_set.members[0].n_pdf_flavours[1] ; ifl++) pdf_set.members[0].val[1][ix][iq][ifl]=1.0;
        pdf_set.members[0].val[1][1][1][1]=5.5;
        pdf_set.members[0].x[1][0]=0;
        pdf_set.members[0].x[1][1]=1;
        pdf_set.members[0].q[1][0]=2.9;
        pdf_set.members[0].q[1][1]=150.0;
        pdf_set.members[0].pdf_flavours[1][0]=-1;
        pdf_set.members[0].pdf_flavours[1][1]=1;

        PdfSet *copy=pdf_set_dup(&pdf_set);

        res = (pdf_set_cmp(copy, &pdf_set));
        return res;
}

int TEST_each_in_pdf() {
        int res;
        int ig, ix, iq, ifl;

        Pdf pdf;
        Info *info=info_add_node_str(NULL, "KEY", "VALUE");
        pdf_initialize(&pdf, info);
        pdf_add_grid(&pdf, 3, 3, 3);
        pdf_add_grid(&pdf, 10, 10, 3);
        pdf_add_grid(&pdf, 100, 5, 3);
        for(ig=0; ig<pdf.ngrids; ig++) {
                for(ix=0; ix<pdf.nx[ig]; ix++) pdf.x[ig][ix]=ix*1./(pdf.nx[ig]);
                for(iq=0; iq<pdf.nq[ig]; iq++) pdf.q[ig][iq]=iq*10.0;
                for(ifl=0; ifl<pdf.n_pdf_flavours[ig]; ifl++) pdf.pdf_flavours[ig][ifl]=-1+ifl;
        }

        EACH_IN_PDF(&pdf, ig, ix, iq, ifl) {
                pdf.val[ig][ix][iq][ifl]=1.0;
        }
        Pdf *orig=pdf_dup(&pdf);
        EACH_IN_PDF(&pdf, ig, ix, iq, ifl) {
                pdf.val[ig][ix][iq][ifl]+=1.0;
        }
        
        EACH_IN_PDF(&pdf, ig, ix, iq, ifl) {
                pdf.val[ig][ix][iq][ifl]/=2.0;
        }

        res=pdf_cmp(orig, &pdf);
        return res;
}
int TEST_each_in_set() {
        int res;
        int im, ig, ix, iq, ifl;

        PdfSet pdf_set;
        Info *info=info_add_node_str(NULL, "KEY", "VALUE");
        Info *minfo=info_add_node_str(NULL, "MKEY", "MVALUE");
        
        pdf_set_initialize(&pdf_set, 2, info);
        Pdf *pdf=&pdf_set.members[0];
        pdf_initialize(pdf, info);
        pdf_add_grid(pdf, 3, 3, 3);
        pdf_add_grid(pdf, 10, 10, 3);
        pdf_add_grid(pdf, 100, 5, 3);
        for(ig=0; ig<pdf->ngrids; ig++) {
                for(ix=0; ix<pdf->nx[ig]; ix++) pdf->x[ig][ix]=ix*1./(pdf->nx[ig]);
                for(iq=0; iq<pdf->nq[ig]; iq++) pdf->q[ig][iq]=iq*10.0;
                for(ifl=0; ifl<pdf->n_pdf_flavours[ig]; ifl++) pdf->pdf_flavours[ig][ifl]=-1+ifl;
        }

        EACH_IN_PDF(pdf, ig, ix, iq, ifl) {
                pdf->val[ig][ix][iq][ifl]=1.0;
        }
        
        pdf_cpy(&pdf_set.members[1], pdf);

        PdfSet *orig=pdf_set_dup(&pdf_set);

        EACH_IN_SET(&pdf_set, im, ig, ix, iq, ifl) pdf_set.members[im].val[ig][ix][iq][ifl]*=2;
        EACH_IN_SET(&pdf_set, im, ig, ix, iq, ifl) pdf_set.members[im].val[ig][ix][iq][ifl]-=1;
        res=pdf_set_cmp(&pdf_set, orig);
        return res;
}
//utility function for int/double convertion to char*
int TEST_n2str();

int main() {
CHECK(TEST_info_add_node_str);
CHECK(TEST_info_add_node_darray);
CHECK(TEST_info_save); 
CHECK(TEST_info_dup);

CHECK(TEST_info_node_where);
CHECK(TEST_info_node_update_str);
CHECK(TEST_info_node_dup);
CHECK(TEST_info_load);

CHECK(TEST_pdf_initialize);
CHECK(TEST_pdf_initialize_as);
CHECK(TEST_pdf_dup);
CHECK(TEST_pdf_cpy);
CHECK(TEST_save_load_lhapdf6_member);
//
CHECK(TEST_pdf_set_initialize);
CHECK(TEST_pdf_set_initialize_as);
CHECK(TEST_pdf_set_dup);
CHECK(TEST_save_load_lhapdf6_set);
CHECK(TEST_each_in_pdf);
CHECK(TEST_each_in_set);

//utility function for int/double convertion to char*
/* TEST_n2str(); */
return 0;
}
