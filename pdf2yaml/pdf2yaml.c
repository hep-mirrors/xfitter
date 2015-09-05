#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include <sys/stat.h>
#include <float.h>
#include <yaml.h>
#include "c2yaml.h"
#include "list.h"
#include "pdf2yaml.h"
//#define DBL_CMP(x, y) (fabs(x-y) < 10 * FLT_EPSILON) 
#define DBL_CMP(x, y) ( fabs(x-y) <= ((fabs(x) > fabs(y)) ? fabs(x) : fabs(y))*25.0*FLT_EPSILON  )
// Pdf info
static Info *parse_yaml_seq(Info *info, yaml_document_t *doc, yaml_node_t *mkey, yaml_node_t *mval);
static Info *parse_yaml_map(Info *info, yaml_document_t *doc, yaml_node_t *node);

static Info *parse_yaml_map(Info *info, yaml_document_t *doc, yaml_node_t *node) {
        yaml_node_t *map_key, *map_value;
        yaml_node_pair_t *mpair;
        mpair=node->data.mapping.pairs.start;
        while(mpair!=node->data.mapping.pairs.end && mpair!=node->data.mapping.pairs.top) {
                map_key=yaml_document_get_node(doc,mpair->key);
                map_value=yaml_document_get_node(doc,mpair->value);
                if(map_value->type==YAML_SCALAR_NODE)
                        info=info_add_node_str(info,
                                        map_key->data.scalar.value,
                                        map_value->data.scalar.value);
                if(map_value->type==YAML_SEQUENCE_NODE) info=parse_yaml_seq(info, doc, map_key, map_value);
                mpair++;
        }
        return info;
}

static Info *parse_yaml_seq( Info *info, yaml_document_t *doc, yaml_node_t *mkey, yaml_node_t *mval) {
        int seq_size=0,i=0;
        yaml_node_item_t *seq_it, *top, *end, *start;
        yaml_node_t *seq_node;
        start= mval->data.sequence.items.start;
        top=mval->data.sequence.items.top;
        end=mval->data.sequence.items.end;
        for(seq_size=0, seq_it=start; seq_it!=end && seq_it!= top; seq_it++) seq_size++;
        double *darray=calloc(seq_size, sizeof(double));
        seq_it=start;
        while(seq_it!=end && seq_it!= top) {
                seq_node=yaml_document_get_node(doc, *seq_it);
                sscanf(seq_node->data.scalar.value, "%lg", darray+(i++));
                seq_it++;
        }
        info=info_add_node_darray(info,
                        mkey->data.scalar.value,
                        darray,seq_size);
        free(darray);
        return info;
}

Info *info_load(FILE *input) {
        Info *info=NULL;
        yaml_parser_t parser;
        yaml_document_t doc;
        yaml_node_t *node, *map_key, *map_value;

        if(!yaml_parser_initialize(&parser)) fputs("Failed to initialize parser!\n", stderr);

        yaml_parser_set_input_file(&parser, input);
        yaml_parser_load(&parser, &doc);
        for(node=doc.nodes.start; node!=doc.nodes.end && node!=doc.nodes.top; node++) {
                if(node->type==YAML_MAPPING_NODE) info=parse_yaml_map(info,&doc,node); 

        }

        yaml_parser_delete(&parser);
        yaml_document_delete(&doc);
        return info;
}


Info *info_add_node_str(Info *info, char* key, char *value) {
        Info_Node *node= malloc(sizeof(Info_Node));
        node->node_type=STRING;
        node->key=malloc(strlen(key)+1);
        strncpy(node->key, key, strlen(key)+1);
        node->value.string=malloc(strlen(value)+1);
        strncpy(node->value.string, value, strlen(value)+1);
        info=list_append(info, node);
        return info;
}

Info *info_add_node_darray(Info *info, char* key, double *darray, int size) {
        Info_Node *node= malloc(sizeof(Info_Node));
        node->node_type=DARRAY;
        node->key=malloc(strlen(key)+1);
        strncpy(node->key, key, strlen(key)+1);
        node->value.darray.vals=malloc(sizeof(double)*size);
        memcpy(node->value.darray.vals, darray, sizeof(double)*size);
        node->value.darray.size=size;
        info=list_append(info, node);
        return info;
}

void info_save(const Info *info, FILE *output) {
        yaml_document_t doc;
        yaml_document_initialize(&doc,NULL,NULL,NULL,1,1);
        yaml_emitter_t emitter;
        Info *it,*r;

        int root = yaml_document_add_mapping(&doc, NULL, YAML_BLOCK_MAPPING_STYLE);

        Info_Node *node;

        if(!yaml_emitter_initialize(&emitter)) 
                fputs("Failed to initialize emitter!\n", stderr);

        yaml_emitter_set_output_file(&emitter,output);

        r=list_reverse(info);
        it=r;
        while(it) {
                node=it->data;
                
                if(node->node_type==STRING)
                yaml_document_append_mapping_pair(&doc, root, s2yaml(&doc,node->key), 
                                s2yaml(&doc,node->value.string));
                if(node->node_type==DARRAY)
                yaml_document_append_mapping_pair(&doc, root, s2yaml(&doc,node->key), 
                                d_array2yaml(&doc,node->value.darray.vals,node->value.darray.size));

                it=it->next;
        }
        list_free(r);
        yaml_emitter_dump(&emitter,&doc);
        yaml_emitter_delete(&emitter);
        yaml_document_delete(&doc);
}

Info_Node *info_node_where(Info *info, char *key) {
        Info_Node *info_node;
        while(info) {
                info_node=info->data;
                if(!strcmp(info_node->key, key)) return info_node;
                info=info->next;
        }
        return NULL;
}

int info_node_update_str(Info_Node *node, const char *new_value) {
        if(node->node_type!=STRING) return 1;
        free(node->value.string);
        node->value.string=strdup(new_value);
        return 0;
}

Info_Node *info_node_dup(Info_Node *node) {
        int size;
        Info_Node *new_node= malloc(sizeof(Info_Node));
        new_node->key=strdup(node->key);
        if(node->node_type==STRING) { 
                new_node->node_type=STRING;
                new_node->value.string=strdup(node->value.string);
        }
        if(node->node_type==DARRAY) {
                new_node->node_type=DARRAY;
                size=node->value.darray.size;
                new_node->value.darray.size=size;
                new_node->value.darray.vals=malloc(sizeof(double)*size);
                memcpy(new_node->value.darray.vals, node->value.darray.vals, sizeof(double)*size);
        }
        return new_node;
}

int info_node_cmp(Info_Node *node1, Info_Node *node2) {
        if(node1==node2) return 1;
        double x,y;
        int res=1;
        int i;
        if(node1->node_type!=node2->node_type) return 0;
        if(node1->node_type==STRING) {
                res*=!strcmp(node1->key, node2->key);
                res*=!strcmp(node1->value.string, node2->value.string); 
        }

        if(node1->node_type==DARRAY) {
                res*=!strcmp(node1->key, node2->key);
                if(node1->value.darray.size!=node2->value.darray.size) return 0;
                for(i=0; i<node1->value.darray.size; i++) {
                        x=node1->value.darray.vals[i];
                        y=node2->value.darray.vals[i];
                        res*=(DBL_CMP(node1->value.darray.vals[i], node2->value.darray.vals[i]));
                        if(!res) {
                                fprintf(stderr, "Error pdf_cmp(info %s): %1.8e %1.8e\n", node1->key,
                                                node1->value.darray.vals[i]-node2->value.darray.vals[i], FLT_EPSILON);
                                return 0;
                        }
                }
        }

        return res;
}

Info *info_dup(const Info *info) {
        Info_Node *node;
        Info *new_info=NULL;
        Info *r=(Info*)list_reverse(info);
        Info *it;
        for(it=r; it; it=it->next) {
                node=it->data;
                node=info_node_dup(node);
                new_info=(Info*)list_append(new_info, node);
        }
        list_free(r);
        return new_info;
}

int info_cmp(Info *info1, Info *info2) {
        while(1) {
                if(info1==info2) return 1;
                if(info1==NULL || info2==NULL) return 0;
                if(!info_node_cmp(info1->data, info2->data)) return 0;
                info1=info1->next;
                info2=info2->next;
        }
        return 0;
}

void info_free(Info *info) {
        Info *it;
        Info_Node *node;
        for(it=info; it; it=it->next) {
                node=it->data;
                free(node->key);
                if(node->node_type==STRING) free(node->value.string);
                if(node->node_type==DARRAY) free(node->value.darray.vals);
                free(node);
        }
        list_free(info);
}

//                     Pdf member
void pdf_add_grid(Pdf *pdf, int nx, int nq, int n_pdf_flavours) {
        int i, j, k;
        int ngrids=pdf->ngrids;
        pdf->nx=realloc(pdf->nx, (ngrids+1)*sizeof(int));
        pdf->nx[ngrids]=nx;
        pdf->x=realloc(pdf->x, (ngrids+1)*sizeof(double*));
        pdf->x[ngrids]=malloc(sizeof(double)*nx);

        pdf->nq=realloc(pdf->nq, (ngrids+1)*sizeof(int));
        pdf->nq[ngrids]=nq;
        pdf->q=realloc(pdf->q, (ngrids+1)*sizeof(double*));
        pdf->q[ngrids]=malloc(sizeof(double)*nq);

        pdf->n_pdf_flavours=realloc(pdf->n_pdf_flavours, (ngrids+1)*sizeof(int));
        pdf->n_pdf_flavours[ngrids]=n_pdf_flavours;
        pdf->pdf_flavours=realloc(pdf->pdf_flavours, (ngrids+1)*sizeof(int*));
        pdf->pdf_flavours[ngrids]=malloc(sizeof(int)*n_pdf_flavours);

        pdf->val=realloc(pdf->val, sizeof(double***)*(ngrids+1));
        pdf->val[ngrids]=malloc(sizeof(double **)*pdf->nx[ngrids]);
        for(i=0; i<pdf->nx[ngrids]; i++) {
                pdf->val[ngrids][i]=malloc(sizeof(double*)*pdf->nq[ngrids]);
                for(j=0; j<pdf->nq[ngrids]; j++) pdf->val[ngrids][i][j]=malloc(sizeof(double)*pdf->n_pdf_flavours[ngrids]);
        }

        pdf->ngrids++;
}

// create empty Pdf member
void pdf_initialize(Pdf *pdf, Info *info) {
        pdf->info=info;
        pdf->ngrids=0;
        pdf->nx=NULL;
        pdf->nq=NULL;
        pdf->n_pdf_flavours=NULL;
        pdf->x=NULL;
        pdf->q=NULL;
        pdf->pdf_flavours=NULL;
        pdf->val=NULL;
}

//initialize new Pdf, copy grid and info from *src and allocate memmory for pdf.val;
void pdf_initialize_as(Pdf *pdf, const Pdf *src) {
        int ig, i, j;
        Info *info=info_dup(src->info);
        pdf_initialize(pdf, info);
        for(ig=0; ig<src->ngrids; ig++) {
                pdf_add_grid(pdf, src->nx[ig], src->nq[ig], src->n_pdf_flavours[ig]);
                memcpy(pdf->x[ig], src->x[ig], sizeof(double)*src->nx[ig]);
                memcpy(pdf->q[ig], src->q[ig], sizeof(double)*src->nq[ig]);
                memcpy(pdf->pdf_flavours[ig], src->pdf_flavours[ig], sizeof(int)*src->n_pdf_flavours[ig]);
        }
}

//initialize new Pdf, and make a deep copy from *src
int pdf_cpy(Pdf *dest, const Pdf *src) {
        int ig, i, j;
        Info *info=info_dup(src->info);
        pdf_initialize(dest, info);
        for(ig=0; ig<src->ngrids; ig++) {
                pdf_add_grid(dest, src->nx[ig], src->nq[ig], src->n_pdf_flavours[ig]);
                memcpy(dest->x[ig], src->x[ig], sizeof(double)*src->nx[ig]);
                memcpy(dest->q[ig], src->q[ig], sizeof(double)*src->nq[ig]);
                memcpy(dest->pdf_flavours[ig], src->pdf_flavours[ig], sizeof(int)*src->n_pdf_flavours[ig]);
                for(i=0; i<src->nx[ig]; i++) 
                        for(j=0; j<src->nq[ig]; j++)
                                memcpy(dest->val[ig][i][j], src->val[ig][i][j], sizeof(double)*src->n_pdf_flavours[ig]);
        }
        

        return 0;
}

Pdf *pdf_dup(Pdf *pdf) {
        int i, j, k;
        Pdf *new_pdf=malloc(sizeof(Pdf));
        Info *info=info_dup(pdf->info);
        pdf_initialize(new_pdf, info);
        for(i=0; i< pdf->ngrids; i++) {
                pdf_add_grid(new_pdf, pdf->nx[i], pdf->nq[i], pdf->n_pdf_flavours[i]);

                memcpy(new_pdf->x[i], pdf->x[i], sizeof(double)*pdf->nx[i]);
                memcpy(new_pdf->q[i], pdf->q[i], sizeof(double)*pdf->nq[i]);
                memcpy(new_pdf->pdf_flavours[i], pdf->pdf_flavours[i], sizeof(int)*pdf->n_pdf_flavours[i]);

                for(j=0; j<pdf->nx[i]; j++) 
                        for(k=0; k<pdf->nq[i]; k++)
                                memcpy(new_pdf->val[i][j][k], pdf->val[i][j][k], sizeof(double)*pdf->n_pdf_flavours[i]);
                
        }
        return new_pdf;
}

int pdf_cmp(Pdf *pdf1, Pdf *pdf2) {
        int ig, ix, iq, ifl;
        int res=1;
        if(pdf1->ngrids!=pdf2->ngrids) {
                fprintf(stderr, "Error pdf_cmp(ngrids): %d %d\n", pdf1->ngrids, pdf2->ngrids);
                return 0;
        }
        for(ig=0; ig<pdf1->ngrids; ig++) {

        if(pdf1->nx[ig]!=pdf2->nx[ig]) {
                fprintf(stderr, "Error pdf_cmp(nx): %d %d\n", pdf1->nx[ig], pdf2->nx[ig] );
                return 0;
        }

        if(pdf1->nq[ig]!=pdf2->nq[ig]) {
                fprintf(stderr, "Error pdf_cmp(nq): %d %d\n", pdf1->nq[ig], pdf2->nq[ig] ); 
                return 0;
        }

        if(pdf1->n_pdf_flavours[ig]!=pdf2->n_pdf_flavours[ig]) { 
                fprintf(stderr, "Error pdf_cmp(n_pdf_flavours): %d %d\n", 
                                pdf1->n_pdf_flavours[ig], pdf2->n_pdf_flavours[ig]);
                return 0;
        }

        if(!info_cmp(pdf1->info, pdf2->info)) {
                fputs("Error pdf_cmp(info)\n", stderr);
                return 0;
        }
        
        for(ix=0; ix< pdf1->nx[ig]; ix++) { 
                res*=(fabs(pdf1->x[ig][ix]-pdf2->x[ig][ix])<FLT_EPSILON);
                if(!res) { 
                        fprintf(stderr, "Error pdf_cmp(x): %g %g\n", pdf1->x[ig][ix], pdf2->x[ig][ix]);
                        return 0;
                }
        }

        for(iq=0; iq< pdf1->nq[ig]; iq++) {
                res*=(fabs(pdf1->q[ig][iq]-pdf2->q[ig][iq])<FLT_EPSILON);
                if(!res) { 
                        fprintf(stderr, "Error pdf_cmp(q): %g %g\n", pdf1->q[ig][iq], pdf2->q[ig][iq]);
                        return 0;
                }
        }

        for(ifl=0; ifl< pdf1->n_pdf_flavours[ig]; ifl++) {
                res*=(pdf1->pdf_flavours[ig][ifl]==pdf2->pdf_flavours[ig][ifl]);
                if(!res) { 
                        fprintf(stderr, "Error pdf_cmp(pdf_flavours): %d %d\n", 
                                        pdf1->pdf_flavours[ig][ifl], pdf2->pdf_flavours[ig][ifl]);
                        return 0;
                }
        }
        if(!res) return 0;

        for(ix=0; ix< pdf1->nx[ig]; ix++) 
        for(iq=0; iq< pdf1->nq[ig]; iq++) 
        for(ifl=0; ifl< pdf1->n_pdf_flavours[ig]; ifl++) {
                res*=(fabs(pdf1->val[ig][ix][iq][ifl]-pdf2->val[ig][ix][iq][ifl])<FLT_EPSILON);
                if(!res) { 
                        fprintf(stderr, "Error pdf_cmp(val): %g %g\n", 
                                        pdf1->val[ig][ix][iq][ifl], pdf2->val[ig][ix][iq][ifl]);
                        return 0;
                }
        }
        }

        return res;
}
// load Pdf member from file. It does not require pdf_initialize for *pdf argument
int load_lhapdf6_member(Pdf *pdf, char *path) {
        double d_tmp;
        int i_tmp, i, j, k;
        size_t size=0, b_len=0;
        char *line=NULL;
        FILE *s_stream=NULL;
        FILE *b_stream=NULL;

        FILE *fp=fopen(path,"r");
        if(!fp) { fprintf(stderr, "cannot open file %s", path); return 1; }

        pdf_initialize(pdf, NULL);

        // load info 
        pdf->info=info_load(fp);

        rewind(fp);
        do {
                free(line);
                line=NULL;
                size=0;
        } while(i_tmp=getline(&line, &size, fp), strcmp(line,"---\n"));

        free(line);
        line=NULL;
        size=0;

        // read grids
        int ngrids=0;
        while(1) {

        // read x line / finalize at the end of file
        if(getline(&line, &size, fp)==-1) { 
                fclose(fp);
                if(pdf->ngrids==0) {fprintf(stderr, "Error in pdf member format, can't read x line\n"); return EXIT_FAILURE;}
                pdf->ngrids=ngrids;
                return 0;
        }
        pdf_add_grid(pdf, 0,0,0);
        s_stream=fmemopen(line,size,"r"); 
        pdf->x[ngrids]=NULL;
        b_stream=open_memstream((char**)&pdf->x[ngrids], &b_len); 
        while(fscanf(s_stream, "%lg", &d_tmp)) { fwrite(&d_tmp, sizeof(double), 1, b_stream); pdf->nx[ngrids]++; }
        fclose(s_stream);
        fclose(b_stream);

        //read q line
        free(line);
        line=NULL;
        size=0;

        if(!getline(&line, &size, fp)) { fprintf(stderr, "Error in pdf member format, cant't read q line\n"); return 1;}
        s_stream=fmemopen(line,size,"r"); 
        pdf->q[ngrids]=NULL;
        b_stream=open_memstream((char**)&pdf->q[ngrids], &b_len); 
        while(fscanf(s_stream, "%lg", &d_tmp)) { fwrite(&d_tmp, sizeof(double), 1, b_stream); pdf->nq[ngrids]++; }
        fclose(s_stream);
        fclose(b_stream);

        //read pdf_flavours line
        free(line);
        line=NULL;
        size=0;

        if(!getline(&line, &size, fp)) { fprintf(stderr, "Error in pdf member format, can't read flavours line\n"); return 1;}
        s_stream=fmemopen(line,size,"r"); 
        pdf->pdf_flavours[ngrids]=NULL;
        b_stream=open_memstream((char**)&pdf->pdf_flavours[ngrids], &b_len); 
        while(fscanf(s_stream, "%d", &i_tmp)) { fwrite(&i_tmp, sizeof(int), 1, b_stream); pdf->n_pdf_flavours[ngrids]++; }
        fclose(s_stream);
        fclose(b_stream);

        free(line);
        line=NULL;
        size=0;

        // read grid
        free(pdf->val[ngrids]);
        pdf->val[ngrids]=malloc(sizeof(double **)*pdf->nx[ngrids]);
        for(i=0; i<pdf->nx[ngrids]; i++) {
                pdf->val[ngrids][i]=malloc(sizeof(double*)*pdf->nq[ngrids]);
                for(j=0; j<pdf->nq[ngrids]; j++) {
                        pdf->val[ngrids][i][j]=malloc(sizeof(double)*pdf->n_pdf_flavours[ngrids]);
                        if(getline(&line, &size, fp)==-1) { fprintf(stderr, "Error in pdf member format, grid row is missed\n"); return EXIT_FAILURE;}
                        s_stream=fmemopen(line,size,"r"); 
                        for(k=0; k < pdf->n_pdf_flavours[ngrids]; k++) {
                                if(!fscanf(s_stream,"%lg", &pdf->val[ngrids][i][j][k])) { 
                                        fprintf(stderr, "Error in pdf member format, grid column is missed\n");
                                        return EXIT_FAILURE;
                                }
                        }
                        free(line);
                        line=NULL;
                        size=0;
                        fclose(s_stream);
                }
        }



        // check grid delimiter ---
        i_tmp=getline(&line, &size, fp);
        if(strcmp(line,"---") && strcmp(line,"---\n")) { 
                fprintf(stderr, "Error in pdf member format, can't find yaml delimiter '---'\n");
                return EXIT_FAILURE;
        }
        free(line);
        line=NULL;
        size=0;

        ngrids++;
        }
        return 1;
}

int save_lhapdf6_member(Pdf *pdf, char *path) {
        int ig,i,j,k;
        FILE *fp=fopen(path,"w");
        info_save(pdf->info, fp);
        fprintf(fp,"---\n");
        for(ig=0; ig<pdf->ngrids; ig++) {
                for(i=0; i<pdf->nx[ig]; i++) fprintf(fp,"%.8e ", pdf->x[ig][i]);
                fprintf(fp,"\n");
                for(i=0; i<pdf->nq[ig]; i++) fprintf(fp,"%.8e ", pdf->q[ig][i]);
                fprintf(fp,"\n");
                for(i=0; i<pdf->n_pdf_flavours[ig]; i++) fprintf(fp,"%i ", pdf->pdf_flavours[ig][i]);
                fprintf(fp,"\n");
                for(i=0; i<pdf->nx[ig]; i++) {
                        for(j=0; j<pdf->nq[ig]; j++) {
                                for(k=0; k < pdf->n_pdf_flavours[ig]; k++) fprintf(fp,"%.8e ", pdf->val[ig][i][j][k]);
                                fprintf(fp,"\n");
                        }
                }
                fprintf(fp,"---\n");
        }
        fclose(fp);
        return 0;
}

void pdf_free(Pdf *pdf) {
        int ig,i,j,k;
        info_free(pdf->info);
        for(ig=0; ig< pdf->ngrids; ig++) {
                free(pdf->x[ig]);
                free(pdf->q[ig]);
                free(pdf->pdf_flavours[ig]);
                for(i=0; i<pdf->nx[ig]; i++) { 
                        for(j=0; j<pdf->nq[ig]; j++) free(pdf->val[ig][i][j]);
                        free(pdf->val[ig][i]);
                }
                free(pdf->val[ig]);
        }
        free(pdf->nx);
        free(pdf->nq);
        free(pdf->n_pdf_flavours);
}

int load_lhapdf6_set(PdfSet *pdf_set, char *path) {
        int i;
        char *path_tmp=strdup(path); 
        char *pdf_name=basename(path_tmp); 

        size_t len=100;
        char *info_path=NULL; 
        FILE *ss=open_memstream(&info_path, &len);
        fprintf(ss, "%s/%s.info", path, pdf_name);
        fflush(ss);
        FILE *info_fp=fopen(info_path, "r");
        if(!info_fp) { 
                fprintf(stderr, "cannot open file %s\n", info_path);
                return 1;
        }
        fclose(ss);
        free(info_path);

        pdf_set->info=info_load(info_fp);
        fclose(info_fp);

        Info_Node *node=info_node_where(pdf_set->info, "NumMembers");
        pdf_set->n_members = atoi(node->value.string);

        pdf_set->members=malloc(sizeof(Pdf)*(pdf_set->n_members)); 
        char *member_path;          //FREE
        for(i=0; i<pdf_set->n_members; i++) { 
                member_path=NULL; 
                ss=open_memstream(&member_path, &len);
                fprintf(ss,"%s/%s_%04d.dat",path,pdf_name,i);
                fflush(ss);
                if(load_lhapdf6_member(&pdf_set->members[i], member_path)) { 
                        fprintf(stderr, "cannot read member %d\n", i);       
                        return 1;
                } 
                fclose(ss);
                free(member_path);
        }

        free(path_tmp);
        return 0;
}

int save_lhapdf6_set(PdfSet *pdf_set, char *path) {
        int i;
        char *path_tmp=strdup(path); 
        char *pdf_name=basename(path_tmp); 
        char *member_path;
        FILE *ss;
        size_t len=100;

        mkdir(path, 0755);

        char *info_path=NULL; 
        ss=open_memstream(&info_path, &len);
        fprintf(ss,"%s/%s.info", path, pdf_name);
        fflush(ss);
        FILE *info_fp=fopen(info_path,"w");
        if(!info_fp) { 
                fprintf(stderr, "Error: can't open file %s\n", info_path);
                return 1;
        }
        info_save(pdf_set->info, info_fp);
        fclose(ss);
        fclose(info_fp);
        free(info_path);

        for(i=0; i< pdf_set->n_members; i++) { 
                member_path=NULL; len=100;
                ss=open_memstream(&member_path, &len);
                fprintf(ss, "%s/%s_%04d.dat", path, pdf_name, i);
                fflush(ss);
                save_lhapdf6_member(&pdf_set->members[i], member_path); 
                fclose(ss);
                free(member_path);
        }
                free(path_tmp);
                return 0;
}

void pdf_set_initialize(PdfSet *pdf_set, int n_members, Info *info) {
        int i;
        pdf_set->n_members=n_members;        
        pdf_set->info=info_dup(info);
        pdf_set->members=malloc(sizeof(Pdf)*pdf_set->n_members);
//        for(i=0; i< pdf_set->n_members; i++) 
//                pdf_initialize(&pdf_set->members[i], minfo);
}

void pdf_set_initialize_as(PdfSet *pdf_set, const PdfSet *src) {
        int i;
        pdf_set->n_members=src->n_members;        
        pdf_set->info=info_dup(src->info);
        pdf_set->members=malloc(sizeof(Pdf)*pdf_set->n_members);
        for(i=0; i< pdf_set->n_members; i++) 
                pdf_initialize_as(&pdf_set->members[i], &src->members[i]);
}

PdfSet *pdf_set_dup(PdfSet *pdf_set) {
        int i;
        Pdf *member;
        PdfSet *copy=malloc(sizeof(PdfSet));
        copy->info=info_dup(pdf_set->info);
        copy->n_members=pdf_set->n_members;
        copy->members=malloc(sizeof(Pdf)*copy->n_members);
        for(i=0; i<copy->n_members; i++) { 
                member=pdf_dup(&pdf_set->members[i]);
                copy->members[i]=*member;
                free(member);
        }
        return copy;
}

int pdf_set_cmp(PdfSet *pdf_set1, PdfSet *pdf_set2) {
        int i;
        if(!info_cmp(pdf_set1->info, pdf_set2->info)) {
                fprintf(stderr, "Error pdf_set_cmp(info)\n");
                return 0;
        }
        if(pdf_set1->n_members!=pdf_set2->n_members) return 0;
        for(i=0; i<pdf_set1->n_members; i++)
                if(!pdf_cmp(&pdf_set1->members[i],&pdf_set2->members[i])) return 0;
        return 1;
}

void pdf_set_free(PdfSet *pdf_set) {
        int i;
        for(i=0; i<pdf_set->n_members; i++) pdf_free(&pdf_set->members[i]);
        info_free(pdf_set->info);
        free(pdf_set->members);
}

// utility function
char *n2str(char *s, double n ){ 
        sprintf(s,"%1.6g",n);
        return s;
}
