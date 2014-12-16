#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <yaml.h>
#include "c2yaml.h"
#include "list.h"
#include "pdf2yaml.h"

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

static Info *parse_yaml_seq(Info *info, yaml_document_t *doc, yaml_node_t *mkey, yaml_node_t *mval) {
        int seq_size=0,i=0;
        yaml_node_item_t *seq_it, *top, *end, *start;
        yaml_node_t *seq_node;
        start= mval->data.sequence.items.start;
        top=mval->data.sequence.items.top;
        end=mval->data.sequence.items.end;
        for(seq_size=0, seq_it=start; seq_it!=end && seq_it!= top; seq_it++) seq_size++;
        double *darray=malloc(sizeof(double)*seq_size);
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

// allocate memory for empty Pdf member
void pdf_initialize(Pdf *pdf, int nx, int nq, int n_pdf_flavours) {
        int i, j, k;
        pdf->nx=nx;
        pdf->nq=nq;
        pdf->n_pdf_flavours=n_pdf_flavours;

        pdf->x=malloc(sizeof(double)*pdf->nx);
        pdf->q=malloc(sizeof(double)*pdf->nq);
        pdf->pdf_flavours=malloc(sizeof(int)*pdf->n_pdf_flavours);

        pdf->val=malloc(sizeof(double **)*pdf->nx);
        for(i=0; i<pdf->nx; i++) {
                pdf->val[i]=malloc(sizeof(double*)*pdf->nq);
                for(j=0; j<pdf->nq; j++) pdf->val[i][j]=malloc(sizeof(double)*pdf->n_pdf_flavours);
        }
}

Pdf *pdf_dup(Pdf *pdf) {
        int i, j;
        Pdf *new_pdf=malloc(sizeof(Pdf));
        pdf_initialize(new_pdf, pdf->nx, pdf->nq, pdf->n_pdf_flavours);
        new_pdf->info=info_dup(pdf->info);
        for(i=0; i<pdf->nx; i++) 
                for(j=0; j<pdf->nq; j++)
                        memcpy(new_pdf->val[i][j], pdf->val[i][j], sizeof(double)*pdf->n_pdf_flavours);
        return new_pdf;
}

// load Pdf member from file. It does not require pdf_initialize for *pdf argument
int load_lhapdf6_member(Pdf *pdf, char *path) {
        double d_tmp;
        int i_tmp, i, j, k;
        size_t s_len=100, b_len=0;
        char *line=NULL;
        FILE *s_stream=NULL;
        FILE *b_stream=NULL;
        FILE *fp=fopen(path,"r");
        if(!fp) { fprintf(stderr, "cannot open file %s", path); return 1; }

        pdf->info=info_load(fp);


        rewind(fp);
        do {
                free(line);
                line=NULL;
                s_len=0;
        } while(i_tmp=getline(&line, &s_len, fp), strcmp(line,"---\n"));

        // read x line
        if(!getline(&line, &s_len, fp)) { fprintf(stderr, "Error in pdf member format\n"); return 1;}
        s_stream=fmemopen(line,s_len+1,"r"); 
        pdf->x=NULL;
        pdf->nx=0;
        b_len=sizeof(double);
        b_stream=open_memstream((char**)&pdf->x, &b_len); 
        while(fscanf(s_stream, "%lg", &d_tmp)) { fwrite(&d_tmp, sizeof(double), 1, b_stream); pdf->nx++; }
        fclose(s_stream);
        fclose(b_stream);

        //read q line
                free(line);
                line=NULL;
                s_len=0;
        if(!getline(&line, &s_len, fp)) { fprintf(stderr, "Error in pdf member format\n"); return 1;}
        s_stream=fmemopen(line,s_len+1,"r"); 
        pdf->q=NULL;
        pdf->nq=0;
        b_len=sizeof(double);
        b_stream=open_memstream((char**)&pdf->q, &b_len); 
        while(fscanf(s_stream, "%lg", &d_tmp)) { fwrite(&d_tmp, sizeof(double), 1, b_stream); pdf->nq++; }
        fclose(s_stream);
        fclose(b_stream);

        //read pdf_flavours line
                free(line);
                line=NULL;
                s_len=0;
        if(!getline(&line, &s_len, fp)) { fprintf(stderr, "Error in pdf member format\n"); return 1;}
        s_stream=fmemopen(line,s_len+1,"r"); 
        pdf->pdf_flavours=NULL;
        pdf->n_pdf_flavours=0;
        b_len=sizeof(int);
        b_stream=open_memstream((char**)&pdf->pdf_flavours, &b_len); 
        while(fscanf(s_stream, "%d", &i_tmp)) { fwrite(&i_tmp, sizeof(int), 1, b_stream); pdf->n_pdf_flavours++; }
        free(line);
        fclose(s_stream);
        fclose(b_stream);

        // read grid
        pdf->val=malloc(sizeof(double **)*pdf->nx);
        for(i=0; i<pdf->nx; i++) {
                pdf->val[i]=malloc(sizeof(double*)*pdf->nq);
                for(j=0; j<pdf->nq; j++) {
                        pdf->val[i][j]=malloc(sizeof(double)*pdf->n_pdf_flavours);
                        for(k=0; k < pdf->n_pdf_flavours; k++) 
                                if(!fscanf(fp,"%lg", &pdf->val[i][j][k])) { 
                                        fprintf(stderr, "Error in pdf member format\n");
                                        return 1;
                                }
                }
        }
        fclose(fp);
        return 0;
}

int save_lhapdf6_member(Pdf *pdf, char *path) {
        int i,j,k;
        FILE *fp=fopen(path,"w");
        info_save(pdf->info, fp);
        fprintf(fp,"---\n");
        for(i=0; i<pdf->nx; i++) fprintf(fp,"%.8e ", pdf->x[i]);
        fprintf(fp,"\n");
        for(i=0; i<pdf->nq; i++) fprintf(fp,"%.8e ", pdf->q[i]);
        fprintf(fp,"\n");
        for(i=0; i<pdf->n_pdf_flavours; i++) fprintf(fp,"%i ", pdf->pdf_flavours[i]);
        fprintf(fp,"\n");
        for(i=0; i<pdf->nx; i++) {
                for(j=0; j<pdf->nq; j++) {
                        for(k=0; k < pdf->n_pdf_flavours; k++) fprintf(fp,"%.8e ", pdf->val[i][j][k]);
                        fprintf(fp,"\n");
                }
        }
        fprintf(fp,"---\n");
        fclose(fp);
        return 0;
}

void pdf_free(Pdf *pdf) {
        int i,j,k;
        info_free(pdf->info);
        free(pdf->x);
        free(pdf->q);
        free(pdf->pdf_flavours);
        for(i=0; i<pdf->nx; i++) { 
                for(j=0; j<pdf->nq; j++) free(pdf->val[i][j]);
                free(pdf->val[i]);
        }
        free(pdf->val);
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
                printf("load: %s\n", member_path);
                if(load_lhapdf6_member(&pdf_set->members[i], member_path)) return 1; //FREE
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
                printf("save: %s\n", member_path);
                save_lhapdf6_member(&pdf_set->members[i], member_path); 
                fclose(ss);
                free(member_path);
        }
                free(path_tmp);
}

void pdf_set_free(PdfSet *pdf_set) {
        int i;
        for(i=0; i<pdf_set->n_members; i++) pdf_free(&pdf_set->members[i]);
        info_free(pdf_set->info);
        free(pdf_set->members);
}

// utility function
char *n2str(char *s, double n ){ 
        sprintf(s,"%g",n);
        return s;
}
