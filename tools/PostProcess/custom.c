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
        if(!strcmp(argv[0],"--help") || argc!=3) help();

        int i, j, ix, iq, ifl, i_tmp; // counters and dummy vars
        char *in_path=argv[1];

        char *in_path_tmp=strdup(in_path); //FREE
        char *pdf_in_name=basename(in_path_tmp); 

        //load info file
        size_t len=0;
        char *info_path=NULL; 
        FILE *ss=open_memstream(&info_path, &len);
        fprintf(ss,"%s/%s.info",in_path,pdf_in_name);
        fflush(ss);
        FILE *info_fp=fopen(info_path, "r");
        fclose(ss);
        
        Info *info=info_load(info_fp);
        fclose(info_fp);

        Info_Node *node=info_node_where(info, "NumMembers");
        int n_members=atoi(node->value.string);
        printf("NumMembers: %d\n", n_members);

        //load members
        Pdf *members=malloc(sizeof(Pdf)*(n_members)); 
        char *member_path;          //FREE
        for(i=0; i<n_members; i++) { 
                member_path=NULL; 
                ss=open_memstream(&member_path, &len);
                fprintf(ss,"%s/%s_%04d.dat",in_path,pdf_in_name,i);
                fflush(ss);
                printf("load: %s\n", member_path);
                load_lhapdf6_member(&members[i], member_path); //FREE
                fclose(ss);
        }

        /***                      CUSTOM CODE                      ***/

        
        /***                    END CUSTOM CODE                    ***/

        //clean up
        for(i=0; i<n_members; i++) pdf_free(&members[i]);
        info_free(info);
        free(in_path_tmp);
        free(members);

        return 0;
}
