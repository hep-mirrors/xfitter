#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "xfitter-process.h"

extern int profile(int argc, char* argv[]);
extern int rotate(int argc, char* argv[]);
extern int symmetrize(int argc, char* argv[]);
extern int custom(int argc, char* argv[]);
extern int reweight(int argc, char* argv[]);
static int help(int argc, char* argv[]);

static const struct command options[]={
        {"symmetrize",symmetrize},
        {"rotate",rotate},
        {"profile",profile},
        {"reweight",reweight},
        {"custom",custom},
        {"help", help},
        {"-h", help},
        {"--help", help},
};

static int help(int argc, char *argv[]) {
        int i;
        char *module_opt[]={"--help"};
        if(!argc) { 
                puts("usage: xfitter-process <module> [<args>]\n\nfor command info use \n\txfitter-process help module");
                puts("\navailable modules:");
                for(i=0; i<sizeof(options)/sizeof(struct command); i++) 
                        if(options[i].function!=help) printf("\t%s\n", options[i].command);
                
        }
        else 
                for(i=0; i<sizeof(options)/sizeof(struct command); i++) {
                        if(!strcmp(options[i].command,argv[0])) {
                                options[i].function(1, module_opt);
                                break;
                        }
                }
      exit(0);
      return(0);
}


int main (int argc, char **argv) {

        int i;
        int result=0;
        argc--;
        argv++;
        if(!argc) { 
                help(0,argv); 
                exit(0);
        }
        for(i=0; i<sizeof(options)/sizeof(struct command); i++) {
                if(!strcmp(options[i].command,argv[0])) {
                        result=options[i].function(--argc, ++argv);
                        break;
                }
        }
        return result;
}
