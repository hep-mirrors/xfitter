#ifndef xfitter_process_h
#define xfitter_process_h

struct command {
        char* command;
        int (*function)(int argc, char* argv[]);
};

#endif
