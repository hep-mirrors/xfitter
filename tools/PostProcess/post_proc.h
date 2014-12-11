#ifndef post_proc_h
#define post_proc_h

struct command {
        char* command;
        int (*function)(int argc, char* argv[]);
};

#endif
