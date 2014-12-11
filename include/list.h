#ifndef List_h
#define List_h

struct List_s {
        struct List_s *next;
        void *data;
};

typedef struct List_s List;

List *list_append(List *node, void *data);

void list_free(List *node);

List *list_reverse(const List *list);
void list_each(List *node, void* (*func)(void*));

#endif
