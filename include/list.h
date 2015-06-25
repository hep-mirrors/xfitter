#ifndef List_h
#define List_h

#ifdef __cplusplus
extern "C" {
#endif

struct List_s {
        struct List_s *next;
        void *data;
};

typedef struct List_s List;

extern List *list_append(List *node, void *data);

extern void list_free(List *node);

extern List *list_reverse(const List *list);

extern void list_each(List *node, void* (*func)(void*));

void list_delete_child(List *parent);

void list_push_child(List *parent, void *data);

List *list_get_parent(List *list, List *child);

#ifdef __cplusplus
}
#endif

#endif
