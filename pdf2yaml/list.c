#include <stdlib.h>
#include "list.h"

List *list_append(List *node, void *data) {
        List *new=malloc(sizeof(struct List_s));
        new->data=data;
        new->next=node;
        return new;
}

void list_free(List *node) {
        List *next;
        while(node!=NULL) {
                next=node->next;
                free(node);
                node=next;
        }
}

void list_each(List *node, void* (*func)(void*)) {
        while(node!=NULL) {
                func(node->data);
                node=node->next;
        }
}

List *list_reverse(const List *list) {
        List* r=NULL;
        while(list) { 
                r=list_append(r, list->data);
                list=list->next;
        }
        return r;
/*         List *this=list; */
/*         List *next=list->next; */
/*         list->next=NULL; */
/*         list=next; */
/*  */
/*         while(list) { */
/*                 next=list->next; */
/*                 list->next=this; */
/*                 this=list; */
/*                 list=next; */
/*         } */
/*         return this; */
}

void list_delete_child(List *parent) {
        List *child=parent->next;
        if(child==NULL) return;
        parent->next=child->next;
        free(child);
}

void list_push_child(List *parent, void *data) {
        List *child=parent->next;
        parent->next=list_append(child, data);
}

List *list_get_parent(List *list, List *child) {
        List *parent;
        for(parent=list; parent->next!=child && parent->next!=NULL; parent=parent->next);
        return parent;
}

