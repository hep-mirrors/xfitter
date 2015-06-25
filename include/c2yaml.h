#ifndef c2yaml_h
#define c2yaml_h

#include <yaml.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef int yaml_id_t;

extern yaml_id_t i2yaml(yaml_document_t *doc,int );
extern yaml_id_t f2yaml(yaml_document_t *doc,float );
extern yaml_id_t d2yaml(yaml_document_t *doc,double );
extern yaml_id_t s2yaml(yaml_document_t *doc,char *);
extern yaml_id_t i_array2yaml(yaml_document_t *doc, const int *arr, int n);
extern yaml_id_t d_array2yaml(yaml_document_t *doc, const double *arr, int n);

#ifdef __cplusplus
}
#endif

#endif
