#include<stdlib.h>
#include<stdio.h>
#include "c2yaml.h"

yaml_id_t i2yaml(yaml_document_t *doc,int i) {
  int id;
  char *txt=alloca(100);
  snprintf(txt,100,"%d",i);
  id=yaml_document_add_scalar(doc, NULL, txt, -1, YAML_PLAIN_SCALAR_STYLE);
  return id;
}

yaml_id_t f2yaml(yaml_document_t *doc,float f){
  int id;
  char *txt=alloca(100);
  snprintf(txt,100,"%g",f);
  id=yaml_document_add_scalar(doc, NULL, txt, -1, YAML_PLAIN_SCALAR_STYLE);
  return id;
}

yaml_id_t d2yaml(yaml_document_t *doc,double d) {
  int id;
  char *txt=alloca(100);
  snprintf(txt,100,"%g",d);
  id=yaml_document_add_scalar(doc, NULL, txt, -1, YAML_PLAIN_SCALAR_STYLE);
  return id;
}

yaml_id_t s2yaml(yaml_document_t *doc,char *txt) {
  int id;
  id=yaml_document_add_scalar(doc, NULL, txt, -1, YAML_PLAIN_SCALAR_STYLE);
  return id;
}

yaml_id_t d_array2yaml(yaml_document_t *doc, const double *arr, int n) {
  int i,id;
  yaml_id_t s=yaml_document_add_sequence(doc, NULL, YAML_FLOW_SEQUENCE_STYLE);
  for(i=0;i<n;i++) yaml_document_append_sequence_item(doc, s, d2yaml(doc, arr[i]));
  return s;
}

yaml_id_t i_array2yaml(yaml_document_t *doc, const int *arr, int n) {
  int i,id;
  yaml_id_t s=yaml_document_add_sequence(doc, NULL, YAML_FLOW_SEQUENCE_STYLE);
  for(i=0;i<n;i++) yaml_document_append_sequence_item(doc, s, i2yaml(doc, arr[i]));
  return s;
}
