#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include "pdf2yaml.h"

static void help() {
  puts("usage: postproc combine outdir dir1 dir2");
  puts("          outdir -  path to output set");
  puts("          dir1 - path to pdf containing grid");
  puts("          dir2 - path to interpolated pdfs");
  exit(0);
}

int combine(int argc, char* argv[0]){
  if(argc!=3) help();
  
  char* out_path=argv[0];
  const char* in_path1=argv[1];
  const char* str = argv[2];
  char in_path2[100];
  int ic;
  for (ic=0;ic<100;ic++) in_path2[ic]='\0';

  PdfSet Set1;
  PdfSet Set2;
  PdfSet* SetOut;
  
  char* colon = NULL;
  char* comma1 = NULL;
  char* comma2 = NULL;
  int n_pdfs, i=0;
  char pdf_number[4];
  for (ic=0;ic<4;ic++) pdf_number[ic]='\0';
  colon = strchr(str, ':');
  
  if(colon == NULL) strcpy(in_path2, str);
  else strncpy(in_path2, str, colon - str);
	
  if(load_lhapdf6_set(&Set1, in_path1)  || load_lhapdf6_set(&Set2, in_path2)) return 1;
  SetOut = pdf_set_dup(&Set1);
  int pdf_numbers[Set2.n_members];

  if (colon == NULL) {
    printf("warning: there are no pdf numbers\n");
    
    n_pdfs=Set2.n_members;
    for (i=0;i<n_pdfs;i++){
      pdf_numbers[i] = i;
    }
  }
  
  else {
    comma2=strchr(str,',');
    comma1 = colon;
    while (comma2!=NULL) {
      strncpy(pdf_number, comma1+1, comma2 - comma1-1);
      pdf_numbers[i]=atoi(pdf_number);
      comma1 = comma2;
      comma2=strchr(comma2+1,',');
      i++;
    }
    strncpy(pdf_number, comma1+1, strlen(comma1));
    pdf_numbers[i]=atoi(pdf_number);
    n_pdfs = i+1;
  }

  Pdf out_pdfs[n_pdfs];

  for (i=0;i<n_pdfs;i++){
    pdf_initialize_as(&out_pdfs[i], &Set1.members[0]);
  }

  char * pdfset_path = basename(in_path2);
  int ig, ix, iq, ifl;
  

  for (i=0; i<n_pdfs; i++) {
    for(ig=0; ig<(out_pdfs[i]).ngrids; ig++) {
      for(ix=0; ix<(out_pdfs[i]).nx[ig]; ix++ ) {
        for(iq=0; iq<(out_pdfs[i]).nq[ig]; iq++ ) {
          interpolation(out_pdfs[i].x[ig][ix], out_pdfs[i].q[ig][iq], 
          pdfset_path, out_pdfs[i].n_pdf_flavours[0],out_pdfs[i].pdf_flavours[0], pdf_numbers[i], out_pdfs[i].val[ig][ix][iq]);
        }	
      }
    }
  }
  save_lhapdf6_set(SetOut, out_path);

  char pdf_out_paths[n_pdfs][256];  
  for(i=0; i<n_pdfs; i++) {
    sprintf(pdf_out_paths[i], "%s%s%s%s%04d%s", \
    out_path, "/",basename(out_path), "_" ,Set1.n_members+i, ".dat");
  }
  
  for (i = 0; i<n_pdfs; i++){
    save_lhapdf6_member(&out_pdfs[i], pdf_out_paths[i]);
    pdf_free(&out_pdfs[i]);  
  }

  pdf_set_free(&Set1);
  pdf_set_free(&Set2);
  pdf_set_free(SetOut);
  
  return 0;

}	
