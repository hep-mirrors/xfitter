#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include "pdf2yaml.h"


static void help() {
  puts("usage: postproc combine outdir dir1 dir2[:pdfs]");
  puts("          outdir -  path to output set");
  puts("          dir1 - path to pdf containing grid");
  puts("          dir2 - path to interpolated pdfs");
  puts("          pdfs - numbers of interpolated pdfs, separated by commas");
  puts("          e.g. 2,7,15");
  puts("          there should be no space between number and comma");
  exit(0);
}

extern double alphas_ipol(double q, char * pdfset_path, int pdf_number);
extern void interpolation(double x, double Q, char* pdfset_path,
    int n_flavours, int *pdf_flavours, int pdf_number, double* values);

void arg_parser(char * colon, int ** pdf_numbers, int * n_pdfs);
void add_grids(Pdf * pdf, Pdf * grid_pdf);
int pdf_update_info(Pdf * pdf, const Info * base_info);


int combine(int argc, char* argv[0]){
  
  if(argc<3) help();
  char* out_path=argv[0];
  const char* in_path1 =argv[1];    
  const char* str[argc];    
  char in_path2[100];
  PdfSet Set1;
  PdfSet* SetOut;
  Info_Node * alphas_qs;
  Info_Node * alphas_vals;

  if(load_lhapdf6_set(&Set1, in_path1)) return 1; 
  Pdf grid_pdf=Set1.members[0];
  SetOut = pdf_set_dup(&Set1);

  int n_pdfs, tot_pdf_number=0, size, ipol;
  int * pdf_numbers;
  int i_added_set, i, k, ic;
  int ig, ix, iq, ifl;
 
  for (i_added_set = 2; i_added_set<argc; i_added_set++) {
    
    char * colon= NULL;
    PdfSet Set2;
    
    str[i_added_set] = argv[i_added_set];
    for (ic=0; ic<100; ic++) in_path2[ic] = '\0';
    
    colon = strchr(str[i_added_set], ':');
    
    if(colon == NULL) strcpy(in_path2, str[i_added_set]);
    else strncpy(in_path2, str[i_added_set], colon - str[i_added_set]);

    if (load_lhapdf6_set(&Set2, in_path2)) return 1;    

    if (colon == NULL) {
      printf("warning: there are no pdf numbers\n");
      n_pdfs=Set2.n_members;
      pdf_numbers=malloc(n_pdfs*sizeof(int));
      for (i=0;i<n_pdfs;i++){
        pdf_numbers[i] = i;
      }
    }
    
    else arg_parser(colon, &pdf_numbers, &n_pdfs); 
    
    Pdf out_pdfs[n_pdfs];

    char * pdfset_path = basename(in_path2);
    
    for (i=0;i<n_pdfs;i++){
      Info * info = info_dup(Set2.members[pdf_numbers[i]].info);
      pdf_initialize(&out_pdfs[i], info);
      add_grids(&out_pdfs[i], &grid_pdf);
      ipol = pdf_update_info(&out_pdfs[i], Set2.info);
      if (!ipol){
        size =0;
        for (ig=0; ig<out_pdfs[i].ngrids; ig++) size += out_pdfs[i].nq[ig]-1;
        double q_vals[size+1];
        double al_vals[size+1];
        k=0;
        for (ig=0; ig<out_pdfs[i].ngrids; ig++){
          for (iq=0; iq<out_pdfs[i].nq[ig]-1; iq++){
            q_vals[k]=out_pdfs[i].q[ig][iq];
            al_vals[k]=alphas_ipol(q_vals[k],pdfset_path, pdf_numbers[i]);
            k++;
          }
        } 
        q_vals[k]=out_pdfs[i].q[ig-1][iq];
        al_vals[k]= alphas_ipol(out_pdfs[i].q[ig-1][iq],pdfset_path, pdf_numbers[i]);
        size++;
        alphas_qs = info_node_where(out_pdfs[i].info, "AlphaS_Qs");
        if (alphas_qs)
          for (k=0; k<size+1; k++) alphas_qs->value.darray.vals[k]=q_vals[k]; 
        else {
          out_pdfs[i].info = info_add_node_darray(out_pdfs[i].info, "AlphaS_Qs", q_vals, size); 
          //printf("K=%d\n", k);
        }
        alphas_vals = info_node_where(out_pdfs[i].info, "AlphaS_Vals");
        if (alphas_vals)
          for (k=0; k<size+1; k++) alphas_vals->value.darray.vals[k]=al_vals[k]; 
        else {
          out_pdfs[i].info = info_add_node_darray(out_pdfs[i].info, "AlphaS_Vals", al_vals, size); 
        }
      } 
    }

    for (i=0; i<n_pdfs; i++) {
      for(ig=0; ig<grid_pdf.ngrids; ig++) {
        for(ix=0; ix<grid_pdf.nx[ig]; ix++) {
          for(iq=0; iq<grid_pdf.nq[ig]; iq++) {
            interpolation(out_pdfs[i].x[ig][ix], out_pdfs[i].q[ig][iq], \
            pdfset_path, out_pdfs[i].n_pdf_flavours[0],out_pdfs[i].pdf_flavours[0], \
            pdf_numbers[i], out_pdfs[i].val[ig][ix][iq]);
          }  
        }
      }
    }
    char pdf_out_paths[n_pdfs][256];  

    for (i = 0; i<n_pdfs; i++){
      Info *base_info0 = Set2.info;
      sprintf(pdf_out_paths[i], "%s%s%s%s%04d%s", \
      out_path, "/",basename(out_path), "_" ,Set1.n_members+tot_pdf_number+i, ".dat");

      save_lhapdf6_member(&out_pdfs[i], pdf_out_paths[i]);

      pdf_free(&out_pdfs[i]);  
    }

    tot_pdf_number+=n_pdfs;  
    pdf_set_free(&Set2);
  }

  char* tot_pdf_number_str[10]; 
  sprintf(tot_pdf_number_str, "%d", Set1.n_members+tot_pdf_number);
  
  info_node_update_str(info_node_where(SetOut->info, "NumMembers"), tot_pdf_number_str);  
  
  save_lhapdf6_set(SetOut, out_path);
  pdf_set_free(&Set1);
  pdf_set_free(SetOut);

  return 0;
}   



void arg_parser(char * colon, int ** pdf_numbers, int * n_pdfs){
  
  char* comma1=NULL;
  char* comma2=NULL;
  int ic,i=0;
  char pdf_number[10];

  for (ic=0;ic<10;ic++) pdf_number[ic]='\0';
  
  comma1 = colon;
  comma2=strchr(colon,',');
  *pdf_numbers=malloc(sizeof(int));
  while (comma2 != NULL) {
    strncpy(pdf_number, comma1+1, comma2 - comma1-1);
    *pdf_numbers = realloc(*pdf_numbers,(i+1)*sizeof(int));
    (*pdf_numbers)[i]=atoi(pdf_number);
    comma1 = comma2;
    comma2=strchr(comma2+1,',');
    i++;
  }
  strncpy(pdf_number, comma1+1, strlen(comma1));
  *n_pdfs = i+1;
  (*pdf_numbers)[*n_pdfs-1]=atoi(pdf_number);
}


void add_grids(Pdf * pdf, Pdf * grid_pdf){
  
  int ig, ix, iq, ifl;
  
  for(ig=0; ig<grid_pdf->ngrids; ig++) {
    pdf_add_grid(pdf, grid_pdf->nx[ig], \
    grid_pdf->nq[ig], grid_pdf->n_pdf_flavours[ig]);
    for (ifl=0; ifl<grid_pdf->n_pdf_flavours[ig]; ifl++) { 
      pdf->pdf_flavours[ig][ifl]=grid_pdf->pdf_flavours[ig][ifl];
    }
    for(iq=0; iq<grid_pdf->nq[ig]; iq++) {
      pdf->q[ig][iq]=grid_pdf->q[ig][iq];
    }
    for(ix=0; ix<grid_pdf->nx[ig]; ix++) {
      pdf->x[ig][ix]=grid_pdf->x[ig][ix];
    }
  } 
}

int pdf_update_info(Pdf * pdf, const Info * base_info){

  const int npars = 6;
  char set_parameters[6][64]={"MCharm", "MBottom", "MTop", \
     "AlphaS_MZ", "AlphaS_OrderQCD", "AlphaS_Type"};
  
  int  j, ipol;
  Info_Node* param; 
  Info_Node* base_param;
  Info_Node* grid_param;
  for(j=0; j<npars; j++){
    base_param = info_node_where(base_info, set_parameters[j]);
    param = info_node_where(pdf->info, set_parameters[j]);
    if (base_param == NULL) {
      continue;
    }
    if (param == NULL){
      pdf->info = info_add_node_str(pdf->info, base_param->key, base_param->value.string);  
      param = info_node_where(pdf->info, set_parameters[j]);
    }
  }
  if (base_param != NULL){ 
    ipol = strcmp(base_param->value.string,"ipol");
  }
  else if (param != NULL){
    ipol = strcmp(param->value.string,"ipol");
  }
  else return 1;
  return ipol;
}
