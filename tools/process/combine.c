#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include "pdf2yaml.h"


static void help() {
  puts("usage: xfitter-process combine pdf_dir_out base_pdf_dir added_pdf_dir1[:pdfs] added_pdf_dir2[:pdfs] ...");
  puts("          pdf_dir_out -  path to output set");
  puts("          base_pdf_dir - path to pdf containing grid");
  puts("          added_pdf_dir1,2,... - paths to interpolated pdfs from 1 or more sets");
  puts("          pdfs - numbers of interpolated pdfs, separated by commas, e.g. 2,7,15");
  puts("          there should be no space between number and comma");
  exit(0);
}

typedef struct Member_List_s{
  int * pdf_numbers;
  int n_pdfs;
  char * in_path;
} Member_List;

extern double alphas_ipol(double q, char * setname, int pdf_number);
extern void interpolation(double x, double Q, char* setname,
    int n_flavours, int *pdf_flavours, int pdf_number, double* values);

Member_List ** arg_parser(int n_added_sets, char ** args);
void add_grids(Pdf * pdf, Pdf * grid_pdf);
int pdf_update_info(Pdf * pdf, const Info * base_info);
void member_lists_free(Member_List *** member_lists_p, int n_lists);


int combine(int argc, char* argv[0]){
  
  if(argc<3) help();

  char* out_path=argv[0];
  int n_added_sets=argc-2, n_pdfs, tot_pdf_number=0, size, ipol;
  int i_added_set, i, k, ic;
  int ig, ix, iq, ifl;
  const char* in_path1 =argv[1];    
  const char* args[n_added_sets];    
  PdfSet Set1;
  PdfSet* SetOut;
  Member_List** member_lists;
  Info_Node * alphas_qs;
  Info_Node * alphas_vals;
  int forcepositive = 1;

  for (i_added_set = 0; i_added_set < n_added_sets; i_added_set++) {
    args[i_added_set] = argv[i_added_set+2];
  }

  member_lists = arg_parser (n_added_sets, args);
  
  if(load_lhapdf6_set(&Set1, in_path1)) return 1; 
  Pdf grid_pdf=Set1.members[0];
  SetOut = pdf_set_dup(&Set1);
  Pdf* out_pdfs[n_added_sets];
  
  for (i_added_set = 0; i_added_set<n_added_sets; i_added_set++) {
    
    PdfSet Set2;
    if (member_lists[i_added_set]->in_path == NULL) puts("NULLL");
    if (load_lhapdf6_set(&Set2, member_lists[i_added_set]->in_path)) return 1;    

    Info_Node* FP_par = info_node_where(Set2.info, "ForcePositive");
    if (FP_par == NULL || FP_par->value.string == "0") {
      forcepositive = 0;
    }

    char * setname = basename(member_lists[i_added_set]->in_path);
    if (member_lists[i_added_set]->n_pdfs == 0) {
      member_lists[i_added_set]->n_pdfs = Set2.n_members;
      member_lists[i_added_set]->pdf_numbers=malloc(Set2.n_members*sizeof(int));
      if( !member_lists[i_added_set]->pdf_numbers ) {
        printf("Memory could not be allocated!");
        exit(1);
      }
      for (i=0; i<Set2.n_members; i++) member_lists[i_added_set]->pdf_numbers[i]=i; 
    }

    n_pdfs = member_lists[i_added_set]->n_pdfs;
    out_pdfs[i_added_set] = malloc(n_pdfs*sizeof(Pdf));
    if( !out_pdfs[i_added_set] ) {
      printf("Memory could not be allocated!");
      exit(1);
    }

    for (i=0;i<n_pdfs;i++){
      Info * info = info_dup(Set2.members[member_lists[i_added_set]->pdf_numbers[i]].info);
      pdf_initialize(&out_pdfs[i_added_set][i], info);
      add_grids(&out_pdfs[i_added_set][i], &grid_pdf);
      ipol = pdf_update_info(&out_pdfs[i_added_set][i], Set2.info);
      if (!ipol){
        size =0;
        for (ig=0; ig<out_pdfs[i_added_set][i].ngrids; ig++) size += out_pdfs[i_added_set][i].nq[ig]-1;
        double q_vals[size+1];
        double al_vals[size+1];
        k=0;
        for (ig=0; ig<out_pdfs[i_added_set][i].ngrids; ig++){
          for (iq=0; iq<out_pdfs[i_added_set][i].nq[ig]-1; iq++){
            q_vals[k]=out_pdfs[i_added_set][i].q[ig][iq];
            al_vals[k]=alphas_ipol(q_vals[k],setname, member_lists[i_added_set]->pdf_numbers[i]);
            k++;
          }
        } 
        q_vals[k]=out_pdfs[i_added_set][i].q[ig-1][iq];
        al_vals[k]= alphas_ipol(out_pdfs[i_added_set][i].q[ig-1][iq],setname, member_lists[i_added_set]->pdf_numbers[i]);
        size++;
        alphas_qs = info_node_where(out_pdfs[i_added_set][i].info, "AlphaS_Qs");
        if (alphas_qs)
          for (k=0; k<size+1; k++) alphas_qs->value.darray.vals[k]=q_vals[k]; 
        else {
          out_pdfs[i_added_set][i].info = info_add_node_darray(out_pdfs[i_added_set][i].info, "AlphaS_Qs", q_vals, size); 
        }
        alphas_vals = info_node_where(out_pdfs[i_added_set][i].info, "AlphaS_Vals");
        if (alphas_vals)
          for (k=0; k<size+1; k++) alphas_vals->value.darray.vals[k]=al_vals[k]; 
        else {
          out_pdfs[i_added_set][i].info = info_add_node_darray(out_pdfs[i_added_set][i].info, "AlphaS_Vals", al_vals, size); 
        }
      } 
    }

    for (i=0; i<n_pdfs; i++) {
      for(ig=0; ig<grid_pdf.ngrids; ig++) {
        for(ix=0; ix<grid_pdf.nx[ig]; ix++) {
          for(iq=0; iq<grid_pdf.nq[ig]; iq++) {
            interpolation(out_pdfs[i_added_set][i].x[ig][ix], out_pdfs[i_added_set][i].q[ig][iq], \
            setname, out_pdfs[i_added_set][i].n_pdf_flavours[0],out_pdfs[i_added_set][i].pdf_flavours[0], \
            member_lists[i_added_set]->pdf_numbers[i], out_pdfs[i_added_set][i].val[ig][ix][iq]);
          }  
        }
      }
    }
    pdf_set_free(&Set2);
  }

  char* tot_pdf_number_str[10]; 
  sprintf(tot_pdf_number_str, "%d", Set1.n_members+tot_pdf_number);
  
  info_node_update_str(info_node_where(SetOut->info, "NumMembers"), tot_pdf_number_str);  
  
  Info_Node* FP_par_orig = info_node_where(Set1.info, "ForcePositive");
  //printf("forcepositive = %i\n", forcepositive);
  if (FP_par_orig != NULL) {
    if (!strcmp(FP_par_orig->value.string, "1") && forcepositive == 0) {
      info_node_update_str (info_node_where (SetOut->info, "ForcePositive"), "0");  
    }
  }

  save_lhapdf6_set(SetOut, out_path);

  tot_pdf_number=0;
  for (i_added_set = 0; i_added_set<n_added_sets; i_added_set++){
    char pdf_out_paths[member_lists[i_added_set]->n_pdfs][256];  
    for (i = 0; i<member_lists[i_added_set]->n_pdfs; i++){
      sprintf(pdf_out_paths[i], "%s%s%s%s%04d%s", \
      out_path, "/",basename(out_path), "_" ,Set1.n_members+tot_pdf_number+i, ".dat");
      save_lhapdf6_member(&out_pdfs[i_added_set][i], pdf_out_paths[i]);

      pdf_free(&out_pdfs[i_added_set][i]);  
    }
    tot_pdf_number+=member_lists[i_added_set]->n_pdfs;  
  }
  
  pdf_set_free(&Set1);
  member_lists_free(&member_lists, n_added_sets);
  pdf_set_free(SetOut);

  return 0;
}   


Member_List ** arg_parser(int n_added_sets, char ** args) {
 
  char* colon=NULL;
  char* comma1=NULL;
  char* comma2=NULL;
  int ic,i=0, i_added_set;
  int length;
  char pdf_number[10];
  Member_List **member_lists;

  member_lists= malloc(n_added_sets*sizeof(Member_List*));
  if( !member_lists ) {
    printf("Memory could not be allocated!");
    exit(1);
  }
  for (i_added_set=0; i_added_set<n_added_sets; i_added_set++){
    i = 0;    
    colon = strchr(args[i_added_set], ':');
    member_lists[i_added_set]=malloc(sizeof(Member_List)); 
    if( !member_lists[i_added_set] ) {
      printf("Memory could not be allocated!");
      exit(1);
    }
    if(colon == NULL)
    {
      printf("warning: there are no pdf numbers for %d added set\n", i_added_set+1);
      member_lists[i_added_set]->in_path = strdup(args[i_added_set]);
      member_lists[i_added_set]->n_pdfs=0;
    }
    
    else {
      length = colon- args[i_added_set];
      member_lists[i_added_set]->in_path = strndup(args[i_added_set], length);
      comma1 = colon;
      comma2=strchr(colon,',');
      member_lists[i_added_set]->pdf_numbers=malloc(sizeof(int));
      if( !member_lists[i_added_set]->pdf_numbers ) {
        printf("Memory could not be allocated!");
        exit(1);
      }
     
      while (comma2 != NULL) {
        strncpy(pdf_number, comma1+1, comma2 - comma1-1);
        member_lists[i_added_set]->pdf_numbers = realloc(member_lists[i_added_set]->pdf_numbers,(i+1)*sizeof(int));
        if( !member_lists[i_added_set]->pdf_numbers ) {
          printf("Reallocation failed!");
          exit(1);
        }
        member_lists[i_added_set]->pdf_numbers[i]=atoi(pdf_number);
        comma1 = comma2;
        comma2=strchr(comma2+1,',');
        i++;
      }
      strncpy(pdf_number, comma1+1, strlen(comma1));
      member_lists[i_added_set]->n_pdfs = i+1;
      member_lists[i_added_set]->pdf_numbers[i]=atoi(pdf_number);
    } 
  }

  return member_lists;
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

void member_lists_free(Member_List *** member_lists_p,int n_lists){
  int i_list;
  for (i_list=0; i_list<n_lists; i_list++){
    free((*member_lists_p)[i_list]->pdf_numbers);
    free((*member_lists_p)[i_list]->in_path);
  }
  free(*member_lists_p);
}
