#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
// fortran common blocks

extern struct { //{{{
    double mz;
    char LHAPDF6OutDir[128], OutDirName[128];
    double grid[500];
    int nx;
    int read_xgrid;
    int nvarl[200];
    int dobands;
    float hf_mass[3];
    int i_fit_order, ipdfset;
    int lead, useGridLHAPDF5, writeLHAPDF6, WriteAlphaSToMemberPDF;
} ccommoninterface_;
//}}}

struct GridQX { //{{{
  int nx,nq2;
  double xmin, xmax, q2min, q2max;
  double *q2,*x;
  enum {QCDNUM_GRID, EXTERNAL_GRID, LHA5_GRID} type;
  //main interface (should include pdf modifications, like lead)
  double (*pdf_ij)(struct GridQX grid, int pid, int ix, int iq2); 
  double (*raw_pdf_ij)(struct GridQX grid, int pid, int ix, int iq2); //pdf from qcdnum grid
};
//}}}

// declarations {{{
void print_lhapdf6_();
void save_info(char *pdf_dir);
void save_data_lhapdf6(int *pdf_set,char *pdf_dir);
void save_data_lhapdf6_(int *pdf_set);
void save_data_lhapdf6_opt_(int *pdf_set);
void print_lhapdf6_();
void print_lhapdf6_opt_();
char* get_flavor_scheme();
char* get_error_type();
char* get_pdf_type(int pdf_set);
double raw_qcdnum_pdf_ij(struct GridQX grid, int pid, int ix, int iq2);
double raw_external_pdf_ij(struct GridQX grid, int pid, int ix, int iq2);
double lead_pdf_ij(struct GridQX grid, int pid, int ix, int iq2);
int get_NumMembers();
void save_alphas_info(FILE* fp, struct GridQX grid);

extern double xfrmix_(int *);
extern double qfrmiq_(int *);
extern double fvalij_(int *,int *,int *,int *,int *);
extern double fvalxq_(int *,int *,double *,double *,int *);
extern double hf_get_alphas_(double *);
extern int getord_(int *);
//}}}

// convert fortran string to C string
char* sfix(char* fstr,int size) { //{{{
  char* cstr=(char*)malloc(sizeof(char)*size);
  memcpy(cstr,fstr,size);
  char* end=cstr+size-1;
  while(*(--end)==' ');
  end++;
  *end='\0';
  cstr=realloc(cstr,strlen(cstr)+1);
  return cstr;
} //}}}

// define grid parameters depending on type
struct GridQX new_grid() { //{{{
  int inull,ix,iq2;
  struct GridQX grid;

  if(ccommoninterface_.read_xgrid) {
    grid.type=EXTERNAL_GRID;
    grpars_(&grid.nx,&grid.xmin,&grid.xmax,&grid.nq2,&grid.q2min,&grid.q2max,&inull);
    grid.nx=ccommoninterface_.nx;
    grid.x=malloc(sizeof(double)*grid.nx);
    grid.q2=malloc(sizeof(double)*grid.nq2);
    for(ix=0;ix<grid.nx;ix++) grid.x[ix]= ccommoninterface_.grid[ix];
    for(iq2=1;iq2<=grid.nq2;iq2++) grid.q2[iq2-1]=qfrmiq_(&iq2);
    grid.raw_pdf_ij=raw_external_pdf_ij;
    grid.pdf_ij=raw_external_pdf_ij;

  } else if(ccommoninterface_.useGridLHAPDF5) {
    grid.type=LHA5_GRID;
    grid.nx=161;
    grid.nq2=161;
    grid.x=malloc(sizeof(double)*grid.nx);
    grid.q2=malloc(sizeof(double)*grid.nq2);
    for(ix=0;ix<grid.nx;ix++) grid.x[ix]= ix<80? pow(10.0,6.0/120.0*ix-6.0) : pow(10.0,2.0/80.0*(ix-80)-2.01);
    for(iq2=0;iq2<grid.nq2;iq2++) grid.q2[iq2]= pow(10.0,(8.30103/160.0*iq2 ));
    grid.raw_pdf_ij=raw_external_pdf_ij;
    grid.pdf_ij=raw_external_pdf_ij;

  } else {
    grid.type=QCDNUM_GRID;
    grpars_(&grid.nx,&grid.xmin,&grid.xmax,&grid.nq2,&grid.q2min,&grid.q2max,&inull);
    grid.x=malloc(sizeof(double)*grid.nx);
    grid.q2=malloc(sizeof(double)*grid.nq2);
    for(ix=1;ix<=grid.nx;ix++) grid.x[ix-1]=xfrmix_(&ix);
    for(iq2=1;iq2<=grid.nq2;iq2++) grid.q2[iq2-1]=qfrmiq_(&iq2);
    grid.raw_pdf_ij=raw_qcdnum_pdf_ij;
    grid.pdf_ij=raw_qcdnum_pdf_ij;
  }
    if(ccommoninterface_.lead) grid.pdf_ij=lead_pdf_ij;
  return grid;
}
//}}}

void delete_grid(struct GridQX grid){ //{{{
  free(grid.q2);
  free(grid.x);
}
//}}}

// pdf in grid point (qcdnum or external grid) 
double raw_qcdnum_pdf_ij(struct GridQX grid, int pid, int ix, int iq2) { //{{{
  int inull;
  ix+=1;
  iq2+=1;
  return fvalij_(&ccommoninterface_.ipdfset,&pid,&ix,&iq2,&inull);
}

double raw_external_pdf_ij(struct GridQX grid, int pid, int ix, int iq2) {
  int inull;
  double x,q2;
  x=grid.x[ix];
  q2=grid.q2[iq2];
  return fvalxq_(&ccommoninterface_.ipdfset,&pid,&x,&q2,&inull);
}
//}}}

double lead_pdf_ij(struct GridQX grid, int pid, int ix, int iq2) { //{{{
  const double A=207.0, Z= 82.0;
  double val1,val2;
  int pid2;
  if(abs(pid)!=1 && abs(pid)!=2) return grid.raw_pdf_ij(grid, pid, ix, iq2);
  pid2= (abs(pid)==1?2:1)*abs(pid)/pid;
  val1=grid.raw_pdf_ij(grid, pid, ix, iq2);
  val2=grid.raw_pdf_ij(grid, pid2, ix, iq2);
  return (Z*val1 + (A-Z)*val2)/A;
}
//}}}

// write central value and info
void print_lhapdf6(char *pdf_dir){ //{{{
                    if(!ccommoninterface_.writeLHAPDF6) return;

  int central_set=0;
  char *outdir=sfix(ccommoninterface_.OutDirName,128);
  char *path=malloc(sizeof(char)*(strlen(outdir)+strlen(pdf_dir)+2));
  sprintf(path,"%s/%s",outdir,pdf_dir);
  mkdir(path,0755);
  save_info(pdf_dir);
  save_data_lhapdf6_(&central_set);
  free(outdir);
  free(path);
} 
// use LHAPDF6OutDir directory, fortran interface
void print_lhapdf6_(){
  char *pdf_dir=sfix(ccommoninterface_.LHAPDF6OutDir,128);
  print_lhapdf6(pdf_dir);
  free(pdf_dir);
} 
// use opt_$LHAPDF6OutDir directory, fortran interface
void print_lhapdf6_opt_(){
  char *pdf_dir=sfix(ccommoninterface_.LHAPDF6OutDir,128);
  char *opt_pdf_dir=malloc((strlen(pdf_dir)+strlen("opt_")+1)*sizeof(char));

  sprintf(opt_pdf_dir,"opt_%s",pdf_dir);
  print_lhapdf6(opt_pdf_dir);
  free(pdf_dir);
  free(opt_pdf_dir);
} //}}}

// store pdf data to herapdf_XXXX.dat yaml file
void save_data_lhapdf6(int *pdf_set,char *pdf_dir){ //{{{
                    if(!ccommoninterface_.writeLHAPDF6) return;

  int iq2,ix,i;
  struct GridQX grid=new_grid();
  int pdg_flavours[]={-5,-4,-3,-2,-1,1,2,3,4,5,21};
  int qcdnum_flavours[]={-5,-4,-3,-2,-1,1,2,3,4,5,0};
  FILE* fp;
//  char file_name[]="output/herapdf/herapdf_XXXX.dat";
  
  char *outdir=sfix(ccommoninterface_.OutDirName,128);
  char *path=malloc(sizeof(char)*(strlen(outdir)+2*strlen(pdf_dir)+3+strlen("_XXXX.dat")));
  sprintf(path,"%s/%s/%s_%04i.dat",outdir,pdf_dir,pdf_dir,*pdf_set);

//  sprintf(file_name,"output/herapdf/herapdf_%04i.dat",*pdf_set);
  if((fp=fopen(path,"w"))==NULL) puts("Cannot open file.");
    fprintf(fp,"PdfType: %s\n", get_pdf_type(*pdf_set));
    fprintf(fp,"Format: lhagrid1\n");
    if(ccommoninterface_.WriteAlphaSToMemberPDF) save_alphas_info(fp,grid);
    fprintf(fp,"---\n");
    for(ix=0;ix<grid.nx;ix++) 
      fprintf(fp,"%e ",grid.x[ix]);
      fprintf(fp,"\n");
    for(iq2=0;iq2<grid.nq2;iq2++) 
      fprintf(fp,"%e ",sqrt(grid.q2[iq2]));
      fprintf(fp,"\n");
    for(i=0;i<sizeof(pdg_flavours)/sizeof(int);i++) 
      fprintf(fp,"%i ",pdg_flavours[i]);
      fprintf(fp,"\n");
    for(ix=0;ix<grid.nx;ix++)
      for(iq2=0;iq2<grid.nq2;iq2++){
          for(i=0;i<sizeof(qcdnum_flavours)/sizeof(int);i++) 
            fprintf(fp,"%e ",grid.pdf_ij(grid,qcdnum_flavours[i],ix,iq2));
          fprintf(fp,"\n");
          }
    fprintf(fp,"---\n");
  fclose(fp);
  delete_grid(grid);
  free(outdir);
  free(path);
} 
// save to LHAPDF6OutDir, fortran interface
void save_data_lhapdf6_(int *pdf_set){ 
  char *pdf_dir=sfix(ccommoninterface_.LHAPDF6OutDir,128);
  save_data_lhapdf6(pdf_set, pdf_dir);
  free(pdf_dir);
}
// use opt_$LHAPDF6OutDir directory, fortran interface
void save_data_lhapdf6_opt_(int *pdf_set){ 
  char *pdf_dir=sfix(ccommoninterface_.LHAPDF6OutDir,128);
  char *opt_pdf_dir=malloc((strlen(pdf_dir)+strlen("opt_")+1)*sizeof(char));

  sprintf(opt_pdf_dir,"opt_%s",pdf_dir);
  save_data_lhapdf6(pdf_set, opt_pdf_dir);
  free(pdf_dir);
  free(opt_pdf_dir);
} //}}}

// store pdf information to herapdf.info yaml file
void save_info(char *pdf_dir) { //{{{
                    if(!ccommoninterface_.writeLHAPDF6) return;

  double xmin,xmax,q2min,q2max,q2,alphas,dnull;
  FILE* fp;
  struct GridQX grid=new_grid();

  char *outdir=sfix(ccommoninterface_.OutDirName,128);
  char *path=malloc(sizeof(char)*(strlen(outdir)+2*strlen(pdf_dir)+3+strlen(".info")));
  sprintf(path,"%s/%s/%s.info",outdir,pdf_dir,pdf_dir);
  
  if((fp=fopen(path,"w"))==NULL) puts("Cannot open file.");
    fprintf(fp,"SetDesc: HERAPDF\n");
    fprintf(fp,"Authors: ...\n");
    fprintf(fp,"Reference: ...\n");
    fprintf(fp,"Format: lhagrid1\n");
    fprintf(fp,"DataVersion: 1\n");
    fprintf(fp,"NumMembers: %i\n",get_NumMembers());
    fprintf(fp,"Flavors: [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21]\n");
    fprintf(fp,"OrderQCD: %i\n", ccommoninterface_.i_fit_order-1); // qcdnum notation LO=1,...; LHAPDF6 LO=0,...
    fprintf(fp,"FlavorScheme: %s\n", get_flavor_scheme());
    fprintf(fp,"ErrorType: %s\n", get_error_type());
    fprintf(fp,"XMin: %g\n", grid.xmin);
    fprintf(fp,"XMax: %g\n",1.0);
    fprintf(fp,"QMin: %g\n", sqrt(grid.q2min));
    fprintf(fp,"QMax: %g\n", sqrt(grid.q2max));
    fprintf(fp,"MZ: %g\n", ccommoninterface_.mz);
    fprintf(fp,"MUp: 0\n");
    fprintf(fp,"MDown: 0\n");
    fprintf(fp,"MStrange: 0\n");
    fprintf(fp,"MCharm: %g\n", ccommoninterface_.hf_mass[0]);
    fprintf(fp,"MBottom: %g\n", ccommoninterface_.hf_mass[1]);
    fprintf(fp,"MTop: %g\n", ccommoninterface_.hf_mass[2]);
    if(!ccommoninterface_.WriteAlphaSToMemberPDF) save_alphas_info(fp,grid);

    fclose(fp);
  delete_grid(grid);
  free(outdir);
  free(path);
} //}}}

void save_alphas_info(FILE* fp, struct GridQX grid) { //{{{
    int iq2,as_order;
    double q2;
    double mz2=ccommoninterface_.mz*ccommoninterface_.mz;

    getord_(&as_order);

    fprintf(fp,"AlphaS_MZ: %g\n", hf_get_alphas_(&mz2));
    fprintf(fp,"AlphaS_OrderQCD: %i\n", as_order-1); // qcdnum notation LO=1,...; LHAPDF6 LO=0,...
    fprintf(fp,"AlphaS_Type: ipol\n");

    fprintf(fp,"AlphaS_Qs: [");
      for(iq2=0;iq2<grid.nq2;iq2++) {
        fprintf(fp,"%g",sqrt(grid.q2[iq2]));
        if(iq2!=grid.nq2-1) fprintf(fp,", ");
      }
      fprintf(fp,"]\n");

    fprintf(fp,"AlphaS_Vals: [");
      for(iq2=0;iq2<grid.nq2;iq2++) {
        q2=grid.q2[iq2];
        fprintf(fp,"%g",hf_get_alphas_(&q2));
        if(iq2!=grid.nq2-1) fprintf(fp,", ");
      }
    fprintf(fp,"]\n");
}
//}}}

char* get_flavor_scheme() { //{{{
  int fixed_scheme;
  double dnull;
  char* fl_scheme;

  getcbt_(&fixed_scheme, &dnull, &dnull, &dnull);
  if(fixed_scheme) fl_scheme="fixed";
    else fl_scheme="variable";
  return fl_scheme;
}
//}}}

char* get_error_type(){ //{{{
  char* error_type;
  if(ccommoninterface_.dobands) 
      error_type="replicas";
    else error_type="hessian";
  return error_type;
}
//}}}

char* get_pdf_type(int pdf_set){ //{{{
  char* pdf_type;
  if(!pdf_set) pdf_type="central";
    else pdf_type="error";
  return pdf_type;
}
//}}}

int get_NumMembers(){ //{{{
  const int MNE=200; // see endmini.inc
  int npar=1, i; // central value
  for(i=0;i<MNE;i++) if(ccommoninterface_.nvarl[i]==1) npar+=2;
  return npar;
}
//}}}
