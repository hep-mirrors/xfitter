#include<string.h>
#include<stdio.h>
#include<malloc.h>
#include<math.h>
// fortran common blocks
extern struct { //{{{
  double mz, mw, mh;
} boson_masses_;  //}}}

extern struct { //{{{
    float q2val[40], starting_scale, strange_frac;
    double chi2maxerror;
    float hf_mass[3], charm_frac;
    int ldebug, dobands, h1qcdfunc;
    int itheory, i_fit_order, iparam, hfscheme;
    int lrand;
    int statype, systype;
    float outxrange[2];
    int outnx, ilenpdf, nchebglu;
    float chebxmin;
    int nchebsea;
    float wmnlen, wmxlen;
    int ichebtypeglu, ichebtypesea, ifsttype, iseedmc, ioffsetchebsea, 
	    ewfit, npolyval;
    int lead;
    int izpopoly, ipolysqr;
    int useapplg;
    int napplgrids;
    int lfitdy, lfastapplgrid, lranddata;
    int idh_mod, ipdfset, icheck_qcdnum;
    int lhapdferrors;
    int nlhapdf_sets, useprevfit;
    int corrsystbyoffset;
    int corsysindex;
    char statscale[32], uncorsysscale[32], corsysscale[32], uncorchi2type[32];
    float corchi2type;
    char hf_scheme[32];
    int scale68;
    int asymerrorsiterations;
    int luseapplgridckm;
} steering_; //}}}

extern struct { //{{{
double grid[500];
int nx;
int read_xgrid;
} ext_xgrid_;
//}}}

struct GridQX { //{{{
  int nx,nq2;
  double xmin, xmax, q2min, q2max;
  double *q2,*x;
  enum {QCDNUM_GRID,EXTERNAL_GRID} type;
  //main interface (should include pdf modifications, like lead)
  double (*pdf_ij)(struct GridQX grid, int pid, int ix, int iq2); 
  double (*raw_pdf_ij)(struct GridQX grid, int pid, int ix, int iq2); //pdf from qcdnum grid
};
//}}}

// declarations // {{{
void print_lhapdf6_();
void save_info(); 
void save_data_lhapdf6_(int *pdf_set);
char* get_flavor_scheme();
char* get_error_type();
char* get_pdf_type(int pdf_set);
double raw_qcdnum_pdf_ij(struct GridQX grid, int pid, int ix, int iq2);
double raw_external_pdf_ij(struct GridQX grid, int pid, int ix, int iq2);
double lead_pdf_ij(struct GridQX grid, int pid, int ix, int iq2);
extern double xfrmix_(int *);
extern double qfrmiq_(int *);
extern double fvalij_(int *,int *,int *,int *,int *);
extern double fvalxq_(int *,int *,double *,double *,int *);
extern double hf_get_alphas_(double *);
//}}}

// convert fortran string to C string
char* sfix(char* fstr,int size) { //{{{
  char* cstr=(char*)malloc(sizeof(char)*size);
  memcpy(cstr,fstr,size);
  cstr[size-1]='\0';
  char* end;
  while( end=rindex(cstr,' ') ) cstr[end-cstr]='\0';
  return cstr;
} //}}}

// define grid parameters depending on type
struct GridQX new_grid() { //{{{
  int inull,ix,iq2;
  struct GridQX grid;
  // qcdnum
  if(ext_xgrid_.read_xgrid) {
    grid.type=EXTERNAL_GRID;
    grpars_(&grid.nx,&grid.xmin,&grid.xmax,&grid.nq2,&grid.q2min,&grid.q2max,&inull);
    grid.nx=ext_xgrid_.nx;
    grid.x=malloc(sizeof(double)*grid.nx);
    grid.q2=malloc(sizeof(double)*grid.nq2);
    for(ix=0;ix<grid.nx;ix++) grid.x[ix]= ext_xgrid_.grid[ix];
    for(iq2=1;iq2<=grid.nq2;iq2++) grid.q2[iq2-1]=qfrmiq_(&iq2);
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
    if(steering_.lead) grid.pdf_ij=lead_pdf_ij;
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
  return fvalij_(&steering_.ipdfset,&pid,&ix,&iq2,&inull);
}

double raw_external_pdf_ij(struct GridQX grid, int pid, int ix, int iq2) {
  int inull;
  double x,q2;
  x=grid.x[ix];
  q2=grid.q2[iq2];
  return fvalxq_(&steering_.ipdfset,&pid,&x,&q2,&inull);
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

// fortran interface to store pdf in LHAPDF6 format
void print_lhapdf6_(){ //{{{
  int central_set=0;
  mkdir("output/herapdf",0755);
  save_info();
  save_data_lhapdf6_(&central_set);
} //}}}

// store pdf data to herapdf_0000.dat yaml file
void save_data_lhapdf6_(int *pdf_set){ //{{{
  int iq2,ix,i;
  struct GridQX grid=new_grid();
  int pdg_flavours[]={-5,-4,-3,-2,-1,1,2,3,4,5,21};
  int qcdnum_flavours[]={-5,-4,-3,-2,-1,1,2,3,4,5,0};
  FILE* fp;
  char file_name[]="output/herapdf/herapdf_XXXX.dat";

  sprintf(file_name,"output/herapdf/herapdf_%04i.dat",*pdf_set);
  if((fp=fopen(file_name,"w"))==NULL) puts("Cannot open file.");
    fprintf(fp,"PdfType: %s\n", get_pdf_type(*pdf_set));
    fprintf(fp,"Format: lhagrid1\n");
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
} //}}}

// store pdf information to herapdf.info yaml file
void save_info() { //{{{

  int iq2,as_order;
  double xmin,xmax,q2min,q2max,q2,alphas,dnull;
  double mz2=boson_masses_.mz*boson_masses_.mz;
  FILE* fp;
  struct GridQX grid=new_grid();
  
  getord_(&as_order);
  if((fp=fopen("output/herapdf/herapdf.info","w"))==NULL) puts("Cannot open file.");
    fprintf(fp,"SetDesc: HERAPDF\n");
    fprintf(fp,"Authors: ...\n");
    fprintf(fp,"Reference: ...\n");
    fprintf(fp,"Format: lhagrid1\n");
    fprintf(fp,"DataVersion: 1\n");
    fprintf(fp,"NumMembers: 1\n");
    fprintf(fp,"Flavors: [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21]\n");
    fprintf(fp,"OrderQCD: %i\n", steering_.i_fit_order-1); // qcdnum notation LO=1,...; LHAPDF6 LO=0,...
    fprintf(fp,"FlavorScheme: %s\n", get_flavor_scheme());
    fprintf(fp,"ErrorType: %s\n", get_error_type());
    fprintf(fp,"XMin: %g\n", grid.xmin);
    fprintf(fp,"XMax: %g\n",1.0);
    fprintf(fp,"QMin: %g\n", sqrt(grid.q2min));
    fprintf(fp,"QMax: %g\n", sqrt(grid.q2max));
    fprintf(fp,"MZ: %g\n", boson_masses_.mz);
    fprintf(fp,"MUp: 0\n");
    fprintf(fp,"MDown: 0\n");
    fprintf(fp,"MStrange: 0\n");
    fprintf(fp,"MCharm: %g\n", steering_.hf_mass[0]);
    fprintf(fp,"MBottom: %g\n", steering_.hf_mass[1]);
    fprintf(fp,"MTop: %g\n", steering_.hf_mass[2]);
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
    fprintf(fp,"]");

    fclose(fp);
} //}}}

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
  if(steering_.dobands) 
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
