#include<string.h>
#include<stdio.h>
#include<malloc.h>
#include<math.h>
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

void save_info();
void save_data();
char* get_flavor_scheme();
char* get_error_type();


// convert fortran string to C string
char* sfix(char* fstr,int size) { //{{{
  char* cstr=(char*)malloc(sizeof(char)*size);
  memcpy(cstr,fstr,size);
  cstr[size-1]='\0';
  char* end;
  while( end=rindex(cstr,' ') ) cstr[end-cstr]='\0';
  return cstr;
} //}}}


// fortran interface to store pdf in LHAPDF6 format
void print_lhapdf6_(){ //{{{
  mkdir("output/herapdf",0755);
  save_info();
  save_data();
} //}}}


// store pdf data to herapdf_0000.dat yaml file
void save_data(){ //{{{
  extern double xfrmix_(int *);
  extern double qfrmiq_(int *);
  extern double fvalij_(int *,int *,int *,int *,int *);
  int nx,nq,null,iq,ix,i;
  double xmin,xmax,qmin,qmax,q2;
  int pdg_flavours[]={-5,-4,-3,-2,-1,1,2,3,4,5,21};
  int qcdnum_flavours[]={-5,-4,-3,-2,-1,1,2,3,4,5,0};
  grpars_(&nx,&xmin,&xmax,&nq,&qmin,&qmax,&null);
  FILE* fp;
  grpars_(&nx,&xmin,&xmax,&nq,&qmin,&qmax,&null);
  if((fp=fopen("output/herapdf/herapdf_0000.dat","w"))==NULL) puts("Cannot open file.");
    fprintf(fp,"PdfType: central\n");
    fprintf(fp,"Format: lhagrid1\n");
    fprintf(fp,"---\n");
    for(ix=1;ix<=nx;ix++) 
      fprintf(fp,"%e ",xfrmix_(&ix)); 
      fprintf(fp,"\n");
    for(iq=1;iq<=nq;iq++) 
      fprintf(fp,"%e ",sqrt(qfrmiq_(&iq))); 
      fprintf(fp,"\n");
    for(i=0;i<sizeof(pdg_flavours)/sizeof(int);i++) 
      fprintf(fp,"%i ",pdg_flavours[i]);
      fprintf(fp,"\n");
    for(ix=1;ix<=nx;ix++)
      for(iq=1;iq<=nq;iq++){
          for(i=0;i<sizeof(qcdnum_flavours)/sizeof(int);i++) 
            fprintf(fp,"%e ",fvalij_(&steering_.ipdfset,&qcdnum_flavours[i],&ix,&iq,&null));
          fprintf(fp,"\n");
          }
    fprintf(fp,"---\n");
  fclose(fp);
} //}}}


// store pdf information to herapdf.info yaml file
void save_info() { //{{{

  extern double hf_get_alphas_(double *);
  extern double qfrmiq_(int *);
  int nx,nq,inull,iq,as_order;
  double xmin,xmax,q2min,q2max,q2,alphas,dnull;
  double mz2=boson_masses_.mz*boson_masses_.mz;
  FILE* fp;
  
  grpars_(&nx,&xmin,&xmax,&nq,&q2min,&q2max,&inull);
  
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
    fprintf(fp,"XMin: %g\n", xmin);
    fprintf(fp,"XMax: %g\n",1.0);
    fprintf(fp,"QMin: %g\n", sqrt(q2min));
    fprintf(fp,"QMax: %g\n", sqrt(q2max));
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
      for(iq=1;iq<=nq;iq++) {
        fprintf(fp,"%g",sqrt(qfrmiq_(&iq)));
        if(iq!=nq) fprintf(fp,", ");
      }
      fprintf(fp,"]\n");

    fprintf(fp,"AlphaS_Vals: [");
      for(iq=1;iq<=nq;iq++) {
        q2=qfrmiq_(&iq);
        fprintf(fp,"%g",hf_get_alphas_(&q2));
        if(iq!=nq) fprintf(fp,", ");
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
