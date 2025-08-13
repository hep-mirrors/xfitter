#include "ABM.h"

extern "C" {
  void initgridconst_();
  void pdffillgrid_();
  
  struct COMMON_masses
  {
    double rmass[150];
    double rmassp[50];
    double rcharge[150];
  };
  extern COMMON_masses masses_;
  
  struct COMMON_gridset
  {
    double delx1,delx2,delxp,dels1[8],dels2[8],xlog1,xlog2,x1,q2ini[8],q2min,q2max,xbmin,xbmax;
    int nxmgrid,nxpgrid,nspgrid,nsmgrid,khalf;
  };
  extern COMMON_gridset gridset_;
  
  struct COMMON_constants_abkm
  {
    double pi,alpha,alphady,rmpr,gfer2,sintc,sintw2,rmw,rmz,rgz,ckm[9],ckm2[9];
  };
  extern COMMON_constants_abkm constants_abkm_;

  struct COMMON_foralpsrenorm
  {
    double q20,q2rep,q2s,q20alphas,alphas0,alpsz,alpss,alpsc,alpsb,alpst,tscale,rscale,fscale,hqscale1,hqscale2;
    int nfeff,kordalps,kfeff,kordhq,kordf2,kordfl,kordf3;
    double alsmz;
  };
  extern COMMON_foralpsrenorm foralpsrenorm_;

  struct COMMON_forpdfset
  {
    int npdftot,kordpdf,kschemepdf,kpdfset,kordkernel;
  };
  extern COMMON_forpdfset forpdfset_;


  struct COMMON_forschemedef
  {
    double ddnnlohq;
    int msbarm,hqnons;
    double bmsnfopt,bmsnnlo,vloop;
  };
  extern COMMON_forschemedef forschemedef_;
}

#include <cstdio>
namespace abm {
  void initgridconst() {
    initgridconst_();
    masses_.rcharge[7] = 0.6666666;
    masses_.rcharge[9] = 0.3333333;
  }

  void pdffillgrid() {
    pdffillgrid_();
  }

  void update_ckm_matrix(const std::vector<const double*>& ckm) {
    for(int i = 0; i < 9; i++) {
      constants_abkm_.ckm[i] = *ckm[i];
    }
  }

  void set_scheme_and_order(const int kschemepdf, const int kordpdf, 
                 const int msbarm, const int flord, 
                 int kordhq/* = -1*/, int kordf2/* = -1*/, int kordfl/* = -1*/, int kordf3/* = -1*/, int kordalps/* = -1*/,
                 const bool hqnons/* = false*/
                 ) {
    forpdfset_.kschemepdf = kschemepdf;
    forschemedef_.msbarm = msbarm;
    forpdfset_.kordpdf = kordpdf;
    // by default set same order for pdf, light, heavy quarks
    foralpsrenorm_.kordhq = (kordhq >= 0) ? kordhq : kordpdf;
    foralpsrenorm_.kordf2 = (kordf2 >= 0) ? kordf2 : kordpdf;
    foralpsrenorm_.kordfl = (kordfl >= 0) ? kordfl : kordpdf + flord; // follow recommendation of Sergey Alekhin for FL
    foralpsrenorm_.kordf3 = (kordf3 >= 0) ? kordf3 : kordpdf;
    foralpsrenorm_.kordalps = (kordalps >= 0) ? kordalps : kordpdf;
    // 10.10.2017 Discussion with Sergey Alekhin:
    // The parameter HQNONS drives the nonsinglet contribution to the charm production.
    // It is infrared unsafe in the NNLO therefore there are pro and contra for including it and it is up to user.
    // In ABMP16 fit it was set to .false.
    // (makes small difference which reaches few % only at highest Q2 of the charm HERA data and is negligible for practical purposes)
    forschemedef_.hqnons = hqnons;
    //abkm_set_input_full_(kschemepdf, kordpdf, kordhq, kordf2, kordfl, kordf3, kordalps, rmass8, rmass10, msbarm, hqscale1, hqscale2, flord);
  }

  void set_hq_masses(const double mc, const double mb) {
    masses_.rmass[7] = mc;
    masses_.rmass[9] = mb;
  }

  void set_hq_scales(const double hqscale1, const double hqscale2) {
    foralpsrenorm_.hqscale1 = hqscale1;
    foralpsrenorm_.hqscale2 = hqscale2;
  }

  void set_xbmin(const double val) {
    gridset_.xbmin = val;
  }

  void set_xbmax(const double val) {
    gridset_.xbmax = val;
  }
}
