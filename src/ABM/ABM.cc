#include "ABM.h"
#include "xfitter_cpp_base.h"

extern "C" {
  // PDFs
  void initgridconst_();
  void pdffillgrid_();
  // structure functions
  double f2qcd_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2);
  double flqcd_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2);
  double f3qcd_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2);
  double f2charm_ffn_(const double& xb, const double& q2, const int& nq);
  double flcharm_ffn_(const double& xb, const double& q2, const int& nq);
  double f2nucharm_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2, const int& nq);
  double ftnucharm_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2, const int& nq);
  double f3nucharm_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2, const int& nq);
  
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

  double calc_point_strfun(const SFproc proc, const SFtype ftype, const SFflav flav, const double q2, const double x, 
    const int order, const int orderDefault, const int orderHQ, const int orderFL, const bool msbar, 
    const int charge, const double sin2thw, const double polar, const double mz, const double* f2c/*=nullptr*/) {
    if (proc == SFproc::nc) {
      return calc_point_strfun_NC(ftype, flav, q2, x, order, orderDefault, orderHQ, orderFL, msbar, charge, sin2thw, polar, mz);
    }
    else if (proc == SFproc::cc) {
      return calc_point_strfun_CC(ftype, flav, q2, x, order, orderDefault, orderHQ, orderFL, msbar, charge, f2c);
    }
    else {
      hf_errlog(19082501, "F: Unsupported process");
      return 0.; // avoid warning
    }
  }

  double calc_point_strfun_NC(const SFtype ftype, const SFflav flav, const double q2, const double x, 
    const int order, const int orderDefault, const int orderHQin, const int orderFLin, const bool msbar, 
    const int charge, const double sin2thw, const double polar, const double mz) {
    int orderALL = (order >= 0) ? order : orderDefault;
    int orderHQ = (order >= 0) ? order : orderHQin;
    int orderFL = (order >= 0) ? order : orderFLin;
    // Take the 3-flavour scheme as a default
    int kschemepdfin = 0;
    abm::set_scheme_and_order(kschemepdfin, orderALL, msbar, orderFL, orderHQ);
    static constexpr int nt = 1; // proton
    switch (flav) {
      case SFflav::l: {
        static constexpr double eleAxial = -0.5;
        const double eleVec = -0.5 + 2 * sin2thw;
        const double facgz = - eleVec - charge * polar * eleAxial;
        const double faczz = eleVec * eleVec + eleAxial * eleAxial + 2 * charge * polar * eleAxial * eleVec;
        const double facgzf3 = -1 * eleAxial - charge * polar * eleVec;
        const double faczzf3 = 2 * eleAxial * eleVec + charge * polar * (eleVec * eleVec + eleAxial * eleAxial);
        const double PZ = 1. / (4 * sin2thw * (1 - sin2thw) * (1 + mz * mz / q2));
        switch (ftype) {
          case SFtype::f2:
            return f2qcd_(3, nt, 22, x, q2) + facgz * PZ * f2qcd_(3, nt, 25, x, q2) + faczz * PZ * PZ * f2qcd_(3, nt, 23, x, q2);
          case SFtype::fl:
            return flqcd_(3, nt, 22, x, q2) + facgz * PZ * flqcd_(3, nt, 25, x, q2) + faczz * PZ * PZ * flqcd_(3, nt, 23, x, q2);
          case SFtype::f3:
            return -1 * charge * (facgzf3 * PZ * f3qcd_(3, nt, 25, x, q2) + faczzf3 * PZ * PZ * f3qcd_(3, nt, 23, x, q2));
        }
      }
      case SFflav::c:
        switch (ftype) {
          case SFtype::f2:
            return f2charm_ffn_(x, q2, 8);
          case SFtype::fl:
            return flcharm_ffn_(x, q2, 8);
          case SFtype::f3:
            return 0.;
        }
      case SFflav::b:
        switch (ftype) {
          case SFtype::f2:
            return f2charm_ffn_(x, q2, 10);
          case SFtype::fl:
            return flcharm_ffn_(x, q2, 10);
          case SFtype::f3:
            return 0.;
        }
    }
    hf_errlog(28022501, "F: Unsupported structure function type or flavour");
    return 0; // avoid warning
  }

  double calc_point_strfun_CC(const SFtype ftype, const SFflav flav, const double q2, const double x, 
      const int order, const int orderDefault, const int orderHQin, const int orderFLin, const bool msbar, 
      const int charge, const double* f2c/*=nullptr*/) {
    if (flav == abm::SFflav::b) {
      return 0;
    }
    int orderALL = (order >= 0) ? order : orderDefault;
    int orderHQ = (order >= 0) ? order : orderHQin;
    int orderFL = (order >= 0) ? order : orderFLin;
    // Take the 3-flavour scheme as a default
    int kschemepdfin = 0;
    abm::set_scheme_and_order(kschemepdfin, orderALL, msbar, orderFL, orderHQ);
    static constexpr int nt = 1; // proton
    static constexpr int ni = 24; // CC
    const int nb = charge > 0 ? 6 : 7;
    switch (flav) {
      case abm::SFflav::l:
        switch (ftype) {
          case abm::SFtype::f2:
            return f2qcd_(nb, nt, ni, x, q2) / 2.;
          case abm::SFtype::fl:
            return flqcd_(nb, nt, ni, x, q2) / 2.;
          case abm::SFtype::f3:
            return f3qcd_(nb, nt, ni, x, q2) / 2.;
        }
      case abm::SFflav::c:
        double f2c_calc = 0.;
        switch (ftype) {
          case abm::SFtype::f2:
            return f2nucharm_(nb, nt, ni, x, q2, 8) / 2.;
          case abm::SFtype::fl:
            if (f2c) {
              f2c_calc = *f2c;
            }
            else {
              f2c_calc = f2nucharm_(nb, nt, ni, x, q2, 8) / 2.;
            }
            return f2c_calc - ftnucharm_(nb, nt, ni, x, q2, 8) / 2.;
          case abm::SFtype::f3:
            return f3nucharm_(nb, nt, ni, x, q2, 8) / 2.;
        }
    }
    hf_errlog(28022501, "F: Unsupported structure function type or flavour");
    return 0; // avoid warning
  }
}
