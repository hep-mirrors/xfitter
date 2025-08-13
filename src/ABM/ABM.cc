#include "ABM.h"
#include "hf_errlog.h"

extern "C" {
  double f2qcd_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2);
  double flqcd_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2);
  double f3qcd_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2);
  double f2charm_ffn_(const double& xb, const double& q2, const int& nq);
  double flcharm_ffn_(const double& xb, const double& q2, const int& nq);
  double f2nucharm_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2, const int& nq);
  double ftnucharm_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2, const int& nq);
  double f3nucharm_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2, const int& nq);
  //void sf_abkm_wrap_(const double& x, const double& q2,
  //                   const double& f2abkm, const double& flabkm, const double& f3abkm,
  //                   const double& f2cabkm, const double& flcabkm, const double& f3cabkm,
  //                   const double& f2babkm, const double& flbabkm, const double& f3babkm,
  //                   const int& ncflag, const double& charge, const double& polar,
  //                   const double& sin2thw, const double& cos2thw, const double& MZ, const int& nt=1);
  //double numufcalflux_(const double& e);
  void abkm_set_input_(const int& kschemepdfin, const int& kordpdfin,
                       const double& rmass8in, const double& rmass10in, const int& msbarmin,
                       const double& hqscale1in, const double& hqscale2in, const int& flord);
  void abkm_set_input_full_(const int& kschemepdfin, const int& kordpdfin,
    const int& kordhqin, const int& kordf2in, const int& kordflin, const int& kordf3in, const int& kordalpsin,
    const double& rmass8in, const double& rmass10in, const int& msbarmin, const double& hqscale1in, const double& hqscale2in, const int& hqnonsin);
  int abkm_change_order_(const int order);
  void initgridconst_();
  void pdffillgrid_();
  
  struct COMMON_masses
  {
    double rmass[150];
    double rmassp[150];
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
}
  
namespace abm {
  void set_input(const int kschemepdfin, const int kordpdfin, 
                 const double rmass8in, const double rmass10in, 
                 const int msbarmin, double hqscale1in, const double hqscale2in, 
                 const int flord, const bool hqnons/* = false*/,
                 int kordpdf/* = -1*/, int kordhq/* = -1*/, int kordf2/* = -1*/, 
                 int kordfl/* = -1*/, int kordf3/* = -1*/, int kordalps/* = -1*/
                ) {
    // by default set same order for pdf, light, heavy quarks
    if (kordpdf == -1) kordpdf = kordpdfin;
    if (kordhq == -1)  kordhq  = kordpdfin;
    if (kordf2 == -1) kordf2 = kordpdfin;
    if (kordfl == -1) kordfl = kordpdfin + flord; // follow recommendation of Sergey Alekhin for FL
    if (kordf3 == -1) kordf3 = kordpdfin;
    if (kordalps == -1) kordalps = kordpdfin;
    // 10.10.2017 Discussion with Sergey Alekhin:
    // The parameter HQNONS drives the nonsinglet contribution to the charm production.
    // It is infrared unsafe in the NNLO therefore there are pro and contra for including it and it is up to user.
    // In ABMP16 fit it was set to .false.
    // (makes small difference which reaches few % only at highest Q2 of the charm HERA data and is negligible for practical purposes)
    abkm_set_input_(kschemepdfin, kordpdfin, rmass8in, rmass10in,msbarmin, hqscale1in, hqscale2in, flord);
  }

  void initgridconst() {
    initgridconst_();
  }

  void pdffillgrid() {
    pdffillgrid_();
  }

  double calc_point_strfun(const SFproc prc, const SFtype ftype, const SFflav flav, const double q2, const double x, const int order, const int charge, const double polar, const double sin2thw, const double mz, const double* f2c_cc) {
    static constexpr int nt = 1; // proton
    int order_old;
    if (order >= 0) {
      order_old = abkm_change_order_(order);
    }
    auto reset_order_and_return = [&](const double val) {
      if (order >= 0) {
        abkm_change_order_(order_old);
      }
      return val;
    };
    switch (prc) {
      case SFproc::nc:
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
                return reset_order_and_return(f2qcd_(3, nt, 22, x, q2) + facgz * PZ * f2qcd_(3, nt, 25, x, q2) + faczz * PZ * PZ * f2qcd_(3, nt, 23, x, q2));
              case SFtype::fl:
                return reset_order_and_return(flqcd_(3, nt, 22, x, q2) + facgz * PZ * flqcd_(3, nt, 25, x, q2) + faczz * PZ * PZ * flqcd_(3, nt, 23, x, q2));
              case SFtype::f3:
                return reset_order_and_return(-1 * charge * (facgzf3 * PZ * f3qcd_(3, nt, 25, x, q2) + faczzf3 * PZ * PZ * f3qcd_(3, nt, 23, x, q2)));
            }
          }
          case SFflav::c:
            switch (ftype) {
              case SFtype::f2:
                return reset_order_and_return(f2charm_ffn_(x, q2, 8));
              case SFtype::fl:
                return reset_order_and_return(flcharm_ffn_(x, q2, 8));
              case SFtype::f3:
                return reset_order_and_return(0.);
            }
          case SFflav::b:
            switch (ftype) {
              case SFtype::f2:
                return reset_order_and_return(f2charm_ffn_(x, q2, 10));
              case SFtype::fl:
                return reset_order_and_return(flcharm_ffn_(x, q2, 10));
              case SFtype::f3:
                return reset_order_and_return(0.);
            }
        }
      case SFproc::cc:
        const int ni = 24; // CC
        const int nb = charge > 0 ? 6 : 7;
        switch (flav) {
          case SFflav::l:
            switch (ftype) {
              case SFtype::f2:
                return reset_order_and_return(f2qcd_(nb, nt, ni, x, q2) / 2.);
              case SFtype::fl:
                return reset_order_and_return(flqcd_(nb, nt, ni, x, q2) / 2.);
              case SFtype::f3:
                return reset_order_and_return(f3qcd_(nb, nt, ni, x, q2) / 2.);
            }
          case SFflav::c:
            switch (ftype) {
              case SFtype::f2:
                return reset_order_and_return(f2nucharm_(nb, nt, ni, x, q2, 8) / 2.);
              case SFtype::fl: {
                double f2c_calc = 0.;
                if (f2c_cc) {
                  f2c_calc = *f2c_cc;
                }
                else {
                  f2c_calc = f2nucharm_(nb, nt, ni, x, q2, 8) / 2.;
                }
                return reset_order_and_return(f2c_calc - ftnucharm_(nb, nt, ni, x, q2, 8) / 2.);
              }
              case SFtype::f3:
                return reset_order_and_return(f3nucharm_(nb, nt, ni, x, q2, 8) / 2.);
            }
          case SFflav::b:
            return reset_order_and_return(0.);
        }
    }
    hf_errlog(28022501, "F: Unsupported structure function process, type or flavour");
  }
}