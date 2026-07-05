double F2(double x, double Q2, double mu2, int NF, int orderAlphas, bool useparam, double EPSABS, double EPSREL, int ITER);
double FL(double x, double Q2, double mu2, int NF, int orderAlphas, bool useparam, double EPSABS, double EPSREL, int ITER);
double F2H(double x, double Q2, double mu2, int NF, int orderAlphas, double m2, int HQidx, double EPSABS, double EPSREL, int ITER);
double FLH(double x, double Q2, double mu2, int NF, int orderAlphas, double m2, int HQidx, double EPSABS, double EPSREL, int ITER);

//extern "C" {
  namespace PDF {
    double alphasQ2(const double Q2) {
      //printf("alphas Q2 = %f\n", Q2);fflush(stdout);
      return alphas_wrapper_(sqrt(Q2));
    }
    double xfxQ2(const int pdgid, const double x, const double Q2) {
      double pdf[13];
      const int id = (pdgid == 21) ? 0 : pdgid;
      //return pdf_ixfxq_wrapper_(id, x, sqrt(Q2));
      pdf_xfxq_wrapper_(x, sqrt(Q2), pdf);
      return pdf[id + 6];
    }
  }
  namespace HQ {
    //extern double t_interp_H2gg3reg_defaultvalue;
    extern std::function<double(double, double)> t_interp_H2gg3reg;
    // {
    //  return 0.5;
    //};
  }
  namespace C2HSCALE {
    void initTables();
  }
  namespace C2HINTERP {
    void initTables();
  }
//}
