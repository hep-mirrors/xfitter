double F2(double x, double Q2, double mu2, int NF, int orderAlphas, bool useparam, double EPSABS, double EPSREL, int ITER);
double FL(double x, double Q2, double mu2, int NF, int orderAlphas, bool useparam, double EPSABS, double EPSREL, int ITER);
double F2H(double x, double Q2, double mu2, int NF, int orderAlphas, double m2, int HQidx, double EPSABS, double EPSREL, int ITER);
double FLH(double x, double Q2, double mu2, int NF, int orderAlphas, double m2, int HQidx, double EPSABS, double EPSREL, int ITER);

//extern "C" {
  namespace PDF {
    double alphasQ2(const double Q2) {
      return alphas_wrapper_(sqrt(Q2));
    }
    double xfxQ2(const int pdgid, const double x, const double Q2) {
      const int id = (pdgid == 21) ? 0 : pdgid;
      return pdf_ixfxq_wrapper_(id, x, sqrt(Q2));
    }
  }
  namespace HQ {
    extern std::function<double(double, double)> t_interp_H2gg3reg;
  }
//}
