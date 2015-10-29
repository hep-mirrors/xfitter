#include <vector>

class fastNLOReader;

extern "C" void evolvepdf_(const double& x, const double& Q, double* xf);

class HoppetInterface {
   public:
      static void InitHoppet(fastNLOReader&);
      static std::vector<double> GetSpl(double, double);
      static std::vector<double> GetXFX(double, double);
      static double EvolveAlphas(double);
      static bool IsInitialized;
      // ---- Alphas vars ---- //
      static double fAlphasMz;
      static double fMz;
      static int fnFlavor;
      static int fnLoop;
      static double QMass[6];
      // ____ //
   private:
      static void StartHoppet();
      static void LHAsub(const double &, const double &, double*);
      static fastNLOReader *fnlo;
};
