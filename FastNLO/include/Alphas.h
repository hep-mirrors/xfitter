#ifndef __ALPHAS_H
#define __ALPHAS_H

class Alphas {

public:
   ~Alphas();

   // initializations
   static void SetMz(double Mz) {
      fMz = Mz;
   };
   static double GetMz() {
      return fMz;
   };
   static void SetAlphasMz(double alphas) {
      fAlphasMz = alphas;
   };
   static double GetAlphasMz() {
      return fAlphasMz;
   };
   static void SetNf(int nf) {
      fNf = nf;
   };
   static int GetNf() {
      return fNf;
   };
   static void SetNLoop(int nLoop) {
      fnLoop = nLoop;
   };
   static int GetNLoop() {
      return fnLoop;
   };
   static void SetFlavorMatchingOn(bool FlavorMatching) {
      bFlavorMatching = FlavorMatching;
   };
   static bool GetFlavorMatchingOn() {
      return bFlavorMatching;
   };
   static void SetFlavorMatchingThresholds(double th1, double th2, double th3, double th4, double th5, double th6);
   static void GetFlavorMatchingThresholds(double& th1, double& th2, double& th3, double& th4, double& th5, double& th6);

   // Getters for Alphas at scale mu
   static double CalcAlphasMu(double mu, double alphasMz = 0, int nLoop = 0, int nFlavors = 0);
   static double CalcAlphasMuFixedNf(double mu, int nf) {                // calculate alpha_s as scale mu for fixed number of flavors nf. Ignore flavor matching thresholds.
      return CalcAlphasMu(mu, fAlphasMz, fnLoop, nf);
   };

   static int CalcNf(double mu);
   static void PrintInfo();

private:
   static Alphas* instance;
   Alphas();

   static double FBeta(double alphasMz , int nLoop , int nf);

public:
   static double fMz;                   // mass of Z0, which is the nominal scale here
   static double fAlphasMz;             // alpha_s at starting scale of Mz
   static int fNf;                      // MAXIMUM number of active flavours. e.g. at low scales mu, number of flavors is calculated with respecting flavor thresholds if FlavorMatching is ON.
   static int fnLoop;                   // n-loop solution of the RGE
   static bool bFlavorMatching;         // switch flaovr matching on or off
   static double fTh[6];                // flavor thresholds (quark masses)

};

#endif // __ALPHAS_H
