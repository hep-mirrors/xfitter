#ifdef __cplusplus
extern "C" {
#endif
  void TRinit(double mc, double mb, int Nf, double Lam);
  void TRsfun(double x, double QQ, double SF[]);
  void TRgetGcontr(double SFg[]);
#ifdef __cplusplus
}
#endif
